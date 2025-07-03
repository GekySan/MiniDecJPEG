#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#define NOMINMAX
#include <windows.h>
#include <shellapi.h>

#include "LocaleInitializer.hpp"

/*
 0  1  5  6 14 15 27 28 |  0  1  2  3  4  5  6  7
 2  4  7 13 16 26 29 42 |  8  9 10 11 12 13 14 15
 3  8 12 17 25 30 41 43 | 16 17 18 19 20 21 22 23
 9 11 18 24 31 40 44 53 | 24 25 26 27 28 29 30 31
10 19 23 32 39 45 52 54 | 32 33 34 35 36 37 38 39
20 22 33 38 46 51 55 60 | 40 41 42 43 44 45 46 47
21 34 37 47 50 56 59 61 | 48 49 50 51 52 53 54 55
35 36 48 49 57 58 62 63 | 56 57 58 59 60 61 62 63
*/
constexpr std::array<int, 64> kZigzagOrder = {
     0,  1,  8, 16,  9,  2,  3, 10,
    17, 24, 32, 25, 18, 11,  4,  5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13,  6,  7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63
};

struct Rgb
{
    uint8_t r, g, b;
};

// Convert YCbCr color space to RGB
// https://en.wikipedia.org/wiki/YCbCr#JPEG_conversion
Rgb YcbcrToRgb(double y, double cb, double cr)
{
    double red = y + 1.402 * (cr - 128.0);
    double green = y - 0.344136 * (cb - 128.0) - 0.714136 * (cr - 128.0);
    double blue = y + 1.772 * (cb - 128.0);

    auto clamp = [](double val)
        {
            return static_cast<uint8_t>(std::max(0.0, std::min(255.0, val)));
        };

    return { clamp(red), clamp(green), clamp(blue) };
}

// Transformée en cosinus discrète inverse (ITDC) sur un bloc de 64 éléments
std::vector<double> PerformIdct(const std::vector<double>& coefficients)
{
    std::vector<double> output(64, 0.0);
    constexpr double kPi = 3.14159265358979323846;

    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 8; ++j)
        {
            double sum = 0.0;
            for (int u = 0; u < 8; ++u)
            {
                for (int v = 0; v < 8; ++v)
                {
                    double coeffU = (u == 0) ? 1.0 / std::sqrt(2.0) : 1.0;
                    double coeffV = (v == 0) ? 1.0 / std::sqrt(2.0) : 1.0;
                    double term = coefficients[v * 8 + u] * coeffU * coeffV * std::cos((2 * j + 1) * u * kPi / 16.0) * std::cos((2 * i + 1) * v * kPi / 16.0);
                    sum += term;
                }
            }
            output[i * 8 + j] = sum / 4.0;
        }
    }
    return output;
}

// Décoder un «nombre codé» en Huffman
int DecodeNumber(int category, int bits)
{
    if (category == 0)
    {
        return 0;
    }
    int limit = 1 << (category - 1);
    if (bits < limit)
    {
        return bits - ((1 << category) - 1);
    }
    else
    {
        return bits;
    }
}

// Lire une valeur 16 bits en big-endian
uint16_t ReadU16BigEndian(std::span<const uint8_t> data)
{
    return (static_cast<uint16_t>(data[0]) << 8) | data[1];
}

class BitStreamReader
{
private:
    std::span<const uint8_t> m_data;
    size_t m_bytePos = 0;
    int m_bitPos = 8;
    uint8_t m_currentByte = 0;

public:
    explicit BitStreamReader(std::span<const uint8_t> data) : m_data(data) {}

    std::optional<int> ReadBit()
    {
        if (m_bytePos >= m_data.size() && m_bitPos >= 8)
        {
            return std::nullopt;
        }

        if (m_bitPos >= 8)
        {
            m_bitPos = 0;
            if (m_bytePos >= m_data.size())
            {
                return std::nullopt;
            }
            m_currentByte = m_data[m_bytePos++];
            if (m_currentByte == 0xFF)
            {
                if (m_bytePos >= m_data.size())
                {
                    return std::nullopt;
                }
                uint8_t next_byte = m_data[m_bytePos];
                if (next_byte == 0x00)
                {
                    m_bytePos++;
                }
                else if (next_byte >= 0xD0 && next_byte <= 0xD7)
                {
                    // Marqueurs RST, ne rien faire
                }
            }
        }
        int bit = (m_currentByte >> (7 - m_bitPos)) & 1;
        m_bitPos++;
        return bit;
    }

    std::optional<int> ReadBits(int count)
    {
        int value = 0;
        for (int i = 0; i < count; ++i)
        {
            auto bitOpt = ReadBit();
            if (!bitOpt)
            {
                return std::nullopt;
            }
            value = (value << 1) | *bitOpt;
        }
        return value;
    }
};

class HuffmanTable
{
private:
    struct Node;
    using NodePtr = std::unique_ptr<Node>;
    struct Node
    {
        std::variant<uint8_t, std::array<NodePtr, 2>> content;
    };

    NodePtr m_root;

    void BuildTree(std::span<const uint8_t> lengths, std::span<const uint8_t> values)
    {
        m_root = std::make_unique<Node>();

        uint32_t code = 0;
        int valueIndex = 0;

        for (int i = 0; i < 16; ++i)
        {
            int numCodesAtLength = lengths[i];
            for (int j = 0; j < numCodesAtLength; ++j)
            {
                if (valueIndex >= values.size())
                {
                    break;
                }

                Node* p_node = m_root.get();
                int codeLength = i + 1;

                for (int bitIndex = codeLength - 1; bitIndex >= 0; --bitIndex)
                {
                    int bit = (code >> bitIndex) & 1;

                    auto* branch = std::get_if<std::array<NodePtr, 2>>(&p_node->content);
                    if (!branch)
                    {
                        p_node->content = std::array<NodePtr, 2>{};
                        branch = &std::get<std::array<NodePtr, 2>>(p_node->content);
                    }

                    if (!(*branch)[bit])
                    {
                        (*branch)[bit] = std::make_unique<Node>();
                    }

                    if (bitIndex == 0)
                    {
                        (*branch)[bit]->content = values[valueIndex];
                    }
                    else
                    {
                        p_node = (*branch)[bit].get();
                    }
                }
                code++;
                valueIndex++;
            }
            code <<= 1;
        }
    }

public:
    HuffmanTable(std::span<const uint8_t> lengths, std::span<const uint8_t> values)
    {
        BuildTree(lengths, values);
    }

    std::optional<uint8_t> Decode(BitStreamReader& stream) const
    {
        const Node* p_node = m_root.get();
        while (true)
        {
            if (auto val_ptr = std::get_if<uint8_t>(&p_node->content))
            {
                return *val_ptr;
            }

            auto bitOpt = stream.ReadBit();
            if (!bitOpt)
            {
                return std::nullopt;
            }
            int bit = *bitOpt;

            const auto* branch = std::get_if<std::array<NodePtr, 2>>(&p_node->content);
            if (!branch || bit >= 2 || !(*branch)[bit])
            {
                return std::nullopt;
            }
            p_node = (*branch)[bit].get();
        }
    }
};

struct DecodedImage
{
    int width = 0;
    int height = 0;
    std::vector<uint8_t> pixel_data;
};

class JpegDecoder
{
private:
    std::vector<uint8_t> m_data;
    std::map<int, std::vector<uint8_t>> m_quantTables;
    std::map<std::string, std::map<int, HuffmanTable>> m_huffmanTables;

    struct Component
    {
        int h_sampling, v_sampling;
        int quant_id;
        int dc_huff_id, ac_huff_id;
    };
    std::map<int, Component> m_components;

    int m_width = 0, m_height = 0;
    int m_maxHSampling = 0, m_maxVSampling = 0;

    int m_ss = 0, m_se = 0, m_ah_al = 0;

    void ParseDqt(std::span<const uint8_t> segmentData)
    {
        std::cout << "  -> Table DQT trouvée\n";
        size_t pos = 0;
        while (pos < segmentData.size())
        {
            int info = segmentData[pos++];
            int tableId = info & 0x0F;
            m_quantTables[tableId].assign(segmentData.begin() + pos, segmentData.begin() + pos + 64);
            pos += 64;
        }
    }

    void ParseDht(std::span<const uint8_t> segmentData)
    {
        std::cout << "  -> Table DHT trouvée\n";
        size_t pos = 0;
        while (pos < segmentData.size())
        {
            int info = segmentData[pos];
            pos++;
            int tableClass = info >> 4;
            int tableId = info & 0x0F;

            auto lengths = segmentData.subspan(pos, 16);
            pos += 16;

            int numValues = std::accumulate(lengths.begin(), lengths.end(), 0);
            auto values = segmentData.subspan(pos, numValues);
            pos += numValues;

            std::string typeStr = (tableClass == 0) ? "DC" : "AC";
            std::cout << "     -> Chargement de la table de Huffman : Classe=" << typeStr << ", ID=" << tableId << "\n";


            m_huffmanTables[typeStr].emplace(tableId, HuffmanTable(lengths, values));
        }
    }

    void ParseSof(std::span<const uint8_t> data)
    {
        m_height = ReadU16BigEndian(data.subspan(1, 2));
        m_width = ReadU16BigEndian(data.subspan(3, 2));
        std::cout << "     Dimensions : " << m_width << "x" << m_height << "\n";

        int numComponents = data[5];
        size_t pos = 6;
        for (int i = 0; i < numComponents; ++i)
        {
            int comp_id = data[pos];
            int sampling = data[pos + 1];
            m_components[comp_id] = {
                .h_sampling = sampling >> 4,
                .v_sampling = sampling & 0x0F,
                .quant_id = data[pos + 2]
            };
            m_maxHSampling = std::max(m_maxHSampling, m_components[comp_id].h_sampling);
            m_maxVSampling = std::max(m_maxVSampling, m_components[comp_id].v_sampling);
            pos += 3;
        }
    }

    void ParseSosHeader(std::span<const uint8_t> data)
    {
        int numComponentsInScan = data[0];
        size_t pos = 1;
        for (int i = 0; i < numComponentsInScan; ++i)
        {
            int comp_id = data[pos];
            if (m_components.count(comp_id))
            {
                int huff_ids = data[pos + 1];
                m_components[comp_id].dc_huff_id = huff_ids >> 4;
                m_components[comp_id].ac_huff_id = huff_ids & 0x0F;
            }
            pos += 2;
        }
        m_ss = data[pos];
        m_se = data[pos + 1];
        m_ah_al = data[pos + 2];
    }

    std::vector<int> DecodeBlock(BitStreamReader& stream, int componentId, int& previousDc)
    {
        std::vector<int> coeffs(64, 0);

        if (m_ss == 0)
        {
            const auto& dcHuffTable = m_huffmanTables.at("DC").at(m_components.at(componentId).dc_huff_id);
            auto categoryOpt = dcHuffTable.Decode(stream);
            if (!categoryOpt)
            {
                return coeffs;
            }

            int category = *categoryOpt;
            if (category > 0)
            {
                auto bitsOpt = stream.ReadBits(category);
                if (!bitsOpt)
                {
                    return coeffs;
                }
                int dcValue = DecodeNumber(category, *bitsOpt);
                coeffs[0] = dcValue + previousDc;
                previousDc = coeffs[0];
            }
            else
            {
                coeffs[0] = previousDc;
            }
        }

        if (m_se > 0)
        {
            if (m_huffmanTables.at("AC").find(m_components.at(componentId).ac_huff_id) == m_huffmanTables.at("AC").end())
            {
                return coeffs;
            }
            const auto& acHuffTable = m_huffmanTables.at("AC").at(m_components.at(componentId).ac_huff_id);

            for (int i = 1; i < 64; )
            {
                auto symbolOpt = acHuffTable.Decode(stream);
                if (!symbolOpt || *symbolOpt == 0x00)
                {
                    break;
                }

                int symbol = *symbolOpt;
                int zeros = symbol >> 4;
                i += zeros;
                int category = symbol & 0x0F;

                if (category > 0 && i < 64)
                {
                    auto bitsOpt = stream.ReadBits(category);
                    if (!bitsOpt)
                    {
                        break;
                    }
                    coeffs[i] = DecodeNumber(category, *bitsOpt);
                }
                i++;
            }
        }
        return coeffs;
    }

    void RenderMcu(int mcuX, int mcuY, const std::map<int, std::vector<std::vector<double>>>& mcuData, DecodedImage& outImage)
    {
        if (!mcuData.count(1) || !mcuData.count(2) || !mcuData.count(3))
        {
            return;
        }

        const auto& yBlocks = mcuData.at(1);
        const auto& cbBlocks = mcuData.at(2);
        const auto& crBlocks = mcuData.at(3);

        int mcuWidth = 8 * m_maxHSampling;
        int mcuHeight = 8 * m_maxVSampling;

        for (int y = 0; y < mcuHeight; ++y)
        {
            for (int x = 0; x < mcuWidth; ++x)
            {
                int pixelX = mcuX + x;
                int pixelY = mcuY + y;
                if (pixelX >= m_width || pixelY >= m_height)
                {
                    continue;
                }

                double yValue = yBlocks[(y / 8) * m_maxHSampling + (x / 8)][(y % 8) * 8 + (x % 8)];

                double cbXRatio = static_cast<double>(m_maxHSampling) / m_components.at(2).h_sampling;
                double cbYRatio = static_cast<double>(m_maxVSampling) / m_components.at(2).v_sampling;
                double cbValue = cbBlocks[0][(static_cast<int>(y / cbYRatio) % 8) * 8 + (static_cast<int>(x / cbXRatio) % 8)];

                double crXRatio = static_cast<double>(m_maxHSampling) / m_components.at(3).h_sampling;
                double crYRatio = static_cast<double>(m_maxVSampling) / m_components.at(3).v_sampling;
                double crValue = crBlocks[0][(static_cast<int>(y / crYRatio) % 8) * 8 + (static_cast<int>(x / crXRatio) % 8)];

                Rgb rgb = YcbcrToRgb(yValue, cbValue, crValue);

                size_t index = (static_cast<size_t>(pixelY) * m_width + pixelX) * 4;
                outImage.pixel_data[index] = rgb.b;
                outImage.pixel_data[index + 1] = rgb.g;
                outImage.pixel_data[index + 2] = rgb.r;
                outImage.pixel_data[index + 3] = 255;
            }
        }
    }

public:
    explicit JpegDecoder(const std::string& filepath)
    {
        std::ifstream file(filepath, std::ios::binary | std::ios::ate);
        if (!file.is_open())
        {
            throw std::runtime_error("Impossible d’ouvrir le fichier : " + filepath);

        }
        std::streamsize size = file.tellg();
        file.seekg(0, std::ios::beg);
        m_data.resize(size);
        if (!file.read(reinterpret_cast<char*>(m_data.data()), size))
        {
            throw std::runtime_error("Erreur lros de la lecture du fichier : " + filepath);
        }
    }

    DecodedImage Decode()
    {
        std::cout << "Démarrage du décodage...\n";
        size_t pos = 0;
        DecodedImage resultImage;

        while (pos < m_data.size())
        {
            if (m_data[pos] != 0xFF)
            {
                pos++;
                continue;
            }
            if (pos + 2 > m_data.size())
            {
                break;
            }
            uint16_t marker = ReadU16BigEndian({ &m_data[pos], 2 });
            pos += 2;

            if (marker == 0xFFD8)
            {
                continue;
            }
            if (marker == 0xFFD9)
            {
                break;
            }
            if (marker >= 0xFFD0 && marker <= 0xFFD9)
            {
                continue;
            }
            if (pos + 2 > m_data.size())
            {
                break;
            }
            uint16_t length = ReadU16BigEndian({ &m_data[pos], 2 });
            if (pos + length > m_data.size())
            {
                break;
            }
            std::span<const uint8_t> segmentData = { &m_data[pos + 2], static_cast<size_t>(length - 2) };

            if (marker == 0xFFDB)
            {
                ParseDqt(segmentData);
            }
            else if (marker == 0xFFC4)
            {
                ParseDht(segmentData);
            }
            else if (marker == 0xFFC0 || marker == 0xFFC2)
            {
                ParseSof(segmentData);
                resultImage.width = m_width;
                resultImage.height = m_height;
                resultImage.pixel_data.resize(static_cast<size_t>(m_width) * m_height * 4);
            }
            else if (marker == 0xFFDA)
            {
                ParseSosHeader(segmentData);
                size_t scanDataStart = pos + length;
                BitStreamReader stream({ &m_data[scanDataStart], m_data.size() - scanDataStart });

                std::map<int, int> previousDcValues = { {1, 0}, {2, 0}, {3, 0} };
                int mcuWidth = 8 * m_maxHSampling;
                int mcuHeight = 8 * m_maxVSampling;

                for (int mcuY = 0; mcuY < m_height; mcuY += mcuHeight)
                {
                    for (int mcuX = 0; mcuX < m_width; mcuX += mcuWidth)
                    {
                        std::map<int, std::vector<std::vector<double>>> mcuData;
                        for (auto const& [componentId, comp] : m_components)
                        {
                            std::vector<std::vector<double>> componentBlocks;
                            for (int v = 0; v < comp.v_sampling; ++v)
                            {
                                for (int h = 0; h < comp.h_sampling; ++h)
                                {
                                    std::vector<int> coeffsZigzag = DecodeBlock(stream, componentId, previousDcValues[componentId]);

                                    const auto& quantTable = m_quantTables.at(comp.quant_id);
                                    std::vector<double> dequantCoeffsZigzag(64);
                                    for (int i = 0; i < 64; ++i)
                                    {
                                        dequantCoeffsZigzag[i] = static_cast<double>(coeffsZigzag[i]) * quantTable[i];
                                    }

                                    std::vector<double> dezigzagCoeffs(64, 0.0);
                                    for (int i = 0; i < 64; ++i)
                                    {
                                        dezigzagCoeffs[kZigzagOrder[i]] = dequantCoeffsZigzag[i];
                                    }

                                    std::vector<double> pixelValuesRaw = PerformIdct(dezigzagCoeffs);

                                    std::vector<double> levelShifted(64);
                                    for (int i = 0; i < 64; ++i)
                                    {
                                        levelShifted[i] = pixelValuesRaw[i] + 128.0;
                                    }
                                    componentBlocks.push_back(levelShifted);
                                }
                            }
                            mcuData[componentId] = componentBlocks;
                        }
                        RenderMcu(mcuX, mcuY, mcuData, resultImage);
                    }
                }
            }
            pos += length;
        }

        std::cout << "Décodage terminé.\n";
        return resultImage;
    }
};

// Windows GUI

HBITMAP g_hBitmap = NULL;

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    switch (uMsg)
    {
    case WM_PAINT:
    {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);
        if (g_hBitmap)
        {
            HDC hdcMem = CreateCompatibleDC(hdc);
            HBITMAP hbmOld = (HBITMAP)SelectObject(hdcMem, g_hBitmap);
            BITMAP bm;
            GetObject(g_hBitmap, sizeof(bm), &bm);
            BitBlt(hdc, 0, 0, bm.bmWidth, bm.bmHeight, hdcMem, 0, 0, SRCCOPY);
            SelectObject(hdcMem, hbmOld);
            DeleteDC(hdcMem);
        }
        EndPaint(hwnd, &ps);
        return 0;
    }
    case WM_DESTROY:
    {
        if (g_hBitmap)
        {
            DeleteObject(g_hBitmap);
        }
        PostQuitMessage(0);
        return 0;
    }
    case WM_CLOSE:
    {
        DestroyWindow(hwnd);
        return 0;
    }
    default:
    {
        return DefWindowProc(hwnd, uMsg, wParam, lParam);
    }
    }
}

void CreateDebugConsole()
{
    if (AllocConsole())
    {
        FILE* pCout, * pCerr, * pCin;

        freopen_s(&pCout, "CONOUT$", "w", stdout);
        freopen_s(&pCerr, "CONOUT$", "w", stderr);
        freopen_s(&pCin, "CONIN$", "r", stdin);

        std::cout.clear();
        std::cerr.clear();
        std::cin.clear();
    }
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE, LPSTR, int nCmdShow)
{
    CreateDebugConsole();

    int argc;
    LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);
    if (argc < 2)
    {
        MessageBox(NULL, L"Veuillez fournir un chemin de fichier JPEG en argument.", L"Erreur d’utilisation", MB_OK | MB_ICONERROR);
        if (argv)
        {
            LocalFree(argv);
        }
        return 1;
    }

    char filepathMb[MAX_PATH];
    WideCharToMultiByte(CP_ACP, 0, argv[1], -1, filepathMb, MAX_PATH, NULL, NULL);
    std::string filepath(filepathMb);

    try
    {
        JpegDecoder decoder(filepath);
        DecodedImage image = decoder.Decode();

        if (image.width == 0 || image.height == 0)
        {
            MessageBoxA(NULL, "Échec du décodage de l’image ou image vide.", "Erreur de décodage", MB_OK | MB_ICONERROR);
            if (argv)
            {
                LocalFree(argv);
            }
            return 1;
        }

        BITMAPINFO bmi = { 0 };
        bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
        bmi.bmiHeader.biWidth = image.width;
        bmi.bmiHeader.biHeight = -image.height;
        bmi.bmiHeader.biPlanes = 1;
        bmi.bmiHeader.biBitCount = 32;
        bmi.bmiHeader.biCompression = BI_RGB;

        HDC hdcScreen = GetDC(NULL);
        g_hBitmap = CreateDIBitmap(hdcScreen, &bmi.bmiHeader, CBM_INIT, image.pixel_data.data(), &bmi, DIB_RGB_COLORS);
        ReleaseDC(NULL, hdcScreen);

        if (!g_hBitmap)
        {
            MessageBoxA(NULL, "Échec de la création du bitmap GDI.", "Erreur", MB_OK | MB_ICONERROR);
            if (argv)
            {
                LocalFree(argv);
            }
            return 1;
        }

        const wchar_t kClassName[] = L"JpegDecoderWindowClass";
        WNDCLASS wc = { 0 };
        wc.lpfnWndProc = WindowProc;
        wc.hInstance = hInstance;
        wc.lpszClassName = kClassName;
        wc.hCursor = LoadCursor(NULL, IDC_ARROW);
        wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
        RegisterClass(&wc);

        std::wstring windowTitle = L"Image decodee : " + std::wstring(argv[1]);
        if (argv)
        {
            LocalFree(argv);
            argv = nullptr;
        }

        HWND hwnd = CreateWindowEx(0, kClassName, windowTitle.c_str(), WS_OVERLAPPEDWINDOW,
            CW_USEDEFAULT, CW_USEDEFAULT, image.width + 16, image.height + 39,
            NULL, NULL, hInstance, NULL);

        if (hwnd == NULL)
        {
            return 0;
        }

        ShowWindow(hwnd, nCmdShow);
        UpdateWindow(hwnd);

        MSG msg = { 0 };
        while (GetMessage(&msg, NULL, 0, 0) > 0)
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
        return (int)msg.wParam;
    }
    catch (const std::exception& e)
    {
        MessageBoxA(NULL, e.what(), "Une erreur inattendue est survenue", MB_OK | MB_ICONERROR);
        if (argv)
        {
            LocalFree(argv);
        }
        return 1;
    }
}