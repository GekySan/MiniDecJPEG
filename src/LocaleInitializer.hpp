#pragma once
#pragma execution_character_set("utf-8")

#include <iostream>
#include <locale>
#include <windows.h>

struct LocaleInitializer {
    LocaleInitializer() {
        try {
            std::locale::global(std::locale("fr_FR.UTF-8"));
        }
        catch (const std::runtime_error& e) {
            std::cerr << "Locale spécifiée non disponible. Utilisation de la locale par défaut." << std::endl;
            std::cerr << "Message d'erreur : " << e.what() << std::endl;
        }
    }
};

static const LocaleInitializer locale_initializer;