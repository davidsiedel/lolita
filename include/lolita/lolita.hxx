//
// Created by dsiedel on 04/05/22.
//

#ifndef LOLITA_LOLITA_HXX
#define LOLITA_LOLITA_HXX
#include <ostream>

//#include "lolita_matrix.hxx"
//#include "lolita_pointers.hxx"

namespace lolita
{

    namespace config
    {

        using character =   char;

        using integer =     int;

        using index =       unsigned short;

        using natural =     unsigned long long;

        using real =        double;

        using boolean =     bool;

    }

    namespace configuration
    {

        using character =   char;

        using integer =     int;

        using index =       unsigned short;

        using natural =     unsigned long long;

        using real =        double;

        using boolean =     bool;

    }

//    namespace lolita = configuration;

    using namespace config;

//    namespace detail
//    {
//
//        template<typename T, typename... U>
//        concept IsAnyOf =   (std::same_as<T, U> || ...);
//
//        template<typename T, typename... U>
//        concept IsSameAs =  (std::same_as<T, U> && ...);
//
//    }

}

#endif //LOLITA_LOLITA_HXX
