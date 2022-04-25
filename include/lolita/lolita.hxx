//
// Created by dsiedel on 13/04/2022.
//

#ifndef FETA_FETA_FETA_TYPES_HXX
#define FETA_FETA_FETA_TYPES_HXX

#include <iostream>
#include <tuple>

namespace lolita
{

    /**
     * @brief
     */
    struct Void{};

    /**
     * @brief
     */
    using Char = char;

    /**
     * @brief
     */
    using Intg = int;

    /**
     * @brief
     */
    using Indx = std::size_t;

    /**
     * @brief
     */
    using Long = unsigned long long;

    /**
     * @brief
     */
    using Real = double;

    /**
     * @brief
     */
    using Strg = std::string;

    /**
     * @brief
     */
    using Bool = bool;

    /**
     * @brief
     */
    using StrgStream = std::stringstream;

    /**
     * @brief
     */
    using OuptStream = std::ostream;

    /**
     * @brief
     */
    using InptStream = std::istream;

    /**
     * @brief
     */
    using FileStream = std::ifstream;

//    template<typename T, auto... A>
//    struct Array;
//
//    template<typename...T>
//    struct Collection;
//
//    template<typename...T>
//    struct Aggregate;
//
//    template<typename K, typename T>
//    struct UnorderedMap;
//
//    template<typename T>
//    struct SharedPointer;
//
//    template<typename T>
//    struct UniquePointer;
//
//    template<typename T, auto... A>
//    struct Matrix;

}

#endif //FETA_FETA_FETA_TYPES_HXX
