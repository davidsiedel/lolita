//
// Created by dsiedel on 09/02/2022.
//

#ifndef FETA__FETA_UNORDERED_MAP_HXX
#define FETA__FETA_UNORDERED_MAP_HXX

//#include <unordered_set>
#include "new/_feta_array.hxx"
//#include <variant>
//#include <ostream>
#include "new/_feta.hxx"
//#include "feta/feta/functions.hxx"
#include <unordered_map>

namespace feta
{

    template<typename K, typename T>
    struct UnorderedMap
    {

        struct Item
        {

            K key;

            T value;

        };

        using Type = T;

        using Data = std::unordered_map<K, T>;

        UnorderedMap() = default;

        UnorderedMap(
                std::initializer_list<T> &&
                init
        )
                :
                data(init)
        {}

        UnorderedMap(
                auto && ...
                init
        )
        {
            Indx const constexpr num_args = sizeof...(init);
            Array<Item, num_args> a = {init...};
            for (Indx i = 0; i < num_args; ++i) {
                data.insert({a(i).key, a(i).value});
            }
        }

        Bool
        operator==(
                UnorderedMap const &
                other
        )
        const
        {
            auto eq_0 = data == other.data;
            return eq_0;
        }

        Bool
        operator!=(
                UnorderedMap const &
                other
        )
        const
        {
            return !(other == * this);
        }

//        T &
//        operator() (
//                K const &
//                i
//        )
//        {
//            return data.operator[](i);
//        }
//
//        T &
//        operator() (
//                K &&
//                i
//        )
//        {
//            return data.operator[](i);
//        }
//
//        T const &
//        operator() (
//                K const &
//                i
//        )
//        const
//        {
//            return data.at(i);
//        }
//
//        T const &
//        operator() (
//                K &&
//                i
//        )
//        const
//        {
//            return data.at(i);
//        }

        T &
        get(
                K const &
                i
        )
        {
            return data.operator[](i);
        }

        T &
        get(
                K &&
                i
        )
        {
            return data.operator[](i);
        }

        T const &
        get(
                K const &
                i
        )
        const
        {
            return data.at(i);
        }

        T const &
        get(
                K &&
                i
        )
        const
        {
            return data.at(i);
        }

        Indx
        getSize()
        const
        {
            return data.size();
        }

        Data data;

    };

}

#endif //FETA__FETA_UNORDERED_MAP_HXX
