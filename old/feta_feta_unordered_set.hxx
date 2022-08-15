//
// Created by dsiedel on 09/02/2022.
//

#ifndef FETA_FETA_FETA_UNORDERED_SET_HXX
#define FETA_FETA_FETA_UNORDERED_SET_HXX

#include <unordered_set>
#include <iterator>
#include "new/_feta.hxx"

namespace feta
{

    template<typename Type, typename HashType>
    struct UnorderedSetHashHelper
    {

        using UnorderedSetDataT = std::unordered_set<Type, HashType>;

    };

    template<typename Type>
    struct UnorderedSetHashHelper<Type, void>
    {

        using UnorderedSetDataT = std::unordered_set<Type>;

    };

    template<typename T, typename HashType = void>
    struct UnorderedSet
    {

        using Type = T;

        using Data = typename UnorderedSetHashHelper<T, HashType>::UnorderedSetDataT;

        UnorderedSet() = default;

        UnorderedSet(
                std::initializer_list<T> &&
                init
        )
                :
                data(init)
        {}

        Bool
        operator==(
                UnorderedSet const &
                other
        )
        const
        {
            auto eq_0 = data == other.data;
            return eq_0;
        }

        Bool
        operator!=(
                UnorderedSet const &
                other
        )
        const
        {
            return !(other == * this);
        }

        T &
        operator() (
                Indx
                i
        )
        {
            auto it = data.begin();
            std::advance(it, i);
            return data.extract(it).value();
        }

        T const &
        operator() (
                Indx
                i
        )
        const
        {
            auto it = data.begin();
            std::advance(it, i);
            return data.extract(it).value();
        }

        constexpr
        Indx
        getSize()
        const
        {
            return data.size();
        }

        Data data;

    };

    template<typename Type>
    using PointerUnorederedSet = UnorderedSet<Type, typename Type::Hash>;

}

#endif //FETA_FETA_FETA_UNORDERED_SET_HXX
