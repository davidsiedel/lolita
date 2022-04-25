//
// Created by dsiedel on 18/04/22.
//

#ifndef LOLITA_LOLITA_POINTERS_HXX
#define LOLITA_LOLITA_POINTERS_HXX

#include <memory>
#include "lolita/lolita.hxx"

namespace lolita
{

    template<typename T>
    struct SharedPointer
    {

        using Type = T;

        using Data = std::shared_ptr<T>;

        struct Hash
        {

            Indx
            operator() (
                    SharedPointer const & other
            )
            const
            {
                return std::hash<Data>()(other.data);
            }

        };

        SharedPointer()
        :
        data(nullptr)
        {}

        explicit
        SharedPointer(
                T const &
                arg
        )
        :
        data(std::make_shared<T>(arg))
        {}

        explicit
        SharedPointer(
                T &&
                arg
        )
        :
        data(std::make_shared<T>(arg))
        {}

        Bool
        operator==(
                SharedPointer const &
                other
        )
        const
        {
            return data == other.data;
        }
//        = default;

        Bool
        operator!=(
                SharedPointer const &
                other
        )
        const
        {
            return !(* this == other);
        }
//        = default;

        T &
        get()
        {
            return * data;
        }

        T const &
        get()
        const
        {
            return * data;
        }

        template<typename U>
        static
        auto
        make(
                U const &
                arg
        )
        {
            auto p = SharedPointer();
            p.data = std::make_shared<U>(arg);
            return p;
        }

        Data data;

    };

    template<typename T>
    struct UniquePointer
    {

        using Type = T;

        using Data = std::unique_ptr<T>;

        struct Hash
        {

            Indx
            operator() (
                    UniquePointer const & other
            )
            const
            {
                return std::hash<Data>()(other.data);
            }

        };

        UniquePointer()
        :
        data(nullptr)
        {}

        explicit
        UniquePointer(
                T const &
                arg
        )
        :
        data(std::make_unique<T>(arg))
        {}

        explicit
        UniquePointer(
                T &&
                arg
        )
        :
        data(std::make_unique<T>(arg))
        {}

        UniquePointer(
                UniquePointer const &
                other
        )
        :
        data(std::make_unique<T>(other.get()))
        {}


        UniquePointer(
                UniquePointer &&
                other
        )
        :
        data(std::move(other.data))
        {}

        UniquePointer &
        operator =(
                UniquePointer const &
                other
        )
        {
            data = std::make_unique<T>(other.get());
            return * this;
        }

        UniquePointer &
        operator =(
                UniquePointer &&
                other
        )
        noexcept
        {
            data = std::move(other.data);
            return * this;
        }

        Bool
        operator==(
                UniquePointer const &
                other
        )
        const = default;

        Bool
        operator!=(
                UniquePointer const &
                other
        )
        const = default;

        T &
        get()
        {
            return * data;
        }

        T const &
        get()
        const
        {
            return * data;
        }

        template<typename U>
        static
        auto
        make(
                U const &
                arg
        )
        {
            auto p = UniquePointer();
            p.data = std::make_unique<U>(arg);
            return p;
        }

        template<typename U>
        static
        auto
        copy(
                UniquePointer<U> const &
                other
        )
        {
            auto p = UniquePointer();
            p.data = std::make_unique<U>(other.get());
            return p;
        }

        Data data;

    };

}

#endif //LOLITA_LOLITA_POINTERS_HXX
