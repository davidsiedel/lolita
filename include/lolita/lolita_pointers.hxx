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
        data(std::make_shared<T>(std::forward<T>(arg)))
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

        Bool
        operator!=(
                SharedPointer const &
                other
        )
        const
        {
            return !(* this == other);
        }

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

        Bool
        exists()
        const
        {
            return data != nullptr;
        }

        template<typename U>
        static
        auto
        make()
        {
            auto p = SharedPointer();
            p.data = std::make_shared<U>(U());
            return p;
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

        template<typename U>
        static
        auto
        make(
                U &&
                arg
        )
        {
            auto p = SharedPointer();
            p.data = std::make_shared<U>(std::forward<T>(arg));
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
        data(std::make_unique<T>(std::forward<T>(arg)))
        {}

        UniquePointer(
                UniquePointer const &
                other
        )
        :
        data(setData(other))
        {}


        UniquePointer(
                UniquePointer &&
                other
        )
        :
        data(setData(std::forward<UniquePointer>(other)))
        {}

        UniquePointer &
        operator =(
                UniquePointer const &
                other
        )
        {
            if (other.data != nullptr) {
                data = std::make_unique<T>(other.get());
            }
            return * this;
        }

        UniquePointer &
        operator =(
                UniquePointer &&
                other
        )
        noexcept
        {
            if (other.data != nullptr) {
                data = std::move(other.data);
            }
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

        Bool
        exists()
        const
        {
            return data != nullptr;
        }

        auto
        setData(
                UniquePointer const &
                other
        )
        {
            return other.data == nullptr ? nullptr : std::make_unique<T>(other.get());
        }

        auto
        setData(
                UniquePointer &&
                other
        )
        {
            return other.data == nullptr ? nullptr : std::move(other.data);
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

    namespace pointer
    {

        template<typename T>
        struct Traits
        {

            auto const static constexpr value = false;

        };

        template<typename T>
        struct Traits<SharedPointer<T>>
        {

            auto const static constexpr value = true;

        };

        template<typename T>
        struct Traits<UniquePointer<T>>
        {

            auto const static constexpr value = true;

        };

        template<typename T>
        void
        make(
                SharedPointer<T> &
                pointer_arg
        )
        {
            pointer_arg.data = std::make_shared<T>(T());
        }

        template<typename T>
        void
        make(
                SharedPointer<T> &
                pointer_arg,
                T &&
                arg
        )
        {
            pointer_arg.data = std::make_shared<T>(std::forward<T>(arg));
        }

        template<typename T>
        void
        make(
                UniquePointer<T> &
                pointer_arg
        )
        {
            pointer_arg.data = std::make_unique<T>(T());
        }

        template<typename T>
        void
        make(
                UniquePointer<T> &
                pointer_arg,
                T &&
                arg
        )
        {
            pointer_arg.data = std::make_unique<T>(std::forward<T>(arg));
        }

    }

    template<typename T>
    concept Pointer = pointer::Traits<T>::value;

}

#endif //LOLITA_LOLITA_POINTERS_HXX
