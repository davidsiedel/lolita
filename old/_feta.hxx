//
// Created by dsiedel on 20/11/2021.
//

#ifndef FETA__FETA_HXX
#define FETA__FETA_HXX

#include <mutex>
#include <thread>
#include <cxxabi.h>
#include <iostream>

namespace feta
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

//    /**
//     * @brief
//     */
//    using Mutx = std::mutex;
//
//    /**
//     * @brief
//     */
//    using Lock = std::scoped_lock<Mutx>;
//
//    /**
//     * @brief
//     */
//    using Thrd = std::jthread;

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

    template<typename T>
    using RawType = std::remove_cv_t<std::remove_reference_t<std::remove_pointer_t<T>>>;

//    template<typename T, typename... U>
//    inline constexpr auto is_in_v = std::disjunction_v<std::is_same<T, U>...>;

//    template<typename T, typename... U>
//    inline constexpr bool are_same_v = std::conjunction_v<std::is_same<T, U>...>;

    template<typename T, typename... U>
    auto const static constexpr is_same_v = std::conjunction_v<std::is_same<T, U>...>;

    template<typename T, typename... U>
    auto const static constexpr is_in_v = std::disjunction_v<std::is_same<T, U>...>;

    template<typename T>
    struct TypeTraits
    {

        using Self = T;

        using RawType = std::remove_cvref_t<std::remove_pointer_t<T>>;

//        template<T>
//        struct TemplatedType {};

        template<typename... U>
        auto const static constexpr is_in = std::disjunction_v<std::is_same<T, U>...>;

        template<typename... U>
        auto const static constexpr is_raw_in = std::disjunction_v<std::is_same<RawType, U>...>;

        template<typename... U>
        auto const static constexpr are_same = std::conjunction_v<std::is_same<T, U>...>;

    };

    template<typename T>
    struct Self
    {

        using Type = T;

        template<T>
        struct TemplatedType
        {};

    };

    namespace detail
    {

        template<typename... L>
        struct UniqueTypePolicy
        {

        private:

            template<typename T, typename... U>
            struct ValueSetter
            {

                auto const static constexpr value = std::conjunction_v<std::is_same<T, U>...>;

            };

        public:

            auto const static constexpr value = ValueSetter<L...>::value;

            using Type = typename std::tuple_element<0, std::tuple<L...>>::type;

        };

        template<>
        struct UniqueTypePolicy<>
        {

            auto const static constexpr value = true;

            using Type = void;

        };

    }

    template<typename... L>
    using UniqueType = typename detail::UniqueTypePolicy<L...>;

//    template<typename... L>
//    using UniqueType2 = typename detail::UniqueTypePolicy<L...>::Type;

    template<typename T, auto... A>
    struct Array;

//    template<typename T>
//    constexpr inline
//    auto
//    diff(
//            T
//            start_arg,
//            T
//            last_arg
//    )
//    {
//        return last_arg - start_arg;
//    }

    template<typename T, typename... U>
    constexpr inline
    auto
    isInn(
            T const &
            val,
            U &&...
            values
    )
    {
        static_assert(is_same_v<T, U...>);
        auto v = Bool(false);
        auto u = Array<T, sizeof...(U)>{values...};
        for (int i = 0; i < sizeof...(U); ++i) {
            if (val == u.get(i)) {
                v = true;
            }
        }
        return v;
    }

    template<typename T>
    constexpr inline
    auto
    range(
            T
            start_arg,
            T
            last_arg
    )
    {
        auto a = Array<T, last_arg - start_arg>();
        for (T i = 0; i < last_arg - start_arg; ++i) {
            a.get(i) = i + start_arg;
        }
        return a;
    }

    template<auto I, typename T>
    constexpr inline
    auto
    range(
            T
            start_arg
    )
    {
        auto a = Array<T, I>();
        for (T i = 0; i < I; ++i) {
            a.get(i) = i + start_arg;
        }
        return a;
    }

    template<typename... T>
    constexpr inline
    auto
    max(
            T &&...
            value_args
    )
    {
        static_assert(UniqueType<T...>::value);
        using type = std::tuple_element_t<0, std::tuple<T...>>;
        auto values = std::array<type, sizeof...(T)>{value_args...};
        auto value = static_cast<type>(0);
        for (auto i: values) {
            if (i > value) {
                value = i;
            }
        }
        return value;
    }

    template<typename T, auto... A>
    constexpr inline
    auto
    getMax(
            Array<T, A...> const &
            array_arg
    )
    {
        auto value = T(0);
        for (auto i: array_arg.data) {
            if (i > value) {
                value = i;
            }
        }
        return value;
    }

    template<typename... T>
    constexpr inline
    auto
    getProduct(
            T &&...
            value_args
    )
    {
        static_assert(UniqueType<T...>::value);
        using type = std::tuple_element_t<0, std::tuple<T...>>;
        auto values = std::array<type, sizeof...(T)>{static_cast<type>(value_args)...};
        auto value = type(1);
        for (Indx i = 0; i < sizeof...(T); ++i) {
            value *= values[i];
        }
        return value;
    }

    template<typename T>
    constexpr inline
    T
    getBinomial(
            T
            n,
            T
            k
    )
    noexcept
    {
        return
                (        k> n  )? 0 :          // out of range
                (k==0 || k==n  )? 1 :          // edge
                (k==1 || k==n-1)? n :          // first
                (     k+k < n  )?              // recursive:
                (getBinomial(n-1,k-1) * n)/k :       //  path to k=1   is faster
                (getBinomial(n-1,k) * n)/(n-k);      //  path to k=n-1 is faster
    }

    template<typename T, typename... U>
    inline static
    void
    print(
            std::ostream &
            out,
            T &&
            arg,
            U &&...
            args
    )
    {
        out << std::forward<T>(arg);
        ((out << " " << std::forward<U>(args)), ...);
    }

    template<typename T, typename... U>
    inline static
    void
    print(
            T &&
            arg,
            U &&...
            args
    )
    {
        std::cout << std::forward<T>(arg);
        ((std::cout << " " << std::forward<U>(args)), ...);
        std::cout << std::endl;
    }

    inline static
    void
    printStructName(
            auto const & arg
    )
    {
        int status;
        std::cout << typeid(arg).name() << std::endl;
//        auto ptr = std::unique_ptr<char, decltype(& std::free)>{
//                abi::__cxa_demangle(typeid(arg).name(), nullptr, nullptr, nullptr),
//                std::free
//        };
//        Strg strg = {ptr.get()};
//        std::cout << strg << std::endl;
//        char * demangled = abi::__cxa_demangle(typeid(arg).name(),0,0,&status);
//        Strg demangled = abi::__cxa_demangle(typeid(arg).name(),0,0,&status);
//        std::cout<<demangled<<"\t"<< quote(A) <<"\n";
//        std::cout << demangled << std::endl;
//        free(demangled);
    }

}

#endif //FETA__FETA_HXX
