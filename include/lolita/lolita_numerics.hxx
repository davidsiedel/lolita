//
// Created by dsiedel on 15/04/2022.
//

#ifndef LOLITA_LOLITA_NUMERICS_HXX
#define LOLITA_LOLITA_NUMERICS_HXX

#include "lolita/lolita.hxx"
#include "lolita/lolita_lolita.hxx"

namespace lolita::numerics
{

    auto const static constexpr pi = Real(3.14);

    template<typename T, typename U>
    static constexpr
    T
    pow(
            T
            x,
            U
            n
    )
    {
        return n == 0 ? T(1) : x * pow(x, n - 1);
    }

    namespace detail
    {

        template<typename T>
        static constexpr
        T
        sqrt(
                T
                x,
                T
                lo,
                T
                hi
        )
        {
            if (lo == hi) {
                return lo;
            }
            const auto mid = T((lo + hi + 1) / 2);
            if (x / mid < mid) {
                return sqrt<T>(x, lo, mid - 1);
            }
            else {
                return sqrt(x, mid, hi);
            }
        }

    }

    template<std::integral T>
    static constexpr
    auto
    sqrt(
            T
            x
    )
    {
        return detail::sqrt<T>(x, T(0), x / T(2) + T(1));
    }

    template<typename T>
    static constexpr
    T
    abs(
            T
            x
    )
    {
        return x < T(0) ? -x : x;
    }

    template<typename T>
    static constexpr
    T
    max(
            T
            x
    )
    {
        return x;
    }

    template<typename T>
    static constexpr
    T
    max(
            T
            x,
            T
            y
    )
    {
        return x < y ? y : x;
    }

    template <typename T, typename... U>
    static constexpr
    T
    max(
            T
            x,
            U...
            y
    )
    {
        return max(x, max(static_cast<T>(y)...));
    }

    template<typename T>
    static constexpr
    T
    prod(
            T
            x
    )
    {
        return x;
    }

    template<typename T>
    static constexpr
    T
    prod(
            T
            x,
            T
            y
    )
    {
        return x * y;
    }

    template<typename T, typename... U>
    static constexpr
    T
    prod(
            T
            x,
            U...
            y
    )
    {
        return prod(x, prod(static_cast<T>(y)...));
    }

    template<typename T>
    static constexpr
    T
    sum(
            T
            x
    )
    {
        return x;
    }

    template<typename T>
    static constexpr
    T
    sum(
            T
            x,
            T
            y
    )
    {
        return x + y;
    }

    template<typename T, typename... U>
    static constexpr
    T
    sum(
            T
            x,
            U...
            y
    )
    {
        return sum(x, sum(static_cast<T>(y)...));
    }

    template<std::integral T>
    static constexpr
    T
    binomial(
            T
            n,
            T
            k
    )
    {
//        return
//                (        k> n  )? 0 :                    // out of range
//                (k==0 || k==n  )? 1 :                    // edge
//                (k==1 || k==n-1)? n :                    // first
//                (     k+k < n  )?                        // recursive:
//                (binomial(n - 1, k - 1) * n) / k :       //  path to k=1   is faster
//                (binomial(n - 1, k) * n) / (n - k);      //  path to k=n-1 is faster
        return (k > n) ? 0 : (k == 0 || k == n ) ? 1 : (k == 1 || k == n-1) ? n : (k + k < n) ? (binomial(n - 1, k - 1) * n) / k : (binomial(n - 1, k) * n) / (n - k);      //  path to k=n-1 is faster
    }

}

#endif //LOLITA_LOLITA_NUMERICS_HXX
