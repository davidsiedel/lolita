#ifndef B7FDA9F4_659F_4665_AE58_B46CDBCCC31A
#define B7FDA9F4_659F_4665_AE58_B46CDBCCC31A

#include <limits>

#include "utility.hxx"

namespace lolita::numerics
{

    static constexpr
    Real
    pi()
    {
        return 3.14;
    }

    template<typename t_T>
    static constexpr
    t_T
    pow(
        t_T x,
        Integer n
    )
    {
        return n == 0 ? t_T(1) : x * pow(x, n - 1);
    }

    template<typename t_T>
    static constexpr
    t_T
    abs(
        t_T x
    )
    {
        return x < t_T(0) ? -x : x;
    }

    template<typename t_T>
    static constexpr
    Boolean
    equal(
        t_T x,
        t_T y
    )
    {
        return abs(x - y) < 1.e-12;
    }

    namespace detail
    {

        static constexpr inline
        Real
        sqrt(
            Real x,
            Real curr,
            Real prev
        )
        {
            return equal(curr, prev) ? curr : sqrt(x, 0.5 * (curr + x / curr), curr);
        }

    }

    template<typename t_T>
    static constexpr
    Real
    sqrt(
        t_T x
    )
    {
        return x >= 0 && x < std::numeric_limits<Real>::infinity() ? detail::sqrt(x, x, 0) : std::numeric_limits<Real>::quiet_NaN();
    }

    template<typename t_T>
    static constexpr
    t_T
    max(
        t_T x
    )
    {
        return x;
    }

    template<typename t_T>
    static constexpr
    t_T
    max(
        t_T x,
        t_T y
    )
    {
        return x < y ? y : x;
    }

    template<typename t_T, typename... t_U>
    static constexpr
    t_T
    max(
        t_T x,
        t_U... y
    )
    {
        return max(x, max(static_cast<t_T>(y)...));
    }

    template<typename t_T>
    static constexpr
    t_T
    prod(
        t_T x
    )
    {
        return x;
    }

    template<typename t_T>
    static constexpr
    t_T
    prod(
        t_T x,
        t_T y
    )
    {
        return x * y;
    }

    template<typename t_T, typename... t_U>
    static constexpr
    t_T
    prod(
        t_T x,
        t_U... y
    )
    {
        return prod(x, prod(static_cast<t_T>(y)...));
    }

    template<typename t_T>
    static constexpr
    t_T
    sum(
        t_T x
    )
    {
        return x;
    }

    template<typename t_T>
    static constexpr
    t_T
    sum(
        t_T x,
        t_T y
    )
    {
        return x + y;
    }

    template<typename t_T, typename... t_U>
    static constexpr
    t_T
    sum(
        t_T x,
        t_U... y
    )
    {
        return sum(x, sum(static_cast<t_T>(y)...));
    }

    template<typename t_T>
    static constexpr
    t_T
    binomial(
        t_T n,
        Integer k
    )
    {
        return (k > n) ? 0 : (k == 0 || k == n) ? 1 : (k == 1 || k == n - 1) ? n : (k + k < n) ? (binomial(n - 1, k - 1) * n) / k : (binomial(n - 1, k) * n) / (n - k);
    }

} // namespace lolita::numerics


#endif /* B7FDA9F4_659F_4665_AE58_B46CDBCCC31A */
