#ifndef C0F85893_A5C0_442B_AA56_F552BEBAD263
#define C0F85893_A5C0_442B_AA56_F552BEBAD263

#define EIGEN_USE_MKL_ALL

//#define EIGEN_USE_LAPACKE
//#define EIGEN_USE_MKL_VML

// #define EIGEN_USE_BLAS
// #define EIGEN_USE_MKL

#include <limits>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/PardisoSupport>

#include "lolita/lolita_utility.hxx"

namespace lolita::numerics
{

    struct Constants
    {

        static constexpr
        Real
        pi()
        {
            return 3.14;
        }

    };

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

}

namespace lolita::matrix
{

    // struct Shape
    // {

    //     constexpr
    //     Shape(
    //         Integer rows,
    //         Integer cols
    //     )
    //     :
    //     rows_(rows),
    //     cols_(cols)
    //     {}

    //     constexpr
    //     Integer
    //     getRows()
    //     const
    //     {
    //         return rows_;
    //     }

    //     constexpr
    //     Integer
    //     getCols()
    //     const
    //     {
    //         return cols_;
    //     }

    //     constexpr
    //     Integer
    //     getSize()
    //     const
    //     {
    //         return rows_ * cols_;
    //     }

    //     Integer rows_;

    //     Integer cols_;

    // };

    using StorageOption = Eigen::StorageOptions;

    static constexpr
    lolita::matrix::StorageOption
    rowMajor()
    {
        return StorageOption::RowMajor;
    }

    static constexpr
    lolita::matrix::StorageOption
    colMajor()
    {
        return StorageOption::ColMajor;
    }

    namespace detail
    {

        static constexpr
        lolita::matrix::StorageOption
        getStorageOption(
            Integer num_rows,
            Integer num_cols
        )
        {
            return num_rows != 1 && num_cols == 1 ? colMajor() : rowMajor();
        }

        static constexpr
        lolita::matrix::StorageOption
        getStorageOption(
            lolita::matrix::StorageOption storage_option,
            Integer num_rows,
            Integer num_cols
        )
        {
            return num_rows != 1 && num_cols == 1 ? colMajor() : num_rows == 1 && num_cols != 1 ? rowMajor() : storage_option;
        }

        template<typename t_T, auto... A>
        struct VectorPolicy;

        template<typename t_T, auto... A>
        struct MatrixPolicy;

        template<typename t_T>
        struct VectorPolicy<t_T>
        {

            using type = Eigen::Matrix<t_T, -1, 1, getStorageOption(-1, 1)>;

        };

        template<typename t_T, auto t_rows>
        struct VectorPolicy<t_T, t_rows>
        {

            using type = Eigen::Matrix<t_T, t_rows, 1, getStorageOption(t_rows, 1)>;

        };

        template<typename t_T>
        struct MatrixPolicy<t_T>
        {

            using type = Eigen::Matrix<t_T, -1, -1, getStorageOption(-1, -1)>;

        };

        template<typename t_T, StorageOption t_storage_option>
        struct MatrixPolicy<t_T, t_storage_option>
        {

            using type = Eigen::Matrix<t_T, -1, -1, getStorageOption(t_storage_option, -1, -1)>;

        };

        template<typename t_T, auto t_rows, auto t_cols>
        struct MatrixPolicy<t_T, t_rows, t_cols>
        {

            using type = Eigen::Matrix<t_T, t_rows, t_cols, getStorageOption(t_rows, t_cols)>;

        };

        template<typename t_T, auto t_rows, auto t_cols, StorageOption t_storage_option>
        struct MatrixPolicy<t_T, t_rows, t_cols, t_storage_option>
        {

            using type = Eigen::Matrix<t_T, t_rows, t_cols, getStorageOption(t_storage_option, t_rows, t_cols)>;

        };

    }

    template<typename t_T, auto... A>
    using Vector = typename detail::VectorPolicy<t_T, A...>::type;

    template<typename t_T, auto... A>
    using Matrix = typename detail::MatrixPolicy<t_T, A...>::type;

    template<typename t_T>
    using Span = Eigen::Map<t_T>;

}

#endif /* C0F85893_A5C0_442B_AA56_F552BEBAD263 */
