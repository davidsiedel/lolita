//
// Created by dsiedel on 06/05/22.
//

#ifndef LOLITA_LOLITA_ALGEBRA_HXX
#define LOLITA_LOLITA_ALGEBRA_HXX

//#define EIGEN_USE_MKL_ALL
//#define EIGEN_USE_LAPACKE
//#define EIGEN_USE_MKL_VML
#define EIGEN_USE_BLAS
#define EIGEN_USE_MKL

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/PardisoSupport>

#include "lolita/lolita.hxx"

namespace lolita::numerics
{

    lolita::real const static constexpr sqrt_2 = 1.41421356237;

    namespace detail
    {

        template<typename T>
        struct IsNaturalConcept
        {

            static constexpr
            lolita::boolean
            isUnsignedShortInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, unsigned short int>;
            }

            static constexpr
            lolita::boolean
            isUnsignedInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, unsigned int>;
            }

            static constexpr
            lolita::boolean
            isUnsignedLongInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, unsigned long int>;
            }

            static constexpr
            lolita::boolean
            isUnsignedLongLongInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, unsigned long long int>;
            }

            lolita::boolean static const constexpr value = isUnsignedShortInt() || isUnsignedInt() || isUnsignedLongInt() || isUnsignedLongLongInt();

        };

        template<typename T>
        struct IsIntegerConcept
        {

            static constexpr
            lolita::boolean
            isShortInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, short int>;
            }

            static constexpr
            lolita::boolean
            isInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, int>;
            }

            static constexpr
            lolita::boolean
            isLongInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, long int>;
            }

            static constexpr
            lolita::boolean
            isLongLongInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, long long int>;
            }

            lolita::boolean static const constexpr value = IsNaturalConcept<T>::value || isShortInt() || isInt() || isLongInt() || isLongLongInt();

        };

        template<typename T>
        struct IsRealConcept
        {

            static constexpr
            lolita::boolean
            isFloat()
            {
                return std::is_same_v<std::remove_cvref_t<T>, float>;
            }

            static constexpr
            lolita::boolean
            isDouble()
            {
                return std::is_same_v<std::remove_cvref_t<T>, double>;
            }

            static constexpr
            lolita::boolean
            isLongDouble()
            {
                return std::is_same_v<std::remove_cvref_t<T>, long double>;
            }

            lolita::boolean static const constexpr value = IsIntegerConcept<T>::value || isFloat() || isDouble() || isLongDouble();

        };

    }

    template<typename T>
    concept NaturalConcept = detail::IsNaturalConcept<T>::value;

    template<typename T>
    concept IntegerConcept = detail::IsIntegerConcept<T>::value;

    template<typename T>
    concept RealConcept = detail::IsRealConcept<T>::value;

//        template<typename T>
//        concept NaturalConcept =   std::same_as<std::remove_cvref_t<T>, unsigned short int>
//                                || std::same_as<std::remove_cvref_t<T>, unsigned int>
//                                || std::same_as<std::remove_cvref_t<T>, unsigned long int>
//                                || std::same_as<std::remove_cvref_t<T>, unsigned long long int>;
//
//        template<typename T>
//        concept IntegerConcept =   NaturalConcept<T>
//                                || std::same_as<std::remove_cvref_t<T>, short int>
//                                || std::same_as<std::remove_cvref_t<T>, int>
//                                || std::same_as<std::remove_cvref_t<T>, long int>
//                                || std::same_as<std::remove_cvref_t<T>, long long int>;
//
//        template<typename T>
//        concept RealConcept =      IntegerConcept<T>
//                                || std::same_as<std::remove_cvref_t<T>, float>
//                                || std::same_as<std::remove_cvref_t<T>, double>
//                                || std::same_as<std::remove_cvref_t<T>, long double>;

    lolita::real const static constexpr pi = 3.14;

    template<lolita::numerics::RealConcept T, lolita::numerics::IntegerConcept U>
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

    template<lolita::numerics::RealConcept T>
    static constexpr
    auto
    sqrt(
            T
            x
    )
    {
        return detail::sqrt<T>(x, T(0), x / T(2) + T(1));
    }

    template<lolita::numerics::RealConcept T>
    static constexpr
    T
    abs(
            T
            x
    )
    {
        return x < T(0) ? -x : x;
    }

    template<lolita::numerics::RealConcept T>
    static constexpr
    T
    max(
            T
            x
    )
    {
        return x;
    }

    template<lolita::numerics::RealConcept T>
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

    template <lolita::numerics::RealConcept T, lolita::numerics::RealConcept... U>
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

    template<lolita::numerics::RealConcept T>
    static constexpr
    T
    prod(
            T
            x
    )
    {
        return x;
    }

    template<lolita::numerics::RealConcept T>
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

    template<lolita::numerics::RealConcept T, lolita::numerics::RealConcept... U>
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

    template<lolita::numerics::RealConcept T>
    static constexpr
    T
    sum(
            T
            x
    )
    {
        return x;
    }

    template<lolita::numerics::RealConcept T>
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

    template<lolita::numerics::RealConcept T, lolita::numerics::RealConcept... U>
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

    template<lolita::numerics::IntegerConcept T, lolita::numerics::IntegerConcept U>
    static constexpr
    T
    binomial(
            T
            n,
            U
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
        return (k > n) ? 0 : (k == 0 || k == n ) ? 1 : (k == 1 || k == n-1) ? n : (k + k < n) ? (binomial(n - 1, k - 1) * n) / k : (binomial(n - 1, k) * n) / (n - k);
    }

}

namespace lolita::matrix
{

    struct Coordinates
    {

        lolita::index row_;

        lolita::index col_;

    };

    struct Shape
    {

        constexpr
        Shape(
                lolita::index rows,
                lolita::index cols
        )
        :
        rows_(rows),
        cols_(cols),
        size_(rows * cols)
        {}

        lolita::index rows_;

        lolita::index cols_;

        lolita::index size_;

    };

    template<typename T>
    concept MatrixConcept = std::derived_from<std::remove_cvref_t<T>, Eigen::DenseBase<std::remove_cvref_t<T>>>;

    using StorageOption = Eigen::StorageOptions;

    lolita::matrix::StorageOption const static constexpr row_major = StorageOption::RowMajor;

    lolita::matrix::StorageOption const static constexpr col_major = StorageOption::ColMajor;

    static constexpr
    lolita::matrix::StorageOption
    storageOption(
            lolita::numerics::IntegerConcept auto
            num_rows,
            lolita::numerics::IntegerConcept auto
            num_cols
    )
    {
        return num_rows != 1 && num_cols == 1 ? col_major : row_major;
    }

    static constexpr
    lolita::matrix::StorageOption
    storageOption(
            lolita::matrix::StorageOption
            storage_option,
            lolita::numerics::IntegerConcept auto
            num_rows,
            lolita::numerics::IntegerConcept auto
            num_cols
    )
    {
        return num_rows != 1 && num_cols == 1 ? col_major : num_rows == 1 && num_cols != 1 ? row_major : storage_option;
    }

    namespace detail
    {

        template<lolita::numerics::RealConcept T, auto... A>
        struct VectorPolicy;

        template<lolita::numerics::RealConcept T, auto... A>
        struct MatrixPolicy;

        template<lolita::numerics::RealConcept T>
        struct VectorPolicy<T>
        {

            using type = Eigen::Matrix<T, -1, 1, storageOption(-1, 1)>;

        };

        template<lolita::numerics::RealConcept T, auto R>
        struct VectorPolicy<T, R>
        {

            using type = Eigen::Matrix<T, R, 1, storageOption(R, 1)>;

        };

        template<lolita::numerics::RealConcept T>
        struct MatrixPolicy<T>
        {

            using type = Eigen::Matrix<T, -1, -1, storageOption(-1, -1)>;

        };

        template<lolita::numerics::RealConcept T, StorageOption O>
        struct MatrixPolicy<T, O>
        {

            using type = Eigen::Matrix<T, -1, -1, storageOption(O, -1, -1)>;

        };

        template<lolita::numerics::RealConcept T, auto R, auto C>
        struct MatrixPolicy<T, R, C>
        {

            using type = Eigen::Matrix<T, R, C, storageOption(R, C)>;

        };

        template<lolita::numerics::RealConcept T, auto R, auto C, StorageOption O>
        struct MatrixPolicy<T, R, C, O>
        {

            using type = Eigen::Matrix<T, R, C, storageOption(O, R, C)>;

        };

    }

    template<lolita::numerics::RealConcept T, auto... A>
    using Vector = typename detail::VectorPolicy<T, A...>::type;

    template<lolita::numerics::RealConcept T, auto... A>
    using Matrix = typename detail::MatrixPolicy<T, A...>::type;

    template<lolita::matrix::MatrixConcept T>
    using Span = Eigen::Map<T>;

    template<lolita::matrix::MatrixConcept T>
    static constexpr
    lolita::integer
    rows()
    {
        return T::CompileTimeTraits::RowsAtCompileTime;
    }

    template<lolita::matrix::MatrixConcept T>
    static constexpr
    lolita::integer
    cols()
    {
        return T::CompileTimeTraits::ColsAtCompileTime;
    }

    template<lolita::matrix::MatrixConcept T>
    static constexpr
    lolita::integer
    size()
    {
        return T::CompileTimeTraits::SizeAtCompileTime;
    }

}

#endif //LOLITA_LOLITA_ALGEBRA_HXX
