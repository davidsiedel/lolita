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

    namespace detail
    {

        template<typename T>
        struct IsNaturalType
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
        struct IsIntegerType
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

            lolita::boolean static const constexpr value = IsNaturalType<T>::value || isShortInt() || isInt() || isLongInt() || isLongLongInt();

        };

        template<typename T>
        struct IsRealType
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

            lolita::boolean static const constexpr value = IsIntegerType<T>::value || isFloat() || isDouble() || isLongDouble();

        };

    }

    template<typename T>
    concept NaturalType = detail::IsNaturalType<T>::value;

    template<typename T>
    concept IntegerType = detail::IsIntegerType<T>::value;

    template<typename T>
    concept RealType = detail::IsRealType<T>::value;

//        template<typename T>
//        concept NaturalType =   std::same_as<std::remove_cvref_t<T>, unsigned short int>
//                                || std::same_as<std::remove_cvref_t<T>, unsigned int>
//                                || std::same_as<std::remove_cvref_t<T>, unsigned long int>
//                                || std::same_as<std::remove_cvref_t<T>, unsigned long long int>;
//
//        template<typename T>
//        concept IntegerType =   NaturalType<T>
//                                || std::same_as<std::remove_cvref_t<T>, short int>
//                                || std::same_as<std::remove_cvref_t<T>, int>
//                                || std::same_as<std::remove_cvref_t<T>, long int>
//                                || std::same_as<std::remove_cvref_t<T>, long long int>;
//
//        template<typename T>
//        concept RealType =      IntegerType<T>
//                                || std::same_as<std::remove_cvref_t<T>, float>
//                                || std::same_as<std::remove_cvref_t<T>, double>
//                                || std::same_as<std::remove_cvref_t<T>, long double>;

    lolita::real const static constexpr pi = 3.14;

    template<lolita::numerics::RealType T, lolita::numerics::IntegerType U>
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

    template<lolita::numerics::RealType T>
    static constexpr
    auto
    sqrt(
            T
            x
    )
    {
        return detail::sqrt<T>(x, T(0), x / T(2) + T(1));
    }

    template<lolita::numerics::RealType T>
    static constexpr
    T
    abs(
            T
            x
    )
    {
        return x < T(0) ? -x : x;
    }

    template<lolita::numerics::RealType T>
    static constexpr
    T
    max(
            T
            x
    )
    {
        return x;
    }

    template<lolita::numerics::RealType T>
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

    template <lolita::numerics::RealType T, lolita::numerics::RealType... U>
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

    template<lolita::numerics::RealType T>
    static constexpr
    T
    prod(
            T
            x
    )
    {
        return x;
    }

    template<lolita::numerics::RealType T>
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

    template<lolita::numerics::RealType T, lolita::numerics::RealType... U>
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

    template<lolita::numerics::RealType T>
    static constexpr
    T
    sum(
            T
            x
    )
    {
        return x;
    }

    template<lolita::numerics::RealType T>
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

    template<lolita::numerics::RealType T, lolita::numerics::RealType... U>
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

    template<lolita::numerics::NaturalType T>
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
        return (k > n) ? 0 : (k == 0 || k == n ) ? 1 : (k == 1 || k == n-1) ? n : (k + k < n) ? (binomial(n - 1, k - 1) * n) / k : (binomial(n - 1, k) * n) / (n - k);
    }

}

namespace lolita::matrix
{

    struct Coordinates
    {

        constexpr
        Coordinates()
                :
                row_(),
                col_()
        {}

        constexpr
        Coordinates(
                lolita::index
                row,
                lolita::index
                col
        )
                :
                row_(row),
                col_(col)
        {}

        lolita::index row_;

        lolita::index col_;

    };

    struct Cardinality
    {

        constexpr
        Cardinality(
                lolita::index
                rows,
                lolita::index
                cols
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
    concept MatrixType = std::derived_from<std::remove_cvref_t<T>, Eigen::DenseBase<std::remove_cvref_t<T>>>;

    using StorageOption = Eigen::StorageOptions;

    lolita::matrix::StorageOption const static constexpr row_major = StorageOption::RowMajor;

    lolita::matrix::StorageOption const static constexpr col_major = StorageOption::ColMajor;

    static constexpr
    lolita::matrix::StorageOption
    storageOption(
            lolita::numerics::IntegerType auto
            num_rows,
            lolita::numerics::IntegerType auto
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
            lolita::numerics::IntegerType auto
            num_rows,
            lolita::numerics::IntegerType auto
            num_cols
    )
    {
        return num_rows != 1 && num_cols == 1 ? col_major : num_rows == 1 && num_cols != 1 ? row_major : storage_option;
    }

    namespace detail
    {

        template<lolita::numerics::RealType T, auto... A>
        struct VectorPolicy;

        template<lolita::numerics::RealType T, auto... A>
        struct MatrixPolicy;

        template<lolita::numerics::RealType T>
        struct VectorPolicy<T>
        {

            using type = Eigen::Matrix<T, -1, 1, storageOption(-1, 1)>;

        };

        template<lolita::numerics::RealType T, auto R>
        struct VectorPolicy<T, R>
        {

            using type = Eigen::Matrix<T, R, 1, storageOption(R, 1)>;

        };

        template<lolita::numerics::RealType T>
        struct MatrixPolicy<T>
        {

            using type = Eigen::Matrix<T, -1, -1, storageOption(-1, -1)>;

        };

        template<lolita::numerics::RealType T, StorageOption O>
        struct MatrixPolicy<T, O>
        {

            using type = Eigen::Matrix<T, -1, -1, storageOption(O, -1, -1)>;

        };

        template<lolita::numerics::RealType T, auto R, auto C>
        struct MatrixPolicy<T, R, C>
        {

            using type = Eigen::Matrix<T, R, C, storageOption(R, C)>;

        };

        template<lolita::numerics::RealType T, auto R, auto C, StorageOption O>
        struct MatrixPolicy<T, R, C, O>
        {

            using type = Eigen::Matrix<T, R, C, storageOption(O, R, C)>;

        };

    }

    template<lolita::numerics::RealType T, auto... A>
    using Vector = typename detail::VectorPolicy<T, A...>::type;

    template<lolita::numerics::RealType T, auto... A>
    using Matrix = typename detail::MatrixPolicy<T, A...>::type;

    template<lolita::matrix::MatrixType T>
    using Span = Eigen::Map<T>;

    template<lolita::matrix::MatrixType T>
    static constexpr
    lolita::integer
    rows()
    {
        return T::CompileTimeTraits::RowsAtCompileTime;
    }

    template<lolita::matrix::MatrixType T>
    static constexpr
    lolita::integer
    cols()
    {
        return T::CompileTimeTraits::ColsAtCompileTime;
    }

    template<lolita::matrix::MatrixType T>
    static constexpr
    lolita::integer
    size()
    {
        return T::CompileTimeTraits::SizeAtCompileTime;
    }

}

#endif //LOLITA_LOLITA_ALGEBRA_HXX
