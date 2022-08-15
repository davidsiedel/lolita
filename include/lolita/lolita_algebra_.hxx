//
// Created by dsiedel on 06/05/22.
//

#ifndef LOLITA_LOLITA_ALGEBRA_HXX
#define LOLITA_LOLITA_ALGEBRA_HXX

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

    Real const static constexpr sqrt_2 = 1.41421356237;

    namespace detail
    {

        template<typename T>
        struct IsNaturalConcept
        {

            static constexpr
            Boolean
            isUnsignedShortInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, unsigned short int>;
            }

            static constexpr
            Boolean
            isUnsignedInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, unsigned int>;
            }

            static constexpr
            Boolean
            isUnsignedLongInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, unsigned long int>;
            }

            static constexpr
            Boolean
            isUnsignedLongLongInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, unsigned long long int>;
            }

            Boolean static const constexpr value = isUnsignedShortInt() || isUnsignedInt() || isUnsignedLongInt() || isUnsignedLongLongInt();

        };

        template<typename T>
        struct IsIntegerConcept
        {

            static constexpr
            Boolean
            isShortInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, short int>;
            }

            static constexpr
            Boolean
            isInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, int>;
            }

            static constexpr
            Boolean
            isLongInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, long int>;
            }

            static constexpr
            Boolean
            isLongLongInt()
            {
                return std::is_same_v<std::remove_cvref_t<T>, long long int>;
            }

            Boolean static const constexpr value = IsNaturalConcept<T>::value || isShortInt() || isInt() || isLongInt() || isLongLongInt();

        };

        template<typename T>
        struct IsRealConcept
        {

            static constexpr
            Boolean
            isFloat()
            {
                return std::is_same_v<std::remove_cvref_t<T>, float>;
            }

            static constexpr
            Boolean
            isDouble()
            {
                return std::is_same_v<std::remove_cvref_t<T>, double>;
            }

            static constexpr
            Boolean
            isLongDouble()
            {
                return std::is_same_v<std::remove_cvref_t<T>, long double>;
            }

            Boolean static const constexpr value = IsIntegerConcept<T>::value || isFloat() || isDouble() || isLongDouble();

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

    Real const static constexpr pi = 3.14;

    template<lolita::numerics::RealConcept T, lolita::numerics::IntegerConcept U>
    static constexpr
    T
    pow(
        T x,
        U n
    )
    {
        return n == 0 ? T(1) : x * pow(x, n - 1);
    }

    template<lolita::numerics::RealConcept T>
    static constexpr
    T
    abs(
        T x
    )
    {
        return x < T(0) ? -x : x;
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

        template<typename T>
        static constexpr
        T
        sqrt2(
            T x,
            T lo,
            T hi
        )
        {
            if (lo == hi)
            {
                return lo;
            }
            const auto mid = T((lo + hi + 1) / 2);
            if (x / mid < mid)
            {
                return sqrt2<T>(x, lo, mid - 1);
            }
            else
            {
                return sqrt2<T>(x, mid, hi);
            }
        }

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

    template<lolita::numerics::RealConcept T>
    static constexpr
    auto
    sqrt2(
        T x
    )
    {
        return detail::sqrt2<T>(x, T(0), x / T(2) + T(1));
    }

    static constexpr inline
    Real
    sqrt(
        Real x
    )
    {
        return x >= 0 && x < std::numeric_limits<Real>::infinity() ? detail::sqrt(x, x, 0) : std::numeric_limits<Real>::quiet_NaN();
    }

    template<lolita::numerics::RealConcept T>
    static constexpr
    T
    max(
        T x
    )
    {
        return x;
    }

    template<lolita::numerics::RealConcept T>
    static constexpr
    T
    max(
        T x,
        T y
    )
    {
        return x < y ? y : x;
    }

    template <lolita::numerics::RealConcept T, lolita::numerics::RealConcept... U>
    static constexpr
    T
    max(
        T x,
        U... y
    )
    {
        return max(x, max(static_cast<T>(y)...));
    }

    template<lolita::numerics::RealConcept T>
    static constexpr
    T
    prod(
        T x
    )
    {
        return x;
    }

    template<lolita::numerics::RealConcept T>
    static constexpr
    T
    prod(
        T x,
        T y
    )
    {
        return x * y;
    }

    template<lolita::numerics::RealConcept T, lolita::numerics::RealConcept... U>
    static constexpr
    T
    prod(
        T x,
        U... y
    )
    {
        return prod(x, prod(static_cast<T>(y)...));
    }

    template<lolita::numerics::RealConcept T>
    static constexpr
    T
    sum(
        T x
    )
    {
        return x;
    }

    template<lolita::numerics::RealConcept T>
    static constexpr
    T
    sum(
        T x,
        T y
    )
    {
        return x + y;
    }

    template<lolita::numerics::RealConcept T, lolita::numerics::RealConcept... U>
    static constexpr
    T
    sum(
        T x,
        U... y
    )
    {
        return sum(x, sum(static_cast<T>(y)...));
    }

    template<lolita::numerics::IntegerConcept T, lolita::numerics::IntegerConcept U>
    static constexpr
    T
    binomial(
        T n,
        U k
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
    
    struct VectorBlock
    {
        
        constexpr
        VectorBlock()
        :
        i_(-1),
        j_(-1)
        {}
        
        constexpr
        VectorBlock(
            Integer i,
            Integer j
        )
        :
        i_(i),
        j_(j)
        {}
        
        Integer i_;
        
        Integer j_;

    };
    
    struct MatrixBlock
    {
        
        constexpr
        MatrixBlock()
        :
        i_(-1),
        j_(-1),
        k_(-1),
        l_(-1)
        {}
        
        constexpr
        MatrixBlock(
            Integer i,
            Integer j,
            Integer k,
            Integer l
        )
        :
        i_(i),
        j_(j),
        k_(k),
        l_(l)
        {}
        
        Integer i_;
        
        Integer j_;
        
        Integer k_;
        
        Integer l_;

    };

    struct Coordinates
    {

        constexpr
        Coordinates()
        :
        row_(-1),
        col_(-1)
        {}

        constexpr
        Coordinates(
            Integer row,
            Integer col
        )
        :
        row_(row),
        col_(col)
        {}

        Integer row_;

        Integer col_;

    };

    struct Shape
    {

        constexpr
        Shape(
            Integer rows,
            Integer cols
        )
        :
        rows_(rows),
        cols_(cols),
        size_(rows * cols)
        {}

        constexpr
        Integer
        rows()
        const
        {
            return rows_;
        }

        constexpr
        Integer
        cols()
        const
        {
            return cols_;
        }

        constexpr
        Integer
        size()
        const
        {
            return size_;
        }

        Integer rows_;

        Integer cols_;

        Integer size_;

    };

    template<typename _T, typename _U, Integer _rows = -1, Integer _cols = -1>
    concept StaticMatrixObject = requires
    {
        requires _T::CompileTimeTraits::RowsAtCompileTime == _rows;
        requires _T::CompileTimeTraits::ColsAtCompileTime == _cols;
        requires std::same_as<_U, typename _T::Scalar>;
        requires std::convertible_to<std::remove_reference_t<_T>, Eigen::Matrix<_U, _rows, _cols>>;
    };

    template<typename T>
    concept MatrixConceptO = std::derived_from<std::remove_cvref_t<T>, Eigen::DenseBase<std::remove_cvref_t<T>>>;

    template<typename _T, typename _U, Integer _rows = -1, Integer _cols = -1>
    concept MatrixConceptE = requires
    {
        requires _T::CompileTimeTraits::RowsAtCompileTime == _rows;
        requires _T::CompileTimeTraits::ColsAtCompileTime == _cols;
        requires std::same_as<_U, typename _T::Scalar>;
        requires std::convertible_to<std::remove_reference_t<_T>, Eigen::Matrix<_U, _rows, _cols>>;
    };

    template<typename T>
    concept MatrixConcept = std::convertible_to<
            std::remove_cvref_t<T>,
            Eigen::Matrix<
                    typename std::remove_cvref_t<T>::Scalar,
                    std::remove_cvref_t<T>::CompileTimeTraits::RowsAtCompileTime,
                    std::remove_cvref_t<T>::CompileTimeTraits::ColsAtCompileTime
            >
    >;

    using StorageOption = Eigen::StorageOptions;

    lolita::matrix::StorageOption const static constexpr row_major = StorageOption::RowMajor;

    lolita::matrix::StorageOption const static constexpr col_major = StorageOption::ColMajor;

    static constexpr
    lolita::matrix::StorageOption
    storageOption(
        lolita::numerics::IntegerConcept auto num_rows,
        lolita::numerics::IntegerConcept auto num_cols
    )
    {
        return num_rows != 1 && num_cols == 1 ? col_major : row_major;
    }

    static constexpr
    lolita::matrix::StorageOption
    storageOption(
        lolita::matrix::StorageOption storage_option,
        lolita::numerics::IntegerConcept auto num_rows,
        lolita::numerics::IntegerConcept auto num_cols
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
    Integer
    rows()
    {
        return T::CompileTimeTraits::RowsAtCompileTime;
    }

    template<lolita::matrix::MatrixConcept T>
    static constexpr
    Integer
    cols()
    {
        return T::CompileTimeTraits::ColsAtCompileTime;
    }

    template<lolita::matrix::MatrixConcept T>
    static constexpr
    Integer
    size()
    {
        return T::CompileTimeTraits::SizeAtCompileTime;
    }

    template<lolita::matrix::MatrixConcept T>
    static constexpr
    Integer
    rows(
            T const & matrix
    )
    {
        return T::CompileTimeTraits::RowsAtCompileTime;
    }

    template<lolita::matrix::MatrixConcept T>
    static constexpr
    Integer
    cols(
            T const & matrix
    )
    {
        return T::CompileTimeTraits::ColsAtCompileTime;
    }

    template<lolita::matrix::MatrixConcept T>
    static constexpr
    Integer
    size(
            T const & matrix
    )
    {
        return T::CompileTimeTraits::SizeAtCompileTime;
    }

    template<lolita::numerics::RealConcept _T, auto... _args>
    static inline
    lolita::matrix::Matrix<_T, _args...>
    Zero()
    {
        return lolita::matrix::Matrix<_T, _args...>::Zero();
    }

    template<typename _Matrix>
    static inline
    _Matrix
    Zero()
    {
        return _Matrix::Zero();
    }

    template<typename T, Integer _rows, Integer _cols>
    concept MatrixConcept2 = MatrixConcept<T> && rows<T>() == _rows && cols<T>() == _cols;

    template<typename T, Integer _rows>
    concept VectorConcept2 = MatrixConcept<T> && rows<T>() == _rows && cols<T>() == 1;

}

#endif //LOLITA_LOLITA_ALGEBRA_HXX
