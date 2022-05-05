//
// Created by dsiedel on 04/05/22.
//

#ifndef LOLITA_LOLITA_HXX
#define LOLITA_LOLITA_HXX

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>
#include <ostream>



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

//#include "lolita_matrix.hxx"
//#include "lolita_pointers.hxx"

namespace lolita
{

    namespace config
    {

        using character =   char;

        using integer =     int;

        using index =       unsigned short;

        using natural =     unsigned long long;

        using real =        double;

        using boolean =     bool;

    }

    namespace configuration
    {

        using character =   char;

        using integer =     int;

        using index =       unsigned short;

        using natural =     unsigned long long;

        using real =        double;

        using boolean =     bool;

    }

//    namespace lolita = configuration;

    using namespace config;

//    namespace detail
//    {
//
//        template<typename T, typename... U>
//        concept IsAnyOf =   (std::same_as<T, U> || ...);
//
//        template<typename T, typename... U>
//        concept IsSameAs =  (std::same_as<T, U> && ...);
//
//    }

    namespace utility
    {

        using Label = std::array<lolita::character, 150>;

        static constexpr
        lolita::utility::Label
        makeLabel(
                std::basic_string_view<lolita::character> &&
                str
        )
        {
            auto label = lolita::utility::Label();
            for (auto i = 0; i < label.size(); ++i) {
                i < str.size() ? label[i] = str[i] : label[i] = '#';
            }
            return label;
        }

        template<typename... U>
        static constexpr
        lolita::utility::Label
        makeLabel(
//                std::basic_string_view<U> &&...
                std::basic_string_view<U>...
                str
        )
        requires(std::same_as<lolita::character, U> && ...)
        {
            auto label = lolita::utility::Label();
            auto j = lolita::index(0);
            auto make = [&] (auto && s) constexpr mutable {
                for (auto i = 0; i < s.size(); ++i) {
                    label[i + j] = s[i];
                }
                j += s.size();
            };
//            (make(std::forward<U>(str)), ...);
            (make(str), ...);
            for (auto i = j; i < label.size(); ++i) {
                label[i] = '#';
            }
            return label;
        }

        static constexpr
        std::basic_string_view<lolita::character>
        readLabel(
                lolita::utility::Label const &
                label
        )
        {
            return std::basic_string_view<lolita::character>(label.data(), std::distance(label.begin(), std::find(label.begin(), label.end(), '#')));
        }

    }

    namespace numerics
    {

        template<typename T>
        concept NaturalType =   std::same_as<std::remove_cvref_t<T>, unsigned short int>
                                || std::same_as<std::remove_cvref_t<T>, unsigned int>
                                || std::same_as<std::remove_cvref_t<T>, unsigned long int>
                                || std::same_as<std::remove_cvref_t<T>, unsigned long long int>;

        template<typename T>
        concept IntegerType =   NaturalType<T>
                                || std::same_as<std::remove_cvref_t<T>, short int>
                                || std::same_as<std::remove_cvref_t<T>, int>
                                || std::same_as<std::remove_cvref_t<T>, long int>
                                || std::same_as<std::remove_cvref_t<T>, long long int>;

        template<typename T>
        concept RealType =      IntegerType<T>
                                || std::same_as<std::remove_cvref_t<T>, float>
                                || std::same_as<std::remove_cvref_t<T>, double>
                                || std::same_as<std::remove_cvref_t<T>, long double>;

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
            return (k > n) ? 0 : (k == 0 || k == n ) ? 1 : (k == 1 || k == n-1) ? n : (k + k < n) ? (binomial(n - 1, k - 1) * n) / k : (binomial(n - 1, k) * n) / (n - k);      //  path to k=n-1 is faster
        }

    }

    namespace matrix
    {

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

//    template<typename T>
//    concept VectorType = MatrixType<T> && matrix::cols<T>() == 1;
//
//    template<typename T>
//    concept DynamicMatrixType = MatrixType<T> && (matrix::cols<T>() == -1 || matrix::rows<T>() == -1);
//
//    template<typename T>
//    concept StaticMatrixType = MatrixType<T> && (matrix::cols<T>() != -1 && matrix::rows<T>() != -1);

    template<lolita::numerics::NaturalType auto... R>
    using Vector1 = matrix::Vector<lolita::real, R...>;

    template<lolita::numerics::NaturalType auto... R>
    using Matrix1 = matrix::Matrix<lolita::real, R...>;

//    using matrix::Vector;
//
//    using matrix::Matrix;

    namespace geometry
    {

        using Point = lolita::matrix::Vector<lolita::real, 3>;

        template<lolita::matrix::MatrixType M>
        static inline
        auto
        getBarycenter(
                M const &
                point_args
        )
        requires(lolita::matrix::rows<M>() > 0)
        {
            auto barycenter = lolita::matrix::Vector<typename M::Scalar , lolita::matrix::rows<M>()>();
            for (auto i = 0; i < lolita::matrix::cols<M>(); ++i) {
                barycenter += point_args.col(i);
            }
            barycenter /= Real(lolita::matrix::cols<M>());
            return barycenter;
        }

        template<lolita::matrix::MatrixType M>
        static inline
        auto
        getNormalVector(
                M const &
                vector_args
        );

        template<lolita::matrix::MatrixType M>
        static inline
        auto
        getNormalVector(
                M const &
                vector_args
        )
        requires(lolita::matrix::rows<M>() == 1 && lolita::matrix::cols<M>() == 1)
        {
            return vector_args;
        }

        template<lolita::matrix::MatrixType M>
        static inline
        auto
        getNormalVector(
                M const &
                vector_args
        )
        requires(lolita::matrix::rows<M>() == 2 && lolita::matrix::cols<M>() == 1)
        {
            auto norm = vector_args.norm();
            return lolita::matrix::Vector<typename M::Scalar, 2>{vector_args(1)/norm, -vector_args(0)/norm};
        }

        template<lolita::matrix::MatrixType M>
        static inline
        auto
        getNormalVector(
                M const &
                vector_args
        )
        requires(lolita::matrix::rows<M>() == 3 && lolita::matrix::cols<M>() == 2)
        {
            auto v0 = vector_args.col(0);
            auto v1 = vector_args.col(1);
            v0 /= v0.norm();
            v1 /= v1.norm();
            return v0.cross(v1);
        }

        template<lolita::matrix::MatrixType M>
        static inline
        auto
        getRotationMatrix(
                M const &
                point_args
        );

        template<lolita::matrix::MatrixType M>
        static inline
        auto
        getRotationMatrix(
                M const &
                point_args
        )
        requires(lolita::matrix::rows<M>() == 1 && lolita::matrix::cols<M>() == 1)
        {
            return point_args;
        }

        template<lolita::matrix::MatrixType M>
        static inline
        auto
        getRotationMatrix(
                M const &
                point_args
        )
        requires(lolita::matrix::rows<M>() == 2 && lolita::matrix::cols<M>() == 2)
        {
            auto rotation_matrix = lolita::matrix::Matrix<typename M::Scalar, 2, 2>();
            auto edge = point_args.col(1) - point_args.col(0);
            edge /= edge.norm();
            rotation_matrix(0, 0) = edge(0);
            rotation_matrix(0, 1) = edge(1);
            rotation_matrix(1, 0) = edge(1);
            rotation_matrix(1, 1) = -edge(0);
            return rotation_matrix;
        }

        template<lolita::matrix::MatrixType M>
        static inline
        auto
        getRotationMatrix(
                M const &
                point_args
        )
        requires(lolita::matrix::rows<M>() == 3 && lolita::matrix::cols<M>() == 3)
        {
            auto rotation_matrix = lolita::matrix::Matrix<typename M::Scalar, 3, 3>();
            auto x_axis = point_args.col(2) - point_args.col(0);
            auto vector = point_args.col(1) - point_args.col(0);
            auto z_axis = x_axis.cross(vector);
            auto y_axis = z_axis.cross(x_axis);
            x_axis /= x_axis.norm();
            y_axis /= y_axis.norm();
            z_axis /= z_axis.norm();
            rotation_matrix.row(0) = x_axis;
            rotation_matrix.row(1) = y_axis;
            rotation_matrix.row(2) = z_axis;
            return rotation_matrix;
        }

        enum struct Frame
        {

            AxiSymmetric,
            Cartesian,

        };

        static constexpr
        std::basic_string_view<lolita::character>
        readFrame(
                lolita::geometry::Frame
                frame
        )
        {
            if (frame == lolita::geometry::Frame::AxiSymmetric) {
                return std::basic_string_view<lolita::character>("AxiSymmetric");
            }
            else {
                return std::basic_string_view<lolita::character>("Cartesian");
            }
        }

        struct Domain
        {

            constexpr
            Domain(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    dim,
                    lolita::geometry::Frame
                    frame
            )
            :
            tag_(lolita::utility::makeLabel(std::forward<std::basic_string_view<lolita::character>>(tag))),
            dim_(dim),
            frame_(frame)
            {}

            constexpr
            lolita::boolean
            operator==(
                    Domain const &
                    other
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    Domain const &
                    other
            )
            const = default;

            constexpr
            lolita::index
            ordIntegration(
                    lolita::index
                    ord
            )
            const
            {
                return frame_ == lolita::geometry::Frame::AxiSymmetric ? 2 * ord + 1 : 2 * ord;
            }

            lolita::utility::Label tag_;

            lolita::index dim_;

            lolita::geometry::Frame frame_;

        };

    }

//    using geometry::Point;

    namespace field
    {

        enum struct Mapping
        {

            Gradient,
            Identity,
            Divergence,
            LargeStrain,
            SmallStrain

        };

        static constexpr
        std::basic_string_view<lolita::character>
        readMapping(
                lolita::field::Mapping
                mapping
        )
        {
            if (mapping == lolita::field::Mapping::Gradient) {
                return std::basic_string_view<lolita::character>("Gradient");
            }
            else if (mapping == lolita::field::Mapping::Identity) {
                return std::basic_string_view<lolita::character>("Identity");
            }
            else if (mapping == lolita::field::Mapping::Divergence) {
                return std::basic_string_view<lolita::character>("Divergence");
            }
            else if (mapping == lolita::field::Mapping::LargeStrain) {
                return std::basic_string_view<lolita::character>("LargeStrain");
            }
            else {
                return std::basic_string_view<lolita::character>("SmallStrain");
            }
        }

        struct Tensor
        {

            constexpr
            Tensor(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    dim
            )
            :
            tag_(lolita::utility::makeLabel(std::forward<std::basic_string_view<lolita::character>>(tag))),
            ord_(dim)
            {}

            constexpr
            lolita::boolean
            operator==(
                    Tensor const &
                    other
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    Tensor const &
                    other
            )
            const = default;

            constexpr
            lolita::matrix::Cardinality
            cardinality(
                    lolita::geometry::Domain const &
                    domain
            )
            const
            {
                if (ord_ == 0) {
                    return lolita::matrix::Cardinality(lolita::numerics::pow(domain.dim_, 0), lolita::numerics::pow(domain.dim_, 0));
                }
                else if (ord_ == 1) {
                    return lolita::matrix::Cardinality(lolita::numerics::pow(domain.dim_, 1), lolita::numerics::pow(domain.dim_, 0));
                }
                else if (ord_ == 2) {
                    return lolita::matrix::Cardinality(lolita::numerics::pow(domain.dim_, 1), lolita::numerics::pow(domain.dim_, 1));
                }
                else if (ord_ == 3) {
                    return lolita::matrix::Cardinality(lolita::numerics::pow(domain.dim_, 2), lolita::numerics::pow(domain.dim_, 1));
                }
                else {
                    return lolita::matrix::Cardinality(lolita::numerics::pow(domain.dim_, 2), lolita::numerics::pow(domain.dim_, 2));
                }
            }

            lolita::utility::Label tag_;

            lolita::index ord_;

        };

        namespace detail
        {

            static constexpr
            lolita::utility::Label
            getOutputFieldLabel(
                    lolita::field::Tensor const &
                    tensor,
                    lolita::field::Mapping
                    mapping
            )
            {
                return lolita::utility::makeLabel(lolita::utility::readLabel(tensor.tag_), lolita::field::readMapping(mapping));
            }

            template<lolita::field::Tensor tensor, lolita::geometry::Domain dom, lolita::field::Mapping map>
            struct MappingPolicy;

            template<lolita::field::Tensor tensor, lolita::geometry::Domain domain, lolita::field::Mapping mapping>
            requires(mapping == lolita::field::Mapping::Identity)
            struct MappingPolicy<tensor, domain, mapping>
            {

                lolita::field::Tensor const static constexpr tensor_output_ = lolita::field::Tensor(getOutputFieldLabel(tensor, mapping), tensor.ord_);;;;;;;;;;
            };

        }

        template<typename... MappingT>
        requires((std::same_as<MappingT, lolita::field::Mapping> && ...) && sizeof...(MappingT) > 0)
        struct Unknown
        {

            constexpr
            Unknown(
                    lolita::field::Tensor const &
                    tensor,
                    MappingT &&...
                    mappings
            )
            :
            tensor_(tensor),
            mappings_({std::forward<lolita::field::Mapping>(mappings)...})
            {}

            constexpr
            Unknown(
                    lolita::field::Tensor &&
                    tensor,
                    MappingT &&...
                    mappings
            )
            :
            tensor_(std::forward<lolita::field::Tensor>(tensor)),
            mappings_({std::forward<lolita::field::Mapping>(mappings)...})
            {}

            constexpr
            Unknown(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    dim,
                    MappingT &&...
                    mappings
            )
            :
            tensor_(std::forward<lolita::field::Tensor>(Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), dim))),
            mappings_({std::forward<lolita::field::Mapping>(mappings)...})
            {}

            constexpr
            lolita::boolean
            operator==(
                    Unknown const &
                    other
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    Unknown const &
                    other
            )
            const = default;

            lolita::field::Tensor tensor_;

            std::array<lolita::field::Mapping, sizeof...(MappingT)> mappings_;

        };

        template<typename... MappingT>
        requires((std::same_as<MappingT, lolita::field::Mapping> && ...) && sizeof...(MappingT) > 0)
        struct DegreeOfFreedom : public Tensor
        {

            constexpr
            DegreeOfFreedom(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    dim,
                    MappingT &&...
                    mappings
            )
            :
            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), dim),
            mappings_({std::forward<lolita::field::Mapping>(mappings)...})
            {}

            std::array<lolita::field::Mapping, sizeof...(MappingT)> mappings_;

        };

        namespace detail
        {

            template<typename T>
            struct IsDegreeOfFreedomType : public std::false_type {};

            template<typename... MappingT>
            struct IsDegreeOfFreedomType<DegreeOfFreedom<MappingT...>> : public std::true_type {};

        }

        template<typename T>
        concept DegreeOfFreedomType = detail::IsDegreeOfFreedomType<T>::value;

        struct MaterialProperty : public Tensor
        {

            constexpr
            MaterialProperty(
                    std::basic_string_view<lolita::character> &&
                    tag
            )
            :
            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), 0)
            {}

        };

        struct InternalVariable : public Tensor
        {

            constexpr
            InternalVariable(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    ord
            )
            :
            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
            {}

        };

        struct ExternalVariable : public Tensor
        {

            constexpr
            ExternalVariable(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    ord
            )
            :
            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
            {}

        };

        struct SomeField : public Tensor
        {

            constexpr
            SomeField(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    ord
            )
            :
            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
            {}

        };

    }

    namespace behaviour
    {



    }

}

//namespace lolita::core
//{
//
//    struct Quadrature
//    {
//
//        enum Type
//        {
//
//            Gauss,
//
//        };
//
//        constexpr
//        Quadrature(
//                lolita::index
//                dim,
//                Quadrature::Type
//                typ
//        )
//        :
//        dim_(dim),
//        typ_(typ)
//        {}
//
//        constexpr
//        lolita::boolean
//        operator==(
//                Quadrature const &
//                other
//        )
//        const = default;
//
//        constexpr
//        lolita::boolean
//        operator!=(
//                Quadrature const &
//                other
//        )
//        const = default;
//
//        lolita::index dim_;
//
//        Quadrature::Type typ_;
//
//    };
//
//    struct Domain
//    {
//
//        enum Type
//        {
//
//            AxiSymmetric,
//            Cartesian,
//
//        };
//
//        constexpr
//        Domain(
//                std::basic_string_view<lolita::character> &&
//                tag,
//                lolita::index
//                dim,
//                Type
//                typ
//        )
//        :
//        tag_(utility::makeLabel(std::forward<std::basic_string_view<lolita::character>>(tag))),
//        dim_(dim),
//        typ_(typ)
//        {}
//
//        constexpr
//        lolita::boolean
//        operator==(
//                Domain const &
//                other
//        )
//        const = default;
//
//        constexpr
//        lolita::boolean
//        operator!=(
//                Domain const &
//                other
//        )
//        const = default;
//
//        constexpr
//        lolita::index
//        ordIntegration(
//                lolita::index
//                ord
//        )
//        const
//        {
//            return typ_ == Domain::AxiSymmetric ? 2 * ord + 1 : 2 * ord;
//        }
//
//        lolita::utility::Label<lolita::character> tag_;
//
//        lolita::index dim_;
//
//        Domain::Type typ_;
//
//    };
//
//    namespace detail
//    {
//
//        struct FieldBase
//        {
//
//            enum Type
//            {
//
//                Material,
//                External,
//                Internal,
//
//            };
//
//            constexpr
//            FieldBase(
//                    std::basic_string_view<character> &&
//                    tag,
//                    index
//                    dim
//            )
//            :
//            tag_(utility::makeLabel(std::forward<std::basic_string_view<character>>(tag))),
//            ord_(dim)
//            {}
//
//            constexpr
//            boolean
//            operator==(
//                    FieldBase const &
//                    other
//            )
//            const = default;
//
//            constexpr
//            boolean
//            operator!=(
//                    FieldBase const &
//                    other
//            )
//            const = default;
//
//            constexpr
//            matrix::Cardinality
//            cardinality(
//                    Domain const &
//                    domain
//            )
//            const
//            {
//                auto r = numerics::max(0, numerics::pow(domain.dim_, ord_ - 1));
//                auto c = numerics::max(0, numerics::pow(domain.dim_, ord_));
//                if (ord_ == 0) {
//                    return matrix::Cardinality(1, 1);
//                } else if (ord_ == 1) {
//                    return matrix::Cardinality(r, c);
//                } else if (ord_ == 2) {
//                    return matrix::Cardinality(c, c);
//                } else if (ord_ == 3) {
//                    return matrix::Cardinality(r, c);
//                } else {
//                    return matrix::Cardinality(numerics::pow(domain.dim_, ord_), numerics::pow(domain.dim_, ord_));
//                }
////                auto dim = domain.dim_;
////                auto pow = numerics::pow(dim, ord_);
////                auto pow1 = numerics::pow(dim, ord_ % 2);
////                auto pow2 = numerics::pow(dim, ord_ - ord_ / 2);
////                auto ren = ord_ / 2;
////                auto rem = ord_ % 2;
//////                if (ord_ == 0) {
//////
//////                }
//////                return ord_ == 0 ? : matrix::Cardinality(0, pow)
//////                return matrix::Cardinality(numerics::pow(dim, rem), numerics::pow(dim, ren));
////                return ord_ % 2 ? matrix::Cardinality(pow2, pow1) : matrix::Cardinality(numerics::pow(dim, ord_ / 2 - 1), numerics::pow(dim, ord_ / 2));
//////                return ord_ == 0 ? matrix::Cardinality(1, 1) : ord_ == 1 : matrix::Cardinality(1, pow) : ord_ == 2 : matrix::Cardinality(pow, pow);
//////                return dim == 1 ? matrix::Cardinality(1, 1) : dim == 2 ? matrix::Cardinality(1, 1)
//            }
//
//            utility::Label<character> tag_;
//
//            index ord_;
//
//        };
//
//    }
//
//    struct Mapping
//    {
//
//        enum Type
//        {
//
//            Gradient,
//            Identity,
//            Divergence,
//            LargeStrain,
//            SmallStrain
//
//        };
//
//        lolita::index const static constexpr mappings_size_ = 3;
//
//        struct MappingData
//        {
//
//            lolita::integer value_;
//
//        };
//
//        std::array<MappingData , mappings_size_> const static constexpr data_ = {
//                +1,
//                +0,
//                -1
//        };
//
//        constexpr
//        Mapping(
//                Mapping::Type
//                typ,
//                lolita::index
//                ord_field,
//                lolita::index
//                dim_field
//        )
//                :
//                typ_(typ),
//                ord_field_(ord_field),
//                dim_field_(dim_field)
//        {}
//
//        constexpr
//        lolita::boolean
//        operator==(
//                Mapping const &
//                other
//        )
//        const = default;
//
//        constexpr
//        lolita::boolean
//        operator!=(
//                Mapping const &
//                other
//        )
//        const = default;
//
//        Mapping::Type typ_;
//
//        lolita::index ord_field_;
//
//        lolita::index dim_field_;
//
//    };
//
//    template<typename... MappingT>
//    requires((std::same_as<MappingT, Mapping::Type> && ...) && sizeof...(MappingT) > 0)
//    struct DegreeOfFreedom : public detail::FieldBase
//    {
//
//        constexpr
//        DegreeOfFreedom(
//                std::basic_string_view<lolita::character> &&
//                tag,
//                lolita::index
//                dim,
//                MappingT &&...
//                mappings
//        )
//        :
//        detail::FieldBase(std::forward<std::basic_string_view<lolita::character>>(tag), dim),
//        mappings_({std::forward<Mapping::Type>(mappings)...})
//        {}
//
//        std::array<Mapping::Type, sizeof...(MappingT)> mappings_;
//
//    };
//
//    namespace detail
//    {
//
//        template<typename T>
//        struct IsDegreeOfFreedomType : public std::false_type {};
//
//        template<typename... MappingT>
//        struct IsDegreeOfFreedomType<DegreeOfFreedom<MappingT...>> : public std::true_type {};
//
//    }
//
//    template<typename T>
//    concept DegreeOfFreedomType = detail::IsDegreeOfFreedomType<T>::value;
//
//    struct MaterialProperty : public detail::FieldBase
//    {
//
//        constexpr
//        MaterialProperty(
//                std::basic_string_view<lolita::character> &&
//                tag
//        )
//        :
//        detail::FieldBase(std::forward<std::basic_string_view<lolita::character>>(tag), 0)
//        {}
//
//    };
//
//    struct InternalVariable : public detail::FieldBase
//    {
//
//        constexpr
//        InternalVariable(
//                std::basic_string_view<lolita::character> &&
//                tag,
//                lolita::index
//                ord
//        )
//        :
//        detail::FieldBase(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
//        {}
//
//    };
//
//    struct ExternalVariable : public detail::FieldBase
//    {
//
//        constexpr
//        ExternalVariable(
//                std::basic_string_view<lolita::character> &&
//                tag,
//                lolita::index
//                ord
//        )
//        :
//        detail::FieldBase(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
//        {}
//
//    };
//
//    struct SomeField : public detail::FieldBase
//    {
//
//        constexpr
//        SomeField(
//                std::basic_string_view<lolita::character> &&
//                tag,
//                lolita::index
//                ord
//        )
//        :
//        detail::FieldBase(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
//        {}
//
//    };
//
//    template<typename... FieldT>
//    requires((std::derived_from<FieldT, detail::FieldBase> && ...) && (DegreeOfFreedomType<FieldT> || ...) && sizeof...(FieldT) > 0)
//    struct BehaviourData
//    {
//
//        constexpr
//        BehaviourData(
//                std::basic_string_view<lolita::character> &&
//                tag,
//                FieldT const &...
//                fields
//        )
//        :
//        tag_(utility::makeLabel(std::forward<std::basic_string_view<lolita::character>>(tag))),
//        fields_(fields...)
//        {}
//
//        utility::Label<lolita::character> tag_;
//
//        std::tuple<FieldT...> fields_;
//
//    };
//
//    template<auto... A>
//    struct BHV
//    {
//
//        constexpr
//        BHV(
//                std::basic_string_view<lolita::character> &&
//                tag
//        )
//        :
//        tag_(utility::makeLabel(std::forward<std::basic_string_view<lolita::character>>(tag)))
//        {}
//
//        utility::Label<lolita::character> tag_;
//
//        std::tuple<decltype(A)...> const static constexpr fields_ = {A...};
//
//    };
//
//    namespace detail
//    {
//
//        struct DiscretizationBase
//        {
//
//            enum Type
//            {
//
//                Material,
//                External,
//                Internal,
//
//            };
//
//            constexpr
//            DiscretizationBase(
//                    std::basic_string_view<lolita::character> &&
//                    tag,
//                    lolita::index
//                    dim
//            )
//            :
//            tag_(utility::makeLabel(std::forward<std::basic_string_view<lolita::character>>(tag))),
//            ord_(dim)
//            {}
//
//            constexpr
//            lolita::boolean
//            operator==(
//                    DiscretizationBase const &
//                    other
//            )
//            const = default;
//
//            constexpr
//            lolita::boolean
//            operator!=(
//                    DiscretizationBase const &
//                    other
//            )
//            const = default;
//
//            utility::Label<lolita::character> tag_;
//
//            lolita::index ord_;
//
//        };
//
//    }
//
//    struct Discretization
//    {
//
//        struct HybridHighOrder
//        {
//
//            constexpr
//            HybridHighOrder(
//                    lolita::index
//                    ord_cell,
//                    lolita::index
//                    ord_face
//            )
//            :
//            ord_cell_(ord_cell),
//            ord_face_(ord_face)
//            {}
//
//            lolita::index ord_cell_;
//
//            lolita::index ord_face_;
//
//        };
//
//
//
//    };
//
//    template<auto A>
//    struct S {};
//
//    template<template<auto> typename T, auto... H>
//    static constexpr
//    std::tuple<T<H>...>
//    expand(
////            auto...
////            items
//    )
//    {
//        return std::tuple<T<H>...>();
//    }
//
//    template<typename... T>
//    struct Holder;
//
//    template<typename Ta, typename Tb>
//    struct Holder<Ta, Tb>
//    {
//
//        Ta a;
//        Tb b;
//
//        template<std::unsigned_integral auto I>
//        constexpr
//        auto
//        get()
//        {
//            if constexpr (I == 0) {
//                return a;
//            }
////            else if constexpr (I == 1) {
////                return b;
////            }
//            else {
//                return b;
//            }
//        }
//
//
//    };
//
//    template<template<auto> typename T, auto H>
//    struct Expander
//    {
//
//
//
//    };
//
//}

#endif //LOLITA_LOLITA_HXX
