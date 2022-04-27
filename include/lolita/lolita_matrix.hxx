//
// Created by dsiedel on 14/04/2022.
//

#ifndef LOLITA_LOLITA_MATRIX_HXX
#define LOLITA_LOLITA_MATRIX_HXX

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

#include "lolita/lolita_containers.hxx"

namespace lolita
{

    template<typename T>
    concept MatrixType = std::derived_from<std::remove_cvref_t<T>, Eigen::DenseBase<std::remove_cvref_t<T>>>;

    namespace matrix
    {

        using StorageOption = Eigen::StorageOptions;

        auto const static constexpr row_major = StorageOption::RowMajor;
        auto const static constexpr col_major = StorageOption::ColMajor;

        static constexpr
        auto
        storageOption(
                auto
                num_rows_arg,
                auto
                num_cols_arg
        )
        {
            return num_rows_arg != 1 && num_cols_arg == 1 ? col_major : row_major;
        }

        static constexpr
        auto
        storageOption(
                StorageOption
                storage_option,
                auto
                num_rows_arg,
                auto
                num_cols_arg
        )
        {
            return num_rows_arg != 1 && num_cols_arg == 1 ? col_major : num_rows_arg == 1 && num_cols_arg != 1 ? row_major : storage_option;
        }

        namespace detail
        {

            template<typename T, auto... A>
            struct VectorPolicy;

            template<typename T, auto... A>
            struct MatrixPolicy;

            template<typename T>
            struct VectorPolicy<T>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, -1, 1, storageOption(-1, 1)>;

            };

            template<typename T, auto R>
            struct VectorPolicy<T, R>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, R, 1, storageOption(R, 1)>;

            };

            template<typename T>
            struct MatrixPolicy<T>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, -1, -1, storageOption(-1, -1)>;

            };

            template<typename T, StorageOption O>
            struct MatrixPolicy<T, O>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, -1, -1, storageOption(O, -1, -1)>;

            };

            template<typename T, auto R, auto C>
            struct MatrixPolicy<T, R, C>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, R, C, storageOption(R, C)>;

            };

            template<typename T, auto R, auto C, StorageOption O>
            struct MatrixPolicy<T, R, C, O>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, R, C, storageOption(O, R, C)>;

            };

        }

        template<typename T, auto... A>
        using Vector = typename detail::VectorPolicy<T, A...>::Type;

        template<typename T, auto... A>
        using Matrix = typename detail::MatrixPolicy<T, A...>::Type;

        template<MatrixType T>
        using MatMap = Eigen::Map<T>;

        template<MatrixType T>
        static constexpr
        auto
        rows()
        {
            return Intg(T::CompileTimeTraits::RowsAtCompileTime);
        }

        template<MatrixType T>
        static constexpr
        auto
        cols()
        {
            return Intg(T::CompileTimeTraits::ColsAtCompileTime);
        }

        template<MatrixType T>
        static constexpr
        auto
        size()
        {
            return Intg(T::CompileTimeTraits::SizeAtCompileTime);
        }

    }

    template<typename T>
    concept VectorType = MatrixType<T> && matrix::cols<T>() == 1;

    template<typename T>
    concept DynamicMatrixType = MatrixType<T> && (matrix::cols<T>() == -1 || matrix::rows<T>() == -1);

    template<typename T>
    concept StaticMatrixType = MatrixType<T> && (matrix::cols<T>() != -1 && matrix::rows<T>() != -1);

    template<typename T, auto... A>
    using Matrix = matrix::Matrix<T, A...>;

    template<typename T, auto... A>
    using Vector = matrix::Vector<T, A...>;

    namespace geometry
    {

        template<MatrixType T>
        static inline
        auto
        getBarycenter(
                T const &
                point_args
        )
        {
            auto barycenter = Vector<Real, matrix::rows<T>()>();
            for (Indx i = 0; i < matrix::cols<T>(); ++i) {
                barycenter += point_args.col(i);
            }
            barycenter /= Real(matrix::cols<T>());
            return barycenter;
        }

        template<MatrixType M>
        static inline
        auto
        getNormalVector(
                M const &
                vector_args
        );

        template<MatrixType M>
        static inline
        auto
        getNormalVector(
                M const &
                vector_args
        )
        requires(matrix::rows<M>() == 1 && matrix::cols<M>() == 1)
        {
            return vector_args;
        }

        template<MatrixType M>
        static inline
        auto
        getNormalVector(
                M const &
                vector_args
        )
        requires(matrix::rows<M>() == 2 && matrix::cols<M>() == 1)
        {
            auto norm = vector_args.norm();
            return Vector<Real, 2>{vector_args(1)/norm, -vector_args(0)/norm};
        }

        template<MatrixType M>
        static inline
        auto
        getNormalVector(
                M const &
                vector_args
        )
        requires(matrix::rows<M>() == 3 && matrix::cols<M>() == 2)
        {
            auto v0 = vector_args.col(0);
            auto v1 = vector_args.col(1);
            v0 /= v0.norm();
            v1 /= v1.norm();
            return v0.cross(v1);
        }


        template<MatrixType M>
        static inline
        auto
        getRotationMatrix(
                M const &
                point_args
        );

        template<MatrixType M>
        static inline
        auto
        getRotationMatrix(
                M const &
                point_args
        )
        requires(matrix::rows<M>() == 1 && matrix::cols<M>() == 1)
        {
            return point_args;
        }

        template<MatrixType M>
        static inline
        auto
        getRotationMatrix(
                M const &
                point_args
        )
        requires(matrix::rows<M>() == 2 && matrix::cols<M>() == 2)
        {
            auto rotation_matrix = Matrix<Real, 2, 2>();
            auto edge = point_args.col(1) - point_args.col(0);
            edge /= edge.norm();
            rotation_matrix(0, 0) = edge(0);
            rotation_matrix(0, 1) = edge(1);
            rotation_matrix(1, 0) = edge(1);
            rotation_matrix(1, 1) = -edge(0);
            return rotation_matrix;
        }

        template<MatrixType M>
        static inline
        auto
        getRotationMatrix(
                M const &
                point_args
        )
        requires(matrix::rows<M>() == 3 && matrix::cols<M>() == 3)
        {
            auto rotation_matrix = Matrix<Real, 3, 3>();
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

//        template<typename M>
//        static inline
//        auto
//        getBarycenter(
//                M const &
//                point_args
//        )
//        {
//            auto const constexpr rows = M::CompileTimeTraits::RowsAtCompileTime;
//            auto const constexpr cols = M::CompileTimeTraits::ColsAtCompileTime;
//            auto barycenter = Vector<Real, rows>();
//            for (Indx i = 0; i < cols; ++i) {
//                barycenter += point_args.col(i);
//            }
//            barycenter /= Real(cols);
//            return barycenter;
//        }
//
//        template<typename M>
//        static inline
//        auto
//        getNormalVector(
//                M const &
//                vector_args
//        )
//        {
//            auto const constexpr rows = M::CompileTimeTraits::RowsAtCompileTime;
//            auto const constexpr cols = M::CompileTimeTraits::ColsAtCompileTime;
//            if constexpr (rows == 1 && cols == 1) {
//                return vector_args;
//            }
//            else if constexpr (rows == 2 && cols == 1) {
//                auto norm = vector_args.norm();
//                return Vector<Real, 2>{vector_args(1)/norm, -vector_args(0)/norm};
//            }
//            else if constexpr (rows == 3 && cols == 2) {
//                auto v0 = vector_args.col(0);
//                auto v1 = vector_args.col(1);
//                v0 /= v0.norm();
//                v1 /= v1.norm();
//                return v0.cross(v1);
//            }
//            else {
//                assert(false);
//            }
//        }
//
//
//        template<typename M>
//        static inline
//        auto
//        getRotationMatrix(
//                M const &
//                point_args
//        )
//        {
//            auto const constexpr rows = M::CompileTimeTraits::RowsAtCompileTime;
//            auto const constexpr cols = M::CompileTimeTraits::ColsAtCompileTime;
//            if constexpr (rows == 1 && cols == 1) {
//                return point_args;
//            }
//            else if constexpr (rows == 2 && cols == 2) {
//                auto rotation_matrix = Matrix<Real, 2, 2>();
//                auto edge = point_args.col(1) - point_args.col(0);
//                edge /= edge.norm();
//                rotation_matrix(0, 0) = edge(0);
//                rotation_matrix(0, 1) = edge(1);
//                rotation_matrix(1, 0) = edge(1);
//                rotation_matrix(1, 1) = -edge(0);
//                return rotation_matrix;
//            }
//            else if constexpr (rows == 3 && cols == 3) {
//                auto rotation_matrix = Matrix<Real, 3, 3>();
//                auto x_axis = point_args.col(2) - point_args.col(0);
//                auto vector = point_args.col(1) - point_args.col(0);
//                auto z_axis = x_axis.cross(vector);
//                auto y_axis = z_axis.cross(x_axis);
//                x_axis /= x_axis.norm();
//                y_axis /= y_axis.norm();
//                z_axis /= z_axis.norm();
//                rotation_matrix.row(0) = x_axis;
//                rotation_matrix.row(1) = y_axis;
//                rotation_matrix.row(2) = z_axis;
//                return rotation_matrix;
//            }
//            else {
//                assert(false);
//            }
//        }

    }

//    namespace matrix2
//    {
//
//        using StorageOption = Eigen::StorageOptions;
//
//        auto const static constexpr row_major = StorageOption::RowMajor;
//        auto const static constexpr col_major = StorageOption::ColMajor;
//
//        namespace detail
//        {
//
//            template<typename T, auto... A>
//            struct VectorPolicy;
//
//            template<typename T, auto... I>
//            struct MatrixPolicy;
//
//            constexpr inline
//            StorageOption
//            getStorageOption(
//                    Intg
//                    num_rows_arg,
//                    Intg
//                    num_cols_arg
//            )
//            {
//                if (num_rows_arg == 1 && num_cols_arg != 1) {
//                    return row_major;
//                }
//                else if (num_rows_arg != 1 && num_cols_arg == 1) {
//                    return col_major;
//                }
//                else {
//                    return row_major;
//                }
//            }
//
//            template<typename T>
//            struct MatrixPolicy<T>
//            {
//
//                using Data = T;
//
//                using Type = Eigen::Matrix<Data, -1, 1, getStorageOption(-1, 1)>;
//
//            };
//
//            template<typename T, auto R>
//            struct MatrixPolicy<T, R>
//            {
//
//                using Data = T;
//
//                using Type = Eigen::Matrix<Data, R, 1, getStorageOption(R, 1)>;
//
//            };
//
//            template<typename T, auto R, auto C>
//            struct MatrixPolicy<T, R, C>
//            {
//
//                using Data = T;
//
//                using Type = Eigen::Matrix<Data, R, C, getStorageOption(R, C)>;
//
//            };
//
//            template<typename T>
//            struct VectorPolicy<T>
//            {
//
//                using Data = T;
//
//                using Type = Eigen::Matrix<Data, -1, 1, getStorageOption(-1, 1)>;
//
//            };
//
//            template<typename T, auto R>
//            struct VectorPolicy<T, R>
//            {
//
//                using Data = T;
//
//                using Type = Eigen::Matrix<Data, R, 1, getStorageOption(R, 1)>;
//
//            };
//
//            template<typename T>
//            using InputType = std::remove_reference_t<std::remove_pointer_t<T>>;
//
//        }
//
//        template<typename T, auto... A>
//        using Matrix = typename detail::MatrixPolicy<detail::InputType<T>, A...>::Type;
//
//        template<typename T, auto... A>
//        using Vector = typename detail::VectorPolicy<detail::InputType<T>, A...>::Type;
//
//        template<typename T, auto... A>
//        using MatrixMap = Eigen::Map<detail::InputType<Eigen::Matrix<T, A...>>>;
//
//    }
//
//

}

#endif //LOLITA_LOLITA_MATRIX_HXX
