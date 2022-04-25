//
// Created by dsiedel on 14/04/2022.
//

#ifndef LOLITA_LOLITA_MATRIX_HXX
#define LOLITA_LOLITA_MATRIX_HXX

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Geometry>

//#include "lolita/lolita_array.hxx"
//#include "lolita/lolita_aggregate.hxx"
#include "lolita/lolita_containers.hxx"

namespace lolita
{

    namespace matrix2
    {

        using StorageOption = Eigen::StorageOptions;

        auto const static constexpr row_major = StorageOption::RowMajor;
        auto const static constexpr col_major = StorageOption::ColMajor;

        namespace detail
        {

            template<typename T, auto... A>
            struct VectorPolicy;

            template<typename T, auto... I>
            struct MatrixPolicy;

            constexpr inline
            StorageOption
            getStorageOption(
                    Intg
                    num_rows_arg,
                    Intg
                    num_cols_arg
            )
            {
                if (num_rows_arg == 1 && num_cols_arg != 1) {
                    return row_major;
                }
                else if (num_rows_arg != 1 && num_cols_arg == 1) {
                    return col_major;
                }
                else {
                    return row_major;
                }
            }

            template<typename T>
            struct MatrixPolicy<T>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, -1, 1, getStorageOption(-1, 1)>;

            };

            template<typename T, auto R>
            struct MatrixPolicy<T, R>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, R, 1, getStorageOption(R, 1)>;

            };

            template<typename T, auto R, auto C>
            struct MatrixPolicy<T, R, C>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, R, C, getStorageOption(R, C)>;

            };

            template<typename T>
            struct VectorPolicy<T>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, -1, 1, getStorageOption(-1, 1)>;

            };

            template<typename T, auto R>
            struct VectorPolicy<T, R>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, R, 1, getStorageOption(R, 1)>;

            };

            template<typename T>
            using InputType = std::remove_reference_t<std::remove_pointer_t<T>>;

        }

        template<typename T, auto... A>
        using Matrix = typename detail::MatrixPolicy<detail::InputType<T>, A...>::Type;

        template<typename T, auto... A>
        using Vector = typename detail::VectorPolicy<detail::InputType<T>, A...>::Type;

        template<typename T, auto... A>
        using MatrixMap = Eigen::Map<detail::InputType<Eigen::Matrix<T, A...>>>;

    }

    namespace matrix
    {

        using StorageOption = Eigen::StorageOptions;

        auto const static constexpr row_major = StorageOption::RowMajor;
        auto const static constexpr col_major = StorageOption::ColMajor;

        namespace detail
        {

            template<typename T, auto... A>
            struct VectorPolicy;

            template<typename T, auto... A>
            struct MatrixPolicy;

            constexpr inline
            StorageOption
            getStorageOption(
                    Intg
                    num_rows_arg,
                    Intg
                    num_cols_arg
            )
            {
                if (num_rows_arg == 1 && num_cols_arg != 1) {
                    return row_major;
                }
                else if (num_rows_arg != 1 && num_cols_arg == 1) {
                    return col_major;
                }
                else {
                    return row_major;
                }
            }


            template<typename T, auto... A>
            struct EigenTraits
            {

                auto const static constexpr rows = Eigen::Matrix<T, A...>::CompileTimeTraits::RowsAtCompileTime;

                auto const static constexpr cols = Eigen::Matrix<T, A...>::CompileTimeTraits::ColsAtCompileTime;

                auto const static constexpr size = Eigen::Matrix<T, A...>::CompileTimeTraits::SizeAtCompileTime;

                auto const static constexpr is_empty = cols == 0 || rows == 0;

                auto const static constexpr is_static_col_vector = cols == 1 && rows > 0;

                auto const static constexpr is_static_row_vector = cols > 0 && rows == 1;

                auto const static constexpr is_dynamic_col_vector = cols == 1 && rows < 0;

                auto const static constexpr is_dynamic_row_vector = cols < 0 && rows == 1;

                auto const static constexpr is_static_vector = is_static_col_vector || is_static_row_vector;

                auto const static constexpr is_dynamic_vector = is_dynamic_col_vector || is_dynamic_row_vector;

                auto const static constexpr is_vector = is_static_vector || is_dynamic_vector;

                auto const static constexpr is_static_matrix = (cols > 1 && rows > 1) || (cols == 1 && rows == 1);

                auto const static constexpr is_dynamic_matrix = cols < 0 && rows < 0;

                auto const static constexpr is_matrix = is_static_matrix || is_dynamic_matrix;

                auto const static constexpr is_static = is_static_matrix || is_static_vector;

                auto const static constexpr is_dynamic = is_dynamic_matrix || is_dynamic_vector;

            };

            template<typename T>
            struct VectorPolicy<T>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, -1, 1, getStorageOption(-1, 1)>;

            };

            template<typename T, auto R>
            struct VectorPolicy<T, R>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, R, 1, getStorageOption(R, 1)>;

            };

            /*
             * EIGEN VIEW
             */

            template<typename T, auto... A>
            struct VectorPolicy<Eigen::Matrix<T, A...>>
            {

                static_assert(EigenTraits<T, A...>::is_vector);

                using Data = Eigen::Matrix<T, A...>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto... A>
            struct VectorPolicy<Eigen::Matrix<T, A...> const>
            {

                using Data = typename VectorPolicy<Eigen::Matrix<T, A...>>::Data;

                using Type = Eigen::Map<Data const>;

            };

//            template<typename T, auto I, auto... A>
//            struct VectorPolicy<Eigen::Matrix<T, A...>, I>
//            {
//
//                static_assert((EigenTraits<T, A...>::is_static && I <= EigenTraits<T, A...>::size) || EigenTraits<T, A...>::is_dynamic);
//
//                using Data = Eigen::Matrix<T, I, getStorageOption(I, 1)>;
//
//                using Type = Eigen::Map<Data>;
//
//            };
//
//            template<typename T, auto I, auto... A>
//            struct VectorPolicy<Eigen::Matrix<T, A...> const, I>
//            {
//
//                using Data = typename VectorPolicy<Eigen::Matrix<T, A...>, I>::Data;
//
//                using Type = Eigen::Map<Data const>;
//
//            };

            /*
             * ARRAY 2D
             */

            template<typename T, auto I, auto J>
            struct VectorPolicy<Array<T, I, J>>
            {

                using Data = Eigen::Matrix<T, I * J, 1, getStorageOption(I * J, 1)>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, auto J>
            struct VectorPolicy<Array<T, I, J> const>
            {

                using Data = typename VectorPolicy<Array<T, I, J>>::Data;

                using Type = Eigen::Map<Data const>;

            };

//            template<typename T, auto I, auto R, auto C>
//            struct VectorPolicy<Array<T, R, C>, I>
//            {
//
//                static_assert(I <= R * C);
//
//                using Data = Eigen::Matrix<T, I, 1, getStorageOption(I, 1)>;
//
//                using Type = Eigen::Map<Data>;
//
//            };
//
//            template<typename T, auto I, auto R, auto C>
//            struct VectorPolicy<Array<T, R, C> const, I>
//            {
//
//                using Data = typename VectorPolicy<Array<T, R, C>, I>::Data;
//
//                using Type = Eigen::Map<Data const>;
//
//            };

            /*
             * ARRAY
             */

            template<typename T>
            struct VectorPolicy<Array<T>>
            {

                using Data = Eigen::Matrix<T, -1, 1, getStorageOption(-1, 1)>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T>
            struct VectorPolicy<Array<T> const>
            {

                using Data = typename VectorPolicy<Array<T>>::Data;

                using Type = Eigen::Map<Data const>;

            };

//            template<typename T, auto I>
//            struct VectorPolicy<Array<T>, I>
//            {
//
//                using Data = Eigen::Matrix<T, I, 1, getStorageOption(I, 1)>;
//
//                using Type = Eigen::Map<Data>;
//
//            };
//
//            template<typename T, auto I>
//            struct VectorPolicy<Array<T> const, I>
//            {
//
//                using Data = typename VectorPolicy<Array<T>, I>::Data;
//
//                using Type = Eigen::Map<Data const>;
//
//            };









            template<typename T>
            struct MatrixPolicy<T>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, -1, -1, getStorageOption(-1, -1)>;

            };

            template<typename T, StorageOption O>
            struct MatrixPolicy<T, O>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, -1, -1, O>;

            };

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
//            template<typename T, auto R, StorageOption O>
//            struct MatrixPolicy<T, R, O>
//            {
//
//                using Data = T;
//
//                using Type = Eigen::Matrix<Data, R, 1, O>;
//
//            };

            template<typename T, auto R, auto C>
            struct MatrixPolicy<T, R, C>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, R, C, getStorageOption(R, C)>;

            };

            template<typename T, auto R, auto C, StorageOption O>
            struct MatrixPolicy<T, R, C, O>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, R, C, O>;

            };

            /*
             * EIGEN VIEW
             */

            template<typename T, auto... A>
            struct MatrixPolicy<Eigen::Matrix<T, A...>>
            {

                using Data = Eigen::Matrix<T, A...>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto... A>
            struct MatrixPolicy<Eigen::Matrix<T, A...> const>
            {

                using Data = typename MatrixPolicy<Eigen::Matrix<T, A...>>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, auto I, auto J, auto... A>
            struct MatrixPolicy<Eigen::Matrix<T, A...>, I, J>
            {

                static_assert(I * J == EigenTraits<T, A...>::size || EigenTraits<T, A...>::size == -1);

                using Data = Eigen::Matrix<T, I, J, getStorageOption(I, J)>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, auto J, auto... A>
            struct MatrixPolicy<Eigen::Matrix<T, A...> const, I, J>
            {

                using Data = typename MatrixPolicy<Eigen::Matrix<T, A...>, I, J>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, StorageOption O, auto... A>
            struct MatrixPolicy<Eigen::Matrix<T, A...>, O>
            {

            private:

                auto const static constexpr eigen_parameter = Aggregate<decltype(A)...>{A...};

            public:

                using Data = Eigen::Matrix<T, eigen_parameter.template get<0>(), eigen_parameter.template get<1>(), O>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, StorageOption O, auto... A>
            struct MatrixPolicy<Eigen::Matrix<T, A...> const, O>
            {

                using Data = typename MatrixPolicy<Eigen::Matrix<T, A...>, O>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, auto I, auto J, StorageOption O, auto... A>
            struct MatrixPolicy<Eigen::Matrix<T, A...>, I, J, O>
            {

                static_assert(I * J == Eigen::Matrix<T, A...>::CompileTimeTraits::SizeAtCompileTime);

                using Data = Eigen::Matrix<T, I, J, O>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, auto J, StorageOption O, auto... A>
            struct MatrixPolicy<Eigen::Matrix<T, A...> const, I, J, O>
            {

                using Data = typename MatrixPolicy<Eigen::Matrix<T, A...>, I, J, O>::Data;

                using Type = Eigen::Map<Data const>;

            };

            /*
             * ARRAY 2
             */

            template<typename T, auto I, auto J>
            struct MatrixPolicy<Array<T, I, J>>
            {

                using Data = Eigen::Matrix<T, I, J, getStorageOption(I, J)>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, auto J>
            struct MatrixPolicy<Array<T, I, J> const>
            {

                using Data = typename MatrixPolicy<Array<T, I, J>>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, auto I, auto J, auto R, auto C>
            struct MatrixPolicy<Array<T, R, C>, I, J>
            {

                static_assert(I * J == R * C);

                using Data = Eigen::Matrix<T, I, J, getStorageOption(I, J)>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, auto J, auto R, auto C>
            struct MatrixPolicy<Array<T, R, C> const, I, J>
            {

                using Data = typename MatrixPolicy<Array<T, R, C>, I, J>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, auto I, auto J, StorageOption O>
            struct MatrixPolicy<Array<T, I, J>, O>
            {

                using Data = Eigen::Matrix<T, I, J, O>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, auto J, StorageOption O>
            struct MatrixPolicy<Array<T, I, J> const, O>
            {

                using Data = typename MatrixPolicy<Array<T, I, J>, O>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, auto I, auto J, auto R, auto C, StorageOption O>
            struct MatrixPolicy<Array<T, R, C>, I, J, O>
            {

                static_assert(I * J == R * C);

                using Data = Eigen::Matrix<T, I, J, O>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, auto J, auto R, auto C, StorageOption O>
            struct MatrixPolicy<Array<T, R, C> const, I, J, O>
            {

                using Data = typename MatrixPolicy<Array<T, R, C>, I, J, O>::Data;

                using Type = Eigen::Map<Data const>;

            };

            /*
             * ARRAY 1
             */

            template<typename T, auto I>
            struct MatrixPolicy<Array<T, I>>
            {

                using Data = Eigen::Matrix<T, I, 1, getStorageOption(I, 1)>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I>
            struct MatrixPolicy<Array<T, I> const>
            {

                using Data = typename MatrixPolicy<Array<T, I>>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, auto I, auto J, auto R>
            struct MatrixPolicy<Array<T, R>, I, J>
            {

                static_assert(I * J == R);

                using Data = Eigen::Matrix<T, I, J, getStorageOption(I, J)>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, auto J, auto R>
            struct MatrixPolicy<Array<T, R> const, I, J>
            {

                using Data = typename MatrixPolicy<Array<T, R>, I, J>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, auto I, StorageOption O>
            struct MatrixPolicy<Array<T, I>, O>
            {

                using Data = Eigen::Matrix<T, I, 1, O>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, StorageOption O>
            struct MatrixPolicy<Array<T, I> const, O>
            {

                using Data = typename MatrixPolicy<Array<T, I>, O>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, auto I, auto J, auto R, StorageOption O>
            struct MatrixPolicy<Array<T, R>, I, J, O>
            {

                static_assert(I * J == R);

                using Data = Eigen::Matrix<T, I, J, O>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, auto J, auto R, StorageOption O>
            struct MatrixPolicy<Array<T, R> const, I, J, O>
            {

                using Data = typename MatrixPolicy<Array<T, R>, I, J, O>::Data;

                using Type = Eigen::Map<Data const>;

            };

            /*
             * ARRAY
             */

            template<typename T>
            struct MatrixPolicy<Array<T>>
            {

                using Data = Eigen::Matrix<T, -1, -1, getStorageOption(-1, -1)>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T>
            struct MatrixPolicy<Array<T> const>
            {

                using Data = typename MatrixPolicy<Array<T>>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, auto I, auto J>
            struct MatrixPolicy<Array<T>, I, J>
            {

                using Data = Eigen::Matrix<T, I, J, getStorageOption(I, J)>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, auto J>
            struct MatrixPolicy<Array<T> const, I, J>
            {

                using Data = typename MatrixPolicy<Array<T>, I, J>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, StorageOption O>
            struct MatrixPolicy<Array<T>, O>
            {

                using Data = Eigen::Matrix<T, -1, -1, O>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, StorageOption O>
            struct MatrixPolicy<Array<T> const, O>
            {

                using Data = typename MatrixPolicy<Array<T>, O>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T, auto I, auto J, StorageOption O>
            struct MatrixPolicy<Array<T>, I, J, O>
            {

                using Data = Eigen::Matrix<T, I, J, O>;

                using Type = Eigen::Map<Data>;

            };

            template<typename T, auto I, auto J, StorageOption O>
            struct MatrixPolicy<Array<T> const, I, J, O>
            {

                using Data = typename MatrixPolicy<Array<T>, I, J, O>::Data;

                using Type = Eigen::Map<Data const>;

            };

            template<typename T>
            using InputType = std::remove_reference_t<std::remove_pointer_t<T>>;

        }

        template<typename T, auto... A>
        using Matrix = typename matrix::detail::MatrixPolicy<detail::InputType<T>, A...>::Type;

        template<typename T, auto... A>
        using Vector = typename matrix::detail::VectorPolicy<detail::InputType<T>, A...>::Type;

    }

    template<typename T, auto... A>
    using Matrix = matrix::Matrix<T, A...>;

    template<typename T, auto... A>
    using Vector = matrix::Vector<T, A...>;

    namespace geometry
    {

        template<typename M>
        static inline
        auto
        getBarycenter(
                M const &
                point_args
        )
        {
            auto const constexpr rows = M::CompileTimeTraits::RowsAtCompileTime;
            auto const constexpr cols = M::CompileTimeTraits::ColsAtCompileTime;
            auto barycenter = Vector<Real, rows>();
            for (Indx i = 0; i < cols; ++i) {
                barycenter += point_args.col(i);
            }
            barycenter /= Real(cols);
            return barycenter;
        }

        template<typename M>
        static inline
        auto
        getNormalVector(
                M const &
                vector_args
        )
        {
            auto const constexpr rows = M::CompileTimeTraits::RowsAtCompileTime;
            auto const constexpr cols = M::CompileTimeTraits::ColsAtCompileTime;
            if constexpr (rows == 1 && cols == 1) {
                return vector_args;
            }
            else if constexpr (rows == 2 && cols == 1) {
                auto norm = vector_args.norm();
                return Vector<Real, 2>{vector_args(1)/norm, -vector_args(0)/norm};
            }
            else if constexpr (rows == 3 && cols == 2) {
                auto v0 = vector_args.col(0);
                auto v1 = vector_args.col(1);
                v0 /= v0.norm();
                v1 /= v1.norm();
                return v0.cross(v1);
            }
            else {
                assert(false);
            }
        }


        template<typename M>
        static inline
        auto
        getRotationMatrix(
                M const &
                point_args
        )
        {
            auto const constexpr rows = M::CompileTimeTraits::RowsAtCompileTime;
            auto const constexpr cols = M::CompileTimeTraits::ColsAtCompileTime;
            if constexpr (rows == 1 && cols == 1) {
                return point_args;
            }
            else if constexpr (rows == 2 && cols == 2) {
                auto rotation_matrix = Matrix<Real, 2, 2>();
                auto edge = point_args.col(1) - point_args.col(0);
                edge /= edge.norm();
                rotation_matrix(0, 0) = edge(0);
                rotation_matrix(0, 1) = edge(1);
                rotation_matrix(1, 0) = edge(1);
                rotation_matrix(1, 1) = -edge(0);
                return rotation_matrix;
            }
            else if constexpr (rows == 3 && cols == 3) {
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
            else {
                assert(false);
            }

        }

    }

}

#endif //LOLITA_LOLITA_MATRIX_HXX
