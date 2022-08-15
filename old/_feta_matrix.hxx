//
// Created by dsiedel on 24/03/2022.
//

#ifndef FETA__FETA__MATRIX_HXX
#define FETA__FETA__MATRIX_HXX

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Geometry>

#include "new/_feta_array.hxx"

namespace feta
{

    namespace matrix
    {

        using StorageOption = Eigen::StorageOptions;

        auto const static constexpr row_major = StorageOption::RowMajor;
        auto const static constexpr col_major = StorageOption::ColMajor;

        namespace detail
        {

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
                } else if (num_rows_arg != 1 && num_cols_arg == 1) {
                    return col_major;
                } else {
                    return row_major;
                }
            }

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

            template<typename T, auto R>
            struct MatrixPolicy<T, R>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, R, 1, getStorageOption(R, 1)>;

            };

            template<typename T, auto R, StorageOption O>
            struct MatrixPolicy<T, R, O>
            {

                using Data = T;

                using Type = Eigen::Matrix<Data, R, 1, O>;

            };

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

                static_assert(I * J == Eigen::Matrix<T, A...>::CompileTimeTraits::SizeAtCompileTime);

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

        }

        template<typename T, auto... A>
        using Matrix = typename matrix::detail::MatrixPolicy<std::remove_reference_t<std::remove_pointer_t<T>>, A...>::Type;

    }

    template<typename T, auto... A>
    using Matrix = matrix::Matrix<T, A...>;

}

#endif //FETA__FETA__MATRIX_HXX
