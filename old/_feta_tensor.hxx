//
// Created by dsiedel on 27/03/2022.
//

#ifndef FETA__FETA_TENSOR_HXX
#define FETA__FETA_TENSOR_HXX

#include "new/_feta_matrix.hxx"

namespace feta
{

    namespace tensor
    {

        template<typename T, auto... A>
        struct Tensor;

        namespace detail
        {

            template<typename U, typename T, auto... A>
            struct TensorAssigmentPolicy;

            template<typename U, typename T, auto I, auto J, auto K>
            struct TensorAssigmentPolicy<Tensor<T, I, J>, Tensor<U, K>>
            {

                static_assert(K == I * J);

                static constexpr
                void
                setData(
                        Tensor<U, K> const &
                        other,
                        Tensor<T, I, J> &
                        tensor
                )
                {
                    tensor.asVector() = other.asVector();
                }

            };

            template<typename T, auto I, auto J>
            struct TensorAssigmentPolicy<Tensor<T, I, J>, T>
            {

                static_assert(I * J == 1);

                static constexpr
                void
                setData(
                        T const &
                        value,
                        Tensor<T, I, J> &
                        tensor
                )
                {
                    tensor.data[0] = value;
                }

            };

            template<typename U, typename T, auto I, auto J>
            struct TensorAssigmentPolicy<Tensor<T, I, J>, U>
            {

                static_assert(std::is_base_of_v<Eigen::MatrixBase<U>, U>);

                static constexpr
                void
                setData(
                        U const &
                        other,
                        Tensor<T, I, J> &
                        tensor
                )
                {
                    auto constexpr rows = Eigen::MatrixBase<U>::RowsAtCompileTime;
                    auto constexpr cols = Eigen::MatrixBase<U>::ColsAtCompileTime;
                    auto constexpr storage_option = Eigen::MatrixBase<U>::IsRowMajor;
                    if constexpr (rows != -1 && cols != -1) {
                        static_assert((rows * cols == I * J) || rows * cols == - I * J);
                        Eigen::Map<Eigen::Matrix<T, rows, cols, storage_option>>(tensor.data.data()) = other;
                    } else {
                        assert(other.size() == I * J);
                        Indx nr = other.rows();
                        Indx nc = other.cols();
                        Eigen::Map<Eigen::Matrix<T, -1, -1, storage_option>>(tensor.data.data(), nr, nc) = other;
                    }
                }

            };

            template<typename U, typename T, auto I, auto R, auto C>
            struct TensorAssigmentPolicy<Tensor<T, I>, Tensor<U, R, C>>
            {

                static_assert(I == R * C);

                static constexpr
                void
                setData(
                        Tensor<U, R, C> const &
                        other,
                        Tensor<T, I> &
                        tensor
                )
                {
                    tensor.asVector() = other.asVector();
                }

            };

            template<typename T, auto I>
            struct TensorAssigmentPolicy<Tensor<T, I>, T>
            {

                static_assert(I == 1);

                static constexpr
                void
                setData(
                        T const &
                        value,
                        Tensor<T, I> &
                        tensor
                )
                {
                    tensor.data[0] = value;
                }

            };

            template<typename U, typename T, auto I>
            struct TensorAssigmentPolicy<Tensor<T, I>, U>
            {

                static_assert(std::is_base_of_v<Eigen::MatrixBase<U>, U>);

                static constexpr
                void
                setData(
                        U const &
                        other,
                        Tensor<T, I> &
                        tensor
                )
                {
                    auto constexpr rows = Eigen::MatrixBase<U>::RowsAtCompileTime;
                    auto constexpr cols = Eigen::MatrixBase<U>::ColsAtCompileTime;
                    auto constexpr storage_option = Eigen::MatrixBase<U>::IsRowMajor;
                    if constexpr (rows != -1 && cols != -1) {
                        static_assert((rows * cols == I) || rows * cols == - I);
                        Eigen::Map<Eigen::Matrix<T, rows, cols, storage_option>>(tensor.data.data()) = other;
                    } else {
                        assert(other.size() == I);
                        Indx nr = other.rows();
                        Indx nc = other.cols();
                        Eigen::Map<Eigen::Matrix<T, -1, -1, storage_option>>(tensor.data.data(), nr, nc) = other;
                    }
                }

            };

        }

        template<typename T, auto I, auto J>
        struct Tensor<T, I, J> : public Array<T, I, J>
        {

            constexpr
            Tensor() = default;

            constexpr
            Tensor(
                    auto &&...
                    values
            )
            {
                static_assert(sizeof...(values) == I * J);
                this->data = {std::forward<T>(values)...};
            }

            constexpr
            Tensor(
                    auto &&
                    other
            )
            {
                setData(other);
            }

            constexpr
            Tensor(
                    auto const &
                    other
            )
            {
                setData(other);
            }

            constexpr
            Tensor &
            operator=(
                    auto &&
                    other
            )
            {
                setData(other);
                return * this;
            }

            constexpr
            Tensor &
            operator=(
                    auto const &
                    other
            )
            {
                setData(other);
                return * this;
            }

            template<auto N, auto M>
            auto
            asMatrix()
            {
                static_assert(N * M == I * J);
                return Matrix<matrix::Matrix<T, N, M>>(this->data.data());
            }

            template<auto N, auto M>
            auto
            asMatrix()
            const
            {
                static_assert(N * M == I * J);
                return Matrix<matrix::Matrix<T, N, M> const>(this->data.data());
            }

            auto
            asMatrix()
            {
                return Matrix<matrix::Matrix<T, I, J>>(this->data.data());
            }

            auto
            asMatrix()
            const
            {
                return Matrix<matrix::Matrix<T, I, J> const>(this->data.data());
            }

            auto
            asVector()
            {
                return Matrix<matrix::Matrix<T, I * J>>(this->data.data());
            }

            auto
            asVector()
            const
            {
                return Matrix<matrix::Matrix<T, I * J> const>(this->data.data());
            }

        private:

            void constexpr
            setData(
                    auto const &
                    other
            )
            {
                detail::TensorAssigmentPolicy<Tensor<T, I, J>, RawType<decltype(other)>>::setData(other, * this);
            }

        };

        template<typename T, auto I>
        struct Tensor<T, I> : public Array<T, I>
        {

            constexpr
            Tensor() = default;

            constexpr
            Tensor(
                    auto &&...
                    values
            )
            {
                static_assert(sizeof...(values) == I);
                this->data = {std::forward<T>(values)...};
            }

            constexpr
            Tensor(
                    auto &&
                    other
            )
            {
                setData(other);
            }

            constexpr
            Tensor(
                    auto const &
                    other
            )
            {
                setData(other);
            }

            constexpr
            Tensor &
            operator=(
                    auto &&
                    other
            )
            {
                setData(other);
                return * this;
            }

            constexpr
            Tensor &
            operator=(
                    auto const &
                    other
            )
            {
                setData(other);
                return * this;
            }

            template<auto N, auto M>
            auto
            asMatrix()
            {
                static_assert(N * M == I);
                return Matrix<matrix::Matrix<T, N, M>>(this->data.data());
            }

            template<auto N, auto M>
            auto
            asMatrix()
            const
            {
                static_assert(N * M == I);
                return Matrix<matrix::Matrix<T, N, M> const>(this->data.data());
            }

            auto
            asVector()
            {
                return Matrix<matrix::Matrix<T, I>>(this->data.data());
            }

            auto
            asVector()
            const
            {
                return Matrix<matrix::Matrix<T, I> const>(this->data.data());
            }

        private:

            void constexpr
            setData(
                    auto const &
                    other
            )
            {
                detail::TensorAssigmentPolicy<Tensor<T, I>, RawType < decltype(other)>>::setData(other, * this);
            }

        };

    }

    template<typename T, auto... A>
    using Tensor = feta::tensor::Tensor<T, A...>;

}

#endif //FETA__FETA_TENSOR_HXX
