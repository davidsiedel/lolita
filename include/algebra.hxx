#ifndef C09B9C41_77CB_4EB9_8DD3_CEAA51EBDEB7
#define C09B9C41_77CB_4EB9_8DD3_CEAA51EBDEB7

#define EIGEN_USE_MKL_ALL

//#define EIGEN_USE_LAPACKE
//#define EIGEN_USE_MKL_VML

// #define EIGEN_USE_BLAS
// #define EIGEN_USE_MKL

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/PardisoSupport>
#include <Eigen/Sparse>

#include "utility.hxx"

namespace lolita
{    

    namespace algebra
    {

        using StorageOption = Eigen::StorageOptions;

        static constexpr
        StorageOption
        rowMajor()
        {
            return StorageOption::RowMajor;
        }

        static constexpr
        StorageOption
        colMajor()
        {
            return StorageOption::ColMajor;
        }

        static constexpr
        StorageOption
        getStorageOption(
            Integer num_rows,
            Integer num_cols
        )
        {
            return num_rows != 1 && num_cols == 1 ? colMajor() : rowMajor();
        }

        static constexpr
        StorageOption
        getStorageOption(
            StorageOption storage_option,
            Integer num_rows,
            Integer num_cols
        )
        {
            return num_rows != 1 && num_cols == 1 ? colMajor() : num_rows == 1 && num_cols != 1 ? rowMajor() : storage_option;
        }

        template<typename t_Scalar, auto... t_args>
        struct VectorTraits;

        template<typename t_Scalar, auto... t_args>
        struct MatrixTraits;

        template<typename t_Scalar>
        struct VectorTraits<t_Scalar>
        {

            using type = Eigen::Matrix<t_Scalar, -1, 1, getStorageOption(-1, 1)>;

        };

        template<typename t_Scalar, auto t_rows>
        struct VectorTraits<t_Scalar, t_rows>
        {

            using type = Eigen::Matrix<t_Scalar, t_rows, 1, getStorageOption(t_rows, 1)>;

        };

        template<typename t_Scalar>
        struct MatrixTraits<t_Scalar>
        {

            using type = Eigen::Matrix<t_Scalar, -1, -1, getStorageOption(-1, -1)>;

        };

        template<typename t_Scalar, StorageOption t_storage_option>
        struct MatrixTraits<t_Scalar, t_storage_option>
        {

            using type = Eigen::Matrix<t_Scalar, -1, -1, getStorageOption(t_storage_option, -1, -1)>;

        };

        template<typename t_Scalar, auto t_rows, auto t_cols>
        struct MatrixTraits<t_Scalar, t_rows, t_cols>
        {

            using type = Eigen::Matrix<t_Scalar, t_rows, t_cols, getStorageOption(t_rows, t_cols)>;

        };

        template<typename t_Scalar, auto t_rows, auto t_cols, StorageOption t_storage_option>
        struct MatrixTraits<t_Scalar, t_rows, t_cols, t_storage_option>
        {

            using type = Eigen::Matrix<t_Scalar, t_rows, t_cols, getStorageOption(t_storage_option, t_rows, t_cols)>;

        };

        template<typename t_Scalar, auto... t_args>
        using Vector = typename VectorTraits<t_Scalar, t_args...>::type;

        template<typename t_Scalar, auto... t_args>
        using Matrix = typename MatrixTraits<t_Scalar, t_args...>::type;
        
        template<typename t_T, typename t_Scalar, Integer t_rows, Integer t_cols>
        struct IsMatrixTraits
        {

            Boolean static constexpr value = std::convertible_to<t_T, Matrix<t_Scalar>> && t_T::RowsAtCompileTime == t_rows && t_T::ColsAtCompileTime == t_cols;

        };

        template<typename t_T, typename t_Scalar>
        struct IsMatrixTraits<t_T, t_Scalar, -1, -1>
        {

            Boolean static constexpr value = std::convertible_to<t_T, Matrix<t_Scalar>>;

        };

        template<typename t_T, typename t_Scalar, Integer t_i>
        struct IsMatrixTraits<t_T, t_Scalar, -1, t_i>
        {

            Boolean static constexpr value = std::convertible_to<t_T, Matrix<t_Scalar>>;

        };

        template<typename t_T, typename t_Scalar, Integer t_i>
        struct IsMatrixTraits<t_T, t_Scalar, t_i, -1>
        {

            Boolean static constexpr value = std::convertible_to<t_T, Matrix<t_Scalar>>;

        };

        template<typename t_T, typename t_Scalar, Integer t_rows = -1, Integer t_cols = -1>
        concept MatrixConcept = IsMatrixTraits<t_T, t_Scalar, t_rows, t_cols>::value;

        template<typename t_T, typename t_Scalar, Integer t_rows>
        struct IsVectorTraits
        {

            Boolean static constexpr value = std::convertible_to<t_T, Matrix<t_Scalar>> && t_T::RowsAtCompileTime == t_rows && t_T::ColsAtCompileTime == 1;

        };

        template<typename t_T, typename t_Scalar>
        struct IsVectorTraits<t_T, t_Scalar, -1>
        {

            Boolean static constexpr value = std::convertible_to<t_T, Matrix<t_Scalar>> && t_T::ColsAtCompileTime == 1;

        };

        template<typename t_T, typename t_Scalar, Integer t_rows = -1>
        concept VectorConcept = IsVectorTraits<t_T, t_Scalar, t_rows>::value;

        template<typename t_Matrix>
        using Span = Eigen::Map<t_Matrix>;

        template<typename t_Matrix>
        using View = Eigen::Map<t_Matrix>;

    } // namespace matrix

    template<typename t_Scalar, auto... t_args>
    using DenseVector = algebra::Vector<t_Scalar, t_args...>;

    template<typename t_Scalar, auto... t_args>
    using DenseMatrix = algebra::Matrix<t_Scalar, t_args...>;

    template<typename t_T, typename t_Scalar, Integer t_rows = -1, Integer t_cols = -1>
    concept DenseMatrixConcept = algebra::MatrixConcept<t_T, t_Scalar, t_rows, t_cols>;

    template<typename t_T, typename t_Scalar, Integer t_rows = -1>
    concept DenseVectorConcept = algebra::VectorConcept<t_T, t_Scalar, t_rows>;

    template<typename t_Scalar, Integer...>
    struct StaticCondensationImplementation;

    struct StaticCondensationInterface
    {

    protected:

        StaticCondensationInterface()
        {}

    public:

        template<typename t_Scalar, Integer... t_system_size>
        auto
        packUp(
            auto &&... args
        )
        {
            return static_cast<StaticCondensationImplementation<t_Scalar, t_system_size...> *>(this)->packUp(std::forward<decltype(args)>(args)...);
        }

        template<typename t_Scalar, Integer... t_system_size>
        auto
        unPack(
            auto &&... args
        )
        {
            return static_cast<StaticCondensationImplementation<t_Scalar, t_system_size...> *>(this)->unPack(std::forward<decltype(args)>(args)...);
        }

    };

    template<typename t_Scalar, Integer t_system_size, Integer t_block_size>
    struct StaticCondensationImplementation<t_Scalar, t_system_size, t_block_size> : StaticCondensationInterface
    {

    private:

        static constexpr
        Integer
        getSize()
        {
            return t_system_size;
        }

        static constexpr
        Integer
        getFTSize()
        {
            return t_system_size - t_block_size;
        }

        static constexpr
        Integer
        getTTSize()
        {
            return t_block_size;
        }

        static
        DenseMatrix<t_Scalar, getTTSize(), getTTSize()>
        make(
            DenseMatrixConcept<t_Scalar> auto && input
        )
        {
            return input.llt().solve(DenseMatrix<t_Scalar, getTTSize(), getTTSize()>::Identity());
        }

        using Base_ = StaticCondensationInterface;

    public:

        StaticCondensationImplementation()
        :
        Base_()
        {}

        StaticCondensationImplementation(
            DenseMatrixConcept<t_Scalar> auto const & matrix,
            DenseVectorConcept<t_Scalar> auto const & vector
        )
        :
        Base_(),
        r_t_(vector.template segment<getTTSize()>(0)),
        k_tf_(matrix.template block<getTTSize(), getFTSize()>(0, getTTSize())),
        k_tt_(make(matrix.template block<getTTSize(), getTTSize()>(0, 0))),
        r_c_(vector.template segment<getFTSize()>(getTTSize()) - matrix.template block<getFTSize(), getTTSize()>(getTTSize(), 0) * k_tt_ * r_t_),
        k_c_(matrix.template block<getFTSize(), getFTSize()>(getTTSize(), getTTSize()) - matrix.template block<getFTSize(), getTTSize()>(getTTSize(), 0) * k_tt_ * k_tf_)
        {}

    private:

        DenseVector<t_Scalar, getTTSize()> r_t_;

        DenseMatrix<t_Scalar, getTTSize(), getFTSize()> k_tf_;
        
        DenseMatrix<t_Scalar, getTTSize(), getTTSize()> k_tt_;

        //

        DenseVector<t_Scalar, getFTSize()> r_c_;

        DenseMatrix<t_Scalar, getFTSize(), getFTSize()> k_c_;
        
    };

    auto static const print_format = Eigen::IOFormat(3, 0, ", ", "\n", "[", "]");

    static inline
    std::basic_string<Character>
    mat2str(
        auto const & mat_in,
        Integer precision = 10
    )
    {
        auto mat = DenseMatrix<Real>(mat_in);
        auto out = std::basic_stringstream<Character>();
        out << std::fixed << std::setprecision(precision);
        if (mat.rows() > 1 && mat.cols() > 1)
        {
            for (auto i = 0; i < mat.rows(); i++)
            {
                out << "[";
                for (auto j = 0; j < mat.cols(); j++)
                {
                    if (mat(i, j) < -1.e-10)
                    {
                        out << std::setprecision(precision) << Real(mat(i, j));
                    }
                    else if (mat(i, j) > 1.e-10)
                    {
                        out << " " << std::setprecision(precision) << Real(mat(i, j));
                    }
                    else
                    {
                        out << " 0.";
                        for (auto k = 0; k < precision; k++)
                        {
                            out << " ";
                        }
                    }
                    if (j < mat.cols() - 1)
                    {
                        out << " ";
                    }
                }
                out << "]\n";
            }
            out << "\n";
        }
        else
        {
            out << "[";
            for (auto j = 0; j < mat.size(); j++)
            {
                if (mat(j) < -1.e-10)
                {
                    out << std::setprecision(precision) << Real(mat(j));
                }
                else if (mat(j) > 1.e-10)
                {
                    out << " " << std::setprecision(precision) << Real(mat(j));
                }
                else
                {
                    out << " 0.";
                    for (auto k = 0; k < precision; k++)
                    {
                        out << " ";
                    }
                }
                if (j < mat.size() - 1)
                {
                    out << " ";
                }
            }
            out << "]\n";
        }
        return out.str();
    }

    // template<typename t_Scalar, Integer...>
    // struct DenseSystem;

    // template<typename t_Scalar>
    // struct DenseSystem<t_Scalar>
    // {

    //     DenseSystem(
    //         DenseVectorConcept<t_Scalar> auto const & vector,
    //         DenseMatrixConcept<t_Scalar> auto const & matrix
    //     )
    //     :
    //     vector_(vector),
    //     matrix_(matrix)
    //     {}

    //     DenseSystem(
    //         DenseVectorConcept<t_Scalar> auto && vector,
    //         DenseMatrixConcept<t_Scalar> auto && matrix
    //     )
    //     :
    //     vector_(std::move(vector)),
    //     matrix_(std::move(matrix))
    //     {}

    //     DenseVector<t_Scalar> const &
    //     getVector()
    //     const
    //     {
    //         return vector_;
    //     }

    //     DenseMatrix<t_Scalar> const &
    //     getMatrix()
    //     const
    //     {
    //         return matrix_;
    //     }

    // private:

    //     DenseVector<t_Scalar> vector_;

    //     DenseMatrix<t_Scalar> matrix_;

    // };

    // template<typename t_Scalar, Integer t_size>
    // struct DenseSystem<t_Scalar, t_size>
    // {

    //     DenseSystem(
    //         DenseVectorConcept<t_Scalar, t_size> auto const & vector,
    //         DenseMatrixConcept<t_Scalar, t_size, t_size> auto const & matrix
    //     )
    //     :
    //     vector_(vector),
    //     matrix_(matrix)
    //     {}

    //     DenseSystem(
    //         DenseVectorConcept<t_Scalar, t_size> auto && vector,
    //         DenseMatrixConcept<t_Scalar, t_size, t_size> auto && matrix
    //     )
    //     :
    //     vector_(std::move(vector)),
    //     matrix_(std::move(matrix))
    //     {}

    //     DenseVector<t_Scalar, t_size> const &
    //     getVector()
    //     const
    //     {
    //         return vector_;
    //     }

    //     DenseMatrix<t_Scalar, t_size, t_size> const &
    //     getMatrix()
    //     const
    //     {
    //         return matrix_;
    //     }

    // private:

    //     DenseVector<t_Scalar, t_size> vector_;

    //     DenseMatrix<t_Scalar, t_size, t_size> matrix_;

    // };

    // template<typename t_Scalar, Integer...>
    // struct ConcreteDenseSystem;
    
    // struct AbstractDenseSystem
    // {

    //     AbstractDenseSystem()
    //     {}

    //     virtual
    //     ~AbstractDenseSystem()
    //     {}

    //     template<typename t_Scalar, Integer... t_args>
    //     auto const &
    //     getMatrix()
    //     const
    //     {
    //         return static_cast<ConcreteDenseSystem<t_Scalar, t_args...> const *>(this)->getMatrix();
    //     }

    //     template<typename t_Scalar, Integer... t_args>
    //     auto const &
    //     getVector()
    //     const
    //     {
    //         return static_cast<ConcreteDenseSystem<t_Scalar, t_args...> const *>(this)->getVector();
    //     }

    // };

    // template<typename t_Scalar>
    // struct ConcreteDenseSystem<t_Scalar> : AbstractDenseSystem, DenseSystem<t_Scalar>
    // {

    //     ConcreteDenseSystem(
    //         DenseVectorConcept<t_Scalar> auto && vector,
    //         DenseMatrixConcept<t_Scalar> auto && matrix
    //     )
    //     :
    //     AbstractDenseSystem(),
    //     DenseSystem<t_Scalar>(std::forward<decltype(vector)>(vector), std::forward<decltype(matrix)>(matrix))
    //     {}

    // };

    // template<typename t_Scalar, Integer t_size>
    // struct ConcreteDenseSystem<t_Scalar, t_size> : AbstractDenseSystem, DenseSystem<t_Scalar, t_size>
    // {

    //     ConcreteDenseSystem(
    //         DenseVectorConcept<t_Scalar, t_size> auto && vector,
    //         DenseMatrixConcept<t_Scalar, t_size, t_size> auto && matrix
    //     )
    //     :
    //     AbstractDenseSystem(),
    //     DenseSystem<t_Scalar, t_size>(std::forward<decltype(vector)>(vector), std::forward<decltype(matrix)>(matrix))
    //     {}
        
    // };

    // template<typename, Integer...>
    // struct ConcreteDenseVector;

    // struct AbstractDenseVector
    // {

    //     AbstractDenseVector()
    //     {}

    //     virtual
    //     ~AbstractDenseVector()
    //     {}

    //     template<typename t_Scalar, Integer... t_args>
    //     auto const &
    //     get()
    //     const
    //     {
    //         return static_cast<ConcreteDenseVector<t_Scalar, t_args...> const *>(this)->get();
    //     }

    // };

    // template<typename t_Scalar, Integer t_num_cols>
    // struct ConcreteDenseVector<t_Scalar, t_num_cols> : AbstractDenseVector
    // {

    //     explicit
    //     ConcreteDenseVector(
    //         DenseVectorConcept<t_Scalar, t_num_cols> auto const & vector
    //     )
    //     :
    //     AbstractDenseVector(),
    //     vector_(vector)
    //     {}
        
    //     explicit
    //     ConcreteDenseVector(
    //         DenseVectorConcept<t_Scalar, t_num_cols> auto && vector
    //     )
    //     :
    //     AbstractDenseVector(),
    //     vector_(std::move(vector))
    //     {}

    //     DenseVector<t_Scalar, t_num_cols> const &
    //     get()
    //     const
    //     {
    //         return vector_;
    //     }

    // private:
        
    //     DenseVector<t_Scalar, t_num_cols> vector_;

    // };

    // template<typename, Integer...>
    // struct ConcreteDenseMatrix;

    // struct AbstractDenseMatrix
    // {

    //     AbstractDenseMatrix()
    //     {}

    //     virtual
    //     ~AbstractDenseMatrix()
    //     {}

    //     template<typename t_Scalar, Integer... t_system_size>
    //     auto const &
    //     get()
    //     const
    //     {
    //         return static_cast<ConcreteDenseMatrix<t_Scalar, t_system_size...> const *>(this)->get();
    //     }

    // };

    // template<typename t_Scalar, Integer t_num_rows, Integer t_num_cols>
    // struct ConcreteDenseMatrix<t_Scalar, t_num_rows, t_num_cols> : AbstractDenseMatrix
    // {

    //     explicit
    //     ConcreteDenseMatrix(
    //         DenseMatrixConcept<t_Scalar, t_num_rows, t_num_cols> auto const & matrix
    //     )
    //     :
    //     AbstractDenseMatrix(),
    //     matrix_(matrix)
    //     {}
        
    //     explicit
    //     ConcreteDenseMatrix(
    //         DenseMatrixConcept<t_Scalar, t_num_rows, t_num_cols> auto && matrix
    //     )
    //     :
    //     AbstractDenseMatrix(),
    //     matrix_(std::move(matrix))
    //     {}

    //     DenseMatrix<t_Scalar, t_num_rows, t_num_cols> const &
    //     get()
    //     const
    //     {
    //         return matrix_;
    //     }

    // private:
        
    //     DenseMatrix<t_Scalar, t_num_rows, t_num_cols> matrix_;
        
    // };

    // // template<typename, Integer...>
    // // struct StaticCondensationImpl;

    // template<typename t_Scalar, Integer...>
    // struct StaticCondensationImplementation;
    
    // struct StaticCondensation
    // {

    //     virtual
    //     ~StaticCondensation()
    //     {}

    //     template<typename t_Scalar, Integer... t_system_size>
    //     auto
    //     packUp(
    //         auto &&... args
    //     )
    //     {
    //         return static_cast<StaticCondensationImpl<t_Scalar, t_system_size...> const *>(this)->packUp(std::forward<decltype(args)>(args)...);
    //     }

    //     template<typename t_Scalar, Integer... t_system_size>
    //     auto
    //     unPack(
    //         auto &&... args
    //     )
    //     {
    //         return static_cast<StaticCondensationImpl<t_Scalar, t_system_size...> const *>(this)->unPack(std::forward<decltype(args)>(args)...);
    //     }

    // };

    
    
    // template<typename t_Scalar, Integer t_system_size, Integer t_block_size>
    // struct StaticCondensationImpl<t_Scalar, t_system_size, t_block_size> : StaticCondensationInterface
    // {

    // private:

    //     using Base_ = StaticCondensationInterface;

    //     using t_SystemIn = DenseSystem<t_Scalar, t_system_size>;

    // public:

    //     StaticCondensationImpl()
    //     :
    //     Base_()
    //     {}

    //     static constexpr
    //     Integer
    //     getSize()
    //     {
    //         return t_system_size;
    //     }

    //     static constexpr
    //     Integer
    //     getFTSize()
    //     {
    //         return t_system_size - t_block_size;
    //     }

    //     static constexpr
    //     Integer
    //     getTTSize()
    //     {
    //         return t_block_size;
    //     }

    //     static
    //     DenseMatrix<t_Scalar, getTTSize(), getTTSize()>
    //     make(
    //         DenseMatrixConcept<t_Scalar, getTTSize(), getTTSize()> auto && input
    //     )
    //     {
    //         return input.llt().solve(DenseMatrix<t_Scalar, getTTSize(), getTTSize()>::Identity());
    //     }

    //     StaticCondensationImpl(
    //         DenseVectorConcept<t_Scalar, getTTSize()> auto const & r_t,
    //         DenseMatrixConcept<t_Scalar, getTTSize(), getTTSize()> auto const & k_tf,
    //         DenseMatrixConcept<t_Scalar, getTTSize(), getFTSize()> auto const & k_tt
    //     )
    //     :
    //     StaticCondensation(),
    //     r_t_(r_t),
    //     k_tf_(k_tf),
    //     k_tt_(k_tt)
    //     {}

    //     StaticCondensationImpl(
    //         DenseVectorConcept<t_Scalar, getTTSize()> auto && r_t,
    //         DenseMatrixConcept<t_Scalar, getTTSize(), getTTSize()> auto && k_tf,
    //         DenseMatrixConcept<t_Scalar, getTTSize(), getFTSize()> auto && k_tt
    //     )
    //     :
    //     StaticCondensation(),
    //     r_t_(std::move(r_t)),
    //     k_tf_(std::move(k_tf)),
    //     k_tt_(std::move(k_tt))
    //     {}

    //     explicit
    //     StaticCondensationImpl(
    //         DenseSystem<t_Scalar, t_system_size> const & system
    //     )
    //     :
    //     StaticCondensation(),
    //     r_t_(system.getVector().template segment<getTTSize()>(0)),
    //     k_tf_(system.getMatrix().template block<getTTSize(), getFTSize()>(0, getTTSize())),
    //     k_tt_(make(system.getMatrix().template block<getTTSize(), getTTSize()>(0, 0)))
    //     {}

    //     explicit
    //     StaticCondensationImpl(
    //         DenseSystem<t_Scalar, t_system_size> && system
    //     )
    //     :
    //     StaticCondensation(),
    //     r_t_(std::move(system).getVector().template segment<getTTSize()>(0)),
    //     k_tf_(std::move(system).getMatrix().template block<getTTSize(), getFTSize()>(0, getTTSize())),
    //     k_tt_(make(std::move(system).getMatrix().template block<getTTSize(), getTTSize()>(0, 0)))
    //     {}

    //     DenseSystem<t_Scalar, getFTSize()>
    //     packUp(
    //         DenseSystem<t_Scalar, t_system_size> && system
    //     )
    //     {
    //         auto const r_f = std::forward<t_SystemIn>(system).getVector().template segment<getFTSize()>(getTTSize());
    //         auto const k_ft = std::forward<t_SystemIn>(system).getMatrix().template block<getFTSize(), getTTSize()>(getTTSize(), 0);
    //         auto const k_ff = std::forward<t_SystemIn>(system).getMatrix().template block<getFTSize(), getFTSize()>(getTTSize(), getTTSize());
    //         return DenseSystem<t_Scalar, getFTSize()>(k_ff - k_ft * k_tt_ * k_tf_, r_f - k_ft * k_tt_ * r_t_);
    //     }

    //     DenseVector<t_Scalar, getTTSize()>
    //     unPack(
    //         DenseVectorConcept<t_Scalar, getFTSize()> auto && vector
    //     )
    //     {
    //         return k_tt_ * (r_t_ - k_tf_ * std::forward<decltype(vector)>(vector));
    //     }

    // private:

    //     DenseVector<t_Scalar, getTTSize()> r_t_;

    //     DenseMatrix<t_Scalar, getTTSize(), getFTSize()> k_tf_;
        
    //     DenseMatrix<t_Scalar, getTTSize(), getTTSize()> k_tt_;

    // };

    // template<typename t_Scalar>
    // struct StaticCondensationImpl<t_Scalar> : StaticCondensation
    // {

    //     StaticCondensationImpl(
    //         Integer size,
    //         Integer size2
    //     )
    //     :
    //     StaticCondensation(),
    //     size_(size),
    //     size2_(size2)
    //     {}

    //     void
    //     make(
    //         DenseMatrixConcept<Real> auto && matrix,
    //         DenseVectorConcept<Real> auto && vector,
    //         Integer size,
    //         Integer size2
    //     )
    //     {
    //         auto const k_tt = std::forward<decltype(matrix)>(matrix).block(0, 0, getTTSize(), getTTSize());
    //         auto const k_ff = std::forward<decltype(matrix)>(matrix).block(getTTSize(), getTTSize(), getFTSize(), getFTSize());
    //         auto const r_f = std::forward<decltype(vector)>(vector).segment(getTTSize(), getFTSize());
    //         k_tt_ = DenseMatrix<t_Scalar>(k_tt.llt().solve(DenseMatrix<t_Scalar>::Identity(getTTSize(), getTTSize())));
    //         k_tf_ = std::forward<decltype(matrix)>(matrix).block(0, getTTSize(), getTTSize(), getFTSize());
    //         k_ft_ = std::forward<decltype(matrix)>(matrix).block(getTTSize(), 0, getFTSize(), getTTSize());
    //         r_t_ = std::forward<decltype(vector)>(vector).segment(0, getTTSize());
    //     }

    //     DenseSystem<t_Scalar>
    //     packUp(
    //         DenseMatrixConcept<Real> auto && matrix,
    //         DenseVectorConcept<Real> auto && vector
    //     )
    //     {
    //         auto const k_ff = std::forward<decltype(matrix)>(matrix).block(getTTSize(), getTTSize(), getFTSize(), getFTSize());
    //         auto const r_f = std::forward<decltype(vector)>(vector).segment(getTTSize(), getFTSize());
    //         return DenseSystem<t_Scalar>(k_ff - k_ft_ * k_tt_ * k_tf_, r_f - k_ft_ * k_tt_ * r_t_);
    //     }

    //     DenseVector<t_Scalar>
    //     unPack(
    //         DenseVectorConcept<t_Scalar> auto && vector
    //     )
    //     {
    //         return k_tt_ * (r_t_ - k_tf_ * std::forward<decltype(vector)>(vector));
    //     }

    // private:

    //     Integer
    //     getFTSize(
    //         Integer size,
    //         Integer size2
    //     )
    //     const
    //     {
    //         return size - size2;
    //     }

    //     Integer
    //     getTTSize(
    //         Integer size,
    //         Integer size2
    //     )
    //     const
    //     {
    //         return size2;
    //     }

    // public:

    //     DenseVector<t_Scalar> r_t_;
        
    //     DenseMatrix<t_Scalar> k_tt_;

    //     DenseMatrix<t_Scalar> k_tf_;

    //     DenseMatrix<t_Scalar> k_ft_;

    // };

} // namespace lolita

#endif /* C09B9C41_77CB_4EB9_8DD3_CEAA51EBDEB7 */
