#ifndef B8C31B4D_CEA8_4010_B9FF_DB2C50BCAE8B
#define B8C31B4D_CEA8_4010_B9FF_DB2C50BCAE8B
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
#include <Eigen/Tensor>

#include "utility.hxx"

namespace lolita
{    

    namespace algebra
    {

        // enum struct Engine
        // {
        //     Eigen,
        // };

        // template<typename Scalar_, Integer rows_, Integer cols_>
        // struct MyTensor
        // {



        //     std::array<Scalar_, rows_ * cols_> data_;

        // };

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

} // namespace lolita

#endif /* C09B9C41_77CB_4EB9_8DD3_CEAA51EBDEB7 */


#endif /* B8C31B4D_CEA8_4010_B9FF_DB2C50BCAE8B */
