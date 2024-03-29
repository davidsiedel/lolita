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
// #include <Eigen/UmfPackSupport>
#include <Eigen/Sparse>

#include "lolita_lolita/utility.hxx"

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

        namespace detail
        {

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

            template<typename t_T, auto... t_args>
            struct VectorTraits;

            template<typename t_T, auto... t_args>
            struct MatrixTraits;

            template<typename t_T>
            struct VectorTraits<t_T>
            {

                using type = Eigen::Matrix<t_T, -1, 1, getStorageOption(-1, 1)>;

            };

            template<typename t_T, auto t_rows>
            struct VectorTraits<t_T, t_rows>
            {

                using type = Eigen::Matrix<t_T, t_rows, 1, getStorageOption(t_rows, 1)>;

            };

            template<typename t_T>
            struct MatrixTraits<t_T>
            {

                using type = Eigen::Matrix<t_T, -1, -1, getStorageOption(-1, -1)>;

            };

            template<typename t_T, StorageOption t_storage_option>
            struct MatrixTraits<t_T, t_storage_option>
            {

                using type = Eigen::Matrix<t_T, -1, -1, getStorageOption(t_storage_option, -1, -1)>;

            };

            template<typename t_T, auto t_rows, auto t_cols>
            struct MatrixTraits<t_T, t_rows, t_cols>
            {

                using type = Eigen::Matrix<t_T, t_rows, t_cols, getStorageOption(t_rows, t_cols)>;

            };

            template<typename t_T, auto t_rows, auto t_cols, StorageOption t_storage_option>
            struct MatrixTraits<t_T, t_rows, t_cols, t_storage_option>
            {

                using type = Eigen::Matrix<t_T, t_rows, t_cols, getStorageOption(t_storage_option, t_rows, t_cols)>;

            };

        } // namespace detail

        template<typename t_T, auto... t_args>
        using Vector = typename detail::VectorTraits<t_T, t_args...>::type;

        template<typename t_T, auto... t_args>
        using Matrix = typename detail::MatrixTraits<t_T, t_args...>::type;

        template<typename t_T>
        using Span = Eigen::Map<t_T>;

        template<typename t_T>
        using View = Eigen::Map<t_T>;

    } // namespace matrix

    template<auto... t_args>
    using RealVector = algebra::Vector<Real, t_args...>;

    template<auto... t_args>
    using RealMatrix = algebra::Matrix<Real, t_args...>;

    template<typename t_T, auto... t_args>
    using Vector = algebra::Vector<t_T, t_args...>;

    template<typename t_T, auto... t_args>
    using Matrix = algebra::Matrix<t_T, t_args...>;

    template<typename t_T>
    concept RealMatrixConcept = std::convertible_to<t_T, Matrix<Real>>;

    template<typename t_T>
    concept RealVectorConcept = std::convertible_to<t_T, Vector<Real>>;

    auto static const print_format = Eigen::IOFormat(3, 0, ", ", "\n", "[", "]");

    static inline
    std::basic_string<Character>
    mat2str(
        auto const & mat_in,
        Integer precision = 10
    )
    {
        auto mat = Matrix<Real>(mat_in);
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
