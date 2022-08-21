#ifndef E11520F9_5899_4E9B_8D95_F11263F73062
#define E11520F9_5899_4E9B_8D95_F11263F73062

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"

namespace lolita
{

    struct Function
    {

        Function(
            Loading const & loading,
            Integer row,
            Integer col
        )
        :
        loading_(loading),
        row_(row),
        col_(col)
        {}

        Function(
            Loading && loading,
            Integer row,
            Integer col
        )
        :
        loading_(std::forward<Loading>(loading)),
        row_(row),
        col_(col)
        {}
        
        inline
        Boolean
        operator==(
            Function const & other
        )
        const = default;
        
        inline
        Boolean
        operator!=(
            Function const & other
        )
        const = default;
        
        inline
        Real
        getImposedValue(
            Point const & point,
            Real const & time
        )
        const
        {
            return loading_(point, time);
        }
        
        inline
        Integer
        getRow()
        const
        {
            return row_;
        }
        
        inline
        Integer
        getCol()
        const
        {
            return col_;
        }

    private:

        Integer row_;

        Integer col_;

        Loading loading_;

    };

    struct Dof
    {

        Dof()
        :
        tag_(0),
        coefficients_(Vector<Real>(0))
        {}
        
        inline
        Boolean
        operator==(
            Dof const & other
        )
        const = default;
        
        inline
        Boolean
        operator!=(
            Dof const & other
        )
        const = default;

        inline
        Natural
        getTag()
        const
        {
            return tag_;
        }

        inline
        void
        setTag(
            Natural tag
        )
        {
            tag_ = tag;
        }

        inline
        Vector<Real> const &
        getCoefficients()
        const
        {
            return coefficients_;
        }

        inline
        Vector<Real> &
        getCoefficients()
        {
            return coefficients_;
        }

        Natural tag_;

        Vector<Real> coefficients_;

    };

    struct System
    {

        using MatrixEntry = Eigen::Triplet<Real>;

        static inline
        std::unique_ptr<System>
        make()
        {
            return std::make_unique<System>();
        }

        System()
        :
        rhs_values_(Vector<Real>(0))
        {}

        inline
        void
        setUnknown(
            std::shared_ptr<Dof> const & dof
        )
        {
            unknowns_.push_back(dof);
        }

        inline
        void
        setBinding(
            std::shared_ptr<Dof> const & dof
        )
        {
            bindings_.push_back(dof);
        }

        inline
        Natural
        getUnknownsSize()
        const
        {
            auto size = Natural(0);
            for (auto const & unknown : unknowns_)
            {
                size += unknown->getCoefficients().size();
            }
            return size;
        }

        inline
        Natural
        getBindingsSize()
        const
        {
            auto size = Natural(0);
            for (auto const & binding : bindings_)
            {
                size += binding->getCoefficients().size();
            }
            return size;
        }

        inline
        Natural
        getSize()
        const
        {
            return rhs_values_.size();
        }

        inline
        void
        initialize()
        {
            for (auto & unknown : unknowns_)
            {
                unknown->setTag(rhs_values_.size());
                rhs_values_.resize(rhs_values_.size() + unknown->getCoefficients().size());
            }
            for (auto & binding : bindings_)
            {
                binding->setTag(rhs_values_.size());
                rhs_values_.resize(rhs_values_.size() + binding->getCoefficients().size());
            }
        }

        void
        addLhsValue(
            Integer i,
            Integer j,
            Real value
        )
        {
            lhs_values_.push_back(MatrixEntry(i, j, value));
        }

        void
        addRhsValue(
            Integer i,
            Real value
        )
        {
            rhs_values_(i) += value;
        }

        Vector<Real>
        getCorrection()
        // const
        {
            auto lhs = Eigen::SparseMatrix<Real>(getSize(), getSize());
            lhs_values_ = std::vector<MatrixEntry>();
            rhs_values_ = Vector<Real>::Zero(getSize());
            // rhs_values_.setZero();
            for(auto i = 0; i < getSize(); i++)
            {
                addLhsValue(i, i, 1);
            }
            lhs.setFromTriplets(lhs_values_.begin(), lhs_values_.end());
            auto solver = Eigen::SparseLU<Eigen::SparseMatrix<Real>>();
            solver.analyzePattern(lhs);
            solver.factorize(lhs);
            auto x = Vector<Real>(solver.solve(rhs_values_));
            std::cout << "solution :" << std::endl;
            std::cout << x << std::endl;
            // Eigen::PardisoLDLT<Eigen::SparseMatrix<T>>  solver;

            // if (params.out_of_core >= 0 && params.out_of_core <= 2)
            //     solver.pardisoParameterArray()[59] = params.out_of_core;

            // if (params.report_factorization_Mflops)
            //     solver.pardisoParameterArray()[18] = -1; //report flops

            // solver.analyzePattern(A);
            // solver.factorize(A);
            // if (solver.info() != Eigen::Success) {
            // std::cerr << "ERROR: Could not factorize the matrix" << std::endl;
            // }

            // x = solver.solve(b);
            // if (solver.info() != Eigen::Success) {
            // std::cerr << "ERROR: Could not solve the linear system" << std::endl;
            // }

            // if (params.report_factorization_Mflops)
            // {
            //     int mflops = solver.pardisoParameterArray()[18];
            //     std::cout << "[PARDISO] Factorization Mflops: " << mflops << std::endl;
            // }
            // auto chol = Eigen::SparseLU<Eigen::SparseMatrix<Real>>(lhs);  // performs a Cholesky factorization of A
            // auto x = Vector<Real>(chol.solve(rhs_values_));         // use the factorization to solve for the given right hand side
            return x;
        }

        std::vector<std::shared_ptr<Dof>> unknowns_;

        std::vector<std::shared_ptr<Dof>> bindings_;

        std::vector<MatrixEntry> lhs_values_;

        Vector<Real> rhs_values_;

    };

    // struct DegreeOfFreedom
    // {

    //     DegreeOfFreedom(
    //         ElementType element_type,
    //         std::basic_string_view<Character> label
    //     )
    //     :
    //     element_type_(element_type),
    //     label_(label)
    //     {}
        
    //     inline
    //     Boolean
    //     operator==(
    //         DegreeOfFreedom const & other
    //     )
    //     const = default;
        
    //     inline
    //     Boolean
    //     operator!=(
    //         DegreeOfFreedom const & other
    //     )
    //     const = default;
        
    //     inline
    //     std::basic_string_view<Character>
    //     getLabel()
    //     const
    //     {
    //         return label_;
    //     }
        
    //     inline
    //     lolita::algebra::Vector<Real> const &
    //     getCoefficients()
    //     const
    //     {
    //         return coefficients_;
    //     }
        
    //     inline
    //     lolita::algebra::Vector<Real> &
    //     getCoefficients()
    //     {
    //         return coefficients_;
    //     }

    //     std::basic_string_view<Character> label_;

    //     ElementType element_type_;

    //     Vector<Real> coefficients_;

    // };

    // struct Load
    // {

    //     Load(
    //         ElementType element_type,
    //         std::basic_string_view<Character> label,
    //         Loading const & loading,
    //         Integer row,
    //         Integer col
    //     )
    //     :
    //     element_type_(element_type),
    //     label_(label),
    //     row_(row),
    //     col_(col),
    //     loading_(loading)
    //     {}

    //     Load(
    //         ElementType element_type,
    //         std::basic_string_view<Character> label,
    //         Loading && loading,
    //         Integer row,
    //         Integer col
    //     )
    //     :
    //     element_type_(element_type),
    //     label_(label),
    //     row_(row),
    //     col_(col),
    //     loading_(std::forward<Loading>(loading))
    //     {}
        
    //     inline
    //     Boolean
    //     operator==(
    //         Load const & other
    //     )
    //     const = default;
        
    //     inline
    //     Boolean
    //     operator!=(
    //         Load const & other
    //     )
    //     const = default;
        
    //     inline
    //     std::basic_string_view<Character>
    //     getLabel()
    //     const
    //     {
    //         return label_;
    //     }
        
    //     inline
    //     Integer
    //     getRow()
    //     const
    //     {
    //         return row_;
    //     }
        
    //     inline
    //     Integer
    //     getCol()
    //     const
    //     {
    //         return col_;
    //     }

    //     std::basic_string_view<Character> label_;

    //     ElementType element_type_;

    //     Integer row_;

    //     Integer col_;

    //     Loading loading_;

    // };

} // namespace lolita

#endif /* E11520F9_5899_4E9B_8D95_F11263F73062 */
