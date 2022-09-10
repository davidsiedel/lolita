#ifndef E11520F9_5899_4E9B_8D95_F11263F73062
#define E11520F9_5899_4E9B_8D95_F11263F73062

#include "lolita.hxx"
#include "core/000_physics_traits.hxx"
#include "core/001_geometry.hxx"

namespace lolita
{
    
    struct DegreeOfFreedom
    {

        inline
        static
        std::shared_ptr<DegreeOfFreedom>
        make(
            std::basic_string<Character> && label
        )
        {
            return std::make_shared<DegreeOfFreedom>(std::forward<std::basic_string<Character>>(label));
        }

        inline
        static
        std::shared_ptr<DegreeOfFreedom>
        make(
            std::basic_string<Character> const & label
        )
        {
            return std::make_shared<DegreeOfFreedom>(label);
        }

    // private:

        DegreeOfFreedom(
            // Field field,
            std::basic_string<Character> && label
        )
        :
        // field_(field),
        label_(std::forward<std::basic_string<Character>>(label)),
        offset_(std::numeric_limits<Natural>::max()),
        size_(0)
        {}

        DegreeOfFreedom(
            // Field field,
            std::basic_string<Character> const & label
        )
        :
        // field_(field),
        label_(label),
        offset_(std::numeric_limits<Natural>::max()),
        size_(0)
        {}

    public:

        inline
        void
        setOffset(
            Integer offset
        )
        {
            offset_ = offset;
        }

        inline
        Natural
        getOffset()
        const
        {
            return offset_;
        }

        inline
        void
        setSize(
            Natural size
        )
        {
            size_ = size;
        }

        inline
        void
        addSize(
            Natural size
        )
        {
            size_ += size;
        }

        inline
        Natural
        getSize()
        const
        {
            return size_;
        }

        inline
        std::basic_string<Character> const &
        getLabel()
        {
            return label_;
        }

        // inline
        // Field
        // getField()
        // const
        // {
        //     return field_;
        // }

    private:

        // Field field_;

        Natural size_;

        Natural offset_;

        std::basic_string<Character> label_;

    };

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
            // auto val = loading_(point, time);
            // std::cout << "val : " << val << std::endl;
            // auto a = 1;
            // return val;
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

    
    struct System
    {

        std::vector<std::shared_ptr<DegreeOfFreedom>> dofss_;

        inline
        void
        addDegreeOfFreedom(
            std::shared_ptr<DegreeOfFreedom> const & degree_of_freedom
        )
        {
            if (std::find(dofss_.begin(), dofss_.end(), degree_of_freedom) == dofss_.end())
            {
                auto offset = Natural(0);
                for (auto const & dof : dofss_)
                {
                    offset += dof->getSize();
                }
                degree_of_freedom->setOffset(offset);
                dofss_.push_back(degree_of_freedom);
            }
        }

        using MatrixEntry = Eigen::Triplet<Real>;

        static inline
        std::unique_ptr<System>
        make()
        {
            return std::make_unique<System>();
        }

        System()
        :
        rhs_values_(Vector<Real>(0)),
        normalization_(1.e-14),
        band_width_(0)
        {}

        inline
        void
        setUnknown(
            std::basic_string<Character> && label,
            Natural size
        )
        {
            unknowns_[label] = size;
        }

        inline
        void
        setBinding(
            std::basic_string<Character> && label,
            Natural size
        )
        {
            bindings_[label] = size;
        }

        inline
        Natural
        getUnknownsSize()
        const
        {
            auto size = Natural(0);
            for (auto const & dof : unknowns_)
            {
                size += dof.second;
            }
            return size;
        }

        inline
        Natural
        getBindingsSize()
        const
        {
            auto size = Natural(0);
            for (auto const & dof : bindings_)
            {
                size += dof.second;
            }
            return size;
        }

        inline
        Natural
        getSize()
        const
        {
            auto size = Natural(0);
            for (auto const & dof : unknowns_)
            {
                size += dof.second;
            }
            for (auto const & dof : bindings_)
            {
                size += dof.second;
            }
            return size;
        }

        inline
        Natural
        getUnknownOffset(
            std::basic_string_view<Character> label
        )
        const
        {
            if (unknowns_.contains(std::string(label)))
            {
                auto offset = Natural(0);
                for (auto const & dof : unknowns_)
                {
                    if (dof.first == label)
                    {
                        break;
                    }
                    offset += dof.second;
                }
                return offset;
            }
            else
            {
                throw std::runtime_error("NO");
            }
        }

        inline
        Natural
        getBindingOffset(
            std::basic_string_view<Character> label
        )
        const
        {
            if (bindings_.contains(std::string(label)))
            {
                auto offset = getUnknownsSize();
                for (auto const & dof : bindings_)
                {
                    if (dof.first == label)
                    {
                        break;
                    }
                    offset += dof.second;
                }
                return offset;
            }
            else
            {
                throw std::runtime_error("NO");
            }
        }

        inline
        algebra::View<Vector<Real> const>
        getUnknownCorrection(
            std::basic_string_view<Character> label
        )
        const
        {
            auto offset = getUnknownOffset(label);
            auto size = unknowns_.at(std::string(label));
            return algebra::View<Vector<Real> const>(correction_values_.data() + offset, size);
        }

        inline
        algebra::View<Vector<Real> const>
        getBindingCorrection(
            std::basic_string_view<Character> label
        )
        const
        {
            auto offset = getBindingOffset(label);
            auto size = bindings_.at(std::string(label));
            return algebra::View<Vector<Real> const>(correction_values_.data() + offset, size);
        }

        inline
        void
        initialize()
        {
            lhs_values_.clear();
            if (rhs_values_.size() == 0)
            {
                rhs_values_ = Vector<Real>::Zero(getSize());
            }
            else
            {
                rhs_values_.setZero();
            }
            normalization_ = 1.e-14;
        }

        inline
        void
        initializeRhs()
        {
            if (rhs_values_.size() == 0)
            {
                rhs_values_ = Vector<Real>::Zero(getSize());
            }
            else
            {
                rhs_values_.setZero();
            }
        }

        inline
        void
        initializeLhs()
        {
            lhs_ = Eigen::SparseMatrix<Real>(getSize(), getSize());
            // lhs_.reserve(Vector<Real>::Constant(getSize(), band_width_));
            lhs_values_.clear();
            // lhs_values_.reserve(40'000 * 18 * 18);
        }

        inline
        void
        initializeNormalization(
            Real value = 1.e-14
        )
        {
            normalization_ = value;
        }

        // inline
        // void
        // addLhsValue(
        //     Integer i,
        //     Integer j,
        //     Real value
        // )
        // {
        //     lhs_values_.push_back(MatrixEntry(i, j, value));
        // }

        // inline
        // void
        // addRhsValue(
        //     Integer i,
        //     Real value
        // )
        // {
        //     rhs_values_(i) += value;
        // }

        inline
        void
        addLhsValue(
            Integer i,
            Integer j,
            Real value
        )
        {
            auto lock = std::scoped_lock<std::mutex>(mutex);
            lhs_values_.push_back(MatrixEntry(i, j, value));
        }

        inline
        void
        addRhsValue(
            Integer i,
            Real value
        )
        {
            auto lock = std::scoped_lock<std::mutex>(mutex);
            rhs_values_(i) += value;
        }

        template<typename t_T>
        void
        addLhsValues(
            Integer i,
            Integer j,
            t_T const & values
        )
        {
            auto lock = std::scoped_lock<std::mutex>(mutex);
            for (auto iii = 0; iii < values.rows(); iii++)
            {
                for (auto jjj = 0; jjj < values.cols(); jjj++)
                {
                    // lhs_.coeffRef(iii + i, jjj + j) += values(iii, jjj);
                    lhs_values_.push_back(MatrixEntry(iii + i, jjj + j, values(iii, jjj)));
                }
            }
        }

        template<typename t_T>
        void
        addRhsValues(
            Integer i,
            t_T const & values
        )
        {
            auto lock = std::scoped_lock<std::mutex>(mutex);
            for (auto iii = 0; iii < values.size(); iii++)
            {
                rhs_values_(iii + i) += values(iii);
            }
        }

        // template<Integer t_rows, Integer t_cols>
        // void
        // addLhsValues(
        //     Integer i,
        //     Integer j,
        //     Matrix<Real, t_rows, t_cols> const & values
        // )
        // {
        //     auto lock = std::scoped_lock<std::mutex>(mutex);
        //     for (auto iii = 0; iii < t_rows; iii++)
        //     {
        //         for (auto jjj = 0; jjj < t_cols; jjj++)
        //         {
        //             lhs_values_.push_back(MatrixEntry(iii + i, jjj + j, values(iii, jjj)));
        //         }
        //     }
        // }

        // template<Integer t_rows>
        // void
        // addRhsValues(
        //     Integer i,
        //     Vector<Real, t_rows> const & values
        // )
        // {
        //     auto lock = std::scoped_lock<std::mutex>(mutex);
        //     for (auto iii = 0; iii < t_rows; iii++)
        //     {
        //         rhs_values_(iii + i) += values(iii);
        //     }
        // }

        inline
        void
        setCorrection()
        {
            // using SOLVER = Eigen::SparseLU<Eigen::SparseMatrix<Real, Eigen::ColMajor>, Eigen::COLAMDOrdering<Eigen::DenseIndex>>;
            // using SOLVER = Eigen::SparseLU<Eigen::SparseMatrix<Real, Eigen::RowMajor>, Eigen::COLAMDOrdering<int> >;
            using SOLVER = Eigen::PardisoLU<Eigen::SparseMatrix<Real>>;
            // using SOLVER = Eigen::PardisoLDLT<Eigen::SparseMatrix<Real>>;
            // auto solver = Eigen::PardisoLU<Eigen::SparseMatrix<Real>>();
            // auto solver = Eigen::PardisoLDLT<Eigen::SparseMatrix<Real>>();
            auto solver = SOLVER();
            // OPT 1
            // auto lhs = Eigen::SparseMatrix<Real>(getSize(), getSize());
            lhs_.setFromTriplets(lhs_values_.begin(), lhs_values_.end());
            // OPT 2
            // lhs_.makeCompressed();
            // std::cout << "lhs : " << "\n";
            // std::cout << mat2str(Matrix<Real>(lhs)) << "\n";
            // std::cout << "rhs : " << "\n";
            // std::cout << mat2str(rhs_values_) << "\n";
            //
            //
            // auto outfile = std::ofstream();
            // outfile << std::fixed << std::setprecision(3);
            // outfile.open("/home/dsiedel/projetcs/lolita/lolita/tests/t0/mat.txt");
            // outfile << "lhs : " << "\n";
            // outfile << Matrix<Real>(lhs).format(print_format) << "\n";
            // outfile << "rhs : " << "\n";
            // outfile << rhs_values_ << "\n";
            //
            //
            // OPT 1
            // solver.analyzePattern(lhs);
            // solver.factorize(lhs);
            // OPT 2
            solver.analyzePattern(lhs_);
            solver.factorize(lhs_);
            if (solver.info() != Eigen::Success)
            {
                std::cerr << "ERROR: Could not factorize the matrix" << std::endl;
            }
            // x = solver.solve(b);
            // auto RHS = Vector<Real>(- rhs_values_);
            correction_values_ = solver.solve(rhs_values_);
            if (solver.info() != Eigen::Success)
            {
                std::cerr << "ERROR: Could not solve the linear system" << std::endl;
            }
            // std::cout << "correction : " << "\n";
            // std::cout << correction_values_ << "\n";
            // lhs_values_ = std::vector<MatrixEntry>();
            // rhs_values_ = Vector<Real>::Zero(getSize());
            // rhs_values_.setZero();
            // for(auto i = 0; i < getSize(); i++)
            // {
            //     addLhsValue(i, i, 1);
            // }
            // lhs.setFromTriplets(lhs_values_.begin(), lhs_values_.end());
            // std::cout << "solution :" << std::endl;
            // std::cout << correction_values_ << std::endl;
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
        }

        inline
        Real
        getNormalization()
        const
        {
            return normalization_;
        }

        inline
        Real
        getResidualInfiniteNorm()
        const
        {
            return rhs_values_.cwiseAbs().maxCoeff();
        }

        inline
        Real
        getResidualEvaluation()
        const
        {
            return getResidualInfiniteNorm() / getNormalization();
        }

        inline
        void
        setNormalization(
            Real value
        )
        {
            if (value > normalization_)
            {
                auto lock = std::scoped_lock<std::mutex>(mutex);
                normalization_ = value;
            }
        }

        inline
        void
        setBandWidth(
            Integer value
        )
        {
            if (value > band_width_)
            {
                auto lock = std::scoped_lock<std::mutex>(mutex);
                band_width_ = value;
            }
        }

        Integer band_width_;

        std::mutex mutex;

        Eigen::SparseMatrix<Real> lhs_;

        Real normalization_;

        std::map<std::basic_string<Character>, Natural> unknowns_;

        std::map<std::basic_string<Character>, Natural> bindings_;

        std::vector<MatrixEntry> lhs_values_;

        Vector<Real> rhs_values_;

        Vector<Real> correction_values_;

    };

    // struct Dof
    // {

    //     Dof()
    //     :
    //     tag_(0),
    //     coefficients_(Vector<Real>(0))
    //     {}
        
    //     inline
    //     Boolean
    //     operator==(
    //         Dof const & other
    //     )
    //     const = default;
        
    //     inline
    //     Boolean
    //     operator!=(
    //         Dof const & other
    //     )
    //     const = default;

    //     inline
    //     Natural
    //     getTag()
    //     const
    //     {
    //         return tag_;
    //     }

    //     inline
    //     void
    //     setTag(
    //         Natural tag
    //     )
    //     {
    //         tag_ = tag;
    //     }

    //     inline
    //     Vector<Real> const &
    //     getCoefficients()
    //     const
    //     {
    //         return coefficients_;
    //     }

    //     inline
    //     Vector<Real> &
    //     getCoefficients()
    //     {
    //         return coefficients_;
    //     }

    //     Natural tag_;

    //     Vector<Real> coefficients_;

    // };
    
    // struct Dof2
    // {

    //     Dof2(
    //         Natural tag,
    //         Integer dim,
    //         std::basic_string_view<Character> label
    //     )
    //     :
    //     tag_(tag),
    //     dim_(dim),
    //     label_(std::make_shared<std::basic_string<Character>>(label))
    //     {}

    //     Natural tag_;

    //     Integer dim_;

    //     std::shared_ptr<std::basic_string<Character>> label_;

    //     std::shared_ptr<Vector<Real>> coefficients_;
        
    //     inline
    //     Boolean
    //     operator==(
    //         Dof2 const & other
    //     )
    //     const = default;
        
    //     inline
    //     Boolean
    //     operator!=(
    //         Dof2 const & other
    //     )
    //     const = default;
        
    //     inline
    //     Boolean
    //     operator==(
    //         std::shared_ptr<Vector<Real>> const & coefficients
    //     )
    //     const
    //     {
    //         return coefficients_ == coefficients;
    //     }
        
    //     inline
    //     Boolean
    //     operator!=(
    //         std::shared_ptr<Vector<Real>> const & coefficients
    //     )
    //     const
    //     {
    //         return !(* this == coefficients);
    //     }

    //     inline
    //     Natural
    //     getTag()
    //     const
    //     {
    //         return tag_;
    //     }

    //     inline
    //     Vector<Real> const &
    //     getCoefficients()
    //     const
    //     {
    //         return * coefficients_;
    //     }

    //     inline
    //     Vector<Real> &
    //     getCoefficients()
    //     {
    //         return * coefficients_;
    //     }

    // };

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
