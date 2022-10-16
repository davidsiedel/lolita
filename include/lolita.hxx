#ifndef A4BCD9B5_985A_4D19_B3E9_7C559F45A353
#define A4BCD9B5_985A_4D19_B3E9_7C559F45A353

#include <MGIS/Behaviour/BehaviourData.hxx>
#include <MGIS/Behaviour/BehaviourData.h>
#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>

#include "utility.hxx"
#include "config.hxx"
#include "numerics.hxx"
#include "algebra.hxx"

namespace lolita
{

    struct Label
    {

    private:

        Integer static constexpr size_ = 100;

        static constexpr
        std::array<Character, size_>
        tag(
            auto &&... str
        )
        {
            auto tag = std::array<Character, size_>();
            auto offset = Integer(0);
            auto make_tag = [&] (
                std::basic_string_view<Character> && str_arg
            )
            constexpr
            {
                for (auto i = 0; i < std::forward<std::basic_string_view<Character>>(str_arg).size(); i++)
                {
                    tag[offset + i] = std::forward<std::basic_string_view<Character>>(str_arg)[i];
                }
                offset += std::forward<std::basic_string_view<Character>>(str_arg).size();
            };
            ((make_tag(std::forward<std::basic_string_view<Character>>(str)), ...));
            return tag;
        }

    public:

        constexpr
        Label(
            auto &&... str
        )
        :
        tag_(tag(std::forward<std::basic_string_view<Character>>(str)...))
        {}

        constexpr
        Boolean
        operator==(
            Label const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Label const & other
        )
        const = default;

        constexpr
        Boolean
        operator==(
            auto && str
        )
        const
        {
            return this->view() == std::forward<decltype(str)>(str);
        }

        constexpr
        Boolean
        operator!=(
            auto && str
        )
        const
        {
            return !(* this == std::forward<decltype(str)>(str));
        }

        constexpr
        Boolean
        is(
            std::basic_string_view<Character> && str
        )
        const
        {
            return * this == std::forward<std::basic_string_view<Character>>(str);
        }

        constexpr
        std::basic_string_view<Character>
        view()
        const
        {
            return std::basic_string_view<Character>(tag_.data(), std::distance(tag_.begin(), std::find(tag_.begin(), tag_.end(), Character())));
        }

        friend inline
        std::ostream &
        operator<<(
            std::ostream & os,
            Label const & label
        )
        {
            os << label.view();
            return os;
        }

        std::array<Character, size_> tag_;

    };
    
    using Point = DenseVector<Real, 3>;

    template<typename t_T>
    concept PointConcept = DenseVectorConcept<t_T, Real, 3>;

    using Loading = std::function<Real(Point const &, Real const &)>;

    struct Strategy
    {

        enum Type
        {
            EigenLU,
            EigenLDLT,
        };

        static constexpr
        Strategy
        eigenLU()
        {
            return Strategy(EigenLU);
        }

        static constexpr
        Strategy
        eigenLU(
            Integer tag
        )
        {
            return Strategy(EigenLU, tag);
        }

    private:

        constexpr
        Strategy(
            Type type
        )
        :
        type_(type),
        tag_(-1)
        {}

        constexpr
        Strategy(
            Type type,
            Integer tag
        )
        :
        type_(type),
        tag_(tag)
        {}

    public:

        constexpr
        Boolean
        operator==(
            Strategy const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Strategy const & other
        )
        const = default;

        Type type_;

        Integer tag_;

    };

    template<Strategy t_s>
    struct LinearSystem
    {

        static
        std::unique_ptr<LinearSystem>
        make_unique()
        {
            return std::make_unique<LinearSystem>();
        }

        LinearSystem()
        :
        rhs_size_(0),
        lhs_size_(0)
        {}

        inline
        Boolean
        operator==(
            LinearSystem const & other
        )
        const = default;

        inline
        Boolean
        operator!=(
            LinearSystem const & other
        )
        const = default;

        std::atomic<Natural> const &
        getSize()
        const
        {
            return rhs_size_;
        }

        std::atomic<Natural> &
        getSize()
        {
            return rhs_size_;
        }

        std::atomic<Natural> const &
        getLhsSize()
        const
        {
            return lhs_size_;
        }

        std::atomic<Natural> const &
        getRhsSize()
        const
        {
            return rhs_size_;
        }

        void
        addRhsSize(
            Natural && size
        )
        {
            rhs_size_ += std::forward<Natural>(size);
        }

        void
        addLhsSize(
            Natural && size
        )
        {
            lhs_size_ += std::forward<Natural>(size);
        }

    private:

        std::atomic<Natural> rhs_size_;

        std::atomic<Natural> lhs_size_;

    };

    struct Domain
    {

        constexpr
        Domain(
            std::basic_string_view<Character> && frame_label,
            Integer dim
        )
        :
        frame_(std::forward<std::basic_string_view<Character>>(frame_label)),
        dim_(dim)
        {}

        constexpr
        Boolean
        operator==(
            Domain const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Domain const & other
        )
        const = default;

        constexpr
        Boolean
        hasDim(
            Integer dim
        )
        const
        {
            return dim_ == dim;
        }
        
        constexpr
        Boolean
        isAxiSymmetric()
        const
        {
            return frame_ == "AxiSymmetric";
        }
        
        constexpr
        Boolean
        isCartesian()
        const
        {
            return frame_ == "Cartesian";
        }

        constexpr
        Integer
        getDim()
        const
        {
            return dim_;
        }

        constexpr
        std::basic_string_view<Character>
        getFrameView()
        const
        {
            return frame_.view();
        }

        Label frame_;

        Integer dim_;

    };
    
    struct Quadrature
    {
        
        constexpr
        Quadrature(
            std::basic_string_view<Character> && rule_label,
            Integer ord
        )
        :
        rule_(std::forward<std::basic_string_view<Character>>(rule_label)),
        ord_(ord)
        {}

        constexpr
        Boolean
        operator==(
            Quadrature const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Quadrature const & other
        )
        const = default;
        
        constexpr
        Label const &
        getRule()
        const
        {
            return rule_;
        }

        constexpr
        Integer
        getOrder()
        const
        {
            return ord_;
        }
        
        constexpr
        Boolean
        isGauss()
        const
        {
            return rule_ == "Gauss";
        }
        
        constexpr
        Boolean
        isGauss(
            Integer ord
        )
        const
        {
            return rule_ == "Gauss" && ord_ == ord;
        }
        
        constexpr
        Boolean
        hasOrd(
            Integer ord
        )
        const
        {
            return ord_ == ord;
        }
        
        Label rule_;
        
        Integer ord_;

    };

    struct Basis
    {
        
        constexpr
        Basis(
            std::basic_string_view<Character> && polynomial,
            Integer ord
        )
        :
        polynomial_(std::forward<std::basic_string_view<Character>>(polynomial)),
        ord_(ord)
        {}

        constexpr
        Boolean
        operator==(
            Basis const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Basis const & other
        )
        const = default;

        constexpr
        Label const &
        getPolynomial()
        const
        {
            return polynomial_;
        }

        constexpr
        std::basic_string_view<Character>
        getPolynomialView()
        const
        {
            return polynomial_.view();
        }
        
        constexpr
        Boolean
        isMonomial()
        const
        {
            return polynomial_ == "Monomial";
        }
        
        constexpr
        Boolean
        isMonomial(
            Integer order
        )
        const
        {
            return polynomial_ == "Monomial" && ord_ == order;
        }

        constexpr
        Integer
        getOrd()
        const
        {
            return ord_;
        }

        constexpr
        Basis
        toOrder(
            Integer order
        )
        const
        {
            return Basis(polynomial_.view(), order);
        }
        
        Label polynomial_;
        
        Integer ord_;

    };

    struct Discretization
    {

        constexpr
        Discretization()
        :
        label_()
        {}

        constexpr explicit
        Discretization(
            std::basic_string_view<Character> && label
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label))
        {}

        constexpr
        Boolean
        operator==(
            Discretization const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Discretization const & other
        )
        const = default;

        constexpr
        Boolean
        isNone()
        const
        {
            return label_ == Label();
        }

        constexpr
        Boolean
        is(
            std::basic_string_view<Character> && label
        )
        const
        {
            return label_ == std::forward<std::basic_string_view<Character>>(label);
        }

        Label label_;

    };
    
    struct HybridDiscontinuousGalerkin : Discretization
    {

        constexpr
        HybridDiscontinuousGalerkin(
            Basis cell_basis,
            Basis face_basis,
            Basis grad_basis
        )
        :
        Discretization("HybridDiscontinuousGalerkin"),
        cell_basis_(cell_basis),
        face_basis_(face_basis),
        grad_basis_(grad_basis)
        {}

        constexpr
        HybridDiscontinuousGalerkin(
            Basis cell_basis,
            Basis face_basis
        )
        :
        Discretization("HybridDiscontinuousGalerkin"),
        cell_basis_(cell_basis),
        face_basis_(face_basis),
        grad_basis_(face_basis)
        {}

        constexpr
        HybridDiscontinuousGalerkin(
            Integer ord_cell_basis,
            Integer ord_face_basis,
            Integer ord_grad_basis
        )
        :
        Discretization("HybridDiscontinuousGalerkin"),
        cell_basis_(Basis("Monomial", ord_cell_basis)),
        face_basis_(Basis("Monomial", ord_face_basis)),
        grad_basis_(Basis("Monomial", ord_grad_basis))
        {}

        constexpr
        HybridDiscontinuousGalerkin(
            Integer ord_cell_basis,
            Integer ord_face_basis
        )
        :
        Discretization("HybridDiscontinuousGalerkin"),
        cell_basis_(Basis("Monomial", ord_cell_basis)),
        face_basis_(Basis("Monomial", ord_face_basis)),
        grad_basis_(Basis("Monomial", ord_face_basis))
        {}

        constexpr
        Boolean
        operator==(
            HybridDiscontinuousGalerkin const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            HybridDiscontinuousGalerkin const & other
        )
        const = default;

        constexpr
        Basis
        getCellBasis()
        const
        {
            return cell_basis_;
        }

        constexpr
        Basis
        getFaceBasis()
        const
        {
            return face_basis_;
        }

        constexpr
        Basis
        getGradBasis()
        const
        {
            return grad_basis_;
        }

        Basis cell_basis_;

        Basis face_basis_;

        Basis grad_basis_;

    };
    
    struct Monomial : Discretization
    {

        constexpr
        Monomial(
            Integer order
        )
        :
        Discretization("Monomial"),
        basis_(Basis("Monomial", order))
        {}

        constexpr
        Boolean
        operator==(
            Monomial const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Monomial const & other
        )
        const = default;

        constexpr
        Basis
        getBasis()
        const
        {
            return basis_;
        }

        Basis basis_;

    };

    template<typename T>
    concept DiscretizationConcept = std::derived_from<std::decay_t<T>, Discretization>;

    struct Field
    {

        constexpr
        Field(
            Integer dim_domain,
            Integer dim_tensor
        )
        :
        label_(),
        dim_domain_(dim_domain),
        dim_tensor_(dim_tensor)
        {}
        
        constexpr
        Field(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        dim_domain_(dim_domain),
        dim_tensor_(dim_tensor)
        {}

        constexpr
        Boolean
        operator==(
            Field const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Field const & other
        )
        const = default;

        constexpr
        Integer
        getDomainDim()
        const
        {
            return dim_domain_;
        }

        constexpr
        Integer
        getTensorDim()
        const
        {
            return dim_tensor_;
        }

        constexpr
        Label const &
        getLabel()
        const
        {
            return label_;
        }

        constexpr
        Boolean
        isTensor(
            Integer dim
        )
        const
        {
            return dim_tensor_ == dim;
        }

        Label label_;

        Integer dim_domain_;

        Integer dim_tensor_;

    };

    template<typename T>
    concept FieldConcept = std::derived_from<std::decay_t<T>, Field>;
    
    struct UnknownField : Field
    {

    private:

        using Base_ = Field;

    public:

        constexpr
        UnknownField(
            Integer dim_domain,
            Integer dim_tensor,
            Basis const & basis
        )
        :
        Base_(dim_domain, dim_tensor),
        basis_(basis)
        {}
        
        constexpr
        UnknownField(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor,
            Basis const & basis
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), dim_domain, dim_tensor),
        basis_(basis)
        {}
        
        constexpr
        UnknownField(
            Field const & field,
            Basis const & basis
        )
        :
        Base_(field),
        basis_(basis)
        {}

        constexpr
        Boolean
        operator==(
            UnknownField const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            UnknownField const & other
        )
        const = default;

        constexpr
        Basis const &
        getBasis()
        const
        {
            return basis_;
        }

        Basis basis_;

    };
    
    struct ImposedField : Field
    {

    private:

        using Base_ = Field;

    public:

        using Base_::Base_;

    };

    template<typename... UnknownFieldT>
    struct LinearOperator : Field
    {

    private:

        using Base_ = Field;

        using UnknownFields_ = utility::Aggregate<UnknownFieldT...>;

    protected:

        constexpr 
        LinearOperator(
            Integer dim_domain,
            Integer dim_tensor,
            UnknownFieldT const &... unknown_fields
        )
        :
        Base_(dim_domain, dim_tensor),
        row_(-1),
        col_(-1),
        unknown_fields_(unknown_fields...)
        {}

        constexpr 
        LinearOperator(
            Integer dim_domain,
            Integer dim_tensor,
            Integer row,
            Integer col,
            UnknownFieldT const &... unknown_fields
        )
        :
        Base_(dim_domain, dim_tensor),
        row_(row),
        col_(col),
        unknown_fields_(unknown_fields...)
        {}

        constexpr 
        LinearOperator(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor,
            UnknownFieldT const &... unknown_fields
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), dim_domain, dim_tensor),
        row_(-1),
        col_(-1),
        unknown_fields_(unknown_fields...)
        {}

        constexpr 
        LinearOperator(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor,
            Integer row,
            Integer col,
            UnknownFieldT const &... unknown_fields
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), dim_domain, dim_tensor),
        row_(row),
        col_(col),
        unknown_fields_(unknown_fields...)
        {}

    public:

        constexpr
        Boolean
        operator==(
            LinearOperator const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            LinearOperator const & other
        )
        const = default;

        constexpr
        Integer
        getRow()
        const
        {
            return row_;
        }

        constexpr
        Integer
        getCol()
        const
        {
            return col_;
        }

        template<Integer i>
        constexpr
        std::tuple_element_t<i, std::tuple<UnknownFieldT...>> const &
        getField()
        const
        {
            return utility::get<i>(unknown_fields_);
        }

        Integer row_;

        Integer col_;

        UnknownFields_ unknown_fields_;

    };

    template<typename CellFieldT>
    struct Gradient : LinearOperator<CellFieldT>
    {

    private:

        using Base_ = LinearOperator<CellFieldT>;

    public:

        explicit constexpr
        Gradient(
            CellFieldT const & cell_field
        )
        :
        Base_(cell_field.getDomainDim(), cell_field.getTensorDim() + 1, cell_field)
        {}

        constexpr 
        Gradient(
            std::basic_string_view<Character> && label,
            CellFieldT const & cell_field
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), cell_field.getDomainDim(), cell_field.getTensorDim() + 1, cell_field)
        {}

        constexpr
        Boolean
        operator==(
            Gradient const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Gradient const & other
        )
        const = default;
        
    };

    template<typename CellFieldT>
    struct Trace : LinearOperator<CellFieldT>
    {

    private:

        using Base_ = LinearOperator<CellFieldT>;

    public:

        explicit constexpr
        Trace(
            CellFieldT const & cell_field
        )
        :
        Base_(cell_field.getDomainDim(), cell_field.getTensorDim(), cell_field)
        {}
        
        constexpr
        Trace(
            std::basic_string_view<Character> && label,
            CellFieldT const & cell_field
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), cell_field.getDomainDim(), cell_field.getTensorDim(), cell_field)
        {}

        constexpr
        Boolean
        operator==(
            Trace const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Trace const & other
        )
        const = default;
        
    };

    template<typename CellFieldT, typename FaceFieldT>
    struct HybridDiscontinuousGalerkinGradient : LinearOperator<CellFieldT, FaceFieldT>
    {

    private:

        using Base_ = LinearOperator<CellFieldT>;

    public:

        constexpr
        HybridDiscontinuousGalerkinGradient(
            CellFieldT const & cell_field,
            CellFieldT const & face_field
        )
        :
        Base_(cell_field.getDomainDim(), cell_field.getTensorDim() + 1, cell_field, face_field)
        {}

        constexpr
        HybridDiscontinuousGalerkinGradient(
            std::basic_string_view<Character> && label,
            CellFieldT const & cell_field,
            CellFieldT const & face_field
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), cell_field.getDomainDim(), cell_field.getTensorDim() + 1, cell_field, face_field)
        {}

        constexpr
        Boolean
        operator==(
            HybridDiscontinuousGalerkinGradient const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            HybridDiscontinuousGalerkinGradient const & other
        )
        const = default;
        
    };

    template<typename CellFieldT, typename FaceFieldT>
    struct HybridDiscontinuousGalerkinStabilization : LinearOperator<CellFieldT, FaceFieldT>
    {

    private:

        using Base_ = LinearOperator<CellFieldT>;

    public:

        constexpr
        HybridDiscontinuousGalerkinStabilization(
            CellFieldT const & cell_field,
            CellFieldT const & face_field
        )
        :
        Base_(face_field.getDomainDim(), cell_field.getTensorDim(), cell_field, face_field)
        {}

        constexpr
        HybridDiscontinuousGalerkinStabilization(
            std::basic_string_view<Character> && label,
            CellFieldT const & cell_field,
            CellFieldT const & face_field
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), face_field.getDomainDim(), cell_field.getTensorDim(), cell_field, face_field)
        {}

        constexpr
        Boolean
        operator==(
            HybridDiscontinuousGalerkinStabilization const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            HybridDiscontinuousGalerkinStabilization const & other
        )
        const = default;
        
    };

    namespace detail
    {

        template<typename T>
        struct IsLinearOperator : std::false_type {};
        
        template<typename... T>
        struct IsLinearOperator<LinearOperator<T...>> : std::true_type {};
        
        template<typename... T>
        struct IsLinearOperator<Gradient<T...>> : std::true_type {};
        
        template<typename... T>
        struct IsLinearOperator<Trace<T...>> : std::true_type {};
        
        template<typename... T>
        struct IsLinearOperator<HybridDiscontinuousGalerkinGradient<T...>> : std::true_type {};
        
        template<typename... T>
        struct IsLinearOperator<HybridDiscontinuousGalerkinStabilization<T...>> : std::true_type {};

    }
    
    template<typename T>
    concept LinearOperatorConcept = detail::IsLinearOperator<std::decay_t<T>>::value;
    
    struct PotentialBase
    {

        constexpr
        PotentialBase(
            Integer dim_domain,
            Quadrature quadrature
        )
        :
        label_(),
        dim_domain_(dim_domain),
        quadrature_(quadrature)
        {}

        constexpr
        PotentialBase(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Quadrature quadrature
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        dim_domain_(dim_domain),
        quadrature_(quadrature)
        {}

        constexpr
        Label const &
        getLabel()
        const
        {
            return label_;
        }

        constexpr
        Integer
        getDomainDim()
        const
        {
            return dim_domain_;
        }

        constexpr
        Quadrature const &
        getQuadrature()
        const
        {
            return quadrature_;
        }

        Label label_;

        Integer dim_domain_;

        Quadrature quadrature_;

    };

    template<typename... LinearOperatorT>
    struct InternalPotential : PotentialBase
    {

    private:

        using Base_ = PotentialBase;

    public:
        
        static constexpr
        Integer
        getNumMappings()
        {
            return sizeof...(LinearOperatorT);
        }

    private:

        using LinearOperators_ = utility::Aggregate<LinearOperatorT...>;

    public:

        constexpr
        InternalPotential(
            Integer dim_domain,
            Quadrature quadrature,
            LinearOperatorT const &... linear_operators
        )
        :
        Base_(dim_domain, quadrature),
        linear_operators_(linear_operators...)
        {}

        constexpr
        InternalPotential(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Quadrature quadrature,
            LinearOperatorT const &... linear_operators
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), dim_domain, quadrature),
        linear_operators_(linear_operators...)
        {}

        constexpr
        Boolean
        operator==(
            InternalPotential const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            InternalPotential const & other
        )
        const = default;

        constexpr
        LinearOperators_ const &
        getLinearOperators()
        const
        {
            return linear_operators_;
        }

        template<Integer i>
        constexpr
        std::tuple_element_t<i, std::tuple<LinearOperatorT...>> const &
        getLinearOperator()
        const
        {
            return utility::get<i>(linear_operators_);
        }

        LinearOperators_ linear_operators_;

    };

    namespace detail
    {

        template<typename T>
        struct IsInternalPotential : std::false_type {};
        
        template<typename... T>
        struct IsInternalPotential<InternalPotential<T...>> : std::true_type {};

    }
    
    template<typename T>
    concept InternalPotentialConcept = detail::IsInternalPotential<std::decay_t<T>>::value;
    
    struct ExternalPotential : PotentialBase
    {

    private:

        using Base_ = PotentialBase;

    public:

        constexpr
        ExternalPotential(
            Integer dim_domain,
            Quadrature quadrature,
            UnknownField const & unknown_field,
            ImposedField const & imposed_field
        )
        :
        Base_(dim_domain, quadrature),
        unknown_field_(unknown_field),
        imposed_field_(imposed_field)
        {}

        constexpr
        ExternalPotential(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Quadrature quadrature,
            UnknownField const & unknown_field,
            ImposedField const & imposed_field
        )
        :
        Base_(std::forward<std::basic_string_view<Character>>(label), dim_domain, quadrature),
        unknown_field_(unknown_field),
        imposed_field_(imposed_field)
        {}

        constexpr
        UnknownField const &
        getUnknownField()
        const
        {
            return unknown_field_;
        }

        constexpr
        ImposedField const &
        getImposedField()
        const
        {
            return imposed_field_;
        }

        UnknownField unknown_field_;

        ImposedField imposed_field_;

    };
    
    template<typename... PotentialT>
    struct InternalEnergy
    {
        
        static constexpr
        Integer
        getNumMappings()
        {
            return sizeof...(PotentialT);
        }

    private:

        using Potentials_ = utility::Aggregate<PotentialT...>;

    public:

        constexpr
        InternalEnergy(
            PotentialT const &... linear_operators
        )
        :
        label_(),
        linear_operators_(linear_operators...)
        {}

        constexpr
        InternalEnergy(
            std::basic_string_view<Character> && label,
            PotentialT const &... linear_operators
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        linear_operators_(linear_operators...)
        {}

        constexpr
        Boolean
        operator==(
            InternalEnergy const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            InternalEnergy const & other
        )
        const = default;

        constexpr
        Potentials_ const &
        getPotentials()
        const
        {
            return linear_operators_;
        }

        template<Integer i>
        constexpr
        std::tuple_element_t<i, std::tuple<PotentialT...>> const &
        getPotential()
        const
        {
            return utility::get<i>(linear_operators_);
        }

        Label label_;

        Potentials_ linear_operators_;

    };

    template<typename... PotentialT>
    struct ExternalEnergy
    {
        
        static constexpr
        Integer
        getNumMappings()
        {
            return sizeof...(PotentialT);
        }

    private:

        using Potentials_ = utility::Aggregate<PotentialT...>;

    public:

        constexpr
        InternalEnergy(
            PotentialT const &... linear_operators
        )
        :
        label_(),
        linear_operators_(linear_operators...)
        {}

        constexpr
        InternalEnergy(
            std::basic_string_view<Character> && label,
            PotentialT const &... linear_operators
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        linear_operators_(linear_operators...)
        {}

        constexpr
        Boolean
        operator==(
            InternalEnergy const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            InternalEnergy const & other
        )
        const = default;

        constexpr
        Potentials_ const &
        getPotentials()
        const
        {
            return linear_operators_;
        }

        template<Integer i>
        constexpr
        std::tuple_element_t<i, std::tuple<PotentialT...>> const &
        getPotential()
        const
        {
            return utility::get<i>(linear_operators_);
        }

        Label label_;

        Potentials_ linear_operators_;
        
    };
    
    template<DiscretizationConcept t_Discretization = Discretization>
    struct DiscreteField : Field
    {

        using FieldDiscretization = t_Discretization;

        constexpr explicit
        DiscreteField(
            Integer dim_domain,
            Integer dim_tensor
        )
        :
        Field(dim_domain, dim_tensor),
        discretization_()
        {}
        
        constexpr
        DiscreteField(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor
        )
        :
        Field(std::forward<std::basic_string_view<Character>>(label), dim_domain, dim_tensor),
        discretization_()
        {}
        
        constexpr
        DiscreteField(
            std::basic_string_view<Character> && label,
            Integer dim_domain,
            Integer dim_tensor,
            FieldDiscretization const & discretization
        )
        :
        Field(std::forward<std::basic_string_view<Character>>(label), dim_domain, dim_tensor),
        discretization_(discretization)
        {}
        
        constexpr explicit
        DiscreteField(
            Field const & field,
            FieldDiscretization const & discretization
        )
        :
        Field(field),
        discretization_(discretization)
        {}

        constexpr
        Boolean
        operator==(
            DiscreteField const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            DiscreteField const & other
        )
        const = default;

        constexpr
        FieldDiscretization const &
        getDiscretization()
        const
        {
            return discretization_;
        }

        FieldDiscretization discretization_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsDiscreteField : std::false_type {};
        
        template<typename t_T>
        struct IsDiscreteField<DiscreteField<t_T>> : std::true_type {};

    }

    template<typename t_T>
    concept DiscreteFieldConcept = detail::IsDiscreteField<std::decay_t<t_T>>::value;
    
    template<DiscreteFieldConcept t_Field>
    struct Mapping
    {

        using Field = t_Field;

        constexpr
        Mapping(
            std::basic_string_view<Character> && transformation,
            Field const & field
        )
        :
        label_(field.getLabel().view(), std::forward<std::basic_string_view<Character>>(transformation)),
        transformation_(std::forward<std::basic_string_view<Character>>(transformation)),
        field_(field),
        row_(-1),
        col_(-1)
        {}

        constexpr
        Mapping(
            std::basic_string_view<Character> && transformation,
            Integer row,
            Integer col,
            Field const & field
        )
        :
        label_(field.getLabel().view(), std::forward<std::basic_string_view<Character>>(transformation)),
        transformation_(std::forward<std::basic_string_view<Character>>(transformation)),
        field_(field),
        row_(row),
        col_(col)
        {}

        constexpr
        Boolean
        operator==(
            Mapping const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Mapping const & other
        )
        const = default;

        constexpr
        Label const &
        getLabel()
        const
        {
            return label_;
        }

        constexpr
        Label const &
        getTransformation()
        const
        {
            return transformation_;
        }

        constexpr
        Integer
        getRow()
        const
        {
            return row_;
        }

        constexpr
        Integer
        getCol()
        const
        {
            return col_;
        }

        constexpr
        Field const &
        getField()
        const
        {
            return field_;
        }

        constexpr
        Boolean
        isGradient()
        const
        {
            return transformation_ == "Gradient";
        }

        constexpr
        Boolean
        isIdentity()
        const
        {
            return transformation_ == "Identity";
        }

        constexpr
        Boolean
        isDivergence()
        const
        {
            return transformation_ == "Divergence";
        }

        constexpr
        Boolean
        isSmallStrain()
        const
        {
            return transformation_ == "SmallStrain";
        }

        constexpr
        Boolean
        isLargeStrain()
        const
        {
            return transformation_ == "LargeStrain";
        }

        Label label_;

        Label transformation_;

        Field field_;

        Integer row_;

        Integer col_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsMapping : std::false_type {};
        
        template<typename t_T>
        struct IsMapping<Mapping<t_T>> : std::true_type {};

    }

    template<typename t_T>
    concept MappingConcept = detail::IsMapping<std::decay_t<t_T>>::value;
    
    template<MappingConcept... t_Strains>
    struct Potential
    {

        using Strains = utility::Aggregate<t_Strains...>;

        using Fields = utility::Aggregate<typename t_Strains::Field...>;
        
        static constexpr
        Integer
        getNumMappings()
        {
            return sizeof...(t_Strains);
        }

        constexpr
        Potential(
            std::basic_string_view<Character> && label,
            Quadrature const & quadrature,
            t_Strains const &... strains
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        quadrature_(quadrature),
        fields_(strains.getField()...),
        strains_(strains...)
        {}

        constexpr
        Boolean
        operator==(
            Potential const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Potential const & other
        )
        const = default;

        constexpr
        Label const &
        getLabel()
        const
        {
            return label_;
        }

        constexpr
        Quadrature const &
        getQuadrature()
        const
        {
            return quadrature_;
        }

        template<Integer t_i>
        constexpr
        utility::aggregate_element_t<t_i, Strains> const &
        getStrain()
        const
        {
            return utility::get<t_i>(strains_);
        }

        constexpr
        Fields const &
        getFields()
        const
        {
            return fields_;
        }

        constexpr
        Strains const &
        getStrains()
        const
        {
            return strains_;
        }

        Label label_;

        Quadrature quadrature_;

        Fields fields_;

        Strains strains_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsPotential : std::false_type {};
        
        template<typename... t_T>
        struct IsPotential<Potential<t_T...>> : std::true_type {};

    }

    template<typename t_T>
    concept PotentialConcept = detail::IsPotential<std::decay_t<t_T>>::value;

    template<PotentialConcept... t_Potentials>
    struct Lagrangian
    {

        using Potentials = utility::Aggregate<t_Potentials...>;
        
        static constexpr
        Integer
        getNumPotentials()
        {
            return sizeof...(t_Potentials);
        }

        constexpr
        Lagrangian(
            std::basic_string_view<Character> && label,
            t_Potentials const &... potentials
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        potentials_(potentials...)
        {}

        constexpr
        Boolean
        operator==(
            Lagrangian const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Lagrangian const & other
        )
        const = default;

        constexpr
        Label const &
        getLabel()
        const
        {
            return label_;
        }

        template<Integer t_i>
        constexpr
        utility::aggregate_element_t<t_i, Potentials> const &
        getPotential()
        const
        {
            return utility::get<t_i>(potentials_);
        }

        constexpr
        Potentials const &
        getPotentials()
        const
        {
            return potentials_;
        }

        Label label_;

        Potentials potentials_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsLagrangian : std::false_type {};
        
        template<typename... t_T>
        struct IsLagrangian<Lagrangian<t_T...>> : std::true_type {};

    }

    template<typename t_T>
    concept LagrangianConcept = detail::IsLagrangian<std::decay_t<t_T>>::value;

} // namespace lolita

#endif /* A4BCD9B5_985A_4D19_B3E9_7C559F45A353 */
