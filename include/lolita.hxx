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

    using Label = utility::Label;
    
    using Point = Vector<Real, 3>;

    template<typename t_T>
    concept PointConcept = VectorConcept<t_T, Real, 3>;

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
        size_(0)
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
            return size_;
        }

        std::atomic<Natural> &
        getSize()
        {
            return size_;
        }

    private:

        std::atomic<Natural> size_;

    };

    struct Domain
    {

        // std::array<Label, 2> static constexpr frame_labels_ = {
        //     "Cartesian",
        //     "AxiSymmetric",
        // };

        // static constexpr
        // Label const &
        // frame(
        //     std::basic_string_view<Character> && frame_label
        // )
        // {
        //     for (auto const & label : frame_labels_)
        //     {
        //         if (label == std::forward<std::basic_string_view<Character>>(frame_label))
        //         {
        //             return label;
        //         }
        //     }
        //     throw std::logic_error("NO !!");
        // }

        // static constexpr
        // std::basic_string_view<Character> const &
        // check(
        //     std::basic_string_view<Character> const & input
        // )
        // {
        //     if constexpr (!utility::isAnyOf(input, "Cartesian", "AxiSymmetric"))
        //     {
        //         throw std::logic_error("NO !!!");
        //     }
        //     return input;
        // }

        constexpr
        Domain(
            std::basic_string_view<Character> && frame_label,
            Integer dim
        )
        :
        // frame_(frame(std::forward<std::basic_string_view<Character>>(frame_label))),
        frame_(std::forward<std::basic_string_view<Character>>(frame_label)),
        dim_(dim)
        {
            // if constexpr (!utility::isAnyOf(std::basic_string_view<Character>(frame_label), "Cartesian", "AxiSymmetric"))
            // {
            //     // throw std::logic_error("NO !!!");
            //     int a = 2;
            // }
        }

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

        // std::array<Label, 1> static constexpr rule_labels_ = {
        //     "Gauss",
        // };

        // static constexpr
        // Label const &
        // rule(
        //     std::basic_string_view<Character> && rule_label
        // )
        // {
        //     for (auto const & label : rule_labels_)
        //     {
        //         if (label == std::forward<std::basic_string_view<Character>>(rule_label))
        //         {
        //             return label;
        //         }
        //     }
        //     throw std::logic_error("NO !!");
        // }
        
        constexpr
        Quadrature(
            std::basic_string_view<Character> && rule_label,
            Integer ord
        )
        :
        // rule_(rule(std::forward<std::basic_string_view<Character>>(rule_label))),
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

        Label label_;

    };
    
    struct HybridDiscontinuousGalerkin : Discretization
    {

        constexpr
        HybridDiscontinuousGalerkin(
            std::basic_string_view<Character> && stabilization,
            Basis cell_basis,
            Basis face_basis,
            Basis grad_basis
        )
        :
        Discretization("HybridDiscontinuousGalerkin"),
        stabilization_(std::forward<std::basic_string_view<Character>>(stabilization)),
        cell_basis_(cell_basis),
        face_basis_(face_basis),
        grad_basis_(grad_basis)
        {}

        constexpr
        HybridDiscontinuousGalerkin(
            std::basic_string_view<Character> && stabilization,
            Basis cell_basis,
            Basis face_basis
        )
        :
        Discretization("HybridDiscontinuousGalerkin"),
        stabilization_(std::forward<std::basic_string_view<Character>>(stabilization)),
        cell_basis_(cell_basis),
        face_basis_(face_basis),
        grad_basis_(face_basis)
        {}

        constexpr
        HybridDiscontinuousGalerkin(
            std::basic_string_view<Character> && stabilization,
            Integer ord_cell_basis,
            Integer ord_face_basis,
            Integer ord_grad_basis
        )
        :
        Discretization("HybridDiscontinuousGalerkin"),
        stabilization_(std::forward<std::basic_string_view<Character>>(stabilization)),
        cell_basis_(Basis("Monomial", ord_cell_basis)),
        face_basis_(Basis("Monomial", ord_face_basis)),
        grad_basis_(Basis("Monomial", ord_grad_basis))
        {}

        constexpr
        HybridDiscontinuousGalerkin(
            std::basic_string_view<Character> && stabilization,
            Integer ord_cell_basis,
            Integer ord_face_basis
        )
        :
        Discretization("HybridDiscontinuousGalerkin"),
        stabilization_(std::forward<std::basic_string_view<Character>>(stabilization)),
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

        constexpr
        Boolean
        isHybridDiscontinuousGalerkin()
        const
        {
            return true;
        }

        constexpr
        Boolean
        isHdg()
        const
        {
            return stabilization_ == "HybridDiscontinuousGalerkin";
        }

        constexpr
        Boolean
        isHho()
        const
        {
            return stabilization_ == "HybridHighOrder";
        }

        Label stabilization_;

        Basis cell_basis_;

        Basis face_basis_;

        Basis grad_basis_;

    };

    // namespace detail
    // {

    //     template<typename t_T>
    //     struct IsHybridDiscontinuousGalerkin : std::false_type {};
        
    //     template<>
    //     struct IsHybridDiscontinuousGalerkin<HybridDiscontinuousGalerkin> : std::true_type {};

    // }

    // template<typename t_T>
    // concept HybridDiscontinuousGalerkinConcept = detail::IsHybridDiscontinuousGalerkin<std::decay_t<t_T>>::value;

    template<typename t_T>
    concept DiscretizationConcept = std::derived_from<t_T, Discretization>;

    template<typename t_T>
    // concept FieldDiscretizationConcept = DiscretizationConcept<t_T> || std::convertible_to<t_T, Basis> || std::convertible_to<t_T, Integer>;
    concept FieldDiscretizationConcept = std::convertible_to<t_T, Discretization> || std::convertible_to<t_T, Basis> || std::convertible_to<t_T, Integer>;
    
    template<FieldDiscretizationConcept t_Discretization = Integer>
    struct Field
    {

        using FieldDiscretization = t_Discretization;
        
        constexpr explicit
        Field(
            Integer dim
        )
        :
        label_(),
        dim_(dim),
        discretization_()
        {}
        
        constexpr
        Field(
            std::basic_string_view<Character> && label,
            Integer dim
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        dim_(dim),
        discretization_()
        {}
        
        constexpr
        Field(
            std::basic_string_view<Character> && label,
            Integer dim,
            FieldDiscretization const & discretization
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        dim_(dim),
        discretization_(discretization)
        {}
        
        constexpr
        Field(
            std::basic_string_view<Character> && label,
            Integer dim,
            FieldDiscretization && discretization
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        dim_(dim),
        discretization_(std::move(discretization))
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
        getDim()
        const
        {
            return dim_;
        }

        constexpr
        Label const &
        getLabel()
        const
        {
            return label_;
        }

        constexpr
        FieldDiscretization const &
        getDiscretization()
        const
        {
            return discretization_;
        }

        constexpr
        Boolean
        hasDiscretization()
        const
        {
            return DiscretizationConcept<FieldDiscretization>;
        }

        constexpr
        Boolean
        hasBasis()
        const
        {
            return std::convertible_to<FieldDiscretization, Basis>;
        }

        constexpr
        Boolean
        hasSize()
        const
        {
            return std::convertible_to<FieldDiscretization, Integer>;
        }

        constexpr
        Boolean
        isTensor(
            Integer dim
        )
        const
        {
            return dim_ == dim;
        }

        Label label_;

        Integer dim_;

        FieldDiscretization discretization_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsField : std::false_type {};
        
        template<typename t_T>
        struct IsField<Field<t_T>> : std::true_type {};

    }

    template<typename t_T>
    concept FieldConcept = detail::IsField<std::decay_t<t_T>>::value;
    
    template<FieldConcept t_Field>
    struct Mapping
    {

        using MappingField = t_Field;

        using MappingFieldDiscretization = typename t_Field::FieldDiscretization;

        constexpr
        Mapping(
            std::basic_string_view<Character> && transformation,
            t_Field const & field
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
            t_Field const & field,
            Integer row,
            Integer col
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
        t_Field const &
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

        t_Field field_;

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

        using PotentialStrains = utility::Aggregate<t_Strains...>;

        using PotentialFields = utility::tuple_unique_t<std::tuple<typename t_Strains::MappingField...>>;

        using PotentialFieldDiscretizations = utility::tuple_unique_t<std::tuple<typename t_Strains::MappingFieldDiscretization...>>;
        
        static constexpr
        Integer
        getNumMappings()
        {
            return sizeof...(t_Strains);
        }

        constexpr
        Potential(
            std::basic_string_view<Character> && label,
            t_Strains const &... strains
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        strains_(strains...)
        {}

        constexpr
        Potential(
            std::basic_string_view<Character> && label,
            t_Strains &&... strains
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        strains_(std::move(strains)...)
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

        template<Integer t_i>
        constexpr
        utility::aggregate_element_t<t_i, PotentialStrains> const &
        getStrain()
        const
        {
            return utility::get<t_i>(strains_);
        }

        constexpr
        PotentialStrains const &
        getStrains()
        const
        {
            return strains_;
        }

        Label label_;

        PotentialStrains strains_;

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

    struct ExternalLoad
    {

        ExternalLoad(
            Integer row,
            Integer col,
            std::function<Real(Point const &, Real const &)> const & function
        )
        :
        row_(row),
        col_(col),
        function_(function)
        {}

        ExternalLoad(
            Integer row,
            Integer col,
            std::function<Real(Point const &, Real const &)> && function
        )
        :
        row_(row),
        col_(col),
        function_(std::move(function))
        {}

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

        inline
        Real
        getValue(
            Point const & point,
            Real const & time
        )
        const
        {
            return function_(point, time);
        }

    private:

        Integer row_;

        Integer col_;

        std::function<Real(Point const &, Real const &)> function_;

    };

} // namespace lolita

#endif /* A4BCD9B5_985A_4D19_B3E9_7C559F45A353 */
