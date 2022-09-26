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

    struct HybridDiscontinuousGalerkin
    {

        constexpr
        HybridDiscontinuousGalerkin(
            std::basic_string_view<Character> && stabilization,
            Basis cell_basis,
            Basis face_basis,
            Basis grad_basis
        )
        :
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

    namespace detail
    {

        template<typename t_T>
        struct IsHybridDiscontinuousGalerkin : std::false_type {};
        
        template<>
        struct IsHybridDiscontinuousGalerkin<HybridDiscontinuousGalerkin> : std::true_type {};

    }

    template<typename t_T>
    concept HybridDiscontinuousGalerkinConcept = detail::IsHybridDiscontinuousGalerkin<std::decay_t<t_T>>::value;

    template<typename t_T>
    concept DiscretizationConcept = HybridDiscontinuousGalerkinConcept<t_T>;
    
    struct Field
    {
        
        constexpr explicit
        Field(
            Integer dim
        )
        :
        label_(),
        dim_(dim)
        {}
        
        constexpr
        Field(
            std::basic_string_view<Character> && label,
            Integer dim
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        dim_(dim)
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
        std::basic_string_view<Character>
        getLabelView()
        const
        {
            return label_.view();
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

    };
    
    struct Mapping
    {

        constexpr
        Mapping(
            std::basic_string_view<Character> && transformation,
            Field const & field
        )
        :
        transformation_(std::forward<std::basic_string_view<Character>>(transformation)),
        field_(field),
        row_(-1),
        col_(-1)
        {}

        constexpr
        Mapping(
            std::basic_string_view<Character> && transformation,
            Field const & field,
            Integer row,
            Integer col
        )
        :
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

        Label transformation_;

        Field field_;

        Integer row_;

        Integer col_;

    };
    
    template<typename... t_Strains>
    struct Potential
    {

        using Strains = std::array<Mapping, sizeof...(t_Strains)>;
        
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
        // strains_(strains...)
        strains_({strains...})
        {}

        constexpr
        Potential(
            std::basic_string_view<Character> && label,
            t_Strains &&... strains
        )
        :
        label_(std::forward<std::basic_string_view<Character>>(label)),
        // strains_(std::move(strains)...)
        strains_({std::move(strains)...})
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
        utility::Label const &
        getLabel()
        const
        {
            return label_;
        }

        // template<Integer t_i>
        // constexpr
        // Mapping
        // getStrain()
        // const
        // {
        //     return utility::get<t_i>(strains_);
        // }

        // template<Integer t_i>
        constexpr
        Mapping
        getStrain(
            Integer i
        )
        const
        {
            return strains_[i];
        }

        constexpr
        Strains const &
        getStrains()
        const
        {
            return strains_;
        }

        Label label_;

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
    
    template<typename... t_Mappings>
    struct GeneralizedStrain
    {

        using Mappings = lolita::utility::Aggregate<t_Mappings...>;
        
        static constexpr
        Integer
        getNumMappings()
        {
            return sizeof...(t_Mappings);
        }

        constexpr
        GeneralizedStrain(
            Integer tag,
            Field field,
            t_Mappings... mappings
        )
        :
        tag_(tag),
        field_(field),
        mappings_(mappings...)
        {}

        constexpr
        Boolean
        operator==(
            GeneralizedStrain const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            GeneralizedStrain const & other
        )
        const = default;

        constexpr
        Integer
        getTag()
        const
        {
            return tag_;
        }

        constexpr
        Field
        getField()
        const
        {
            return field_;
        }

        template<Integer t_i>
        constexpr
        Mapping
        getMapping()
        const
        {
            return utility::get<t_i>(mappings_);
        }

        Integer tag_;

        Field field_;

        Mappings mappings_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsGeneralizedStrain : std::false_type {};
        
        template<typename... t_T>
        struct IsGeneralizedStrain<GeneralizedStrain<t_T...>> : std::true_type {};

    }

    template<typename t_T>
    concept GeneralizedStrainConcept = detail::IsGeneralizedStrain<t_T>::value;

    template<GeneralizedStrainConcept... t_GeneralizedStrains>
    struct Behavior
    {

        using GeneralizedStrains = lolita::utility::Aggregate<t_GeneralizedStrains...>;

        static constexpr
        Integer
        getNumGeneralizedStrains()
        {
            return sizeof...(t_GeneralizedStrains);
        }

        constexpr
        Behavior(
            Integer tag,
            t_GeneralizedStrains... generalized_strains
        )
        :
        tag_(tag),
        generalized_strains_(generalized_strains...)
        {}

        constexpr
        Boolean
        operator==(
            Behavior const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Behavior const & other
        )
        const = default;

        constexpr
        Integer
        getTag()
        const
        {
            return tag_;
        }

        template<Integer t_i>
        constexpr
        std::tuple_element_t<t_i, std::tuple<t_GeneralizedStrains...>> const &
        getGeneralizedStrain()
        const
        {
            return utility::get<t_i>(generalized_strains_);
        }

        Integer tag_;

        GeneralizedStrains generalized_strains_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsBehavior : std::false_type {};
        
        template<typename... t_T>
        struct IsBehavior<Behavior<t_T...>> : std::true_type {};

    }

    template<typename t_T>
    concept BehaviorConcept = detail::IsBehavior<t_T>::value;

    template<GeneralizedStrainConcept t_GeneralizedStrain, BehaviorConcept t_Behavior, typename t_Discretization>
    struct FiniteElementMethod
    {

        constexpr
        FiniteElementMethod(
            t_GeneralizedStrain generalized_strain,
            t_Behavior behavior,
            t_Discretization discretization,
            Quadrature quadrature
        )
        :
        generalized_strain_(generalized_strain),
        behavior_(behavior),
        discretization_(discretization),
        quadrature_(quadrature)
        {}

        constexpr
        Boolean
        operator==(
            FiniteElementMethod const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            FiniteElementMethod const & other
        )
        const = default;

        constexpr
        Boolean
        isHdg()
        const
        {
            return HybridDiscontinuousGalerkinConcept<t_Discretization>;
        }

        constexpr
        Field
        getField()
        const
        {
            return generalized_strain_.field_;
        }

        template<Integer t_i>
        constexpr
        Mapping
        getMapping()
        const
        {
            return utility::get<t_i>(generalized_strain_.mappings_);
        }

        constexpr
        t_GeneralizedStrain const &
        getGeneralizedStrain()
        const
        {
            return generalized_strain_;
        }

        constexpr
        t_Behavior const &
        getBehavior()
        const
        {
            return behavior_;
        }

        constexpr
        t_Discretization const &
        getDiscretization()
        const
        {
            return discretization_;
        }

        constexpr
        Quadrature const &
        getQuadrature()
        const
        {
            return quadrature_;
        }

        t_GeneralizedStrain generalized_strain_;

        t_Behavior behavior_;

        t_Discretization discretization_;

        Quadrature quadrature_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsFiniteElementMethod : std::false_type {};
        
        template<typename... t_T>
        struct IsFiniteElementMethod<FiniteElementMethod<t_T...>> : std::true_type {};

    }

    struct ElementType
    {

        enum Type
        {

            // Points,
            // Curves,
            // Facets,
            // Solids,
            Cells,
            Faces,
            Edges,
            Nodes,

        };

        static constexpr
        ElementType
        cells(
            Domain domain
        )
        {
            return ElementType(domain, Type::Cells);
        }

        static constexpr
        ElementType
        faces(
            Domain domain
        )
        {
            return ElementType(domain, Type::Faces);
        }

        constexpr
        ElementType(
            Domain domain,
            Type type
        )
        :
        domain_(domain),
        type_(type)
        {}
        
        constexpr
        Integer
        getDim()
        const
        {
            switch (type_)
            {
                case Type::Cells: return domain_.getDim() - 0;
                case Type::Faces: return domain_.getDim() - 1;
                case Type::Edges: return 1;
                case Type::Nodes: return 0;
                default : return -1;
            }
        }

        Domain domain_;

        Type type_;

    };
    
    template<typename t_T>
    concept FiniteElementMethodConcept = detail::IsFiniteElementMethod<std::decay_t<t_T>>::value;

    
    
    // struct FieldDescription
    // {

    //     struct DataBase
    //     {

    //         DataBase()
    //         {}
        
    //         inline
    //         void
    //         addReal(
    //             std::basic_string<Character> && label
    //         )
    //         {
    //             if (std::find(reals_.begin(), reals_.end(), std::forward<std::basic_string<Character>>(label)) == reals_.end())
    //             {
    //                 reals_.push_back(std::forward<std::basic_string<Character>>(label));
    //             }
    //         }
            
    //         inline
    //         Integer
    //         getRealIndex(
    //             std::basic_string<Character> && label
    //         )
    //         const
    //         {
    //             return std::distance(reals_.begin(), std::find(reals_.begin(), reals_.end(), std::forward<std::basic_string<Character>>(label)));
    //         }
            
    //         inline
    //         Integer
    //         getNumReals()
    //         const
    //         {
    //             return reals_.size();
    //         }
            
    //         inline
    //         void
    //         addMatrix(
    //             std::basic_string<Character> && label
    //         )
    //         {
    //             if (std::find(matrices_.begin(), matrices_.end(), std::forward<std::basic_string<Character>>(label)) == matrices_.end())
    //             {
    //                 matrices_.push_back(std::forward<std::basic_string<Character>>(label));
    //             }
    //         }
            
    //         inline
    //         Integer
    //         getMatrixIndex(
    //             std::basic_string<Character> && label
    //         )
    //         const
    //         {
    //             return std::distance(matrices_.begin(), std::find(matrices_.begin(), matrices_.end(), std::forward<std::basic_string<Character>>(label)));
    //         }
            
    //         inline
    //         Integer
    //         getNumMatrices()
    //         const
    //         {
    //             return matrices_.size();
    //         }
            
    //         inline
    //         void
    //         addVector(
    //             std::basic_string<Character> && label
    //         )
    //         {
    //             if (std::find(vectors_.begin(), vectors_.end(), std::forward<std::basic_string<Character>>(label)) == vectors_.end())
    //             {
    //                 vectors_.push_back(std::forward<std::basic_string<Character>>(label));
    //             }
    //         }
            
    //         inline
    //         Integer
    //         getVectorIndex(
    //             std::basic_string<Character> && label
    //         )
    //         const
    //         {
    //             return std::distance(vectors_.begin(), std::find(vectors_.begin(), vectors_.end(), std::forward<std::basic_string<Character>>(label)));
    //         }
            
    //         inline
    //         Integer
    //         getNumVectors()
    //         const
    //         {
    //             return matrices_.size();
    //         }

    //     private:
        
    //         std::vector<std::basic_string<Character>> reals_;
            
    //         std::vector<std::basic_string<Character>> vectors_;
            
    //         std::vector<std::basic_string<Character>> matrices_;

    //     };

    //     struct IntegrationPointData : DataBase
    //     {

    //         inline
    //         Boolean
    //         operator==(
    //             IntegrationPointData const & other
    //         )
    //         const = default;

    //         inline
    //         Boolean
    //         operator!=(
    //             IntegrationPointData const & other
    //         )
    //         const = default;
        
    //     };

    //     struct ElementData : DataBase
    //     {

    //         template<auto t_discretization>
    //         static
    //         std::shared_ptr<ElementData>
    //         make()
    //         {
    //             auto element_data = std::make_shared<ElementData>();
    //             element_data->addReal("StabilizationParamater");
    //             element_data->addMatrix("StabilizationOperator");
    //             element_data->addMatrix("KTT");
    //             element_data->addVector("RT");
    //             return element_data;
    //         }

    //         inline
    //         Boolean
    //         operator==(
    //             ElementData const & other
    //         )
    //         const = default;

    //         inline
    //         Boolean
    //         operator!=(
    //             ElementData const & other
    //         )
    //         const = default;
        
    //     };

    //     struct MeshData : DataBase
    //     {
            
    //         struct ExternalLoad
    //         {
                
    //             ExternalLoad(
    //                 std::function<Real(Point const &, Real const &)> const & function,
    //                 Integer row,
    //                 Integer col
    //             )
    //             :
    //             function_(function),
    //             row_(row),
    //             col_(col)
    //             {}
                
    //             ExternalLoad(
    //                 std::function<Real(Point const &, Real const &)> && function,
    //                 Integer row,
    //                 Integer col
    //             )
    //             :
    //             function_(std::move(function)),
    //             row_(row),
    //             col_(col)
    //             {}

    //             inline
    //             Boolean
    //             operator==(
    //                 ExternalLoad const & other
    //             )
    //             const = default;

    //             inline
    //             Boolean
    //             operator!=(
    //                 ExternalLoad const & other
    //             )
    //             const = default;

    //             inline
    //             Integer
    //             getRow()
    //             const
    //             {
    //                 return row_;
    //             }

    //             inline
    //             Integer
    //             getCol()
    //             const
    //             {
    //                 return row_;
    //             }

    //             inline
    //             Real
    //             getValue(
    //                 Point const & point,
    //                 Real const & time
    //             )
    //             const
    //             {
    //                 return function_(point, time);
    //             }

    //         private:

    //             Integer row_;

    //             Integer col_;

    //             std::function<Real(Point const &, Real const &)> function_;
            
    //         };

    //         inline
    //         Boolean
    //         operator==(
    //             MeshData const & other
    //         )
    //         const = default;

    //         inline
    //         Boolean
    //         operator!=(
    //             MeshData const & other
    //         )
    //         const = default;

    //         std::vector<ExternalLoad> loads_;

    //     };
        
    //     explicit
    //     FieldDescription(
    //         std::basic_string<Character> const & label
    //     )
    //     :
    //     label_(label)
    //     {}
        
    //     explicit
    //     FieldDescription(
    //         std::basic_string<Character> && label
    //     )
    //     :
    //     label_(std::move(label))
    //     {}

    //     inline
    //     Boolean
    //     operator==(
    //         FieldDescription const & other
    //     )
    //     const = default;

    //     inline
    //     Boolean
    //     operator!=(
    //         FieldDescription const & other
    //     )
    //     const = default;
        
    //     inline
    //     std::basic_string<Character> const &
    //     getLabel()
    //     const
    //     {
    //         return label_;
    //     }

    //     IntegrationPointData const &
    //     getIntegrationPointData()
    //     const
    //     {
    //         return integration_point_data_;
    //     }

    //     IntegrationPointData &
    //     getIntegrationPointData()
    //     {
    //         return integration_point_data_;
    //     }

    //     ElementData const &
    //     getElementData()
    //     const
    //     {
    //         return element_data_;
    //     }

    //     ElementData &
    //     getElementData()
    //     {
    //         return element_data_;
    //     }

    //     MeshData const &
    //     getMeshData()
    //     const
    //     {
    //         return mesh_data_;
    //     }

    //     MeshData &
    //     getMeshData()
    //     {
    //         return mesh_data_;
    //     }

    // private:

    //     MeshData mesh_data_;

    //     ElementData element_data_;

    //     IntegrationPointData integration_point_data_;
    
    //     std::basic_string<Character> label_;

    // };

} // namespace lolita

#endif /* A4BCD9B5_985A_4D19_B3E9_7C559F45A353 */
