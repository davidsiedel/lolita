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
    
    using Point = Vector<Real, 3>;

    template<typename t_T>
    concept PointConcept = VectorConcept<t_T, Real, 3>;

    using Loading = std::function<Real(Point const &, Real const &)>;

    struct Domain
    {

        enum Frame
        {

            Cartesian,
            AxiSymmetric

        };

        constexpr
        Domain(
            Integer dim,
            Frame frame
        )
        :
        dim_(dim),
        frame_(frame)
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

        // constexpr
        // Domain
        // getSubDomain()
        // const
        // {
        //     return Domain(dim_ - 1, frame_);
        // }

        constexpr
        Boolean
        hasDim(
            Integer dim
        )
        const
        {
            return dim_ == dim;
        }

        // static constexpr
        // Domain
        // axiSymmetric(
        //     Integer dim
        // )
        // {
        //     return Domain(dim, Frame::AxiSymmetric);
        // }
        
        constexpr
        Boolean
        isAxiSymmetric()
        const
        {
            return frame_ == Frame::AxiSymmetric;
        }

        // static constexpr
        // Domain
        // cartesian(
        //     Integer dim
        // )
        // {
        //     return Domain(dim, Frame::Cartesian);
        // }
        
        constexpr
        Boolean
        isCartesian()
        const
        {
            return frame_ == Frame::Cartesian;
        }

        constexpr
        Integer
        getDim()
        const
        {
            return dim_;
        }

        // constexpr
        // Frame
        // getFrame()
        // const
        // {
        //     return frame_;
        // }

        Integer dim_;

        Frame frame_;

    };

    static constexpr
    Domain
    domain(
        std::basic_string_view<Character> label,
        Integer dim
    )
    {
        if (label == "Cartesian")
        {
            return Domain(dim, Domain::Cartesian);
        }
        else if (label == "AxiSymmetric")
        {
            return Domain(dim, Domain::AxiSymmetric);
        }
        else
        {
            throw std::logic_error("NO !!!");
        }
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

    struct Quadrature
    {
        
        enum Rule
        {

            Gauss,

        };
        
        constexpr
        Quadrature(
            Integer ord,
            Rule rule
        )
        :
        ord_(ord),
        rule_(rule)
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

        static constexpr
        Quadrature
        gauss(
            Integer ord
        )
        {
            return Quadrature(ord, Gauss);
        }
        
        constexpr
        Boolean
        isGauss()
        const
        {
            return rule_ == Gauss;
        }
        
        constexpr
        Boolean
        isGauss(
            Integer ord
        )
        const
        {
            return rule_ == Gauss && ord_ == ord;
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
        
        Integer ord_;
        
        Rule rule_;

    };

    static constexpr
    Quadrature
    quadrature(
        std::basic_string_view<Character> label,
        Integer order
    )
    {
        if (label == "Gauss")
        {
            return Quadrature::gauss(order);
        }
        else
        {
            throw std::logic_error("NO !!!");
        }
    }

    struct Basis
    {

        enum Polynomial
        {

            Lagrange,
            Monomial,

        };

        static constexpr
        Basis
        monomial(
            Integer order
        )
        {
            return Basis(order, Basis::Monomial);
        }
        
        constexpr
        Basis(
            Integer ord,
            Polynomial polynomial
        )
        :
        ord_(ord),
        polynomial_(polynomial)
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
        Polynomial
        getPolynomial()
        const
        {
            return polynomial_;
        }
        
        constexpr
        Boolean
        isMonomial()
        const
        {
            return polynomial_ == Monomial;
        }
        
        constexpr
        Boolean
        isMonomial(
            Integer order
        )
        const
        {
            return polynomial_ == Monomial && ord_ == order;
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
            return Basis(order, polynomial_);
        }
        
        Integer ord_;
        
        Polynomial polynomial_;

    };

    struct Field
    {

        static constexpr
        Field
        scalar()
        {
            return Field(0);
        }

        static constexpr
        Field
        scalar(
            std::basic_string_view<Character> && label
        )
        {
            return Field(std::forward<std::basic_string_view<Character>>(label), 0);
        }

        static constexpr
        Field
        vector()
        {
            return Field(1);
        }

        static constexpr
        Field
        vector(
            std::basic_string_view<Character> && label
        )
        {
            return Field(std::forward<std::basic_string_view<Character>>(label), 1);
        }

        static constexpr
        Field
        tensor(
            Integer dim
        )
        {
            return Field(dim);
        }

        static constexpr
        Field
        tensor(
            std::basic_string_view<Character> && label,
            Integer dim
        )
        {
            return Field(std::forward<std::basic_string_view<Character>>(label), dim);
        }
        
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
        utility::Label const &
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

        utility::Label label_;

        Integer dim_;

    };

    static constexpr
    Field
    field(
        std::basic_string_view<Character> && label,
        std::basic_string_view<Character> && label2
    )
    {
        if (std::forward<std::basic_string_view<Character>>(label) == "Scalar")
        {
            return Field(std::forward<std::basic_string_view<Character>>(label2), 0);
        }
        else if (std::forward<std::basic_string_view<Character>>(label) == "Vector")
        {
            return Field(std::forward<std::basic_string_view<Character>>(label2), 1);
        }
        else if (std::forward<std::basic_string_view<Character>>(label) == "Tensor")
        {
            return Field(std::forward<std::basic_string_view<Character>>(label2), 2);
        }
        else
        {
            throw std::logic_error("NO !!!");
        }
    }
    
    struct Mapping
    {

        enum Transformation
        {

            Gradient,
            Identity,
            Divergence,
            SmallStrain,
            LargeStrain,

        };

        static constexpr
        Mapping
        gradient(
            Field const & field
        )
        {
            return Mapping(field, Gradient);
        }

        static constexpr
        Mapping
        identity(
            Field const & field
        )
        {
            return Mapping(field, Identity);
        }

        static constexpr
        Mapping
        divergence(
            Field const & field
        )
        {
            return Mapping(field, Divergence);
        }

        static constexpr
        Mapping
        smallStrain(
            Field const & field
        )
        {
            return Mapping(field, SmallStrain);
        }

        static constexpr
        Mapping
        largeStrain(
            Field const & field
        )
        {
            return Mapping(field, LargeStrain);
        }

        constexpr
        Mapping(
            Field const & field,
            Transformation transformation
        )
        :
        field_(field),
        transformation_(transformation)
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
            return transformation_ == Gradient;
        }

        constexpr
        Boolean
        isIdentity()
        const
        {
            return transformation_ == Identity;
        }

        constexpr
        Boolean
        isDivergence()
        const
        {
            return transformation_ == Divergence;
        }

        constexpr
        Boolean
        isSmallStrain()
        const
        {
            return transformation_ == SmallStrain;
        }

        constexpr
        Boolean
        isLargeStrain()
        const
        {
            return transformation_ == LargeStrain;
        }

        Field field_;

        Transformation transformation_;

    };
    
    static constexpr
    Mapping
    mapping(
        std::basic_string_view<Character> label,
        Field const & field
    )
    {
        if (label == "Gradient")
        {
            return Mapping(field, Mapping::Gradient);
        }
        else if (label == "Identity")
        {
            return Mapping(field, Mapping::Identity);
        }
        else if (label == "Divergence")
        {
            return Mapping(field, Mapping::Divergence);
        }
        else if (label == "SmallStrain")
        {
            return Mapping(field, Mapping::SmallStrain);
        }
        else if (label == "LargeStrain")
        {
            return Mapping(field, Mapping::LargeStrain);
        }
        else
        {
            throw std::logic_error("NO");
        }
    }
    
    template<typename... t_Strains>
    struct Potential
    {

        using Strains = utility::Aggregate<t_Strains...>;
        
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
        utility::Label const &
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

        template<Integer t_i>
        constexpr
        Mapping
        getStrain()
        const
        {
            return utility::get<t_i>(strains_);
        }
        
        // constexpr
        // Integer
        // getNumFields()
        // const
        // {
        //     auto arr = std::array<Field, getNumMappings()>();
                        
        //     auto gnum = Integer(0);
        //     auto set = [&] <Integer t_i = 0> (
        //         auto & self
        //     )
        //     constexpr mutable
        //     {
        //         auto num = Integer(0);
        //         auto hlr = Boolean(true);
        //         auto const & field_i = getStrain<t_i>().getField();
        //         auto set2 = [&] <Integer t_j = t_i> (
        //             auto & self2
        //         )
        //         constexpr mutable
        //         {
        //             auto const & field_j = getStrain<t_j>().getField();
        //             if (hlr)
        //             {
        //                 if (field_j != field_i)
        //                 {
        //                     num += 1;
        //                     hlr = false;
        //                 }
        //             }
        //             if constexpr (t_j < getNumMappings() - 1)
        //             {
        //                 self2.template operator()<t_j + 1>(self2);
        //             }
        //         };
        //         if constexpr (t_i < getNumMappings() - 1)
        //         {
        //             self.template operator()<t_i + 1>(self);
        //         }
        //     };
        //     set(set);
        // }

        utility::Label label_;

        Strains strains_;

    };

    // struct TestConst
    // {

    //     constexpr
    //     TestConst(
    //         std::basic_string_view<Character> && label
    //     )
    //     :
    //     tag_(std::forward<std::basic_string_view<Character>>(label))
    //     {}

    //     utility::Label tag_;

    // };

    // template<TestConst value>
    // struct TestStruct
    // {};

    namespace detail
    {

        template<typename t_T>
        struct IsPotential : std::false_type {};
        
        template<typename... t_T>
        struct IsPotential<Potential<t_T...>> : std::true_type {};

    }

    template<typename t_T>
    concept PotentialConcept = detail::IsPotential<std::decay_t<t_T>>::value;

    template<typename... t_Strains>
    static constexpr
    Potential<t_Strains...>
    potential(
        std::basic_string_view<Character> && label,
        t_Strains const &... strains
    )
    {
        return Potential<t_Strains...>(std::forward<std::basic_string_view<Character>>(label), strains...);
    }

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

        Integer
        getRow()
        const
        {
            return row_;
        }

        Integer
        getCol()
        const
        {
            return col_;
        }

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

    // /**
    //  * @brief Linked to a mesh set, and a field. can be one for faces, one for cells, one for nodes, etc.
    //  * Not necessarily linked to a Dof, can be in a placeholder where no dofs are present.
    //  * But the counterpart in the element contains the Dofs, which is done by the discretization part.
    //  * 
    //  */
    // struct GeneralData
    // {

    //     explicit
    //     GeneralData(
    //         Field const & field
    //     )
    //     :
    //     tag_(field.getTag())
    //     {}

    //     Character
    //     getTag()
    //     const
    //     {
    //         return tag_;
    //     }

    //     inline
    //     void
    //     addScalarField(
    //         std::basic_string<Character> && label
    //     )
    //     {
    //         for (auto const & scalar_field : scalar_fields_)
    //         {
    //             if (scalar_field == std::forward<std::basic_string<Character>>(label))
    //             {
    //                 return;
    //             }
    //         }
    //         scalar_fields_.push_back(std::forward<std::basic_string<Character>>(label));
    //         // auto itr_field = std::find(scalar_fields_.begin(), scalar_fields_.end(), std::forward<std::basic_string<Character>>(label));
    //         // if (itr_field == scalar_fields_.end())
    //         // {
    //         //     scalar_fields_.push_back(std::forward<std::basic_string<Character>>(label));
    //         // }
    //     }

    //     inline
    //     void
    //     addVectorField(
    //         std::basic_string<Character> && label
    //     )
    //     {
    //         for (auto const & vector_field : vector_fields_)
    //         {
    //             if (vector_field == std::forward<std::basic_string<Character>>(label))
    //             {
    //                 return;
    //             }
    //         }
    //         vector_fields_.push_back(std::forward<std::basic_string<Character>>(label));
    //         // auto itr_field = std::find(vector_fields_.begin(), vector_fields_.end(), std::forward<std::basic_string<Character>>(label));
    //         // if (itr_field == vector_fields_.end())
    //         // {
    //         //     vector_fields_.push_back(std::forward<std::basic_string<Character>>(label));
    //         // }
    //     }

    //     inline
    //     void
    //     addMatrixField(
    //         std::basic_string<Character> && label
    //     )
    //     {
    //         for (auto const & matrix_field : matrix_fields_)
    //         {
    //             if (matrix_field == std::forward<std::basic_string<Character>>(label))
    //             {
    //                 return;
    //             }
    //         }
    //         matrix_fields_.push_back(std::forward<std::basic_string<Character>>(label));
    //         // auto itr_field = std::find(matrix_fields_.begin(), matrix_fields_.end(), std::forward<std::basic_string<Character>>(label));
    //         // if (itr_field == matrix_fields_.end())
    //         // {
    //         //     matrix_fields_.push_back(std::forward<std::basic_string<Character>>(label));
    //         // }
    //     }

    //     inline
    //     Integer
    //     getScalarTag(
    //         std::basic_string<Character> && label
    //     )
    //     const
    //     {
    //         auto tag = Integer(0);
    //         for (auto const & scalar_field : scalar_fields_)
    //         {
    //             if (scalar_field == std::forward<std::basic_string<Character>>(label))
    //             {
    //                 return tag;
    //             }
    //             tag ++;
    //         }
    //         throw std::runtime_error("NO!");
    //         // auto itr_field = std::find(scalar_fields_.begin(), scalar_fields_.end(), std::forward<std::basic_string<Character>>(label));
    //         // if (itr_field != scalar_fields_.end())
    //         // {
    //         //     return std::distance(scalar_fields_.begin(), itr_field);
    //         // }
    //         // else
    //         // {
    //         //     throw std::runtime_error("NO!");
    //         // }
    //     }

    //     inline
    //     Integer
    //     getVectorTag(
    //         std::basic_string<Character> && label
    //     )
    //     const
    //     {
    //         auto tag = Integer(0);
    //         for (auto const & vector_field : vector_fields_)
    //         {
    //             if (vector_field == std::forward<std::basic_string<Character>>(label))
    //             {
    //                 return tag;
    //             }
    //             tag ++;
    //         }
    //         throw std::runtime_error("NO!");
    //         // auto itr_field = std::find(vector_fields_.begin(), vector_fields_.end(), std::forward<std::basic_string<Character>>(label));
    //         // if (itr_field != vector_fields_.end())
    //         // {
    //         //     return std::distance(vector_fields_.begin(), itr_field);
    //         // }
    //         // else
    //         // {
    //         //     throw std::runtime_error("NO!");
    //         // }
    //     }

    //     inline
    //     Integer
    //     getMatrixTag(
    //         std::basic_string<Character> && label
    //     )
    //     const
    //     {
    //         auto tag = Integer(0);
    //         for (auto const & matrix_field : matrix_fields_)
    //         {
    //             if (matrix_field == std::forward<std::basic_string<Character>>(label))
    //             {
    //                 return tag;
    //             }
    //             tag ++;
    //         }
    //         throw std::runtime_error("NO!");
    //         // auto itr_field = std::find(matrix_fields_.begin(), matrix_fields_.end(), std::forward<std::basic_string<Character>>(label));
    //         // if (itr_field != matrix_fields_.end())
    //         // {
    //         //     return std::distance(matrix_fields_.begin(), itr_field);
    //         // }
    //         // else
    //         // {
    //         //     throw std::runtime_error("NO!");
    //         // }
    //     }

    //     inline
    //     Integer
    //     getNumScalarFields()
    //     const
    //     {
    //         return scalar_fields_.size();
    //     }

    //     inline
    //     Integer
    //     getNumVectorFields()
    //     const
    //     {
    //         return vector_fields_.size();
    //     }

    //     inline
    //     Integer
    //     getNumMatrixFields()
    //     const
    //     {
    //         return matrix_fields_.size();
    //     }

    //     inline
    //     void
    //     addLoadField(
    //         Integer row,
    //         Integer col,
    //         std::function<Real(Point const &, Real const &)> && function
    //     )
    //     {
    //         for (auto & load : loads_)
    //         {
    //             if (load.getRow() == row && load.getCol() == col)
    //             {
    //                 load = ExternalLoad(row, col, std::forward<std::function<Real(Point const &, Real const &)>>(function));
    //                 return;
    //             }
    //         }
    //         loads_.push_back(ExternalLoad(row, col, std::forward<std::function<Real(Point const &, Real const &)>>(function)));
    //         // auto find_load = [&] (
    //         //     ExternalLoad const & load
    //         // )
    //         // {
    //         //     return load.getRow() == row && load.getCol() == col;
    //         // };
    //         // auto itr_load = std::find_if(loads_.begin(), loads_.end(), find_load);
    //         // if (itr_load == loads_.end())
    //         // {
    //         //     loads_.push_back(ExternalLoad(row, col, std::forward<std::function<Real(Point const &, Real const &)>>(function)));
    //         // }
    //         // else
    //         // {
    //         //     * itr_load = ExternalLoad(row, col, std::forward<std::function<Real(Point const &, Real const &)>>(function));
    //         // }
    //     }

    // private:

    //     Character tag_;

    //     std::vector<std::basic_string<Character>> scalar_fields_;

    //     std::vector<std::basic_string<Character>> vector_fields_;

    //     std::vector<std::basic_string<Character>> matrix_fields_;

    //     std::vector<ExternalLoad> loads_;

    // };

    // static inline
    // std::shared_ptr<GeneralData>
    // data(
    //     Field const & field
    // )
    // {
    //     return std::make_shared<GeneralData>(field);
    // }


    // struct PotentialData
    // {

    //     explicit
    //     PotentialData(
    //         auto const & pot
    //     )
    //     :
    //     tag_(pot.getTag())
    //     {}

    //     inline
    //     void
    //     addBehavior(
    //         auto const &... args
    //     )
    //     {
    //         bhv_ = mgis::behaviour::load(args...);
    //     }

    // private:

    //     Character tag_;

    //     mgis::behaviour::Behaviour bhv_;

    // };

    // struct Dummy
    // {
    //     Dummy(GeneralData const & data) : data_(data) {}
    //     GeneralData const & data_;
    // };
    
    struct HybridDiscontinuousGalerkin
    {

        enum Stabilization
        {

            Hdg,
            HybridHighOrder,

        };

        static constexpr
        HybridDiscontinuousGalerkin
        hybridDiscontinuousGalerkin(
            Integer ord_cell,
            Integer ord_face
        )
        {
            return HybridDiscontinuousGalerkin(Basis::monomial(ord_cell), Basis::monomial(ord_face), Basis::monomial(ord_face), Hdg);
        }

        static constexpr
        HybridDiscontinuousGalerkin
        hybridHighOrder(
            Integer ord_cell,
            Integer ord_face
        )
        {
            return HybridDiscontinuousGalerkin(Basis::monomial(ord_cell), Basis::monomial(ord_face), Basis::monomial(ord_face), HybridHighOrder);
        }

        constexpr
        HybridDiscontinuousGalerkin(
            Basis cell_basis,
            Basis face_basis,
            Basis grad_basis,
            Stabilization stabilization
        )
        :
        cell_basis_(cell_basis),
        face_basis_(face_basis),
        grad_basis_(grad_basis),
        stabilization_(stabilization)
        {}

        constexpr
        HybridDiscontinuousGalerkin(
            Basis cell_basis,
            Basis face_basis,
            Stabilization stabilization
        )
        :
        cell_basis_(cell_basis),
        face_basis_(face_basis),
        grad_basis_(face_basis),
        stabilization_(stabilization)
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
            return stabilization_ == Stabilization::Hdg;
        }

        constexpr
        Boolean
        isHho()
        const
        {
            return stabilization_ == Stabilization::HybridHighOrder;
        }

        Basis cell_basis_;

        Basis face_basis_;

        Basis grad_basis_;

        Stabilization stabilization_;

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

    template<typename t_T>
    concept FiniteElementMethodConcept = detail::IsFiniteElementMethod<std::decay_t<t_T>>::value;

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
