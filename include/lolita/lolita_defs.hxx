#ifndef AD31B868_1D56_48DB_A2CC_F898A668702D
#define AD31B868_1D56_48DB_A2CC_F898A668702D

#include <MGIS/Behaviour/BehaviourData.hxx>
#include <MGIS/Behaviour/BehaviourData.h>
#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"

namespace lolita2
{

    struct MeshFileFormat
    {

        enum Format
        {

            Gmsh,

        };

        constexpr
        MeshFileFormat(
            Format format
        )
        :
        format_(format)
        {}

        constexpr
        lolita::boolean
        operator==(
            MeshFileFormat const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            MeshFileFormat const & other
        )
        const = default;

        constexpr
        lolita::boolean
        isGmsh()
        const
        {
            return format_ == Format::Gmsh;
        }

        Format format_;

    };
    
    using Point = lolita::matrix::Vector<lolita::real, 3>;

    using Loading = std::function<lolita::real(lolita2::Point const &, lolita::real const &)>;

    struct Load
    {

        auto static constexpr zero = [] (auto const &, auto const &) constexpr { return lolita::real(0); };

        Load()
        :
        loading_(std::make_shared<Loading>(zero))
        {}

        Load(
            Loading const & loading
        )
        :
        loading_(std::make_shared<Loading>(loading))
        {}

        Load(
            Loading && loading
        )
        :
        loading_(std::make_shared<Loading>(std::forward<Loading>(loading)))
        {}
        
        lolita::real
        getImposedValue(
            lolita2::Point const & point,
            lolita::real const & time
        )
        const
        {
            return loading_->operator ()(point, time);
        }

        std::shared_ptr<Loading> loading_;

    };

    struct Label
    {

        using Tag = std::array<lolita::character, 50>;

    private:

        static constexpr
        Label::Tag
        setTag(
            std::basic_string_view<lolita::character> str
        )
        {
            auto tag = Label::Tag();
            auto count = lolita::integer(0);
            for (auto c : str) {
                tag[count] = c;
                count ++;
            }
            return tag;
        }

    public:

        constexpr
        Label(
            std::basic_string_view<lolita::character> str
        )
        :
        tag_(setTag(str))
        {}

        constexpr
        lolita::boolean
        operator==(
            Label const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Label const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator==(
            std::basic_string_view<lolita::character> str
        )
        const
        {
            return str == this->view();
        }

        constexpr
        lolita::boolean
        operator!=(
            std::basic_string_view<lolita::character> str
        )
        const
        {
            return !(* this == str);
        }

        constexpr
        std::basic_string_view<lolita::character>
        view()
        const
        {
            return std::basic_string_view<lolita::character>(tag_.data(), std::distance(tag_.begin(), std::find(tag_.begin(), tag_.end(), lolita::character())));
        }

        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            Label const & label
        )
        {
            os << label.view();
            return os;
        }

        Label::Tag tag_;

    };

    struct Domain
    {

        enum Frame
        {

            Cartesian,
            AxiSymmetric

        };

        constexpr
        Domain(
            lolita::integer dim,
            Frame frame
        )
        :
        dim_(dim),
        frame_(frame)
        {}

        constexpr
        lolita::boolean
        operator==(
            Domain const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Domain const & other
        )
        const = default;

        constexpr
        Domain
        getSubDomain()
        const
        {
            return Domain(dim_ - 1, frame_);
        }

        constexpr
        lolita::boolean
        hasDim(
            lolita::integer dim
        )
        const
        {
            return dim_ == dim;
        }
        
        constexpr
        lolita::boolean
        isAxiSymmetric()
        const
        {
            return frame_ == Frame::AxiSymmetric;
        }
        
        constexpr
        lolita::boolean
        isCartesian()
        const
        {
            return frame_ == Frame::Cartesian;
        }

        lolita::integer dim_;

        Frame frame_;

    };

    struct Quadrature
    {
        
        enum Rule
        {

            Gauss,

        };
        
        constexpr
        Quadrature(
            Quadrature::Rule rule,
            lolita::integer ord
        )
        :
        rule_(rule),
        ord_(ord)
        {}

        constexpr
        lolita::boolean
        operator==(
            Quadrature const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Quadrature const & other
        )
        const = default;

        static constexpr
        Quadrature
        gauss(
            lolita::integer ord
        )
        {
            return Quadrature(Rule::Gauss, ord);
        }
        
        constexpr
        lolita::boolean
        isGauss()
        const
        {
            return rule_ == Rule::Gauss;
        }
        
        Quadrature::Rule rule_;
        
        lolita::integer ord_;

    };

    struct Basis
    {

        enum Polynomial
        {

            Lagrange,
            Monomial,

        };
        
        constexpr
        Basis(
            Polynomial polynomial,
            lolita::integer ord
        )
        :
        polynomial_(polynomial),
        ord_(ord)
        {}

        constexpr
        lolita::boolean
        operator==(
            Basis const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
            Basis const & other
        )
        const = default;

        static constexpr
        Basis
        monomial(
            lolita::integer ord
        )
        {
            return Basis(Polynomial::Monomial, ord);
        }
        
        constexpr
        lolita::boolean
        isMonomial()
        const
        {
            return polynomial_ == Polynomial::Monomial;
        }
        
        Polynomial polynomial_;
        
        lolita::integer ord_;

    };

    struct Mapping
    {

        enum Type
        {

            Gradient,
            Identity,
            Divergence,
            SmallStrain,
            LargeStrain,

        };

        constexpr
        Mapping(
            Type type
        )
        :
        type_(type)
        {}

        constexpr
        lolita::boolean
        isGradient()
        const
        {
            return type_ == Type::Gradient;
        }

        constexpr
        lolita::boolean
        isIdentity()
        const
        {
            return type_ == Type::Identity;
        }

        constexpr
        lolita::boolean
        isDivergence()
        const
        {
            return type_ == Type::Divergence;
        }

        constexpr
        lolita::boolean
        isSmallStrain()
        const
        {
            return type_ == Type::SmallStrain;
        }

        constexpr
        lolita::boolean
        isLargeStrain()
        const
        {
            return type_ == Type::LargeStrain;
        }

        Type type_;

    };

    struct Field
    {
        
        constexpr
        Field(
            lolita::integer dim
        )
        :
        dim_(dim)
        {}

        constexpr
        lolita::boolean
        isTensor(
            lolita::integer dim
        )
        const
        {
            return dim_ == dim;
        }

        lolita::integer dim_;

    };

    template<typename... t_Mappings>
    struct Unknown
    {

        constexpr
        Unknown(
            Field field,
            t_Mappings... mappings
        )
        :
        field_(field),
        mappings_({mappings...})
        {}

        Field field_;

        std::array<Mapping, sizeof...(t_Mappings)> mappings_;

    };

    template<typename UnknownT>
    struct HybridDiscontinuousGalerkin
    {

        constexpr
        HybridDiscontinuousGalerkin(
            UnknownT unknown,
            Basis cell_basis,
            Basis face_basis
        )
        :
        unknown_(unknown),
        cell_basis_(cell_basis),
        face_basis_(face_basis)
        {}

        UnknownT unknown_;

        Basis cell_basis_;

        Basis face_basis_;

    };

    namespace detail
    {

        template<typename t_T>
        struct IsHybridDiscontinuousGalerkin : std::true_type {};

        template<typename t_T>
        struct IsHybridDiscontinuousGalerkin<HybridDiscontinuousGalerkin<t_T>> : std::false_type {};

    }

    template<typename t_T>
    concept HybridDiscontinuousGalerkinConcept = detail::IsHybridDiscontinuousGalerkin<t_T>::value;
    
    // template<auto... t_finite_elements>
    // struct ElementGroup
    // {

    // private:
    
    //     using t_FiniteElements = std::tuple<std::remove_cvref_t<decltype(t_finite_elements)>...>;

    // public:
    
    //     template<template<auto, auto, auto> typename t_T, auto t_element, auto t_domain>
    //     using ElementPointers = std::tuple<std::shared_ptr<t_T<t_element, t_domain, t_finite_elements>>...>;
        
    //     template<template<auto, auto, auto> typename t_T, auto t_element, auto t_domain>
    //     using Elements = std::tuple<t_T<t_element, t_domain, t_finite_elements>...>;
    
    // };

    // namespace detail
    // {

    //     template<typename t_T>
    //     struct IsElementGroup : std::true_type {};

    //     template<auto... t_args>
    //     struct IsElementGroup<ElementGroup<t_args...>> : std::false_type {};

    // }

    // template<typename t_T>
    // concept ElementGroupConcept = detail::IsElementGroup<t_T>::value;



}


#endif /* AD31B868_1D56_48DB_A2CC_F898A668702D */
