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

    struct MeshData
    {

        enum Format
        {

            Gmsh,

        };

        constexpr
        MeshData(
            Format format
        )
        :
        format_(format)
        {}

        constexpr
        lolita::boolean
        operator==(
                MeshData const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                MeshData const & other
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

        auto const static constexpr zero = [] (auto const &, auto const &) constexpr { return lolita::real(0); };

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
                lolita2::Point const &
                point,
                lolita::real const &
                time
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
            auto count = lolita::index(0);
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
                std::basic_string_view<lolita::character> const & str
        )
        const
        {
            return str == this->view();
        }

        constexpr
        lolita::boolean
        operator!=(
                std::basic_string_view<lolita::character> const & str
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
    
    struct Unknown
    {
        
        constexpr
        Unknown(
            std::basic_string_view<lolita::character> label
        )
        :
        label_(label)
        {}

        constexpr
        lolita::boolean
        operator==(
                Unknown const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                Unknown const & other
        )
        const = default;

        Label label_;

    };

    struct Fieldd
    {
        
        enum Frame
        {

            Cell,
            Face,
            Edge,
            Node,

        };
        
        constexpr
        Fieldd(
            std::basic_string_view<lolita::character> label,
            Frame frame,
            lolita::integer ord
        )
        :
        label_(label),
        frame_(frame),
        ord_(ord)
        {}

        constexpr
        lolita::boolean
        operator==(
                Fieldd const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                Fieldd const & other
        )
        const = default;
        
        constexpr
        lolita::boolean
        isFaceField()
        const
        {
            return frame_ == Frame::Face;
        }
        
        constexpr
        lolita::boolean
        isCellField()
        const
        {
            return frame_ == Frame::Cell;
        }
        
        constexpr
        lolita::boolean
        isTensor(
                lolita::integer ord
        )
        const
        {
            return ord_ == ord;
        }

        Label label_;

        Frame frame_;
        
        lolita::integer ord_;

    };
    
    namespace detail
    {

        template<typename t_T>
        struct IsUnknown : std::false_type {};

        template<>
        struct IsUnknown<Unknown> : std::true_type {};

    }

    template<typename t_T>
    concept UnknownConcept = detail::IsUnknown<std::remove_cvref_t<t_T>>::value;

    struct ElementaryUnknown
    {
        
        constexpr
        ElementaryUnknown(
            Basis const & basis,
            Fieldd const & field
        )
        :
        basis_(basis),
        field_(field)
        {}

        constexpr
        lolita::boolean
        operator==(
                ElementaryUnknown const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                ElementaryUnknown const & other
        )
        const = default;

        Basis basis_;

        Fieldd field_;

    };

    template<ElementaryUnknown... t_elementary_unknowns>
    struct ElementaryField
    {



    };

    struct FieldG
    {
        
        enum Basis
        {

            Monomial,
            Lagrange,

        };
        
        enum Frame
        {

            Cell,
            Face,
            Edge,
            Node,

        };
        
        constexpr
        FieldG(
            lolita::integer ord_field,
            lolita::integer ord_basis,
            FieldG::Basis basis,
            FieldG::Frame frame
        )
        :
        ord_field_(ord_field),
        ord_basis_(ord_basis),
        basis_(basis),
        frame_(frame)
        {}

        constexpr
        lolita::boolean
        operator==(
                FieldG const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                FieldG const & other
        )
        const = default;
        
        constexpr
        lolita::boolean
        isMonomial()
        const
        {
            return basis_ == FieldG::Basis::Monomial;
        }
        
        constexpr
        lolita::boolean
        isFaceField()
        const
        {
            return frame_ == FieldG::Frame::Face;
        }
        
        constexpr
        lolita::boolean
        isCellField()
        const
        {
            return frame_ == FieldG::Frame::Cell;
        }
        
        constexpr
        lolita::boolean
        isTensor(
                lolita::integer ord
        )
        const
        {
            return ord_field_ == ord;
        }
        
        lolita::integer ord_field_;
        
        lolita::integer ord_basis_;
        
        FieldG::Basis basis_;
        
        FieldG::Frame frame_;

    };

    template<auto...>
    struct Field;

    namespace detail
    {

        template<typename t_T>
        struct IsField : std::false_type {};

        template<auto... t_args>
        struct IsField<Field<t_args...>> : std::true_type {};

    }

    template<typename t_T>
    concept FieldConcept = detail::IsField<std::remove_cvref_t<t_T>>::value;

    template<UnknownConcept auto t_unknown>
    struct Field<t_unknown> : FieldG
    {

        Unknown static constexpr unknown_ = t_unknown;
        
        constexpr
        Field(
            std::basic_string_view<lolita::character> label,
            lolita::integer ord_field,
            lolita::integer ord_basis,
            FieldG::Basis basis,
            FieldG::Frame frame
        )
        :
        label_(label),
        FieldG(ord_field, ord_basis, basis, frame)
        {}
        
        constexpr
        Field(
            std::basic_string_view<lolita::character> label,
            FieldG const & field
        )
        :
        label_(label),
        FieldG(field)
        {}

        constexpr
        lolita::boolean
        operator==(
                Field const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                Field const & other
        )
        const = default;

        Label label_;

    };

    template<FieldConcept auto... t_fields>
    struct Field<t_fields...> : FieldG
    {
        
        constexpr
        Field(
            std::basic_string_view<lolita::character> label,
            lolita::integer ord_field,
            lolita::integer ord_basis,
            FieldG::Basis basis,
            FieldG::Frame frame
        )
        :
        label_(label),
        FieldG(ord_field, ord_basis, basis, frame)
        {}
        
        constexpr
        Field(
            std::basic_string_view<lolita::character> label,
            FieldG const & field
        )
        :
        label_(label),
        FieldG(field)
        {}

        constexpr
        lolita::boolean
        operator==(
                Field const & other
        )
        const = default;

        constexpr
        lolita::boolean
        operator!=(
                Field const & other
        )
        const = default;

        Label label_;

    };

    

}


#endif /* AD31B868_1D56_48DB_A2CC_F898A668702D */
