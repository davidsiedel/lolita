#ifndef D16FC3C7_A3D0_4425_9E72_D7BB9F9534D8
#define D16FC3C7_A3D0_4425_9E72_D7BB9F9534D8

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_n_0.hxx"

namespace lolita2::geometry
{
    
    template<Element t_element, Quadrature t_quadrature>
    struct ElementQuadratureRuleTraits;
    
    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.isNode() || !t_element.isNode())
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {
        
        lolita::index const static constexpr dim_ = 1;

        std::array<std::array<lolita::real, 3>, dim_> const static constexpr reference_points_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };
        
        std::array<lolita::real, dim_> const static constexpr reference_weights_ = {
                +1.0000000000000000
        };

    };
    
    template<Element t_element, Quadrature t_quadrature>
    struct ElementQuadratureTraits : public ElementQuadratureRuleTraits<t_element, t_quadrature>
    {
        lolita::finite_element::Quadrature const static constexpr quadrature_ = t_quadrature;

        static constexpr
        Quadrature
        getQuadrature()
        {
            return t_quadrature;
        }
        
        static constexpr
        lolita::integer
        getDim()
        {
            return ElementQuadratureRuleTraits<t_element, t_quadrature>::dim_;
        }

    };

    template<Element t_element, Domain t_domain>
    struct FiniteElement;

    template<Element t_element, Domain t_domain>
    struct FiniteElementConnectivity
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element, t_domain>;
        
        template<Element t__element, Domain t__domain>
        using t_ElementPointer = std::shared_ptr<FiniteElement<t__element, t__domain>>;

    public:
    
        using Components = typename t_ElementTraits::template InnerConnectivity<t_ElementPointer>;
        
        using Neighbours = typename t_ElementTraits::template OuterConnectivity<t_ElementPointer>;

        // FiniteElementConnectivity(
        //     lolita::natural tag
        // )
        // :
        // tag_(tag)
        // {}
        
        lolita::natural tag_;
        
        Neighbours neighbours_;
        
        Components components_;
        
        std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;
        
        std::shared_ptr<lolita::domain::Point> coordinates_;
        
        lolita::boolean
        operator==(
                FiniteElementConnectivity const & other
        )
        const = default;
        
        lolita::boolean
        operator!=(
                FiniteElementConnectivity const & other
        )
        const = default;
        
        std::basic_string<lolita::character>
        hash()
        const
        requires(t_element.isNode())
        {
            return std::to_string(this->tag_);
        }
        
        std::basic_string<lolita::character>
        hash()
        const
        {
            std::basic_stringstream<lolita::character> hash;
            auto const & nodes = getComponents<t_element.dim_ - 1, 0>();
            for (auto const & node : nodes) {
                hash << node->hash();
            }
            return hash.str();
        }
        
        template<lolita::index t_i, lolita::index t_j>
        static constexpr
        lolita::index
        getComponentNodeConnection(
                lolita::index i,
                lolita::index j
        )
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element, t_domain>::node_connectivity_))[i][j];
        }
        
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> &
        getNeighbours()
        {
            return std::get<t_j>(std::get<t_i>(neighbours_));
        }
        
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> const &
        getNeighbours()
        const
        {
            return std::get<t_j>(std::get<t_i>(neighbours_));
        }
        
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> &
        getComponents()
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(components_));
        }
        
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> const &
        getComponents()
        const
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(components_));
        }
        
        template<lolita::index t_i, lolita::index t_j>
        lolita::index
        getComponentIndex(
                lolita::index i
        )
        const
        requires(!t_element.isNode())
        {
            auto constexpr t_component = t_ElementTraits::template getComponent<t_i, t_j>();
            using t_NeighbourTraits = ElementTraits<t_component, t_domain>;
            auto constexpr t_coordinates = t_NeighbourTraits::template getNeighbourCoordinates<t_element>();
            auto const & items = getComponents<t_i, t_j>()[i]->template getNeighbours<t_coordinates.dim_, t_coordinates.tag_>();
            auto is_equal = [&] (t_ElementPointer<t_element, t_domain> const & ptr_element) {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }
        
        template<lolita::index t_i, lolita::index t_j>
        lolita::integer
        getComponentOrientation(
                lolita::index i
        )
        const
        requires(!t_element.isNode())
        {
            return getComponentIndex<t_i, t_j>(i) == 0 ? 1 : -1;
        }

    };

}

#endif /* D16FC3C7_A3D0_4425_9E72_D7BB9F9534D8 */
