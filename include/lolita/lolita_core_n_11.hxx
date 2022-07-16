#ifndef CC79BDC5_49DB_4A81_8A93_18ABD6551AF1
#define CC79BDC5_49DB_4A81_8A93_18ABD6551AF1

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_00.hxx"
#include "lolita/lolita_core_n_0.hxx"
#include "lolita/lolita_core_n_111.hxx"

namespace lolita2::geometry
{

    template<template<Element, Domain, auto...> typename t_T, Element t_element, Domain t_domain, auto... t_args>
    struct FiniteElementConnectivity
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element, t_domain>;
        
        template<Element t__element, Domain t__domain>
        using t_ElementPointer = std::shared_ptr<t_T<t__element, t__domain, t_args...>>;

    public:
    
        using Components = typename t_ElementTraits::template InnerConnectivity<t_ElementPointer>;
        
        using Neighbours = typename t_ElementTraits::template OuterConnectivity<t_ElementPointer>;
        
        lolita::natural tag_;
        
        Neighbours neighbours_;
        
        Components components_;

        std::vector<std::shared_ptr<MeshDomain>> domains_;
        
        std::shared_ptr<Point> coordinates_;
        
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
            for (auto const & node : nodes)
            {
                hash << node->hash();
            }
            return hash.str();
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita::integer
        getComponentNodeConnection(
            lolita::integer i,
            lolita::integer j
        )
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element, t_domain>::node_connectivity_))[i][j];
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> &
        getNeighbours()
        {
            return std::get<t_j>(std::get<t_i>(neighbours_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> const &
        getNeighbours()
        const
        {
            return std::get<t_j>(std::get<t_i>(neighbours_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> &
        getComponents()
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(components_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> const &
        getComponents()
        const
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(components_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        lolita::integer
        getComponentIndex(
            lolita::integer i
        )
        const
        requires(!t_element.isNode())
        {
            auto constexpr t_component = t_ElementTraits::template getComponent<t_i, t_j>();
            using t_NeighbourTraits = ElementTraits<t_component, t_domain>;
            auto constexpr t_coordinates = t_NeighbourTraits::template getNeighbourCoordinates<t_element>();
            auto const & items = getComponents<t_i, t_j>()[i]->template getNeighbours<t_coordinates.dim_, t_coordinates.tag_>();
            auto is_equal = [&] (t_ElementPointer<t_element, t_domain> const & ptr_element)
            {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        lolita::integer
        getComponentOrientation(
            lolita::integer i
        )
        const
        requires(!t_element.isNode())
        {
            return getComponentIndex<t_i, t_j>(i) == 0 ? 1 : -1;
        }

    };

    template<Element t_element, Domain t_domain, auto t_arg>
    struct Fefinal : FiniteElementConnectivity<Fefinal, t_element, t_domain, t_arg>
    {

        void
        make()
        {

        }

    };

    template<Element t_element, Domain t_domain, ElementGroupConcept auto... t_element_group>
    struct FiniteElementHolder : FiniteElementConnectivity<FiniteElementHolder, t_element, t_domain, t_element_group...>
    {

        using Elements = std::tuple<std::shared_ptr<Fefinal<t_element, t_domain, t_element_group>>...>;

        Elements elements_;
    
        template<lolita::integer t_i>
        std::tuple_element_t<t_i, Elements> const &
        getElement()
        const
        {
            return std::get<t_i>(elements_);
        }
        
        template<lolita::integer t_i>
        std::tuple_element_t<t_i, Elements> &
        getElement()
        {
            return std::get<t_i>(elements_);
        }

    };

    template<Domain t_domain, auto... t_args>
    struct FiniteElementSet : ElementSet<FiniteElementHolder, t_domain, t_args...>
    {

        FiniteElementSet(
            std::basic_string_view<lolita::character> str
        )
        :
        ElementSet<FiniteElementHolder, t_domain, t_args...>(MeshFileParser<FiniteElementHolder, t_domain, t_args...>::makeIt(str))
        {}
        
        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            FiniteElementSet const & mesh
        )
        {
            auto print_element_components = [&] <Element t_element, lolita::integer t_i = 0, lolita::integer t_j = 0> (
                    auto const & elem_arg,
                    auto & self
            )
            mutable
            {
                if constexpr (!t_element.isNode())
                {
                    auto const constexpr _component = ElementTraits<t_element, t_domain>::template getComponent<t_i, t_j>();
                    for (auto const & c_ : elem_arg->template getComponents<t_i, t_j>())
                    {
                        os << "layer : " << t_i << " type : " << t_j << " <-- " << _component << " " << c_->hash() << std::endl;
                    }
                    if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumComponents<t_i>() - 1)
                    {
                        self.template operator()<t_element, t_i, t_j + 1>(elem_arg, self);
                    }
                    else if constexpr (t_i < ElementTraits<t_element, t_domain>::getNumComponents() - 1)
                    {
                        self.template operator()<t_element, t_i + 1, 0>(elem_arg, self);
                    }
                }
            };
            auto print_element_neighbours = [&] <Element t_element, lolita::integer t_i = 0, lolita::integer t_j = 0> (
                    auto const & elem_arg,
                    auto & self
            )
            mutable
            {
                auto const constexpr _neighbour = ElementTraits<t_element, t_domain>::template getNeighbour<t_i, t_j>();
                for (auto const & c_ : elem_arg->template getNeighbours<t_i, t_j>())
                {
                    if constexpr (!t_element.isNode() && t_i == 0)
                    {
                        os << "layer : " << t_i << " type : " << t_j << " <-> " << _neighbour << " " << c_->hash() << std::endl;
                    }
                    else
                    {
                        os << "layer : " << t_i << " type : " << t_j << " --> " << _neighbour << " " << c_->hash() << std::endl;
                    }
                }
                if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumNeighbours<t_i>() - 1)
                {
                    self.template operator()<t_element, t_i, t_j + 1>(elem_arg, self);
                }
                else if constexpr (t_i < ElementTraits<t_element, t_domain>::getNumNeighbours() - 1)
                {
                    self.template operator()<t_element, t_i + 1, 0>(elem_arg, self);
                }
            };
            auto print_elements2 = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                if constexpr (t_i == 0 && t_j == 0)
                {
                    os << "*** Elements : " << std::endl;
                }
                for (auto const & element : mesh.template getElements<t_i, t_j>())
                {
                    os << "* Element : " << t_element << " " << element.second->hash() << std::endl;
                    os << "domains : ";
                    for (auto const & domain : element.second->domains_)
                    {
                        os << domain->tag_ << " ";
                    }
                    os << std::endl;
                    print_element_components.template operator()<t_element>(element.second, print_element_components);
                    print_element_neighbours.template operator()<t_element>(element.second, print_element_neighbours);
                }                
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            print_elements2(print_elements2);
            return os;
        }

    };

}

#endif /* CC79BDC5_49DB_4A81_8A93_18ABD6551AF1 */
