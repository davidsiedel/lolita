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

        // std::basic_string_view<lolita::character> hash_;
        
        Neighbours neighbours_;
        
        Components components_;
        
        // std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;

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

        static
        std::basic_string<lolita::character>
        getElementHash(
            std::array<lolita::integer, t_element.num_nodes_> node_tags
        )
        requires(!t_element.isNode())
        {
            auto element_hash = std::basic_stringstream<lolita::character>();
            std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
            for (auto node_tag : node_tags)
            {
                element_hash << node_tag;
            }
            return element_hash.str();
        }

        static
        std::basic_string<lolita::character>
        getElementHash(
            lolita::natural node_tag
        )
        requires(t_element.isNode())
        {
            return std::to_string(node_tag);
        }
        
        template<lolita::integer t_i = 0, lolita::integer t_j = 0>
        static
        void
        makeElement(
            auto & list,
            std::shared_ptr<FiniteElementHolder> & ptr_element,
            std::array<lolita::integer, t_element.num_nodes_> const & node_tags
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = ElementTraits<t_element, t_domain>::template getComponent<t_i, t_j>();
            auto const constexpr _is_initialized = t_i == 0 && t_j == 0;
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto const constexpr t_node_coordinates = ElementTraits<t_element, t_domain>::template getComponentCoordinates<Element::node()>();
            auto const constexpr _component_coordinates = DomainTraits<t_domain>::template getElementCoordinates<_component>();
            auto const constexpr _neighbour_coordinates = ElementTraits<_component, t_domain>::template getNeighbourCoordinates<t_element>();
            auto & components = list.template getElements<_component_coordinates.dim_, _component_coordinates.tag_>();
            auto & element_component_array = ptr_element->template getComponents<t_i, t_j>();
            for (auto i = 0; i < element_component_array.size(); ++i)
            {
                auto component_hash = std::basic_string<lolita::character>();
                if constexpr(!_component.isNode())
                {
                    auto component_node_tags = std::array<lolita::integer, _component.num_nodes_>();
                    for (lolita::integer j = 0; j < _component.num_nodes_; ++j)
                    {
                        auto const k = FiniteElementHolder::template getComponentNodeConnection<t_i, t_j>(i, j);
                        component_node_tags[j] = node_tags[k];
                    }
                    component_hash = lolita2::geometry::FiniteElementHolder<_component, t_domain, t_element_group...>::getElementHash(component_node_tags);
                    if (!components.contains(component_hash))
                    {
                        using t_Component = lolita2::geometry::FiniteElementHolder<_component, t_domain, t_element_group...>;
                        auto ptr_component = std::make_shared<t_Component>(t_Component());
                        lolita2::geometry::FiniteElementHolder<_component, t_domain, t_element_group...>::template makeElement<0, 0>(ptr_component, component_node_tags);
                    }
                }
                else
                {
                    component_hash = std::to_string(node_tags[FiniteElementHolder::template getComponentNodeConnection<t_i, t_j>(i, 0)]);
                }
                element_component_array[i] = components[component_hash];
                components[component_hash]->template getNeighbours<_neighbour_coordinates.dim_, _neighbour_coordinates.tag_>().push_back(ptr_element);
            }
            if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumComponents<t_i>() - 1)
            {
                makeElement<t_i, t_j + 1u>(ptr_element, node_tags);
            }
            else if constexpr (t_i < ElementTraits<t_element, t_domain>::getNumComponents() - 1)
            {
                makeElement<t_i + 1u, 0u>(ptr_element, node_tags);
            }
            if constexpr (_is_initialized)
            {
                auto const & nodes = ptr_element->template getComponents<t_node_coordinates.dim_, t_node_coordinates.tag_>();
                auto domains = std::set<std::shared_ptr<std::basic_string<lolita::character>>>();
                for (auto const & domain : nodes[0]->domains_)
                {
                    auto has_domain = true;
                    for (lolita::integer j = 1; j < t_element.num_nodes_; ++j)
                    {
                        has_domain = std::find(nodes[j]->domains_.begin(), nodes[j]->domains_.end(), domain) != nodes[j]->domains_.end();
                    }
                    if (has_domain)
                    {
                        domains.insert(domain);
                    }
                }
                ptr_element->tag_ = list.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>().size();
                ptr_element->domains_.assign(domains.begin(), domains.end());
                list.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>()[getElementHash(node_tags)] = ptr_element;
            }
        }

        // static
        // void
        // makeElement(
        //     auto & list,
        //     std::shared_ptr<FiniteElementHolder> & ptr_element,
        //     lolita::natural const & tag,
        //     Point const & coordinates,
        //     std::vector<std::shared_ptr<std::basic_string<lolita::character>>> const & domains
        // )
        // requires(t_element.isNode())
        // {
        //     auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
        //     auto hklm = std::forward<std::shared_ptr<FiniteElementHolder>>(ptr_element);
        //     hklm->tag_ = std::forward<lolita::natural>(tag);
        //     hklm->domains_ = std::forward<std::vector<std::shared_ptr<std::basic_string<lolita::character>>>>(domains);
        //     hklm->coordinates_ = std::make_shared<Point>(std::forward<Point>(coordinates));
        //     list.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>()[std::to_string(tag)] = hklm;
        // }

        template<lolita::integer t_i = 1>
        static
        void
        hello(
            auto & list,
            std::shared_ptr<FiniteElementHolder> ptr_element,
            lolita::natural const & tag,
            Point const & coordinates,
            std::vector<std::shared_ptr<std::basic_string<lolita::character>>> const & domains
        )
        {
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            list.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
            std::cout << "yo !" << std::endl;
            // auto hklm = std::forward<std::shared_ptr<FiniteElementHolder>>(ptr_element);
            if constexpr (t_i > 0)
            {
                hello<t_i - 1>(list, ptr_element, tag, coordinates, domains);
            }
            // auto tagg = std::forward<lolita::natural>(tag);
            // auto coordinatess = std::forward<Point>(coordinates);
            // auto domainss = std::forward<std::vector<std::shared_ptr<std::basic_string<lolita::character>>>>(domains);
        }

    };

    template<template<Element, Domain, auto...> typename t_FiniteElement, Domain t_domain, auto... t_args>
    struct MeshSetBase
    {

    private:

        template<template<Element, Domain, auto...> typename t__FiniteElement, auto... t__args>
        struct FiniteElementCollectionTraits
        {

            template<Element t_element, Domain t__domain>
            using t_FiniteElements = std::map<std::basic_string<lolita::character>, std::shared_ptr<t__FiniteElement<t_element, t__domain, t__args...>>>;

            using t_FiniteElementCollection = ElementCollection<t_FiniteElements, t_domain>;

        };

        template<template<Element, Domain, auto...> typename t__FiniteElement, auto... t__args>
        using t_FiniteElements = typename FiniteElementCollectionTraits<t__FiniteElement, t__args...>::t_FiniteElementCollection;

    public:

        using FiniteElements = typename FiniteElementCollectionTraits<t_FiniteElement, t_args...>::t_FiniteElementCollection;

        using Domains = std::map<std::basic_string<lolita::character>, std::shared_ptr<std::basic_string<lolita::character>>>;
        
        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            MeshSetBase const & mesh
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
                for (auto const & element : mesh.elements_.template getElements<t_i, t_j>())
                {
                    os << "* Element : " << t_element << " " << element.second->hash() << std::endl;
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

        FiniteElements elements_;

        Domains domains_;

    };

    template<Domain t_domain, auto... t_args>
    using ICILA = MeshFileParser<FiniteElementHolder, t_domain, t_args...>;

}

#endif /* CC79BDC5_49DB_4A81_8A93_18ABD6551AF1 */
