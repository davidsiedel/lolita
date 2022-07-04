#ifndef C31C92CF_B3D3_4D98_9746_7FD9EEFC72FC
#define C31C92CF_B3D3_4D98_9746_7FD9EEFC72FC

#include <execution>
#include "lolita/lolita_core_n_1.hxx"

namespace lolita2::geometry
{
    
    template<Domain _domain>
    struct Mesh;

    template<Element t_element, Domain t_domain>
    struct FiniteElement : FiniteElementConnectivity<t_element, t_domain>
    {

        // FiniteElement()
        // {}

        FiniteElement()
        requires(t_element.isNode())
        {}

        FiniteElement()
        requires(!t_element.isNode())
        {}

        template<Element ___element = t_element>
        static
        void
        makeElement(
                std::array<lolita::index, ___element.num_nodes_> node_tags_,
                Mesh<t_domain> & mesh_data
        )
        requires(!___element.isNode())
        {
            auto get_element_hash = [&] <Element __element> (
                    std::array<lolita::index, __element.num_nodes_> node_tags
            )
            {
                auto element_hash = std::basic_stringstream<lolita::character>();
                if (node_tags.size() > 0) {
                    std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
                }
                for (int i = 0; i < node_tags.size(); ++i) {
                    element_hash << node_tags[i];
                }
                return element_hash.str();
            };
            auto init_elem = [&] <Element __element> (
                    std::shared_ptr<FiniteElement<__element, t_domain>> & ptr_element,
                    std::array<lolita::index, __element.num_nodes_> node_tags
            )
            {
                using __ElementDescription = ElementTraits<__element, t_domain>;
                using __MeshDescription = DomainTraits<t_domain>;
                auto const constexpr _element_coordinates = __MeshDescription::template getElementCoordinates<__element>();
                auto const constexpr _node_coordinates = __ElementDescription::template getComponentCoordinates<Element::node()>();
                auto & elements = mesh_data.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>();
                auto const & nodes = ptr_element->template getComponents<_node_coordinates.dim_, _node_coordinates.tag_>();
                auto domains = std::unordered_set<std::shared_ptr<std::basic_string<lolita::character>>>();
                for (auto const & domain : nodes[0]->domains_) {
                    auto has_domain = true;
                    for (int j = 1; j < __element.num_nodes_; ++j) {
                        has_domain = std::find(nodes[j]->domains_.begin(), nodes[j]->domains_.end(), domain) != nodes[j]->domains_.end();
                    }
                    if (has_domain) {
                        domains.insert(domain);
                    }
                }
                ptr_element->tag_ = elements.size();
                ptr_element->domains_.assign(domains.begin(), domains.end());
                auto element_hash = get_element_hash.template operator ()<__element>(node_tags);
                elements.insert({element_hash, ptr_element});
            };
            auto make_element = [&] <Element __element = t_element, lolita::index _i = 0, lolita::index _j = 0> (
                    std::shared_ptr<FiniteElement<__element, t_domain>> & ptr_element,
                    std::array<lolita::index, __element.num_nodes_> node_tags,
                    auto & make_element_imp
            )
            mutable
            {
                if constexpr (!__element.isNode()) {
                    using __ElementDescription = ElementTraits<__element, t_domain>;
                    using __ComponentDescription = ElementTraits<__ElementDescription::template getComponent<_i, _j>(), t_domain>;
                    using __MeshDescription = DomainTraits<t_domain>;
                    auto const constexpr _is_initialized = _i == 0 && _j == 0;
                    auto const constexpr _component = __ElementDescription::template getComponent<_i, _j>();
                    auto const constexpr _component_coordinates = __MeshDescription::template getElementCoordinates<_component>();
                    auto const constexpr _neighbour_coordinates = __ComponentDescription::template getNeighbourCoordinates<__element>();
                    using __Component = FiniteElement<_component, t_domain>;
                    using __Self = FiniteElement<__element, t_domain>;
                    auto & components = mesh_data.elements_.template getElements<_component_coordinates.dim_, _component_coordinates.tag_>();
                    auto & element_component_array = ptr_element->template getComponents<_i, _j>();
                    for (auto i = 0; i < element_component_array.size(); ++i) {
                        auto component_hash = std::basic_string<lolita::character>();
                        if constexpr(!_component.isNode()) {
                            auto component_node_tags = std::array<lolita::index, _component.num_nodes_>();
                            for (int j = 0; j < _component.num_nodes_; ++j) {
                                auto const k = __Self::template getComponentNodeConnection<_i, _j>(i, j);
                                component_node_tags[j] = node_tags[k];
                            }
                            component_hash = get_element_hash.template operator ()<_component>(component_node_tags);
                            if (!components.contains(component_hash)) {
                                auto ptr_component = std::make_shared<__Component>(__Component());
                                make_element_imp.template operator ()<_component, 0, 0>(ptr_component, component_node_tags, make_element_imp);
                            }
                        }
                        else {
                            component_hash = std::to_string(node_tags[__Self::template getComponentNodeConnection<_i, _j>(i, 0)]);
                        }
                        element_component_array[i] = components[component_hash];
                        components[component_hash]->template getNeighbours<_neighbour_coordinates.dim_, _neighbour_coordinates.tag_>().push_back(ptr_element);
                    }
                    if constexpr (_j < __ElementDescription::template getNumComponents<_i>() - 1) {
                        make_element_imp.template operator()<__element, _i, _j + 1u>(ptr_element, node_tags, make_element_imp);
                    }
                    else if constexpr (_i < __ElementDescription::getNumComponents() - 1) {
                        make_element_imp.template operator()<__element, _i + 1u, 0u>(ptr_element, node_tags, make_element_imp);
                    }
                    if constexpr (_is_initialized) {
                        init_elem.template operator ()<__element>(ptr_element, node_tags);
                    }
                }
            };
            auto ptr_element = std::make_shared<FiniteElement>(FiniteElement());
            make_element(ptr_element, node_tags_, make_element);
        }

        template<Element ___element = t_element>
        static
        void
        makeElement(
                lolita::natural const & tag_,
                lolita::domain::Point const & coordinates_,
                std::vector<std::shared_ptr<std::basic_string<lolita::character>>> const & domains_,
                // ElementInitializationData<t_element> const & initialization_data,
                Mesh<t_domain> & mesh_data
        )
        requires(___element.isNode())
        {
            
            auto make_element = [&] (
                    std::shared_ptr<FiniteElement<t_element, t_domain>> & ptr_element
            )
            mutable
            {
                auto const constexpr _element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
                // ptr_element->tag_ = initialization_data.tag_;
                // ptr_element->domains_ = initialization_data.domains_;
                // ptr_element->coordinates_ = std::make_shared<lolita::domain::Point>(initialization_data.coordinates_);
                // auto elem_hash = std::to_string(initialization_data.tag_);
                // mesh_data.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>().insert({elem_hash, ptr_element});
                ptr_element->tag_ = tag_;
                ptr_element->domains_ = domains_;
                ptr_element->coordinates_ = std::make_shared<lolita::domain::Point>(coordinates_);
                auto elem_hash = std::to_string(tag_);
                mesh_data.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>().insert({elem_hash, ptr_element});
            };
            auto ptr_element = std::make_shared<FiniteElement>(FiniteElement());
            make_element(ptr_element);
        }

        template<Element __element = t_element>
        void
        initialize(
                Mesh<t_domain> & mesh_data
        );

        template<Element __element = t_element>
        void
        initialize(
                Mesh<t_domain> & mesh_data
        )
        requires(!__element.isNode())
        {
            /*
             *
             */
            auto initialize_element = [&] <lolita::index _k = 0u> (
                    auto & initialize_element_imp
            )
                    mutable
            {
                /*
                 *
                 */
                auto initialize_components = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_components_imp
                )
                        mutable
                {
                    for (int i = 0; i < ElementTraits<t_element, t_domain>::template getNumComponents<__i, __j>(); ++i) {
                        auto & rhs = this->template getComponents<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getComponents<__i, __j>()[i];
                        lhs = rhs;
                    }
                    if constexpr (__j < ElementTraits<t_element, t_domain>::template getNumComponents<__i>() - 1) {
                        initialize_components_imp.template operator()<__i, __j + 1u>(initialize_components_imp);
                    }
                    else if constexpr (__i < ElementTraits<t_element, t_domain>::getNumComponents() - 1) {
                        initialize_components_imp.template operator()<__i + 1u, 0u>(initialize_components_imp);
                    }
                };
                /*
                 *
                 */
                auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_neighbours_imp
                )
                        mutable
                {
                    for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                        auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                        lhs.push_back(rhs);
                    }
                    if constexpr (__j < ElementTraits<t_element, t_domain>::template getNumNeighbours<__i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < ElementTraits<t_element, t_domain>::getNumNeighbours() - 1) {
                        initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                    }
                };
                /*
                 *
                 */
                initialize_components(initialize_components);
                initialize_neighbours(initialize_neighbours);
                this->template getElement<_k>()->tag_ = this->tag_;
                this->template getElement<_k>()->domains_ = this->domains_;
                this->template getElement<_k>()->initialize(mesh_data);
                // if constexpr (_k < std::tuple_size_v<t_Elements> - 1) {
                //     initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                // }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

        template<Element __element = t_element>
        void
        initialize(
                Mesh<t_domain> & mesh_data
        )
        requires(__element.isNode())
        {
            /*
              *
              */
            auto initialize_element = [&] <lolita::index _k = 0u> (
                    auto & initialize_element_imp
            )
                    constexpr mutable
            {
                /*
                 *
                 */
                auto initialize_coordinates = [&] ()
                {
                    this->template getElement<_k>()->components_ = this->components_;
                    this->template getElement<_k>()->coordinates_ = this->coordinates_;
                };
                /*
                 *
                 */
                auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_neighbours_imp
                )
                        mutable
                {
                    for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                        auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                        lhs.push_back(rhs);
                    }
                    if constexpr (__j < ElementTraits<t_element, t_domain>::template getNumNeighbours<__i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < ElementTraits<t_element, t_domain>::getNumNeighbours() - 1) {
                        initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                    }
                };
                /*
                 *
                 */
                initialize_coordinates();
                initialize_neighbours(initialize_neighbours);
                this->template getElement<_k>()->tag_ = this->tag_;
                this->template getElement<_k>()->domains_ = this->domains_;
                this->template getElement<_k>()->initialize(mesh_data);
                // if constexpr (_k < std::tuple_size_v<t_Elements> - 1) {
                //     initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                // }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

    };

}


#endif /* C31C92CF_B3D3_4D98_9746_7FD9EEFC72FC */
