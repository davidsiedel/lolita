//
// Created by dsiedel on 04/06/22.
//

#ifndef LOLITA_LOLITA_CORE_5_HXX
#define LOLITA_LOLITA_CORE_5_HXX

#include <execution>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_1.hxx"
#include "lolita/lolita_core_2.hxx"
#include "lolita/lolita_core_4.hxx"
#include "lolita/lolita_core_5.hxx"

namespace lolita::core2::finite_element
{

    /**
     * @brief
     * @tparam _T
     * @tparam _element
     * @tparam _domain
     * @tparam _finite_element
     */
    template<
            template<lolita::core2::geometry::Element, lolita::domain::Domain, auto...> typename _T,
            lolita::core2::geometry::Element _element,
            lolita::domain::Domain _domain,
            auto... _finite_element
    >
    struct FiniteElementGeometry
    {

    private:

        /**
         * @brief
         */
        using _ElementDescription = lolita::core2::geometry::ElementGeometryTraits<_element, _domain>;

        /**
         * @brief
         */
        template<lolita::core2::geometry::Element __element, lolita::domain::Domain __domain, auto... __finite_element>
        using _ElementPointer = std::shared_ptr<_T<__element, __domain, __finite_element...>>;

        /**
         * @brief
         */
        template<lolita::core2::geometry::Element>
        struct _ConnexionPolicy;

        /**
         * @brief
         * @tparam __element
         */
        template<lolita::core2::geometry::Element __element>
        requires(!__element.isPoint())
        struct _ConnexionPolicy<__element>
        {

            /**
             * @brief
             */
            using _Components = typename lolita::core2::geometry::ElementGeometryTraits<__element, _domain>::template InnerConnectivity<_ElementPointer, _finite_element...>;

            /**
             * @brief
             */
            using _Neighbours = typename lolita::core2::geometry::ElementGeometryTraits<__element, _domain>::template OuterConnectivity<_ElementPointer, _finite_element...>;

        };

        /**
         * @brief
         * @tparam __element
         */
        template<lolita::core2::geometry::Element __element>
        requires(__element.isPoint())
        struct _ConnexionPolicy<__element>
        {

            /**
             * @brief
             */
            using _Components = std::shared_ptr<lolita::domain::Point>;

            /**
             * @brief
             */
            using _Neighbours = typename lolita::core2::geometry::ElementGeometryTraits<__element, _domain>::template OuterConnectivity<_ElementPointer, _finite_element...>;

        };

    public:

        /**
         * @brief The element object
         */
        lolita::core2::geometry::Element const static constexpr element_ = _element;

        /**
         * @brief Element components as shared pointers
         */
        using Components = typename _ConnexionPolicy<_element>::_Components;

        /**
         * @brief Element neighbours as shared pointers
         */
        using Neighbours = typename _ConnexionPolicy<_element>::_Neighbours;

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        lolita::boolean
        operator==(
                FiniteElementGeometry const & other
        )
        const = default;

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        lolita::boolean
        operator!=(
                FiniteElementGeometry const & other
        )
        const = default;

        /**
         * @brief
         * @return
         */
        std::basic_string<lolita::character>
        hash()
        const
        {
            if constexpr (_element.isPoint()) {
                return std::to_string(this->tag_);
            }
            else {
                std::basic_stringstream<lolita::character> hash;
                auto const & nodes = getComponents<_element.dim_ - 1, 0>();
                for (auto const & node : nodes) {
                    hash << node->hash();
                }
                return hash.str();
            }
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @param i
         * @param j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        static constexpr
        lolita::index
        getComponentNodeConnection(
                lolita::index i,
                lolita::index j
        )
        requires(!_element.isPoint())
        {
            return std::get<_j>(std::get<_i>(lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::node_connectivity_))[i][j];
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Neighbours>> &
        getNeighbours()
        {
            return std::get<_j>(std::get<_i>(neighbours_));
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Neighbours>> const &
        getNeighbours()
        const
        {
            return std::get<_j>(std::get<_i>(neighbours_));
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Components>> &
        getComponents()
        requires(!_element.isPoint())
        {
            return std::get<_j>(std::get<_i>(components_));
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Components>> const &
        getComponents()
        const
        requires(!_element.isPoint())
        {
            return std::get<_j>(std::get<_i>(components_));
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @param i
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        lolita::index
        getComponentIndex(
                lolita::index i
        )
        const
        requires(!_element.isPoint())
        {
//            auto const constexpr _component = lolita::core2::component<_element, _domain, _i, _j>();
//            auto const constexpr _position = lolita::core2::neighbourPosition<_component, _domain, _element>();
            auto const constexpr _component = _ElementDescription::template getComponent<_i, _j>();
            auto const constexpr _position = lolita::core2::geometry::ElementGeometryTraits<_component, _domain>::template getNeighbourCoordinates<_element>();
            auto const & items = getComponents<_i, _j>()[i]->template getNeighbours<_position.dim_, _position.tag_>();
            auto is_equal = [&] (
                    std::shared_ptr<FiniteElementGeometry> const & ptr_element
            )
            {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @param i
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        lolita::integer
        getComponentOrientation(
                lolita::index i
        )
        const
        requires(!_element.isPoint())
        {
            return getComponentIndex<_i, _j>(i) == 0 ? 1 : -1;
        }

        /**
         * @brief
         */
        Neighbours neighbours_;

        /**
         * @brief
         */
        Components components_;

        /**
         * @brief
         */
        std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;

        /**
         * @brief
         */
        lolita::natural tag_;

    };

    template<lolita::core2::geometry::Element _element, lolita::domain::Domain _domain, auto _finite_element>
    struct FEObject;

    template<lolita::core2::geometry::Element _element, lolita::domain::Domain _domain, auto _element_group>
    struct FEObject : public lolita::core2::finite_element::FiniteElementGeometry<FEObject, _element, _domain, _element_group>
    {

        lolita::boolean const static constexpr is_finite_element_ = false;

    private:

        /**
         * @brief
         */
        struct ElementPolicy
        {

            /**
             * @brief
             */
            using ElementPointers = typename std::remove_cvref_t<decltype(_element_group)>::template ElementPointers<FEObject, _element, _domain>;

            /**
             * @brief
             */
            using Elements = typename std::remove_cvref_t<decltype(_element_group)>::template Elements<FEObject, _element, _domain>;

        };

        /**
         * @brief
         */
        using _ElementPointers = typename ElementPolicy::ElementPointers;

        /**
         * @brief
         */
        using _Elements = typename ElementPolicy::Elements;

    public:

        /**
         * @brief
         * @tparam _i
         * @return
         */
        template<lolita::index _i>
        std::tuple_element_t<_i, _ElementPointers> const &
        getElement()
        const
        {
            return std::get<_i>(elements_);
        }

        /**
         * @brief
         * @tparam _i
         * @return
         */
        template<lolita::index _i>
        std::tuple_element_t<_i, _ElementPointers> &
        getElement()
        {
            return std::get<_i>(elements_);
        }

        /**
         * @brief
         */
        void
        make()
        {
            auto make_elements = [&] <lolita::index _k = 0> (auto & make_elements_imp) mutable {
                using _Element = std::tuple_element_t<_k, _Elements>;
                this->template getElement<_k>() = std::make_shared<_Element>(_Element());
                if constexpr (!_Element::is_finite_element_) {
                    this->template getElement<_k>()->make();
                }
                if constexpr (_k < std::tuple_size_v<_Elements> - 1) {
                    make_elements_imp.template operator()<_k + 1u>(make_elements_imp);
                }
            };
            make_elements(make_elements);
        }

        template<lolita::core2::geometry::Element __element = _element, auto _arg>
        static
        void
        makeElement(
                lolita::core2::mesh::ElementInitializationData<_element> const & initialization_data,
                lolita::core2::mesh::Mesh<_domain, _arg> & mesh_data
        );

        template<lolita::core2::geometry::Element ___element = _element, auto _arg>
        static
        void
        makeElement(
                lolita::core2::mesh::ElementInitializationData<_element> const & initialization_data,
                lolita::core2::mesh::Mesh<_domain, _arg> & mesh_data
        )
        requires(!___element.isPoint())
        {
            /*
             *
             */
            auto get_element_hash = [&] <lolita::core2::geometry::Element __element> (
                    std::array<lolita::index, __element.num_nodes_> node_tags
            )
            {
                auto element_hash = std::basic_stringstream<lolita::character>();
                if constexpr (node_tags.size() > 0) {
                    std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
                }
                for (int i = 0; i < node_tags.size(); ++i) {
                    element_hash << node_tags[i];
                }
                return element_hash.str();
            };
            /*
             * make element
             */
            auto make_element = [&] <lolita::core2::geometry::Element __element = _element, lolita::index _i = 0, lolita::index _j = 0> (
                    std::shared_ptr<lolita::core2::finite_element::FEObject<__element, _domain, _element_group>> & ptr_element,
                    lolita::core2::mesh::ElementInitializationData<__element> const & element_initialization_data,
                    auto & make_element_imp
            )
                    mutable
            {
                if constexpr (!__element.isPoint()) {
                    using __ElementDescription = lolita::core2::geometry::ElementGeometryTraits<__element, _domain>;
                    using __ComponentDescription = lolita::core2::geometry::ElementGeometryTraits<__ElementDescription::template getComponent<_i, _j>(), _domain>;
                    using __MeshDescription = lolita::core2::geometry::DomainGeometryTraits<_domain>;
                    auto const constexpr _is_initialized = _i == 0 && _j == 0;
                    auto const constexpr _component = __ElementDescription::template getComponent<_i, _j>();
                    auto const constexpr _component_coordinates = __MeshDescription::template getElementCoordinates<_component>();
                    auto const constexpr _neighbour_coordinates = __ComponentDescription::template getNeighbourCoordinates<__element>();
                    auto const constexpr _element_coordinates = __MeshDescription::template getElementCoordinates<__element>();
                    auto const constexpr _node_coordinates = __ElementDescription::template getComponentCoordinates<lolita::core2::geometry::Element::Node()>();
                    using __Component = lolita::core2::finite_element::FEObject<_component, _domain, _element_group>;
                    using __Self = lolita::core2::finite_element::FEObject<__element, _domain, _element_group>;
                    auto & components = mesh_data.elements_.template getElements<_component_coordinates.dim_, _component_coordinates.tag_>();
                    auto & element_component_array = ptr_element->template getComponents<_i, _j>();
                    for (auto i = 0; i < element_component_array.size(); ++i) {
                        auto component_hash = std::basic_string<lolita::character>();
                        if constexpr(!_component.isPoint()) {
                            auto component_initialization_data = lolita::core2::mesh::ElementInitializationData<_component>();
                            for (int j = 0; j < _component.num_nodes_; ++j) {
                                auto const k = __Self::template getComponentNodeConnection<_i, _j>(i, j);
                                component_initialization_data.node_tags_[j] = element_initialization_data.node_tags_[k];
                            }
                            component_initialization_data.tag_ = components.size();
                            component_hash = get_element_hash.template operator ()<_component>(component_initialization_data.node_tags_);
                            if (!components.contains(component_hash)) {
                                auto ptr_component = std::make_shared<__Component>(__Component());
                                make_element_imp.template operator ()<_component, 0, 0>(ptr_component, component_initialization_data, make_element_imp);
                            }
                        }
                        else {
                            component_hash = std::to_string(element_initialization_data.node_tags_[__Self::template getComponentNodeConnection<_i, _j>(i, 0)]);
                        }
                        element_component_array[i] = components[component_hash];
                        components[component_hash]->template getNeighbours<_neighbour_coordinates.dim_, _neighbour_coordinates.tag_>().push_back(ptr_element);
                    }
                    if constexpr (_j < __ElementDescription::template getNumComponents<_i>() - 1) {
                        make_element_imp.template operator()<__element, _i, _j + 1u>(ptr_element, element_initialization_data, make_element_imp);
                    }
                    else if constexpr (_i < __ElementDescription::getNumComponents() - 1) {
                        make_element_imp.template operator()<__element, _i + 1u, 0u>(ptr_element, element_initialization_data, make_element_imp);
                    }
                    if constexpr (_is_initialized) {
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
                        ptr_element->tag_ = element_initialization_data.tag_;
                        ptr_element->domains_.assign(domains.begin(), domains.end());
                        ptr_element->make();
                        auto element_hash = get_element_hash.template operator ()<__element>(element_initialization_data.node_tags_);
                        elements.insert({element_hash, ptr_element});
                    }
                }
            };
            /*
             *
             */
            auto ptr_element = std::make_shared<FEObject>(FEObject());
            make_element(ptr_element, initialization_data, make_element);
        }

        template<lolita::core2::geometry::Element ___element = _element, auto _arg>
        static
        void
        makeElement(
                lolita::core2::mesh::ElementInitializationData<_element> const & initialization_data,
                lolita::core2::mesh::Mesh<_domain, _arg> & mesh_data
        )
        requires(___element.isPoint())
        {
            /*
             *
             */
            auto make_element = [&] (
                    std::shared_ptr<lolita::core2::finite_element::FEObject<_element, _domain, _element_group>> & ptr_element
            )
                    mutable
            {
                auto const constexpr _element_coordinates = lolita::core2::geometry::DomainGeometryTraits<_domain>::template getElementCoordinates<_element>();
                ptr_element->tag_ = initialization_data.tag_;
                ptr_element->domains_ = initialization_data.domains_;
                ptr_element->components_ = std::make_shared<lolita::domain::Point>(initialization_data.coordinates_);
                ptr_element->make();
                auto elem_hash = std::to_string(initialization_data.tag_);
                mesh_data.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>().insert({elem_hash, ptr_element});
            };
            /*
             *
             */
            auto ptr_element = std::make_shared<FEObject>(FEObject());
            make_element(ptr_element);
        }

        template<lolita::core2::geometry::Element __element = _element, auto _arg>
        void
        initialize(
                lolita::core2::mesh::Mesh<_domain, _arg> & mesh_data
        );

        template<lolita::core2::geometry::Element __element = _element, auto _arg>
        void
        initialize(
                lolita::core2::mesh::Mesh<_domain, _arg> & mesh_data
        )
        requires(!__element.isPoint())
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
                    for (int i = 0; i < lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::template getNumComponents<__i, __j>(); ++i) {
                        auto & rhs = this->template getComponents<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getComponents<__i, __j>()[i];
                        lhs = rhs;
                    }
                    if constexpr (__j < lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::template getNumComponents<__i>() - 1) {
                        initialize_components_imp.template operator()<__i, __j + 1u>(initialize_components_imp);
                    }
                    else if constexpr (__i < lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::getNumComponents() - 1) {
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
                    if constexpr (__j < lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::template getNumNeighbours<__i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::getNumNeighbours() - 1) {
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
                if constexpr (_k < std::tuple_size_v<_Elements> - 1) {
                    initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

        template<lolita::core2::geometry::Element __element = _element, auto _arg>
        void
        initialize(
                lolita::core2::mesh::Mesh<_domain, _arg> & mesh_data
        )
        requires(__element.isPoint())
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
                    if constexpr (__j < lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::template getNumNeighbours<__i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::getNumNeighbours() - 1) {
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
                if constexpr (_k < std::tuple_size_v<_Elements> - 1) {
                    initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

        _ElementPointers elements_;

    };

    template<lolita::core2::geometry::Element _element, lolita::domain::Domain _domain, lolita::finite_element::FiniteElementConcept auto _element_group>
    struct FEObject<_element, _domain, _element_group> : lolita::core2::finite_element::FiniteElementGeometry<FEObject, _element, _domain, _element_group>
//    : public lolita::core2::finite_element::unknown::FiniteElementUnknowns<_element, _domain, _element_group>
    {

        lolita::boolean const static constexpr is_finite_element_ = true;

        template<auto... __args>
        void
        initialize(
                lolita::core2::mesh::Mesh<_domain, __args...> & mesh
        )
        {
            std::cout << "making element : " << this->hash() << std::endl;
//            this->setBehaviour(mesh);
//            this->setLoads(mesh);
//            this->setUnknowns(mesh);
        }

    };

}

#endif //LOLITA_LOLITA_CORE_5_HXX
