//
// Created by dsiedel on 14/06/22.
//

#ifndef LOLITA_LOLITA_CORE_5_FINAL_HXX
#define LOLITA_LOLITA_CORE_5_FINAL_HXX

#include <execution>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_1.hxx"
#include "lolita/lolita_core_2.hxx"
#include "lolita/lolita_core_4.hxx"
#include "lolita/lolita_core_5_000_connectivity.hxx"
#include "lolita/lolita_core_5_001_base.hxx"
#include "lolita/lolita_core_5_002_basis.hxx"
#include "lolita/lolita_core_5_003_unknown.hxx"
#include "lolita/lolita_core_5_004_module.hxx"

namespace lolita::core2::finite_element
{

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementTraits
    {

    private:

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        getOrdQuadrature()
        {
            auto ord_quadrature = t_finite_element.ord_quadrature_;
            return t_domain.frame_.isAxiSymmetric() ? 2 * ord_quadrature + 1 : 2 * ord_quadrature;
        }


    public:

        /**
         * @brief
         */
        lolita::finite_element::FiniteElementConcept auto const static constexpr finite_element_ = t_finite_element;

        /**
         * @brief
         */
        using ElementDescription = lolita::core2::geometry::ElementTraits<t_element, t_domain>;

        /**
         * @brief
         */
        using Quadrature = lolita::core2::geometry::ElementQuadratureTraits<t_element, t_finite_element.quadrature_, getOrdQuadrature()>;

        /**
         * @brief
         */
        using Field = lolita::core2::field::TensorPolicy<t_finite_element.unknown_.tensor_, t_domain.dim_>;

        /**
         * @brief
         * @tparam t_method
         * @return
         */
        template<lolita::finite_element::FiniteElementMethod t_method>
        static constexpr
        lolita::boolean
        hasMethod()
        {
            return finite_element_.discretization_ == t_method;
        }

//        /**
//         * @brief
//         * @tparam t_unknown
//         * @return
//         */
//        template<lolita_fld::Field... t_field>
//        static constexpr
//        lolita::integer
//        getDimUnknowns()
//        {
//            auto num_element_unknowns = FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>::template getDimUnknowns<t_unknown>();
//            auto get_num_components_unknowns = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
//                auto const constexpr t_component = ElementDescription::template getComponent<t_i, t_j>();
//                auto num_component_unknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>::template getDimUnknowns<t_unknown>();
//                auto num_components = ElementDescription::template getNumComponents<t_i, t_j>();
//                num_element_unknowns += num_component_unknowns * num_components;
//                if constexpr (t_j < ElementDescription::template getNumComponents<t_i>() - 1) {
//                    self.template operator()<t_i, t_j + 1>(self);
//                }
//                else if constexpr (t_i < ElementDescription::template getNumComponents<>() - 1) {
//                    self.template operator()<t_i + 1, 0>(self);
//                }
//            };
//            if constexpr (!t_element.isPoint()) {
//                get_num_components_unknowns(get_num_components_unknowns);
//            }
//            return num_element_unknowns;
//        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        getDimUnknowns()
        {
            using t_ElementFieldUnknowns = FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>;
            auto num_structural_unknowns = t_ElementFieldUnknowns::template getDimUnknowns<unknown::Unknown::Structural()>();
            auto num_subsidiary_unknowns = t_ElementFieldUnknowns::template getDimUnknowns<unknown::Unknown::Subsidiary()>();
            auto get_num_components_unknowns = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                auto constexpr t_component = ElementDescription::template getComponent<t_i, t_j>();
                using t_NeighbourFieldUnknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>;
                auto num_structural_component_unknowns = t_NeighbourFieldUnknowns::template getDimUnknowns<unknown::Unknown::Structural()>();
                auto num_subsidiary_component_unknowns = t_NeighbourFieldUnknowns::template getDimUnknowns<unknown::Unknown::Subsidiary()>();
                auto num_components = ElementDescription::template getNumComponents<t_i, t_j>();
                num_structural_unknowns += num_structural_component_unknowns * num_components;
                num_subsidiary_unknowns += num_subsidiary_component_unknowns * num_components;
                if constexpr (t_j < ElementDescription::template getNumComponents<t_i>() - 1) {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < ElementDescription::template getNumComponents<>() - 1) {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            if constexpr (!t_element.isPoint()) {
                get_num_components_unknowns(get_num_components_unknowns);
            }
            return num_structural_unknowns + num_subsidiary_unknowns;
        }

        /**
         * @brief
         * @tparam t_unknown
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown>
        static constexpr
        lolita::integer
        getDimUnknowns()
        {
            auto num_element_unknowns = FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>::template getDimUnknowns<t_unknown>();
            auto get_num_components_unknowns = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                auto const constexpr t_component = ElementDescription::template getComponent<t_i, t_j>();
                auto num_component_unknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>::template getDimUnknowns<t_unknown>();
                auto num_components = ElementDescription::template getNumComponents<t_i, t_j>();
                num_element_unknowns += num_component_unknowns * num_components;
                if constexpr (t_j < ElementDescription::template getNumComponents<t_i>() - 1) {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < ElementDescription::template getNumComponents<>() - 1) {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            if constexpr (!t_element.isPoint()) {
                get_num_components_unknowns(get_num_components_unknowns);
            }
            return num_element_unknowns;
        }

        /**
         * @brief
         * @tparam t_unknown
         * @tparam t_i
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown, lolita::integer t_i>
        static constexpr
        lolita::integer
        getDimUnknowns()
        {
            if constexpr (t_i == 0) {
                return FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>::template getDimUnknowns<t_unknown>();
            }
            else {
                auto num_unknowns = 0;
                auto get_num_components_unknowns = [&] <lolita::integer t_j = 0> (auto & self) constexpr mutable {
                    auto const constexpr t_component = ElementDescription::template getComponent<t_i - 1, t_j>();
                    auto num_component_unknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>::template getDimUnknowns<t_unknown>();
                    auto num_components = ElementDescription::template getNumComponents<t_i - 1, t_j>();
                    num_unknowns += num_component_unknowns * num_components;
                    if constexpr (t_j < ElementDescription::template getNumComponents<t_i - 1>() - 1) {
                        self.template operator()<t_j + 1>(self);
                    }
                };
                if constexpr (!t_element.isPoint()) {
                    get_num_components_unknowns(get_num_components_unknowns);
                }
                return num_unknowns;
            }
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            using t_ElementFieldUnknowns = FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>;
            auto num_structural_unknowns = t_ElementFieldUnknowns::template getNumUnknowns<unknown::Unknown::Structural()>();
            auto num_subsidiary_unknowns = t_ElementFieldUnknowns::template getNumUnknowns<unknown::Unknown::Subsidiary()>();
            auto get_num_components_unknowns = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                auto constexpr t_component = ElementDescription::template getComponent<t_i, t_j>();
                using t_NeighbourFieldUnknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>;
                auto num_structural_component_unknowns = t_NeighbourFieldUnknowns::template getNumUnknowns<unknown::Unknown::Structural()>();
                auto num_subsidiary_component_unknowns = t_NeighbourFieldUnknowns::template getNumUnknowns<unknown::Unknown::Subsidiary()>();
                auto num_components = ElementDescription::template getNumComponents<t_i, t_j>();
                num_structural_unknowns += num_structural_component_unknowns * num_components;
                num_subsidiary_unknowns += num_subsidiary_component_unknowns * num_components;
                if constexpr (t_j < ElementDescription::template getNumComponents<t_i>() - 1) {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < ElementDescription::template getNumComponents<>() - 1) {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            if constexpr (!t_element.isPoint()) {
                get_num_components_unknowns(get_num_components_unknowns);
            }
            return num_structural_unknowns + num_subsidiary_unknowns;
        }

        /**
         * @brief
         * @tparam t_unknown
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown>
        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            auto num_element_unknowns = FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>::template getNumUnknowns<t_unknown>();
            auto get_num_components_unknowns = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                auto const constexpr t_component = ElementDescription::template getComponent<t_i, t_j>();
                auto num_component_unknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>::template getNumUnknowns<t_unknown>();
                auto num_components = ElementDescription::template getNumComponents<t_i, t_j>();
                num_element_unknowns += num_component_unknowns * num_components;
                if constexpr (t_j < ElementDescription::template getNumComponents<t_i>() - 1) {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < ElementDescription::template getNumComponents<>() - 1) {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            if constexpr (!t_element.isPoint()) {
                get_num_components_unknowns(get_num_components_unknowns);
            }
            return num_element_unknowns;
        }

        /**
         * @brief
         * @tparam t_unknown
         * @tparam t_i
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown, lolita::integer t_i>
        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            if constexpr (t_i == 0) {
                return FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>::template getNumUnknowns<t_unknown>();
            }
            else {
                auto num_unknowns = 0;
                auto get_num_components_unknowns = [&] <lolita::integer t_j = 0> (auto & self) constexpr mutable {
                    auto const constexpr t_component = ElementDescription::template getComponent<t_i - 1, t_j>();
                    auto num_component_unknowns = FiniteElementFieldUnknowns<t_component, t_domain, t_finite_element>::template getNumUnknowns<t_unknown>();
                    auto num_components = ElementDescription::template getNumComponents<t_i - 1, t_j>();
                    num_unknowns += num_component_unknowns * num_components;
                    if constexpr (t_j < ElementDescription::template getNumComponents<t_i - 1>() - 1) {
                        self.template operator()<t_j + 1>(self);
                    }
                };
                if constexpr (!t_element.isPoint()) {
                    get_num_components_unknowns(get_num_components_unknowns);
                }
                return num_unknowns;
            }
        }

        /**
         * @brief
         * @tparam t_mapping
         * @return
         */
        template<lolita::field::Mapping t_mapping>
        static constexpr
        lolita::integer
        getMappingSize()
        {
            return lolita::core2::field::MappingPolicy<t_finite_element.unknown_.tensor_, t_domain, t_mapping>::shape_.size_;
        }

        /**
         * @brief
         * @tparam t_mapping
         * @return
         */
        template<lolita::field::Mapping t_mapping>
        static constexpr
        lolita::matrix::VectorBlock
        getMappingBlock()
        {
            auto mapping_row = lolita::integer(0);
            for (auto mapping : t_finite_element.unknown_.mappings_) {
                if (mapping == t_mapping) {
                    return lolita::matrix::VectorBlock(mapping_row, mapping_row + getMappingSize<t_mapping>());
                }
                mapping_row += getMappingSize<t_mapping>();
            }
            return lolita::matrix::VectorBlock();
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        getGeneralizedStrainNumRows()
        {
            auto mapping_size = lolita::integer(0);
            auto set_dim_mapping = [&] <lolita::integer t_i = 0> (auto & self) constexpr mutable {
                mapping_size += getMappingSize<t_finite_element.unknown_.mappings_[t_i]>();
                if constexpr (t_i < t_finite_element.unknown_.mappings_.size() - 1) {
                    self.template operator()<t_i + 1>(self);
                }
            };
            set_dim_mapping(set_dim_mapping);
            return mapping_size;
        }

    };

//    template<std::array<lolita::character, 3> t_str>
//    static constexpr
//    lolita::integer
//    makeTest()
//    {
//        if constexpr (t_str == std::array<lolita::character, 3>{"AAA"}) {
//            return 1;
//        }
//        else {
//            return -1;
//        }
//    }

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_element_group
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, auto t_element_group>
    struct FiniteElement : FiniteElementConnectivity<t_element, t_domain, t_element_group>
    {

        /**
         * @brief
         */
        lolita::boolean const static constexpr is_finite_element_ = false;

    private:

        /**
         * @brief
         */
        struct ElementGroupTraits
        {

            /**
             * @brief
             */
            using ElementPointers = typename std::remove_cvref_t<decltype(t_element_group)>::template ElementPointers<FiniteElement, t_element, t_domain>;

            /**
             * @brief
             */
            using Elements = typename std::remove_cvref_t<decltype(t_element_group)>::template Elements<FiniteElement, t_element, t_domain>;

        };

        /**
         * @brief
         */
        using t_ElementPointers = typename ElementGroupTraits::ElementPointers;

        /**
         * @brief
         */
        using t_Elements = typename ElementGroupTraits::Elements;

    public:

        /**
         * @brief
         * @tparam t_i
         * @return
         */
        template<lolita::index t_i>
        std::tuple_element_t<t_i, t_ElementPointers> const &
        getElement()
        const
        {
            return std::get<t_i>(elements_);
        }

        /**
         * @brief
         * @tparam t_i
         * @return
         */
        template<lolita::index t_i>
        std::tuple_element_t<t_i, t_ElementPointers> &
        getElement()
        {
            return std::get<t_i>(elements_);
        }

        /**
         * @brief
         */
        void
        make()
        {
            auto make_elements = [&] <lolita::index t_k = 0> (auto & make_elements_imp) mutable {
                using t_Element = std::tuple_element_t<t_k, t_Elements>;
                this->template getElement<t_k>() = std::make_shared<t_Element>(t_Element());
                if constexpr (!t_Element::is_finite_element_) {
                    this->template getElement<t_k>()->make();
                }
                if constexpr (t_k < std::tuple_size_v<t_Elements> - 1) {
                    make_elements_imp.template operator()<t_k + 1u>(make_elements_imp);
                }
            };
            make_elements(make_elements);
        }

        template<lolita::core2::geometry::Element t__element = t_element, auto t_arg>
        static
        void
        makeElement(
                lolita::core2::mesh::ElementInitializationData<t_element> const & initialization_data,
                lolita::core2::mesh::Mesh<t_domain, t_arg> & mesh_data
        );

        template<lolita::core2::geometry::Element ___element = t_element, auto t_arg>
        static
        void
        makeElement(
                lolita::core2::mesh::ElementInitializationData<t_element> const & initialization_data,
                lolita::core2::mesh::Mesh<t_domain, t_arg> & mesh_data
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
            auto make_element = [&] <lolita::core2::geometry::Element __element = t_element, lolita::index _i = 0, lolita::index _j = 0> (
                    std::shared_ptr<lolita::core2::finite_element::FiniteElement<__element, t_domain, t_element_group>> & ptr_element,
                    lolita::core2::mesh::ElementInitializationData<__element> const & element_initialization_data,
                    auto & make_element_imp
            )
                    mutable
            {
                if constexpr (!__element.isPoint()) {
                    using __ElementDescription = lolita::core2::geometry::ElementTraits<__element, t_domain>;
                    using __ComponentDescription = lolita::core2::geometry::ElementTraits<__ElementDescription::template getComponent<_i, _j>(), t_domain>;
                    using __MeshDescription = lolita::core2::geometry::DomainTraits<t_domain>;
                    auto const constexpr _is_initialized = _i == 0 && _j == 0;
                    auto const constexpr _component = __ElementDescription::template getComponent<_i, _j>();
                    auto const constexpr _component_coordinates = __MeshDescription::template getElementCoordinates<_component>();
                    auto const constexpr _neighbour_coordinates = __ComponentDescription::template getNeighbourCoordinates<__element>();
                    auto const constexpr _element_coordinates = __MeshDescription::template getElementCoordinates<__element>();
                    auto const constexpr _node_coordinates = __ElementDescription::template getComponentCoordinates<lolita::core2::geometry::Element::Node()>();
                    using __Component = lolita::core2::finite_element::FiniteElement<_component, t_domain, t_element_group>;
                    using __Self = lolita::core2::finite_element::FiniteElement<__element, t_domain, t_element_group>;
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
            auto ptr_element = std::make_shared<FiniteElement>(FiniteElement());
            make_element(ptr_element, initialization_data, make_element);
        }

        template<lolita::core2::geometry::Element ___element = t_element, auto t_arg>
        static
        void
        makeElement(
                lolita::core2::mesh::ElementInitializationData<t_element> const & initialization_data,
                lolita::core2::mesh::Mesh<t_domain, t_arg> & mesh_data
        )
        requires(___element.isPoint())
        {
            /*
             *
             */
            auto make_element = [&] (
                    std::shared_ptr<lolita::core2::finite_element::FiniteElement<t_element, t_domain, t_element_group>> & ptr_element
            )
                    mutable
            {
                auto const constexpr _element_coordinates = lolita::core2::geometry::DomainTraits<t_domain>::template getElementCoordinates<t_element>();
                ptr_element->tag_ = initialization_data.tag_;
                ptr_element->domains_ = initialization_data.domains_;
                ptr_element->coordinates_ = std::make_shared<lolita::domain::Point>(initialization_data.coordinates_);
//                ptr_element->coordinates_ = std::make_unique<lolita::domain::Point>(initialization_data.coordinates_);
                ptr_element->make();
                auto elem_hash = std::to_string(initialization_data.tag_);
                mesh_data.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>().insert({elem_hash, ptr_element});
            };
            /*
             *
             */
            auto ptr_element = std::make_shared<FiniteElement>(FiniteElement());
            make_element(ptr_element);
        }

        template<lolita::core2::geometry::Element __element = t_element, auto t_arg>
        void
        initialize(
                lolita::core2::mesh::Mesh<t_domain, t_arg> & mesh_data
        );

        template<lolita::core2::geometry::Element __element = t_element, auto t_arg>
        void
        initialize(
                lolita::core2::mesh::Mesh<t_domain, t_arg> & mesh_data
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
                    for (int i = 0; i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<__i, __j>(); ++i) {
                        auto & rhs = this->template getComponents<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getComponents<__i, __j>()[i];
                        lhs = rhs;
                    }
                    if constexpr (__j < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<__i>() - 1) {
                        initialize_components_imp.template operator()<__i, __j + 1u>(initialize_components_imp);
                    }
                    else if constexpr (__i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::getNumComponents() - 1) {
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
                    if constexpr (__j < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumNeighbours<__i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::getNumNeighbours() - 1) {
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
                if constexpr (_k < std::tuple_size_v<t_Elements> - 1) {
                    initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

        template<lolita::core2::geometry::Element __element = t_element, auto _arg>
        void
        initialize(
                lolita::core2::mesh::Mesh<t_domain, _arg> & mesh_data
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
                    if constexpr (__j < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumNeighbours<__i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::getNumNeighbours() - 1) {
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
                if constexpr (_k < std::tuple_size_v<t_Elements> - 1) {
                    initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

        t_ElementPointers elements_;

    };

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElement<t_element, t_domain, t_finite_element> : FiniteElementCell_n<t_element, t_domain, t_finite_element>
    {

        lolita::boolean const static constexpr is_finite_element_ = true;

        template<auto t_element_group>
        void
        initialize(
                lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
        )
        {
            std::cout << "making element : " << this->hash() << std::endl;
            this->setBehaviour(mesh);
            this->setLoads(mesh);
            this->setUnknowns(mesh);
            this->getCurrentCoordinates();
            this->template getUnknowns<unknown::Unknown::Subsidiary()>();
            this->template getUnknowns<unknown::Unknown::Structural()>();
            this->setIntegrationPoint();
//            if constexpr (t_element.isSub(t_domain, 0)) {
//                this->setGeneralizedGradients();
//            }
            this->action();
        }

        void
        action()
        const
        {
            auto pt = lolita::domain::Point();
            this->template getBasisEvaluation<lolita::core2::finite_element::basis::Basis::Monomial(), 1>(pt);
        }

    };

}

#endif //LOLITA_LOLITA_CORE_5_FINAL_HXX
