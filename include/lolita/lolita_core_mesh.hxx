//
// Created by dsiedel on 10/05/22.
//

#ifndef LOLITA_LOLITA_CORE_MESH_HXX
#define LOLITA_LOLITA_CORE_MESH_HXX

#include <execution>
#include <ostream>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_element.hxx"
#include "lolita/lolita_core_finite_element.hxx"

namespace lolita::core::mesh
{

    using Element = core::element::Element;

//    template<Element E, lolita::geometry::Domain D, auto Med>
//    using FiniteElementE = core::element::FiniteElementF<E, D, Med>;

    template<Element E, lolita::geometry::Domain D, lolita::finite_element::FiniteElementConcept auto... _finite_element>
    using FiniteElementFinal = core::element::FiniteElementFinal<E, D, _finite_element...>;

    template<lolita::mesh::MeshFormatType, lolita::geometry::Domain D, lolita::finite_element::FiniteElementConcept auto... _finite_element>
    struct MeshModule;

    struct MeshModuleBase
    {

        MeshModuleBase(
                lolita::utility::File const &
                mesh_file
        )
        :
        file_(mesh_file)
        {}

        lolita::utility::File const & file_;

    };

    template<lolita::mesh::MeshFormatType Mft, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto... _finite_element>
    struct MeshBase
    {

    private:

        template<Element __element, lolita::geometry::Domain __domain, lolita::finite_element::FiniteElementConcept auto... __finite_element>
        using ElementPointerMap = std::unordered_map<std::basic_string<lolita::character>, std::shared_ptr<FiniteElementFinal<__element, __domain, __finite_element...>>>;

        using Module = typename MeshModule<Mft, _domain, _finite_element...>::Module;

        using Implementation = typename MeshModule<Mft, _domain, _finite_element...>::Implementation;

        struct DegreeOfFreedomIndex
        {

            constexpr
            DegreeOfFreedomIndex()
            :
            unknown_index(0),
            binding_index(0)
            {}

            constexpr
            DegreeOfFreedomIndex(
                    lolita::natural &&
                    unknown_index_arg,
                    lolita::natural &&
                    binding_index_arg
            )
            :
            unknown_index(std::forward<lolita::natural>(unknown_index_arg)),
            binding_index(std::forward<lolita::natural>(binding_index_arg))
            {}

            lolita::natural unknown_index;

            lolita::natural binding_index;

        };

//        using DegreeOfFreedomIndices = std::array<DegreeOfFreedomIndex, sizeof...(_finite_element)>;
//
//        using Loads = std::vector<lolita::finite_element::Load>;
//
//        using Behaviours = std::vector<lolita::behaviour::MgisBehaviour>;

//        using ElementPointersSets = std::unordered_map<std::basic_string<lolita::character>, ElementPointers>;

        std::shared_ptr<lolita::finite_element::LoadComponent> const null_load_;

    public:

        using ElementPointers = lolita::core::element::ElementCollection<ElementPointerMap, _domain, _finite_element...>;

        std::array<DegreeOfFreedomIndex, sizeof...(_finite_element)> unknown_indices;

        Module module_;

        ElementPointers elements_;

        std::unordered_map<std::basic_string<lolita::character>, ElementPointers> element_sets_;

        std::vector<lolita::finite_element::Load> loads_;

        std::vector<lolita::behaviour::MgisBehaviour> behaviours_;

        MeshBase(
                lolita::utility::File const &
                mesh_file,
                std::vector<lolita::finite_element::Load> const &
                loads_arg,
                std::vector<lolita::behaviour::MgisBehaviour> const &
                behaviours_arg
        )
        :
        module_(mesh_file),
        null_load_(std::make_shared<lolita::finite_element::LoadComponent>(lolita::finite_element::LoadComponent())),
        loads_(loads_arg),
        behaviours_(behaviours_arg)
        {
            makeNodes();
            makeElementSets();
            makeNodeSets();
            makeCells();
            makeNeighbourhood();
        }

        void
        makeElementSets()
        {
            static_cast<Implementation *>(this)->makeElementSets();
        }

        void
        makeNodeSets()
        {
            static_cast<Implementation *>(this)->makeNodeSets();
        }

        void
        makeNodes()
        {
            static_cast<Implementation *>(this)->makeNodes();
        }

        template<Element _element>
        std::shared_ptr<FiniteElementFinal<_element, _domain, _finite_element...>>
        makeElementT()
        const
        {
            using _Element = FiniteElementFinal<_element, _domain, _finite_element...>;
            return std::make_shared<_Element>(_Element());
        }

        void
        makeNode(
                lolita::index
                node_tag,
                lolita::geometry::Point const &
                node_coordinates
        )
        {
            /*
             *
             */
            auto set_finite_elements = [&] <lolita::index _k = 0u> (
                    std::shared_ptr<FiniteElementFinal<element::pnt_00, _domain, _finite_element...>> &
                    ptr_node,
                    auto &
                    self
            )
            constexpr mutable
            {
                ptr_node->template makeElement<_k>();
                if constexpr (_k < sizeof...(_finite_element) - 1) {
                    self.template operator()<_k + 1u>(ptr_node, self);
                }
            };
            /*
             *
             */
            auto set_finite_elements_components = [&] <lolita::index _k = 0u> (
                    std::shared_ptr<FiniteElementFinal<element::pnt_00, _domain, _finite_element...>> &
                    ptr_node,
                    auto &
                    self
            )
            constexpr mutable
            {
                ptr_node->template getElement<_k>()->coordinates_ = ptr_node->coordinates_;
                ptr_node->template getElement<_k>()->tag_ = ptr_node->tag_;
                if constexpr (_k < sizeof...(_finite_element) - 1) {
                    self.template operator()<_k + 1u>(ptr_node, self);
                }
            };
            /*
             *
             */
            auto set_node = [&] (
                    std::shared_ptr<FiniteElementFinal<element::pnt_00, _domain, _finite_element...>> &
                    ptr_node
            )
            {
                auto hash = std::to_string(node_tag);
                auto & nodes = elements_.template getElements<0, 0>();
                ptr_node->tag_ = nodes.size();
                ptr_node->coordinates_ = std::make_shared<lolita::geometry::Point>(node_coordinates);
                set_finite_elements(ptr_node, set_finite_elements);
                set_finite_elements_components(ptr_node, set_finite_elements_components);
                nodes.insert({hash, ptr_node});
            };
            /*
             *
             */
            auto ptr_node = makeElementT<element::pnt_00>();
            set_node(ptr_node);
        }

        template<Element EE>
        void
        makeCell(
                std::array<lolita::index, EE.num_nodes_> const &
                cell_node_tags
        )
        {
            /*
             *
             */
            auto get_element_domains = [&] <Element __element> (
                    std::array<lolita::index, __element.num_nodes_> const &
                    node_tags
            )
            {
                if constexpr (!element::PointConcept<__element>) {
                    std::vector<std::basic_string<lolita::character>> domain_names;
                    for (auto const & [set_name, element_tags]: element_sets_) {
//                    auto const & node_hashes = std::get<0>(std::get<0>(element_tags));
                        auto const & node_hashes = element_tags.template getElements<0, 0>();
                        lolita::boolean belongs = true;
                        for (auto node_tag: node_tags) {
                            auto node_hash = std::to_string(node_tag);
                            if (!node_hashes.contains(node_hash)) {
                                belongs = false;
                                break;
                            }
                        }
                        if (belongs) {
                            domain_names.push_back(set_name);
                        }
                    }
                    return domain_names;
                }
            };
            /*
             *
             */
            auto get_component_node_tags = [] <Element __element, auto _i, auto _j> (
                    std::array<lolita::index, __element.num_nodes_> const &
                    node_tags,
                    lolita::index
                    index
            )
            {
                if constexpr (!element::PointConcept<__element>) {
                    auto component_node_tags = std::array<lolita::index, element::component<__element, _domain, _i, _j>().num_nodes_>();
                    for (lolita::index j = 0; j < element::component<__element, _domain, _i, _j>().num_nodes_; ++j) {
                        auto const k = FiniteElementFinal<__element, _domain, _finite_element...>::template getComponentNodeConnection<_i, _j>(index, j);
                        component_node_tags[j] = node_tags[k];
                    }
                    return component_node_tags;
                }
            };
            /*
             *
             */
            auto get_element_hash = [] <Element __element> (
                    std::array<lolita::index, __element.num_nodes_>
                    node_tags
            )
            {
                std::basic_stringstream<lolita::character> element_hash;
                if constexpr (node_tags.size() > 0) {
                    std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
                }
                for (lolita::index i = 0; i < node_tags.size(); ++i) {
                    element_hash << node_tags[i];
                }
                return element_hash.str();
            };
            /*
             *
             */
            auto set_finite_elements = [&] <Element __element, lolita::index _k = 0u> (
                    std::shared_ptr<FiniteElementFinal<__element, _domain, _finite_element...>> &
                    ptr_element,
                    auto &
                    self
            )
            mutable
            {
                ptr_element->template makeElement<_k>();
                if constexpr (_k < sizeof...(_finite_element) - 1) {
                    self.template operator()<__element, _k + 1u>(ptr_element, self);
                }
            };
            /*
             *
             */
            auto set_finite_elements_components = [&] <Element __element, lolita::index _k = 0u, lolita::index _i = 0u, lolita::index _j = 0u> (
                    std::shared_ptr<FiniteElementFinal<__element, _domain, _finite_element...>> &
                    ptr_element,
                    auto &
                    self
            )
            mutable
            {
                for (int i = 0; i < element::numComponents<__element, _domain, _i, _j>(); ++i) {
                    auto & rhs = ptr_element->template getComponents<_i, _j>()[i]->template getElement<_k>();
                    auto & lhs = ptr_element->template getElement<_k>()->template getComponents<_i, _j>()[i];
                    lhs = rhs;
                }
                if constexpr (_j < element::numComponents<__element, _domain, _i>() - 1) {
                    self.template operator()<__element, _k, _i, _j + 1u>(ptr_element, self);
                }
                else if constexpr (_i < element::numComponents<__element, _domain>() - 1) {
                    self.template operator()<__element, _k, _i + 1u, 0u>(ptr_element, self);
                }
                else if constexpr (_k < sizeof...(_finite_element) - 1) {
                    self.template operator()<__element, _k + 1u, 0u, 0u>(ptr_element, self);
                }
            };
            /*
             *
             */
            auto set_element_n2 = [&] <Element __element> (
                    std::array<lolita::index, __element.num_nodes_> const &
                    element_node_tags,
                    std::shared_ptr<FiniteElementFinal<__element, _domain, _finite_element...>> &
                    ptr_element
            )
            {
                auto const constexpr _crd_elm = element::elementPosition<_domain, __element>();
                auto & elements = elements_.template getElements<_crd_elm[0], _crd_elm[1]>();
                auto element_hash = get_element_hash.template operator ()<__element>(element_node_tags);
                auto element_domains = get_element_domains.template operator()<__element>(element_node_tags);
                ptr_element->tag_ = elements.size();
                set_finite_elements.template operator()<__element>(ptr_element, set_finite_elements);
                set_finite_elements_components.template operator()<__element>(ptr_element, set_finite_elements_components);
                elements.insert({element_hash, ptr_element});
            };
            /*
             *
             */
            auto set_element_n = [&] <Element __element, lolita::index _i = 0, lolita::index _j = 0> (
                    std::array<lolita::index, __element.num_nodes_> const &
                    element_node_tags,
                    std::shared_ptr<FiniteElementFinal<__element, _domain, _finite_element...>> &
                    ptr_element,
                    auto &
                    self
            )
            mutable
            {
                if constexpr (!element::PointConcept<__element>) {
                    auto const constexpr _component = element::component<__element, _domain, _i, _j>();
                    auto const constexpr _component_coordinates = element::elementPosition<_domain, _component>();
                    auto & components = elements_.template getElements<_component_coordinates[0], _component_coordinates[1]>();
                    auto & element_component_array = ptr_element->template getComponents<_i, _j>();
                    for (auto i = 0; i < element_component_array.size(); ++i) {
                        auto component_node_tags = get_component_node_tags.template operator ()<__element, _i, _j>(element_node_tags, i);
                        auto component_hash = get_element_hash.template operator ()<_component>(component_node_tags);
                        if (!components.contains(component_hash)) {
                            auto ptr_component = makeElementT<_component>();
                            self.template operator ()<_component>(component_node_tags, ptr_component, self);
                        }
                        element_component_array[i] = components[component_hash];
                        auto const constexpr _element_coordinates = element::neighbourPosition<_component, _domain, __element>();
                        components[component_hash]->template getNeighbours<_element_coordinates[0], _element_coordinates[1]>().push_back(ptr_element);
                    }
                    if constexpr (_j < element::numComponents<__element, _domain, _i>() - 1) {
                        self.template operator()<__element, _i, _j + 1>(element_node_tags, ptr_element, self);
                    }
                    else if constexpr (_i < element::numComponents<__element, _domain>() - 1) {
                        self.template operator()<__element, _i + 1, 0>(element_node_tags, ptr_element, self);
                    }
                    if constexpr (_i == 0 && _j == 0) {
                        set_element_n2.template operator()<__element>(element_node_tags, ptr_element);
                    }
                }
            };
            /*
             *
             */
            auto ptr_cell = makeElementT<EE>();
            set_element_n.template operator ()<EE>(cell_node_tags, ptr_cell, set_element_n);
        }

        void
        makeCells(
        )
        {
            auto make_cells = [&] <auto I = 0> (
                    auto &
                    self
            )
            mutable
            {
                static_cast<Implementation *>(this)->template makeCells<element::element<_domain, _domain.dim_, I>()>();
                if constexpr (I < element::numElements<_domain, _domain.dim_>() - 1) {
                    self.template operator ()<I + 1>(self);
                }
            };
            make_cells(make_cells);
        }

        void
        makeNeighbourhood()
        {
            /*
             *
             */
            auto make_element_neighbourhood = [&] <Element E, lolita::index K = 0> (auto & self)
            mutable
            {
                if constexpr (!element::PointConcept<E>) {
                    auto const constexpr poss = element::elementPosition<_domain, E>();
                    auto & element_map = elements_.template getElements<poss[0], poss[1]>();
                    for (auto & element_map_item: element_map) {
                        auto & element = element_map_item.second;
                        auto const & element_nds = element->template getComponents<poss[0] - 1, 0>();
                        for (int i = 0; i < element_nds.size(); ++i) {
                            auto const & nde = element_nds[i];
                            auto const & ngs = nde->template getNeighbours<poss[0] - 1, K>();
                            for (auto const & neighbour : ngs) {
                                if (((neighbour->tag_ != element->tag_) && K == poss[1]) || (K != poss[1])) {
                                    auto & element_ngs = element->template getNeighbours<0, K>();
                                    lolita::boolean found = false;
                                    for (auto const & ngb: element_ngs) {
                                        if (ngb->tag_ == neighbour->tag_) {
                                            found = true;
                                            break;
                                        }
                                    }
                                    if (!found) {
                                        element_ngs.push_back(neighbour);
                                    }
                                }
                            }
                        }
                    }
                    if constexpr (K < element::numElements<_domain, poss[0]>() - 1) {
                        self.template operator()<E, K + 1>(self);
                    }
                }
            };
            /*
             *
             */
            auto make_neighbourhood = [&] <lolita::index I = 1, lolita::index J = 0> (
                    auto & self
            )
            mutable
            {
                make_element_neighbourhood.template operator()<element::element<_domain, I, J>()>(make_element_neighbourhood);
                if constexpr (J < element::numElements<_domain, I>() - 1) {
                    self.template operator()<I, J + 1>(self);
                }
                else if constexpr (I < element::numElements<_domain>() - 2) {
                    self.template operator()<I + 1, 0>(self);
                }
            };
            make_neighbourhood(make_neighbourhood);
        }

        friend
        std::ostream &
        operator<<(
                std::ostream & os,
                MeshBase const & base
        )
        {
            /*
             *
             */
            auto print_element_components = [&] <Element EE, auto I = 0, auto J = 0> (
                    auto const & elem_arg,
                    auto & self
            )
            mutable
            {
                if constexpr (!element::PointConcept<EE>) {
                    for (auto const & c_ : elem_arg->template getComponents<I, J>()) {
                        os << "layer : " << I << " type : " << J << " <-- " << c_->hash() << std::endl;
                    }
                    if constexpr (J < element::numComponents<EE, _domain, I>() - 1) {
                        self.template operator()<EE, I, J + 1>(elem_arg, self);
                    }
                    else if constexpr (I < element::numComponents<EE, _domain>() - 1) {
                        self.template operator()<EE, I + 1, 0>(elem_arg, self);
                    }
                }
            };
            /*
             *
             */
            auto print_element_neighbours = [&] <Element EE, auto I = 0, auto J = 0> (
                    auto const & elem_arg,
                    auto & self
            )
            mutable
            {
                for (auto const & c_ : elem_arg->template getNeighbours<I, J>()) {
                    if constexpr (!element::PointConcept<EE> && I == 0) {
                        os << "layer : " << I << " type : " << J << " <-> " << c_->hash() << std::endl;
                    }
                    else {
                        os << "layer : " << I << " type : " << J << " --> " << c_->hash() << std::endl;
                    }
                }
                if constexpr (J < element::numNeighbours<EE, _domain, I>() - 1) {
                    self.template operator()<EE, I, J + 1>(elem_arg, self);
                }
                else if constexpr (I < element::numNeighbours<EE, _domain>() - 1) {
                    self.template operator()<EE, I + 1, 0>(elem_arg, self);
                }
            };
            /*
             *
             */
            auto print_set = [&] <lolita::index I = 0, lolita::index J = 0> (
                    auto const & set_arg,
                    auto & self
            )
            mutable
            {
                for (auto const & item : set_arg.template getElements<I, J>()) {
                    os << "layer : " << I << " type : " << J << " <-- " << item.second->hash() << std::endl;
                }
                if constexpr (J < element::numElements<_domain, I>() - 1) {
                    self.template operator()<I, J + 1>(set_arg, self);
                }
                else if constexpr (I < element::numElements<_domain>() - 2) {
                    self.template operator()<I + 1, 0>(set_arg, self);
                }
            };
            /*
             *
             */
            auto print_sets = [&] ()
            mutable
            {
                for (auto const & set : base.element_sets_) {
                    os << "*** lolita::geometry::Domain : " << set.first << std::endl;
                    print_set(set.second, print_set);
                }
            };
            /*
             *
             */
            auto print_elements = [&] <lolita::index I = 0, lolita::index J = 0> (
                    auto & self
            )
            mutable
            {
                auto constexpr elt = element::element<_domain, I, J>();
                if constexpr (I == 0 && J == 0) {
                    os << "*** Elements : " << std::endl;
                }
                for (auto const & element : base.elements_.template getElements<I, J>()) {
                    os << "* Element : " << element.second->hash() << std::endl;
                    print_element_components.template operator()<elt>(element.second, print_element_components);
                    print_element_neighbours.template operator()<elt>(element.second, print_element_neighbours);
                }
                if constexpr (J < element::numElements<_domain, I>() - 1) {
                    self.template operator()<I, J + 1>(self);
                }
                else if constexpr (I < element::numElements<_domain>() - 2) {
                    self.template operator()<I + 1, 0>(self);
                }
            };
            /*
             *
             */
            print_elements(print_elements);
            print_sets();
            return os;
        }

    };

    template<lolita::geometry::Domain D, lolita::finite_element::FiniteElementConcept auto... _finite_element>
    struct MeshModule<lolita::mesh::MeshFormatType::Gmsh, D, _finite_element...>
    {

        auto const static constexpr element_table = std::array<Element, 10>{
            element::seg_02,
            element::tri_03,
            element::qua_04,
            element::pnt_00,
        };

        static
        Element
        getElementType(
                lolita::index
                tag_arg
        )
        {
            if (tag_arg == 15) {
                return element::pnt_00;
            }
            else if (tag_arg == 1) {
                return element::seg_02;
            }
            else if (tag_arg == 2) {
                return element::tri_03;
            }
            else if (tag_arg == 3) {
                return element::qua_04;
            }
            else {
                assert(false);
            }
        }

        struct Module : MeshModuleBase
        {

            explicit
            Module(
                    lolita::utility::File const &
                    mesh_file
            )
            :
            MeshModuleBase(mesh_file)
            {
                setGeometricalEntities();
                setPhysicalEntities();
                setPhysicalGroups();
            }

            lolita::boolean
            operator==(
                    Module const &
                    other
            )
            const = default;

            lolita::boolean
            operator!=(
                    Module const &
                    other
            )
            const = default;

            struct PhysicalEntity
            {

                lolita::boolean
                operator==(
                        PhysicalEntity const &
                        other
                )
                const = default;

                lolita::boolean
                operator!=(
                        PhysicalEntity const &
                        other
                )
                const = default;

                lolita::index tag_;

                lolita::index dim_;

                std::basic_string<lolita::character> name_;

            };

            struct GeometricalEntity
            {

                lolita::boolean
                operator==(
                        GeometricalEntity const &
                        other
                )
                const = default;

                lolita::boolean
                operator!=(
                        GeometricalEntity const &
                        other
                )
                const = default;

                lolita::index tag_;

                lolita::index dim_;

                std::vector<lolita::index> physical_entities_tags_;

                std::vector<lolita::index> bounding_entities_tags_;

            };

            struct PhysicalGroup
            {

                lolita::boolean
                operator==(
                        PhysicalGroup const &
                        other
                )
                const = default;

                lolita::boolean
                operator!=(
                        PhysicalGroup const &
                        other
                )
                const = default;

                std::basic_string<lolita::character> name_;

                std::array<std::vector<lolita::index>, 4> geometrical_entities_tags_;

            };

            using PhysicalGroups = std::vector<PhysicalGroup>;

            using GeometricalEntities = std::array<std::vector<GeometricalEntity>, 4>;

            using PhysicalEntities = std::vector<PhysicalEntity>;

        private:

            void
            setGeometricalEntities()
            {
//                lolita::index line_start = array::index(this->file_, "$Entities");
                lolita::index line_start = std::distance(this->file_.lines_.begin(), std::find(this->file_.lines_.begin(), this->file_.lines_.end(), "$Entities"));
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(this->file_.lines_[line_start + offset]);
                lolita::index num_points;
                lolita::index num_curves;
                lolita::index num_surfaces;
                lolita::index num_volumes;
                line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
                auto num_domains = std::array<lolita::index, 4>{num_points, num_curves, num_surfaces, num_volumes};
                offset += 1;
                for (lolita::index i = 0; i < 4; ++i) {
                    for (lolita::index j = 0; j < num_domains[i]; ++j) {
                        line_stream = std::basic_stringstream<lolita::character>(this->file_.lines_[line_start + offset]);
                        lolita::index tag;
                        line_stream >> tag;
                        if (i == 0) {
                            for (lolita::index k = 0; k < 3; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        } else {
                            for (lolita::index k = 0; k < 6; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        }
                        lolita::index num_physical_entities;
                        line_stream >> num_physical_entities;
                        std::vector<lolita::index> physical_entities_tags;
                        for (lolita::index k = 0; k < num_physical_entities; ++k) {
                            lolita::integer physical_entity_tag;
                            line_stream >> physical_entity_tag;
                            physical_entities_tags.push_back(physical_entity_tag);
                        }
                        std::vector<lolita::index> bounding_entities_tags = {};
                        if (i > 0) {
                            lolita::index num_bounding_entities;
                            line_stream >> num_bounding_entities;
                            for (lolita::index k = 0; k < num_bounding_entities; ++k) {
                                lolita::integer bounding_entity_tag;
                                line_stream >> bounding_entity_tag;
                                bounding_entities_tags.push_back(std::abs(bounding_entity_tag));
                            }
                        }
                        geometrical_entities[i].push_back(GeometricalEntity{
                                tag,
                                i,
                                physical_entities_tags,
                                bounding_entities_tags
                        });
                        offset += 1;
                    }
                }
            }

            void
            setPhysicalEntities()
            {
//                lolita::index line_start = array::index(this->file_.lines_, "$PhysicalNames");
                lolita::index line_start = std::distance(this->file_.lines_.begin(), std::find(this->file_.lines_.begin(), this->file_.lines_.end(), "$PhysicalNames"));
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(this->file_.lines_[line_start + offset]);
                lolita::index num_physical_names;
                line_stream >> num_physical_names;
                offset += 1;
                for (lolita::index i = 0; i < num_physical_names; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(this->file_.lines_[line_start + offset]);
                    lolita::index dim;
                    lolita::index tag;
                    std::basic_string<lolita::character> name;
                    line_stream >> dim >> tag >> name;
                    lolita::utility::removeCharacter(name, '"');
                    physical_entities.push_back(PhysicalEntity{tag, dim, name});
                    offset += 1;
                }
            }

            std::array<std::unordered_set<lolita::index>, 4>
            getSubGeometricalEntities(
                    lolita::index d,
                    lolita::index t,
                    std::array<std::unordered_set<lolita::index>, 4> & a
            )
            {
                a[d].insert(t);
                for (lolita::index i = 0; i < geometrical_entities[d][t - 1].bounding_entities_tags_.size(); ++i) {
                    lolita::index const d2 = d - 1;
                    lolita::index const t2 = geometrical_entities[d][t - 1].bounding_entities_tags_[i];
                    a = getSubGeometricalEntities(d2, t2, a);
                }
                return a;
            }

            void
            setPhysicalGroups()
            {
                auto get_subs = [&] (
                        auto d,
                        auto t,
                        std::array<std::unordered_set<lolita::index>, 4> & a,
                        auto & self
                )
                mutable
                {
                    a[d].insert(t);
                    for (lolita::index i = 0; i < geometrical_entities[d][t - 1].bounding_entities_tags.size(); ++i) {
                        a = self(d - 1, geometrical_entities[d][t - 1].bounding_entities_tags[i], a, self);
                    }
                    return a;
                };
                /*
                 *
                 */
                for (lolita::index i = 0; i < physical_entities.size(); ++i) {
                    std::array<std::unordered_set<lolita::index>, 4> tags;
                    for (lolita::index j = 0; j < 4; ++j) {
                        for (lolita::index k = 0; k < geometrical_entities[j].size(); ++k) {
                            for (lolita::index l = 0; l < geometrical_entities[j][k].physical_entities_tags_.size(); ++l) {
                                lolita::index const & t0 = geometrical_entities[j][k].physical_entities_tags_[l];
                                lolita::index const & t1 = physical_entities[i].tag_;
                                if (t0 == t1) {
//                                    tags = get_subs(j, k + 1, tags, get_subs);
                                    tags = getSubGeometricalEntities(j, k + 1, tags);
                                }
                            }
                        }
                    }
                    std::array<std::vector<lolita::index>, 4> group_tags;
                    for (lolita::index j = 0; j < 4; ++j) {
                        group_tags[j].assign(tags[j].begin(), tags[j].end());
                    }
                    physical_groups.push_back(PhysicalGroup{physical_entities[i].name_, group_tags});
                }
                //return groups;
            }

        public:

            GeometricalEntities geometrical_entities;
            PhysicalEntities physical_entities;
            PhysicalGroups physical_groups;

        };

        struct Implementation : public MeshBase<lolita::mesh::MeshFormatType::Gmsh, D, _finite_element...>
        {

            //using Base = MeshBase<lolita::mesh::MeshFormatType::Gmsh, D, Med>;

            void
            makeNodes()
            {
//                auto & nodes = this->elements_.template get<0>().template get<0>();
//                auto & nodes = std::get<0>(std::get<0>(this->elements_));
                auto const & file_content = this->module_.file_.lines_;
                auto & nodes = this->elements_.template getElements<0, 0>();
                lolita::index line_start = std::distance(file_content.begin(), std::find(file_content.begin(), file_content.end(), "$Nodes"));
//                lolita::index line_start = array::index(this->module_.file_, "$Nodes");
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_content[line_start + offset]);
                lolita::index num_entity_block;
                lolita::index num_nodes;
                lolita::index min_node_tag;
                lolita::index max_node_tag;
                line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                offset += 1;
                for (lolita::index i = 0; i < num_entity_block; ++i) {
//                    line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                    line_stream = std::basic_stringstream<lolita::character>(file_content[line_start + offset]);
                    lolita::index entity_dim;
                    lolita::index entity_tag;
                    lolita::index parametric;
                    lolita::index num_nodes_in_block;
                    line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                    offset += 1;
                    for (lolita::index j = 0; j < num_nodes_in_block; ++j) {
//                        line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                        line_stream = std::basic_stringstream<lolita::character>(file_content[line_start + offset]);
                        lolita::natural node_tag;
                        line_stream >> node_tag;
//                        line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset + num_nodes_in_block));
                        line_stream = std::basic_stringstream<lolita::character>(file_content[line_start + offset + num_nodes_in_block]);
                        lolita::geometry::Point coordinates = lolita::geometry::Point::Zero();
                        for (int k = 0; k < D.dim_; ++k) {
                            line_stream >> coordinates(k);
                        }
//                        this->setElement(node_tag, coordinates);
//                        this->makeNode(node_tag, coordinates);
                        this->makeNode(node_tag, coordinates);
                        offset += 1;
                    }
                    offset += num_nodes_in_block;
                }
            }

            void
            makeElementSets()
            {
                auto const & physical_groups = this->module_.physical_groups;
                for (int i = 0; i < physical_groups.size(); ++i) {
                    this->element_sets_[physical_groups[i].name_];
                }
            }

            void
            makeNodeSets()
            {
                auto const & file_content = this->module_.file_.lines_;
                auto const & physical_groups = this->module_.physical_groups;
                for (int m = 0; m < physical_groups.size(); ++m) {
                    auto & node_set_name = physical_groups[m].name_;
//                    auto & node_set = std::get<0>(std::get<0>(this->element_sets_[node_set_name]));
                    auto & node_set = this->element_sets_[node_set_name].template getElements<0, 0>();
//                    auto & nodes = std::get<0>(std::get<0>(this->elements_));
                    auto & nodes = this->elements_.template getElements<0, 0>();
                    for (lolita::index i = 0; i < 4; ++i) {
                        for (lolita::index j = 0; j < physical_groups[m].geometrical_entities_tags_[i].size(); ++j) {
                            lolita::index tag = physical_groups[m].geometrical_entities_tags_[i][j];
//                            lolita::index line_start = array::index(this->module_.file_, "$Nodes");
                            lolita::index line_start = std::distance(file_content.begin(), std::find(file_content.begin(), file_content.end(), "$Nodes"));
                            lolita::index offset = 1;
//                            std::basic_string<lolita::character>Stream line_stream(this->module_.file_.get(line_start + offset));
                            std::basic_stringstream<lolita::character> line_stream(file_content[line_start + offset]);
                            lolita::index num_entity_block;
                            lolita::index num_nodes;
                            lolita::index min_node_tag;
                            lolita::index max_node_tag;
                            line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                            offset += 1;
                            for (lolita::index k = 0; k < num_entity_block; ++k) {
//                                line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                                line_stream = std::basic_stringstream<lolita::character>(file_content[line_start + offset]);
                                lolita::index entity_dim;
                                lolita::index entity_tag;
                                lolita::index parametric;
                                lolita::index num_nodes_in_block;
                                line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                                offset += 1;
                                if (entity_dim == i && entity_tag == tag) {
                                    for (lolita::index l = 0; l < num_nodes_in_block; ++l) {
//                                        line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                                        line_stream = std::basic_stringstream<lolita::character>(file_content[line_start + offset]);
                                        std::basic_string<lolita::character> node_hash = line_stream.str();
                                        node_set.insert({node_hash, nodes[node_hash]});
                                        offset += 1;
                                    }
                                } else {
                                    offset += num_nodes_in_block;
                                }
                                offset += num_nodes_in_block;
                            }
                        }
                    }
                }
            }

            template<Element _cell>
            void
            makeCells()
            requires(!element::PointConcept <_cell>)
            {
//                lolita::index line_start = array::index(this->module_.file_, "$Elements");
                auto const & file_content = this->module_.file_.lines_;
                lolita::index line_start = std::distance(file_content.begin(), std::find(file_content.begin(), file_content.end(), "$Elements"));
                lolita::index offset = 1;
//                std::basic_string<lolita::character>Stream line_stream(this->module_.file_.get(line_start + offset));
                std::basic_stringstream<lolita::character> line_stream(file_content[line_start + offset]);
                lolita::index num_entity_blocks;
                lolita::index num_elements;
                lolita::index min_element_tag;
                lolita::index max_element_tag;
                line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
                offset += 1;
                for (lolita::index i = 0; i < num_entity_blocks; ++i) {
//                    line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                    line_stream = std::basic_stringstream<lolita::character>(file_content[line_start + offset]);
                    lolita::index entity_dim = 0;
                    lolita::index entity_tag = 0;
                    lolita::index element_type_tag = 0;
                    lolita::index num_elements_in_block = 0;
                    line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                    offset += 1;
                    Element const element_type_arg = getElementType(element_type_tag);
                    if (entity_dim == D.dim_ and getElementType(element_type_tag) == _cell) {
                        for (lolita::index j = 0; j < num_elements_in_block; ++j) {
//                            line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                            line_stream = std::basic_stringstream<lolita::character>(file_content[line_start + offset]);
                            lolita::index tag;
                            line_stream >> tag;
                            std::array<lolita::index, _cell.num_nodes_> node_tags;
                            for (lolita::index k = 0; k < _cell.num_nodes_; ++k) {
                                line_stream >> node_tags[k];
                            }
//                            this->template makeElement<EE>(node_tags);
                            this->template makeCell<_cell>(node_tags);
                            offset += 1;
                        }
                    } else {
                        offset += num_elements_in_block;
                    }
                }
            }

        };

    };

//        template<Element EE>
//        void
//        makeElement(
//                std::ranges::range auto const &
//                cell_node_tags_arg
//        )
//        {
//            /*
//             *
//             */
//            auto get_element_domains = [&] (
//                    std::ranges::range auto const &
//                    node_tags_arg
//            )
//            {
//                std::vector<std::basic_string<lolita::character>> domain_names;
//                for (auto const & [set_name, element_tags] : element_sets_) {
//                    auto const & node_hashes = std::get<0>(std::get<0>(element_tags));
//                    lolita::boolean belongs = true;
//                    for (auto const & node_tag: node_tags_arg) {
//                        auto node_hash = std::to_string(node_tag);
//                        if (!node_hashes.contains(node_hash)) {
//                            belongs = false;
//                            break;
//                        }
//                    }
//                    if (belongs) {
//                        domain_names.push_back(set_name);
//                    }
//                }
//                return domain_names;
//            };
//            /*
//             *
//             */
//            auto get_component_node_tags = [] <Element E, auto I, auto J> (
//                    std::ranges::range auto const &
//                    node_tags_arg,
//                    lolita::index
//                    index_arg
//            )
//            {
//                auto component_node_tags = std::array<lolita::index, element::neighbour<E, 0, I, J>().num_nodes>();
//                for (lolita::index j = 0; j < element::neighbour<E, 0, I, J>().num_nodes; ++j) {
//                    auto const & element_node_connectivity = FiniteElementE<E, D, Med>::node_connectivity;
//                    auto const & bound_node_connectivity = std::get<J>(std::get<I>(element_node_connectivity));//.template get<I>().template get<J>();
////                    auto const k = bound_node_connectivity.get(index_arg).get(j);
//                    auto const k = bound_node_connectivity[index_arg][j];
//                    component_node_tags[j] = node_tags_arg[k];
//                }
//                return component_node_tags;
//            };
//            /*
//             *
//             */
//            auto get_element_hash = [] (
//                    auto
//                    node_tags_arg
//            )
//            {
//                std::basic_stringstream<lolita::character> element_hash;
//                if constexpr (node_tags_arg.size() > 0) {
//                    std::sort(std::execution::par_unseq, node_tags_arg.begin(), node_tags_arg.end());
//                }
//                for (lolita::index i = 0; i < node_tags_arg.size(); ++i) {
//                    element_hash << node_tags_arg[i];
//                }
//                return element_hash.str();
//            };
//            /*
//             *
//             */
//            auto set_finite_elements = [&] <auto K = 0> (
//                    auto &
//                    ptr_element_arg,
//                    auto &
//                    self
//            )
//            mutable
//            {
////                pointer::make(ptr_element_arg.get().template get<K>());
//                if constexpr (K < Med.size() - 1) {
//                    self.template operator()<K + 1>(ptr_element_arg, self);
//                }
//            };
//            /*
//             *
//             */
//            auto set_finite_elements_components = [&] <Element E, auto K = 0, auto I = 0, auto J = 0> (
//                    auto &
//                    ptr_element_arg,
//                    auto &
//                    set_components_imp
//            )
//            mutable
//            {
//                auto & elt = ptr_element_arg.get();
//                for (int i = 0; i < element::numNeighbours<E, 0, I, J>(); ++i) {
//                    auto & rhs = elt.components.template get<I>().template get<J>().get(i);
//                    auto & lhs = elt.template get<K>().get().components.template get<I>().template get<J>().get(i);
//                    lhs.ptr = rhs.ptr.get().template get<K>();
////                    lhs.orientation = rhs.orientation;
//                    lhs.index = rhs.index;
//                }
//                if constexpr (J < element::numNeighbours<E, 0, I>() - 1) {
//                    set_components_imp.template operator()<E, K, I, J + 1>(ptr_element_arg, set_components_imp);
//                }
//                else if constexpr (I < element::numNeighbours<E, 0>() - 1) {
//                    set_components_imp.template operator()<E, K, I + 1, 0>(ptr_element_arg, set_components_imp);
//                }
//                else if constexpr (K < Med.size() - 1) {
//                    set_components_imp.template operator()<E, K + 1, 0, 0>(ptr_element_arg, set_components_imp);
//                }
//            };
//            /*
//             *
//             */
//            auto set_finite_element_unknowns = [&] <auto K = 0> (
//                    auto &
//                    ptr_element_arg,
//                    auto &
//                    self
//            )
//            mutable
//            {
//                auto & elt = ptr_element_arg.get();
//                elt.template get<K>().get().setUnknowns(unknown_indices[K].unknown_index);
//                if constexpr (K < Med.size() - 1) {
//                    self.template operator()<K + 1>(ptr_element_arg, self);
//                }
//            };
//            /*
//             *
//             */
//            auto set_finite_element_materials = [&] <auto K = 0> (
//                    auto &
//                    ptr_element_arg,
//                    auto const &
//                    element_domains_arg,
//                    auto &
//                    self
//            )
//            mutable
//            {
//                auto & elt = ptr_element_arg.get().template get<K>().get();
//                for (auto const & bhv : behaviours_)
//                {
////                    auto same_finite_element = Med.template get<K>().field().tag == bhv.finite_element_tag;
////                    if (same_finite_element) {
////                        auto found_domain = false;
////                        for (auto const & element_domain : element_domains_arg.data) {
////                            for (auto const & behaviour_domain : bhv.domains.data) {
////                                if (element_domain == behaviour_domain) {
////                                    found_domain = true;
////                                    print("setting material :", bhv.finite_element_tag, behaviour_domain);
////                                    elt.setMaterial(bhv.ptr_behaviour.get());
////                                    break;
////                                }
////                            }
////                        }
////                    }
//                }
//                if constexpr (K < Med.size() - 1) {
//                    self.template operator()<K + 1>(ptr_element_arg, element_domains_arg, self);
//                }
//            };
//            /*
//             *
//             */
//            auto set_finite_element_loads = [&] <auto K = 0> (
//                    auto &
//                    ptr_element_arg,
//                    auto const &
//                    element_domains_arg,
//                    auto &
//                    self
//            )
//            mutable
//            {
//                auto constexpr field = Field(Med.template get<K>().field().ord_field, D.dim_);
//                auto & elt = ptr_element_arg.get().template get<K>().get();
//                for (int i = 0; i < field.rows(); i++)
//                {
//                    for (int j = 0; j < field.cols(); j++)
//                    {
//                        auto found_load = false;
//                        for (auto const & load : loads_)
//                        {
//                            auto same_finite_element = Med.template get<K>().field().tag == load.unknown_tag_;
//                            auto same_direction = load.components_.row_ == i && load.components_.col_ == j;
//                            if (same_direction && same_finite_element)
//                            {
//                                for (auto const & element_domain : element_domains_arg.data)
//                                {
//                                    if (element_domain == load.domain_tag_)
//                                    {
////                                        print("Adding load ! : ", i, j);
//                                        elt.setLoad(unknown_indices[K].binding_index, load.load_, i, j);
//                                        found_load = true;
//                                        break;
//                                    }
//                                }
//                            }
//                        }
//                        if (!found_load)
//                        {
////                            print("Zero load ! : ", i, j);
//                            elt.setLoad(unknown_indices[K].binding_index, null_load_, i, j);
//                        }
//                    }
//                }
//                if constexpr (K < Med.size() - 1) {
//                    self.template operator()<K + 1>(ptr_element_arg, element_domains_arg, self);
//                }
//            };
//            /*
//             *
//             */
//            auto initialize_finite_elements = [&] <auto K = 0> (
//                    auto &
//                    ptr_element_arg,
//                    auto &
//                    self
//            )
//            mutable
//            {
//                ptr_element_arg.get().template get<K>().get().initialize();
//                if constexpr (K < Med.size() - 1) {
//                    self.template operator()<K + 1>(ptr_element_arg, self);
//                }
//            };
//            /*
//             *
//             */
//            auto set_element = [&] <Element E, auto I = 0, auto J = 0> (
//                    std::array<lolita::index, E.num_nodes> const &
//                    element_node_tags_arg,
//                    std::shared_ptr<FiniteElementE<E, D, Med>> &
//                    ptr_element_arg,
//                    auto &
//                    set_element_imp
//            )
//            mutable
//            {
//                auto const constexpr cmp = element::neighbour<E, 0, I, J>();
//                auto const constexpr bd = cmp.dim;
//                auto const constexpr bt = element::elementIndex<cmp>();
//                auto const constexpr ed = E.dim;
//                auto const constexpr et = element::elementIndex<E>();
////                auto & components = elements_.template get<cmp.dim>().template get<bt>();
////                auto & elements = elements_.template get<E.dim>().template get<et>();
//                auto & components = std::get<bt>(std::get<cmp.dim>(elements_));//.template get<cmp.dim>().template get<bt>();
//                auto & elements = std::get<et>(std::get<E.dim>(elements_));//elements_.template get<E.dim>().template get<et>();
////                auto & element_component_array = ptr_element_arg.get().components.template get<I>().template get<J>();
//                auto & element_component_array = ptr_element_arg->template getNeighbours<0, I, J>();
//                auto element_hash = get_element_hash(element_node_tags_arg);
//                for (auto i = 0; i < element::numNeighbours<E, 0, I, J>(); ++i) {
//                    auto component_node_tags = get_component_node_tags.template operator ()<E, I, J>(element_node_tags_arg, i);
//                    auto component_hash = get_element_hash(component_node_tags);
//                    //
////                    if (!components.data.contains(component_hash)) {
////                        auto ptr_component = std::shared_ptr<FiniteElementE<cmp, D, Med>>(FiniteElementE<cmp, D, Med>());
////                        set_element_imp.template operator ()<cmp>(component_node_tags, ptr_component, set_element_imp);
////                    }
////                    element_component_array.get(i) = {components.get(component_hash), components.get(component_hash).get().count};
////                    components.get(component_hash).get().count ++;
////                    auto & component_neighbours = components.get(component_hash).get().neighbours;
////                    if constexpr (cmp == element::pnt_00) {
////                        lolita::index const constexpr nd = ed - 1;
////                        component_neighbours.template get<nd>().template get<et>().data.push_back({ptr_element_arg});
////                    }
////                    else {
////                        lolita::index const constexpr nd = ed - bd;
////                        component_neighbours.template get<nd>().template get<et>().data.push_back({ptr_element_arg});
////                    }
//                    //
//                    if constexpr (cmp == element::pnt_00) {
//                        //element_component_array.get(i) = {components.get(component_hash), 1};
////                        element_component_array.get(i) = {components.get(component_hash), components.get(component_hash).get().count};
////                        element_component_array[i] = {components.get(component_hash), components.get(component_hash).get().count};
//                        element_component_array[i] = components[component_hash];
////                        components.get(component_hash).get().count ++;
//                    }
//                    else {
//                        if (components.contains(component_hash)) {
//                            //element_component_array.get(i) = {components.get(component_hash), -1};
////                            element_component_array.get(i) = {components.get(component_hash), components.get(component_hash).get().count};
////                            components.get(component_hash).get().count ++;
//                            element_component_array[i] = components[component_hash];
//                        }
//                        else {
//                            auto ptr_component = std::make_shared<FiniteElementE<cmp, D, Med>>(FiniteElementE<cmp, D, Med>());
//                            set_element_imp.template operator ()<cmp>(component_node_tags, ptr_component, set_element_imp);
//                            //element_component_array.get(i) = {components.get(component_hash), 1};
////                            element_component_array.get(i) = {components.get(component_hash), components.get(component_hash).get().count};
////                            components.get(component_hash).get().count ++;
//                            element_component_array[i] = components[component_hash];
//                        }
//                    }
////                    auto & component_neighbours = components[component_hash]->neighbours;
//                    if constexpr (cmp == element::pnt_00) {
//                        lolita::index const constexpr nd = ed - 1;
////                        component_neighbours.template get<nd>().template get<et>().data.push_back({ptr_element_arg, component_neighbours.template get<nd>().template get<et>().size()});
////                        std::get<et>(std::get<nd>(component_neighbours)).push_back(ptr_element_arg);
//                        components[component_hash]->template getNeighbours<1, nd, et>().push_back(ptr_element_arg);
//                    }
//                    else {
//                        lolita::index const constexpr nd = ed - bd;
////                        component_neighbours.template get<nd>().template get<et>().data.push_back({ptr_element_arg, component_neighbours.template get<nd>().template get<et>().size()});
////                        std::get<et>(std::get<nd>(component_neighbours)).push_back(ptr_element_arg);
//                        components[component_hash]->template getNeighbours<1, nd, et>().push_back(ptr_element_arg);
//                    }
//                }
//                if constexpr (J < element::numNeighbours<E, 0, I>() - 1) {
//                    set_element_imp.template operator()<E, I, J + 1>(element_node_tags_arg, ptr_element_arg, set_element_imp);
//                }
//                else if constexpr (I < element::numNeighbours<E, 0>() - 1) {
//                    set_element_imp.template operator()<E, I + 1, 0>(element_node_tags_arg, ptr_element_arg, set_element_imp);
//                }
//                if constexpr (I == 0 && J == 0) {
////                    ptr_element_arg.get().tag = elements.size();
//                    ptr_element_arg->tag = elements.size();
//                    auto element_domains = get_element_domains(element_node_tags_arg);
////                    for (auto const & element_domain : element_domains.data) {
////                        auto & element_set = element_sets_.get(element_domain).template get<E.dim>().template get<et>();
////                        element_set.data.insert({element_hash, ptr_element_arg});
////                    }
////                    set_finite_elements(ptr_element_arg, set_finite_elements);
////                    set_finite_elements_components.template operator ()<E>(ptr_element_arg, set_finite_elements_components);
////                    set_finite_element_unknowns(ptr_element_arg, set_finite_element_unknowns);
////                    set_finite_element_loads(ptr_element_arg, element_domains, set_finite_element_loads);
////                    set_finite_element_materials(ptr_element_arg, element_domains, set_finite_element_materials);
////                    initialize_finite_elements(ptr_element_arg, initialize_finite_elements);
//                    elements.insert({element_hash, ptr_element_arg});
//                }
//            };
//            /*
//             *
//             */
//            auto set_element_n = [&] <Element E, auto I = 0, auto J = 0> (
//                    std::array<lolita::index, E.num_nodes> const &
//                    element_node_tags_arg,
//                    std::shared_ptr<FiniteElementE<E, D, Med>> &
//                    ptr_element_arg,
//                    auto &
//                    self
//            )
//            mutable
//            {
//                auto const constexpr cmp = element::neighbour<E, 0, I, J>();
//                auto const constexpr bd = cmp.dim;
//                auto const constexpr bt = element::elementIndex<cmp>();
//                auto const constexpr ed = E.dim;
//                auto const constexpr et = element::elementIndex<E>();
//                auto & components = std::get<bt>(std::get<cmp.dim>(elements_));//.template get<cmp.dim>().template get<bt>();
//                auto & elements = std::get<et>(std::get<E.dim>(elements_));//elements_.template get<E.dim>().template get<et>();
//                auto & element_component_array = ptr_element_arg->template getNeighbours<0, I, J>();
//                auto element_hash = get_element_hash(element_node_tags_arg);
//                for (auto i = 0; i < element::numNeighbours<E, 0, I, J>(); ++i) {
//                    auto component_node_tags = get_component_node_tags.template operator ()<E, I, J>(element_node_tags_arg, i);
//                    auto component_hash = get_element_hash(component_node_tags);
//                    if (!components.contains(component_hash)) {
//                        auto ptr_component = std::make_shared<FiniteElementE<cmp, D, Med>>(FiniteElementE<cmp, D, Med>());
//                        self.template operator ()<cmp>(component_node_tags, ptr_component, self);
//                    }
//                    element_component_array[i] = components[component_hash];
//                    auto const constexpr nhj = element::PointConcept<cmp> ? ed - 1 : ed - bd;
//                    components[component_hash]->template getNeighbours<1, nhj, et>().push_back(ptr_element_arg);
//                }
//                if constexpr (J < element::numNeighbours<E, 0, I>() - 1) {
//                    self.template operator()<E, I, J + 1>(element_node_tags_arg, ptr_element_arg, self);
//                }
//                else if constexpr (I < element::numNeighbours<E, 0>() - 1) {
//                    self.template operator()<E, I + 1, 0>(element_node_tags_arg, ptr_element_arg, self);
//                }
//                if constexpr (I == 0 && J == 0) {
//                    ptr_element_arg->tag = elements.size();
//                    auto element_domains = get_element_domains(element_node_tags_arg);
//                    elements.insert({element_hash, ptr_element_arg});
//                }
//            };
//            /*
//             *
//             */
//            auto ptr_cell = std::make_shared<FiniteElementE<EE, D, Med>>(FiniteElementE<EE, D, Med>());
////            set_element.template operator ()<EE>(cell_node_tags_arg, ptr_cell, set_element);
//            set_element_n.template operator ()<EE>(cell_node_tags_arg, ptr_cell, set_element_n);
//        }

}

#endif //LOLITA_LOLITA_CORE_MESH_HXX
