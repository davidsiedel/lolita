//
// Created by dsiedel on 17/04/22.
//

#include <execution>
#include <ostream>

#include "lolita/lolita_core_finite_element.hxx"
#include "lolita/lolita_file.hxx"

#ifndef LOLITA_LOLITA_CORE_MESH_HXX
#define LOLITA_LOLITA_CORE_MESH_HXX

namespace lolita::core::mesh
{

    using Element = core::element::Element;

    //using DegreeOfFreedomIndex = core::element::DegreeOfFreedomIndex;

    template<Element E, Domain D, auto Med>
    using FiniteElementE = core::element::FiniteElement<E, D, Med>;

    template<MeshFormatType, Domain D, auto Med>
    struct MeshModule;

    struct MeshModuleBase
    {

        MeshModuleBase(
                file::File const &
                mesh_file
        )
        :
        file_(mesh_file)
        {}

        file::File const & file_;

    };

    template<MeshFormatType Mft, Domain D, auto Med>
    struct MeshBase
    {

    private:

        template<Element E>
        using ElementPointerMap = UnorderedMap<Strg, SharedPointer<FiniteElementE<E, D, Med>>>;

        using Module = typename MeshModule<Mft, D, Med>::Module;

        using Implementation = typename MeshModule<Mft, D, Med>::Implementation;

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
                    auto
                    unknown_index_arg,
                    auto
                    binding_index_arg
            )
            :
            unknown_index(unknown_index_arg),
            binding_index(binding_index_arg)
            {}

            Long unknown_index;

            Long binding_index;

        };

        using DegreeOfFreedomIndices = Array<DegreeOfFreedomIndex, Med.size()>;

        using Loads = Array<Load<D>>;

        using Behaviours = Array<Behaviour>;

        SharedPointer<LoadComponent<D>> const null_load_;

    public:

        using ElementPointers = element::Elements<D.dim, ElementPointerMap>;

        using ElementPointersSets = UnorderedMap<Strg, ElementPointers>;

        DegreeOfFreedomIndices unknown_indices;

        Module module_;

        ElementPointers elements_;

        ElementPointersSets element_sets_;

        Loads loads_;

        Behaviours behaviours_;

        MeshBase(
                file::File const &
                mesh_file,
                Loads const &
                loads_arg,
                Behaviours const &
                behaviours_arg
        )
        :
        module_(mesh_file),
        null_load_(SharedPointer<LoadComponent<D>>(LoadComponent<D>())),
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

        void
        makeNode(
                auto
                node_tag_arg,
                auto const &
                node_coordinates_arg
        )
        {
            auto hash = std::to_string(node_tag_arg);
            auto ptr_node = SharedPointer<FiniteElementE<element::pnt_00, D, Med>>(FiniteElementE<element::pnt_00, D, Med>());
            auto ptr_coordinates = SharedPointer<Vector<Real, D.dim>>(node_coordinates_arg);
            /*
             *
             */
            auto make_node = [&] <auto K = 0> (auto & self) {
                pointer::make(ptr_node.get().template get<K>());
                ptr_node.get().template get<K>().get().coordinates = ptr_coordinates;
                if constexpr (K < FiniteElementE<element::pnt_00, D, Med>::size() - 1) {
                    self.template operator()<K + 1>(self);
                }
            };
            /*
             *
             */
            make_node(make_node);
            auto & nodes = elements_.template get<0>().template get<0>();
            ptr_node.get().tag = nodes.size();
            ptr_node.get().coordinates = ptr_coordinates;
            nodes.data.insert({hash, ptr_node});
        }

        template<Element EE>
        void
        makeElement(
                auto const &
                cell_node_tags_arg
        )
        {
            /*
             *
             */
            auto get_element_domains = [&] (
                    auto const &
                    node_tags_arg
            )
            {
                Array<Strg> domain_names;
                for (auto const & [set_name, element_tags] : element_sets_.data) {
                    auto const & node_hashes = element_tags.template get<0>().template get<0>();
                    Bool belongs = true;
                    for (auto const & node_tag: node_tags_arg.data) {
                        auto node_hash = std::to_string(node_tag);
                        if (!node_hashes.data.contains(node_hash)) {
                            belongs = false;
                            break;
                        }
                    }
                    if (belongs) {
                        domain_names.data.push_back(set_name);
                    }
                }
                return domain_names;
            };
            /*
             *
             */
            auto get_component_node_tags = [] <Element E, auto I, auto J> (
                    auto const &
                    node_tags_arg,
                    auto
                    index_arg
            )
            {
                auto component_node_tags = Array<Indx, element::component<E, I, J>().num_nodes>();
                for (Indx j = 0; j < element::component<E, I, J>().num_nodes; ++j) {
                    auto const & element_node_connectivity = FiniteElementE<E, D, Med>::node_connectivity;
                    auto const & bound_node_connectivity = element_node_connectivity.template get<I>().template get<J>();
                    auto const k = bound_node_connectivity.get(index_arg).get(j);
                    component_node_tags.get(j) = node_tags_arg.get(k);
                }
                return component_node_tags;
            };
            /*
             *
             */
            auto get_element_hash = [] (
                    auto
                    node_tags_arg
            )
            {
                StrgStream element_hash;
                if constexpr (node_tags_arg.size() > 0) {
                    std::sort(std::execution::par_unseq, node_tags_arg.data.begin(), node_tags_arg.data.end());
                }
                for (Indx i = 0; i < node_tags_arg.size(); ++i) {
                    element_hash << node_tags_arg.get(i);
                }
                return element_hash.str();
            };
            /*
             *
             */
            auto set_finite_elements = [&] <auto K = 0> (
                    auto &
                    ptr_element_arg,
                    auto &
                    self
            )
            mutable
            {
                pointer::make(ptr_element_arg.get().template get<K>());
                if constexpr (K < Med.size() - 1) {
                    self.template operator()<K + 1>(ptr_element_arg, self);
                }
            };
            /*
             *
             */
            auto set_finite_elements_components = [&] <Element E, auto K = 0, auto I = 0, auto J = 0> (
                    auto &
                    ptr_element_arg,
                    auto &
                    set_components_imp
            )
            mutable
            {
                auto & elt = ptr_element_arg.get();
                for (int i = 0; i < element::numComponents<E, I, J>(); ++i) {
                    auto & rhs = elt.components.template get<I>().template get<J>().get(i);
                    auto & lhs = elt.template get<K>().get().components.template get<I>().template get<J>().get(i);
                    lhs.ptr = rhs.ptr.get().template get<K>();
                    lhs.orientation = rhs.orientation;
                }
                if constexpr (J < element::numComponents<E, I>() - 1) {
                    set_components_imp.template operator()<E, K, I, J + 1>(ptr_element_arg, set_components_imp);
                }
                else if constexpr (I < element::numComponents<E>() - 1) {
                    set_components_imp.template operator()<E, K, I + 1, 0>(ptr_element_arg, set_components_imp);
                }
                else if constexpr (K < Med.size() - 1) {
                    set_components_imp.template operator()<E, K + 1, 0, 0>(ptr_element_arg, set_components_imp);
                }
            };
            /*
             *
             */
            auto set_finite_element_unknowns = [&] <auto K = 0> (
                    auto &
                    ptr_element_arg,
                    auto &
                    self
            )
            mutable
            {
                auto & elt = ptr_element_arg.get();
                elt.template get<K>().get().setUnknowns(unknown_indices.get(K).unknown_index);
                if constexpr (K < Med.size() - 1) {
                    self.template operator()<K + 1>(ptr_element_arg, self);
                }
            };
            /*
             *
             */
            auto set_finite_element_materials = [&] <auto K = 0> (
                    auto &
                    ptr_element_arg,
                    auto const &
                    element_domains_arg,
                    auto &
                    self
            )
            mutable
            {
                auto & elt = ptr_element_arg.get().template get<K>().get();
                for (auto const & bhv : behaviours_.data)
                {
                    auto same_finite_element = Med.template get<K>().field().tag == bhv.finite_element_tag;
                    if (same_finite_element) {
                        auto found_domain = false;
                        for (auto const & element_domain : element_domains_arg.data) {
                            for (auto const & behaviour_domain : bhv.domains.data) {
                                if (element_domain == behaviour_domain) {
                                    found_domain = true;
                                    print("setting material :", bhv.finite_element_tag, behaviour_domain);
                                    elt.setMaterial(bhv.ptr_behaviour.get());
                                    break;
                                }
                            }
                        }
                    }
                }
                if constexpr (K < Med.size() - 1) {
                    self.template operator()<K + 1>(ptr_element_arg, element_domains_arg, self);
                }
            };
            /*
             *
             */
            auto set_finite_element_loads = [&] <auto K = 0> (
                    auto &
                    ptr_element_arg,
                    auto const &
                    element_domains_arg,
                    auto &
                    self
            )
            mutable
            {
                auto constexpr field = Field(Med.template get<K>().field().ord_field, D.dim);
                auto & elt = ptr_element_arg.get().template get<K>().get();
                for (int i = 0; i < field.rows(); i++)
                {
                    for (int j = 0; j < field.cols(); j++)
                    {
                        auto found_load = false;
                        for (auto const & load : loads_.data)
                        {
                            auto same_finite_element = Med.template get<K>().field().tag == load.finite_element_label;
                            auto same_direction = load.components.i == i && load.components.j == j;
                            if (same_direction && same_finite_element)
                            {
                                for (auto const & element_domain : element_domains_arg.data)
                                {
                                    if (element_domain == load.domain_label)
                                    {
                                        print("Adding load ! : ", i, j);
                                        elt.setLoad(unknown_indices.get(K).binding_index, load.load, i, j);
                                        found_load = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if (!found_load)
                        {
                            print("Zero load ! : ", i, j);
                            elt.setLoad(unknown_indices.get(K).binding_index, null_load_, i, j);
                        }
                    }
                }
                if constexpr (K < Med.size() - 1) {
                    self.template operator()<K + 1>(ptr_element_arg, element_domains_arg, self);
                }
            };
            /*
             *
             */
            auto initialize_finite_elements = [&] <auto K = 0> (
                    auto &
                    ptr_element_arg,
                    auto &
                    self
            )
            mutable
            {
                ptr_element_arg.get().template get<K>().get().initialize();
                if constexpr (K < Med.size() - 1) {
                    self.template operator()<K + 1>(ptr_element_arg, self);
                }
            };
            /*
             *
             */
            auto set_element = [&] <Element E, auto I = 0, auto J = 0> (
                    auto &&
                    element_node_tags_arg,
                    auto &
                    ptr_element_arg,
                    auto &
                    set_element_imp
            )
            mutable
            {
                auto const constexpr cmp = element::component<E, I, J>();
                auto const constexpr bd = cmp.dim;
                auto const constexpr bt = element::elementIndex<cmp>();
                auto const constexpr ed = E.dim;
                auto const constexpr et = element::elementIndex<E>();
                auto & components = elements_.template get<cmp.dim>().template get<bt>();
                auto & elements = elements_.template get<E.dim>().template get<et>();
                auto & element_component_array = ptr_element_arg.get().components.template get<I>().template get<J>();
                auto element_hash = get_element_hash(element_node_tags_arg);
                for (auto i = 0; i < element::numComponents<E, I, J>(); ++i) {
                    auto component_node_tags = get_component_node_tags.template operator ()<E, I, J>(element_node_tags_arg, i);
                    auto component_hash = get_element_hash(component_node_tags);
                    if constexpr (cmp == element::pnt_00) {
                        element_component_array.get(i) = {components.get(component_hash), 1};
                    }
                    else {
                        if (components.data.contains(component_hash)) {
                            element_component_array.get(i) = {components.get(component_hash), -1};
                        }
                        else {
                            auto ptr_component = SharedPointer<FiniteElementE<cmp, D, Med>>(FiniteElementE<cmp, D, Med>());
                            set_element_imp.template operator ()<cmp>(component_node_tags, ptr_component, set_element_imp);
                            element_component_array.get(i) = {components.get(component_hash), 1};
                        }
                    }
                    auto & component_neighbours = components.get(component_hash).get().neighbours;
                    if constexpr (cmp == element::pnt_00) {
                        Indx const constexpr nd = ed - 1;
                        component_neighbours.template get<nd>().template get<et>().data.push_back({ptr_element_arg});
                    }
                    else {
                        Indx const constexpr nd = ed - bd;
                        component_neighbours.template get<nd>().template get<et>().data.push_back({ptr_element_arg});
                    }
                }
                if constexpr (J < element::numComponents<E, I>() - 1) {
                    set_element_imp.template operator()<E, I, J + 1>(element_node_tags_arg, ptr_element_arg, set_element_imp);
                }
                else if constexpr (I < element::numComponents<E>() - 1) {
                    set_element_imp.template operator()<E, I + 1, 0>(element_node_tags_arg, ptr_element_arg, set_element_imp);
                }
                if constexpr (I == 0 && J == 0) {
                    ptr_element_arg.get().tag = elements.size();
                    auto element_domains = get_element_domains(element_node_tags_arg);
                    for (auto const & element_domain : element_domains.data) {
                        auto & element_set = element_sets_.get(element_domain).template get<E.dim>().template get<et>();
                        element_set.data.insert({element_hash, ptr_element_arg});
                    }
                    set_finite_elements(ptr_element_arg, set_finite_elements);
                    set_finite_elements_components.template operator ()<E>(ptr_element_arg, set_finite_elements_components);
                    set_finite_element_unknowns(ptr_element_arg, set_finite_element_unknowns);
                    set_finite_element_loads(ptr_element_arg, element_domains, set_finite_element_loads);
                    set_finite_element_materials(ptr_element_arg, element_domains, set_finite_element_materials);
                    initialize_finite_elements(ptr_element_arg, initialize_finite_elements);
                    elements.data.insert({element_hash, ptr_element_arg});
                }
            };
            /*
             *
             */
            auto ptr_cell = SharedPointer<FiniteElementE<EE, D, Med>>(FiniteElementE<EE, D, Med>());
            set_element.template operator ()<EE>(cell_node_tags_arg, ptr_cell, set_element);
        }

        void
        makeCells(
        )
        {
            auto make_cells = [&] <auto I = 0> (
                    auto &
                    make_cells_imp
            )
            mutable
            {
                static_cast<Implementation *>(this)->template makeCells<element::element<D.dim, I>()>();
                if constexpr (I < element::numElements<D.dim>() - 1) {
                    make_cells_imp.template operator ()<I + 1>(make_cells_imp);
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
            auto make_element_neighbourhood = [&] <Element E, auto K = 0> (auto & self)
            mutable
            {
                auto const constexpr ed = E.dim;
                auto const constexpr et = element::elementIndex<E>();
                auto & element_map = elements_.template get<ed>().template get<et>();
                for (auto & element_map_item: element_map.data) {
                    auto & element = element_map_item.second;
                    auto const & element_nds = element.get().components.template get<ed - 1>().template get<0>();
                    for (int i = 0; i < element_nds.size(); ++i) {
                        auto const & nde = element_nds.get(i).ptr;
                        auto const & ngs = nde.get().neighbours.template get<ed - 1>().template get<K>();
                        for (auto const & neighbour: ngs.data) {
                            if (((neighbour.ptr.get().tag != element.get().tag) && K == et) || (K != et)) {
                                auto & element_ngs = element.get().neighbours.template get<0>().template get<K>();
                                Bool found = false;
                                for (auto const & ngb: element_ngs.data) {
                                    if (ngb.ptr.get().tag == neighbour.ptr.get().tag) {
                                        found = true;
                                        break;
                                    }
                                }
                                if (!found) {
                                    element_ngs.data.push_back({neighbour.ptr});
                                }
                            }
                        }
                    }
                }
                if constexpr (K < element::numElements<ed>() - 1) {
                    self.template operator()<E, K + 1>(self);
                }
            };
            /*
             *
             */
            auto make_neighbourhood = [&] <auto I = 1, auto J = 0> (auto & self)
            mutable
            {
                make_element_neighbourhood.template operator()<element::element<I, J>()>(make_element_neighbourhood);
                if constexpr (J < element::numElements<I>() - 1) {
                    self.template operator()<I, J + 1>(self);
                }
                else if constexpr (I < D.dim) {
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
                if constexpr (!element::Point<EE>) {
                    for (auto const & c_ : elem_arg.get().components.template get<I>().template get<J>().data) {
                        os << "layer : " << I << " type : " << J << " <-- " << c_.ptr.get().hash() << std::endl;
                    }
                    if constexpr (J < element::numComponents<EE, I>() - 1) {
                        self.template operator()<EE, I, J + 1>(elem_arg, self);
                    }
                    else if constexpr (I < element::numComponents<EE>() - 1) {
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
                for (auto const & c_ : elem_arg.get().neighbours.template get<I>().template get<J>().data) {
                    if constexpr (!element::Point<EE> && I == 0) {
                        os << "layer : " << I << " type : " << J << " <-> " << c_.ptr.get().hash() << std::endl;
                    }
                    else {
                        os << "layer : " << I << " type : " << J << " --> " << c_.ptr.get().hash() << std::endl;
                    }
                }
                if constexpr (J < element::numNeighbours<EE, D.dim, I>() - 1) {
                    self.template operator()<EE, I, J + 1>(elem_arg, self);
                }
                else if constexpr (I < element::numNeighbours<EE, D.dim>() - 1) {
                    self.template operator()<EE, I + 1, 0>(elem_arg, self);
                }
            };
            /*
             *
             */
            auto print_set = [&] <auto I = 0, auto J = 0> (
                    auto const & set_arg,
                    auto & self
            )
            mutable
            {
                for (auto const & item : set_arg.template get<I>().template get<J>().data) {
                    os << "layer : " << I << " type : " << J << " <-- " << item.second.get().hash() << std::endl;
                }
                if constexpr (J < element::numElements<I>() - 1) {
                    self.template operator()<I, J + 1>(set_arg, self);
                }
                else if constexpr (I < D.dim) {
                    self.template operator()<I + 1, 0>(set_arg, self);
                }
            };
            /*
             *
             */
            auto print_sets = [&] ()
            mutable
            {
                for (auto const & set : base.element_sets_.data) {
                    os << "*** Domain : " << set.first << std::endl;
                    print_set(set.second, print_set);
                }
            };
            /*
             *
             */
            auto print_elements = [&] <auto I = 0, auto J = 0> (
                    auto & self
            )
            mutable
            {
                auto constexpr elt = element::element<I, J>();
                if constexpr (I == 0 && J == 0) {
                    os << "*** Elements : " << std::endl;
                }
                for (auto const & element : base.elements_.template get<I>().template get<J>().data) {
                    os << "* Element : " << element.second.get().hash() << std::endl;
                    print_element_components.template operator()<elt>(element.second, print_element_components);
                    print_element_neighbours.template operator()<elt>(element.second, print_element_neighbours);
                }
                if constexpr (J < element::numElements<I>() - 1) {
                    self.template operator()<I, J + 1>(self);
                }
                else if constexpr (I < D.dim) {
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

    template<Domain D, auto Med>
    struct MeshModule<MeshFormatType::Gmsh, D, Med>
    {

        auto const static constexpr element_table = Array<Element, 10>{
            element::seg_02,
            element::tri_03,
            element::qua_04,
            element::pnt_00,
        };

        static
        Element
        getElementType(
                Indx
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
                    file::File const &
                    mesh_file
            )
            :
            MeshModuleBase(mesh_file)
            {
                setGeometricalEntities();
                setPhysicalEntities();
                setPhysicalGroups();
            }

            Bool
            operator==(
                    Module const &
                    other
            )
            const = default;

            Bool
            operator!=(
                    Module const &
                    other
            )
            const = default;

            struct PhysicalEntity
            {

                Bool
                operator==(
                        PhysicalEntity const &
                        other
                )
                const = default;

                Bool
                operator!=(
                        PhysicalEntity const &
                        other
                )
                const = default;

                Indx tag;

                Indx dim;

                Strg name;

            };

            struct GeometricalEntity
            {

                Bool
                operator==(
                        GeometricalEntity const &
                        other
                )
                const = default;

                Bool
                operator!=(
                        GeometricalEntity const &
                        other
                )
                const = default;

                Indx tag;

                Indx dim;

                Array<Indx> physical_entities_tags;

                Array<Indx> bounding_entities_tags;

            };

            struct PhysicalGroup
            {

                Bool
                operator==(
                        PhysicalGroup const &
                        other
                )
                const = default;

                Bool
                operator!=(
                        PhysicalGroup const &
                        other
                )
                const = default;

                Strg name;

                Array<Array<Indx>, 4> geometrical_entities_tags;

            };

            using PhysicalGroups = Array<PhysicalGroup>;

            using GeometricalEntities = Array<Array<GeometricalEntity>, 4>;

            using PhysicalEntities = Array<PhysicalEntity>;

        private:

            void
            setGeometricalEntities()
            {
                Indx line_start = array::index(this->file_, "$Entities");
                Indx offset = 1;
                StrgStream line_stream(this->file_.get(line_start + offset));
                Indx num_points;
                Indx num_curves;
                Indx num_surfaces;
                Indx num_volumes;
                line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
                auto num_domains = Array<Indx, 4>{num_points, num_curves, num_surfaces, num_volumes};
                offset += 1;
                for (Indx i = 0; i < 4; ++i) {
                    for (Indx j = 0; j < num_domains.get(i); ++j) {
                        line_stream = StrgStream(this->file_.get(line_start + offset));
                        Indx tag;
                        line_stream >> tag;
                        if (i == 0) {
                            for (Indx k = 0; k < 3; ++k) {
                                Real a;
                                line_stream >> a;
                            }
                        } else {
                            for (Indx k = 0; k < 6; ++k) {
                                Real a;
                                line_stream >> a;
                            }
                        }
                        Indx num_physical_entities;
                        line_stream >> num_physical_entities;
                        Array<Indx> physical_entities_tags;
                        for (Indx k = 0; k < num_physical_entities; ++k) {
                            Intg physical_entity_tag;
                            line_stream >> physical_entity_tag;
                            physical_entities_tags.data.push_back(physical_entity_tag);
                        }
                        Array<Indx> bounding_entities_tags = {};
                        if (i > 0) {
                            Indx num_bounding_entities;
                            line_stream >> num_bounding_entities;
                            for (Indx k = 0; k < num_bounding_entities; ++k) {
                                Intg bounding_entity_tag;
                                line_stream >> bounding_entity_tag;
                                bounding_entities_tags.data.push_back(std::abs(bounding_entity_tag));
                            }
                        }
                        geometrical_entities.get(i).data.push_back(GeometricalEntity{
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
                Indx line_start = array::index(this->file_, "$PhysicalNames");
                Indx offset = 1;
                StrgStream line_stream(this->file_.get(line_start + offset));
                Indx num_physical_names;
                line_stream >> num_physical_names;
                offset += 1;
                for (Indx i = 0; i < num_physical_names; ++i) {
                    line_stream = StrgStream(this->file_.get(line_start + offset));
                    Indx dim;
                    Indx tag;
                    Strg name;
                    line_stream >> dim >> tag >> name;
                    file::removeCharacter(name, '"');
                    physical_entities.data.push_back(PhysicalEntity{tag, dim, name});
                    offset += 1;
                }
            }

            Array<UnorderedSet<Indx>, 4>
            getSubGeometricalEntities(
                    Indx d,
                    Indx t,
                    Array<UnorderedSet<Indx>, 4> & a
            )
            {
                a.get(d).data.insert(t);
                for (Indx i = 0; i < geometrical_entities.get(d).get(t - 1).bounding_entities_tags.size(); ++i) {
                    Indx const d2 = d - 1;
                    Indx const t2 = geometrical_entities.get(d).get(t - 1).bounding_entities_tags.get(i);
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
                        Array<UnorderedSet<Indx>, 4> & a,
                        auto & self
                )
                mutable
                {
                    a.get(d).data.insert(t);
                    for (Indx i = 0; i < geometrical_entities.get(d).get(t - 1).bounding_entities_tags.size(); ++i) {
                        a = self(d - 1, geometrical_entities.get(d).get(t - 1).bounding_entities_tags.get(i), a, self);
                    }
                    return a;
                };
                /*
                 *
                 */
                for (Indx i = 0; i < physical_entities.size(); ++i) {
                    Array<UnorderedSet<Indx>, 4> tags;
                    for (Indx j = 0; j < 4; ++j) {
                        for (Indx k = 0; k < geometrical_entities.get(j).size(); ++k) {
                            for (Indx l = 0; l < geometrical_entities.get(j).get(k).physical_entities_tags.size(); ++l) {
                                Indx const & t0 = geometrical_entities.get(j).get(k).physical_entities_tags.get(l);
                                Indx const & t1 = physical_entities.get(i).tag;
                                if (t0 == t1) {
//                                    tags = get_subs(j, k + 1, tags, get_subs);
                                    tags = getSubGeometricalEntities(j, k + 1, tags);
                                }
                            }
                        }
                    }
                    Array<Array<Indx>, 4> group_tags;
                    for (Indx j = 0; j < 4; ++j) {
                        group_tags.get(j).data.assign(tags.get(j).data.begin(), tags.get(j).data.end());
                    }
                    physical_groups.data.push_back(PhysicalGroup{physical_entities.get(i).name, group_tags});
                }
                //return groups;
            }

        public:

            GeometricalEntities geometrical_entities;
            PhysicalEntities physical_entities;
            PhysicalGroups physical_groups;

        };

        struct Implementation : public MeshBase<MeshFormatType::Gmsh, D, Med>
        {

            //using Base = MeshBase<MeshFormatType::Gmsh, D, Med>;

            void
            makeNodes()
            {
                auto & nodes = this->elements_.template get<0>().template get<0>();
                Indx line_start = array::index(this->module_.file_, "$Nodes");
                Indx offset = 1;
                StrgStream line_stream(this->module_.file_.get(line_start + offset));
                Indx num_entity_block;
                Indx num_nodes;
                Indx min_node_tag;
                Indx max_node_tag;
                line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                offset += 1;
                for (Indx i = 0; i < num_entity_block; ++i) {
                    line_stream = StrgStream(this->module_.file_.get(line_start + offset));
                    Indx entity_dim;
                    Indx entity_tag;
                    Indx parametric;
                    Indx num_nodes_in_block;
                    line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                    offset += 1;
                    for (Indx j = 0; j < num_nodes_in_block; ++j) {
                        line_stream = StrgStream(this->module_.file_.get(line_start + offset));
                        Long node_tag;
                        line_stream >> node_tag;
                        line_stream = StrgStream(this->module_.file_.get(line_start + offset + num_nodes_in_block));
                        Vector<Real, D.dim> coordinates;
                        for (int k = 0; k < D.dim; ++k) {
                            line_stream >> coordinates(k);
                        }
//                        this->setElement(node_tag, coordinates);
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
                    this->element_sets_.get(physical_groups.get(i).name);
                }
            }

            void
            makeNodeSets()
            {
                auto const & physical_groups = this->module_.physical_groups;
                for (int m = 0; m < physical_groups.size(); ++m) {
                    auto & node_set_name = physical_groups.get(m).name;
                    auto & node_set = this->element_sets_.get(node_set_name).template get<0>().template get<0>();
                    auto & nodes = this->elements_.template get<0>().template get<0>();
                    for (Indx i = 0; i < 4; ++i) {
                        for (Indx j = 0; j < physical_groups.get(m).geometrical_entities_tags.get(i).size(); ++j) {
                            Indx tag = physical_groups.get(m).geometrical_entities_tags.get(i).get(j);
                            Indx line_start = array::index(this->module_.file_, "$Nodes");
                            Indx offset = 1;
                            StrgStream line_stream(this->module_.file_.get(line_start + offset));
                            Indx num_entity_block;
                            Indx num_nodes;
                            Indx min_node_tag;
                            Indx max_node_tag;
                            line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                            offset += 1;
                            for (Indx k = 0; k < num_entity_block; ++k) {
                                line_stream = StrgStream(this->module_.file_.get(line_start + offset));
                                Indx entity_dim;
                                Indx entity_tag;
                                Indx parametric;
                                Indx num_nodes_in_block;
                                line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                                offset += 1;
                                if (entity_dim == i && entity_tag == tag) {
                                    for (Indx l = 0; l < num_nodes_in_block; ++l) {
                                        line_stream = StrgStream(this->module_.file_.get(line_start + offset));
                                        Strg node_hash = line_stream.str();
                                        node_set.data.insert({node_hash, nodes.get(node_hash)});
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

            template<Element EE>
            void
            makeCells()
            {
                Indx line_start = array::index(this->module_.file_, "$Elements");
                Indx offset = 1;
                StrgStream line_stream(this->module_.file_.get(line_start + offset));
                Indx num_entity_blocks;
                Indx num_elements;
                Indx min_element_tag;
                Indx max_element_tag;
                line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
                offset += 1;
                for (Indx i = 0; i < num_entity_blocks; ++i) {
                    line_stream = StrgStream(this->module_.file_.get(line_start + offset));
                    Indx entity_dim = 0;
                    Indx entity_tag = 0;
                    Indx element_type_tag = 0;
                    Indx num_elements_in_block = 0;
                    line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                    offset += 1;
                    Element const element_type_arg = getElementType(element_type_tag);
                    if (entity_dim == D.dim and getElementType(element_type_tag) == EE) {
                        for (Indx j = 0; j < num_elements_in_block; ++j) {
                            line_stream = StrgStream(this->module_.file_.get(line_start + offset));
                            Indx tag;
                            line_stream >> tag;
                            Array<Indx, EE.num_nodes> node_tags;
                            for (Indx k = 0; k < EE.num_nodes; ++k) {
                                line_stream >> node_tags.get(k);
                            }
                            this->template makeElement<EE>(node_tags);
                            offset += 1;
                        }
                    } else {
                        offset += num_elements_in_block;
                    }
                }
            }

        };

    };

}

#endif //LOLITA_LOLITA_CORE_MESH_HXX
