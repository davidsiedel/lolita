//
// Created by dsiedel on 17/04/22.
//

#include <execution>

#include "lolita/lolita_core.hxx"
//#include "lolita/lolita_unordered_map.hxx"
//#include "lolita/lolita_unordered_set.hxx"
#include "lolita/lolita_core_finite_element.hxx"
#include "lolita/lolita_file.hxx"

#ifndef LOLITA_LOLITA_CORE_MESH_HXX
#define LOLITA_LOLITA_CORE_MESH_HXX

namespace lolita::core::mesh
{

    using Element = core::element::Element;

    template<Element E, Domain D, auto Med>
    using FiniteElementE = core::element::FiniteElement<E, D, Med>;

//    template<Element E, Domain D, auto Med>
//    using ElementPointerMap = UnorderedMap<Strg, SharedPointer<FiniteElementE<E, D, Med>>>;

    template<MeshFormatType, Domain D, auto Med>
    struct MeshModule;

    template<MeshFormatType Mft, Domain D, auto Med>
    struct MeshBase
    {

    protected:

        template<Element E>
        using ElementPointerMap = UnorderedMap<Strg, SharedPointer<FiniteElementE<E, D, Med>>>;

        using Module = typename MeshModule<Mft, D, Med>::Module;

        using Implementation = typename MeshModule<Mft, D, Med>::Implementation;

    public:

        using ElementPointers = element::Elements<D.dim, ElementPointerMap>;

        using ElementPointersSets = UnorderedMap<Strg, ElementPointers>;

        using DomainPointers = Array<SharedPointer<MeshInteriorDomain<Med>>>;

        file::File const & file_;

        Module module_;

        ElementPointers elements_;

        ElementPointersSets element_sets_;

        DomainPointers mesh_domains_;

        MeshBase(
                file::File const &
                mesh_file,
                Array<MeshInteriorDomain<Med>> &&
                domains
        )
        :
        file_(mesh_file),
        module_(mesh_file)
        {
            //Module const auxiliary_data(mesh_file);
            this->setDomainPointers(domains);
            this->makeNodes();
            this->makeElementSets();
            for (auto const & item : this->element_sets_.data) {
                this->makeNodeSet(item.first);
            }
            this->makeCells();
            this->makeNeighbourhood();
        }

        void
        makeElementSets()
        {
            static_cast<Implementation *>(this)->makeElementSets();
        }

        void
        makeNodeSet(
                Strg element_set_name
        )
        {
            static_cast<Implementation *>(this)->makeNodeSet(element_set_name);
        }

        void
        makeNodes()
        {
            static_cast<Implementation *>(this)->makeNodes();
        }

        /*
         *
         */

        void
        setDomainPointers(
                Array<MeshInteriorDomain<Med>> const &
                domains
        )
        {
            for (auto && domain: domains.data) {
                mesh_domains_.data.push_back(SharedPointer<MeshInteriorDomain<Med>>(domain));
            }
        }

    protected:

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
        makeCell(
                auto const &
                cell_node_tags_arg
        )
        {
            /*
             *
             */
            auto get_element_domains = [&] <Element E> (
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
                auto bound_node_tags = Array<Indx, element::component<E, I, J>().num_nodes>();
                for (Indx j = 0; j < element::component<E, I, J>().num_nodes; ++j) {
                    auto const & element_node_connectivity = FiniteElementE<E, D, Med>::node_connectivity;
                    auto const & bound_node_connectivity = element_node_connectivity.template get<I>().template get<J>();
                    auto const k = bound_node_connectivity.get(index_arg).get(j);
                    bound_node_tags.get(j) = node_tags_arg.get(k);
                }
                return bound_node_tags;
            };
            /*
             *
             */
            auto get_element_hash = [] <Element E> (
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
            auto set_components = [&] <Element E, auto K = 0, auto I = 0, auto J = 0> (
                    auto &
                    ptr_element_arg,
                    auto &
                    set_components_imp
            )
            mutable
            {
                if constexpr (I == 0 && J == 0) {
                    pointer::make(ptr_element_arg.get().template get<K>());
                }
                auto & elt = ptr_element_arg.get();
                for (int i = 0; i < element::numComponents<E, I, J>(); ++i) {
                    auto & rhs = elt.components.template get<I>().template get<J>().get(i);
                    auto & lhs = elt.template get<K>().get().components.template get<I>().template get<J>().get(i);
                    lhs.ptr = rhs.ptr.get().template get<K>();
                    lhs.orientation = rhs.orientation;
                    //print("hello", K, I, J, i);
                }
                if constexpr (J < element::numComponents<E, I>() - 1) {
                    set_components_imp.template operator()<E, K, I, J + 1>(ptr_element_arg, set_components_imp);
                }
                else if constexpr (I < element::numComponents<E>() - 1) {
                    set_components_imp.template operator()<E, K, I + 1, 0>(ptr_element_arg, set_components_imp);
                }
                else if constexpr (K < elt.size() - 1) {
                    set_components_imp.template operator()<E, K + 1, 0, 0>(ptr_element_arg, set_components_imp);
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
                auto element_hash = get_element_hash.template operator ()<E>(element_node_tags_arg);
                for (auto i = 0; i < element::numComponents<E, I, J>(); ++i) {
                    //auto component_node_tags = getComponentNodeTags<E, I, J>(element_node_tags_arg, i);
                    auto component_node_tags = get_component_node_tags.template operator ()<E, I, J>(element_node_tags_arg, i);
                    auto component_hash = get_element_hash.template operator ()<cmp>(component_node_tags);
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
                    set_components.template operator ()<E>(ptr_element_arg, set_components);
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
        makeCells()
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

        struct Module
        {

            explicit
            Module(
                    file::File const &
                    mesh_file
            )
            :
            file_(mesh_file)
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

            file::File const & file_;

            using PhysicalGroups = Array<PhysicalGroup>;

            using GeometricalEntities = Array<Array<GeometricalEntity>, 4>;

            using PhysicalEntities = Array<PhysicalEntity>;

        private:

            void
            setGeometricalEntities()
            {
                //GeometricalEntities _geometrical_entities;
                Indx line_start = array::index(file_, "$Entities");
                Indx offset = 1;
                StrgStream line_stream(file_.get(line_start + offset));
                Indx num_points;
                Indx num_curves;
                Indx num_surfaces;
                Indx num_volumes;
                line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
                auto num_domains = Array<Indx, 4>{num_points, num_curves, num_surfaces, num_volumes};
                offset += 1;
                for (Indx i = 0; i < 4; ++i) {
                    for (Indx j = 0; j < num_domains.get(i); ++j) {
                        line_stream = StrgStream(file_.get(line_start + offset));
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
                //return geometrical_entities;
            }

            void
            setPhysicalEntities()
            {
                //PhysicalEntities _physical_entities;
                Indx line_start = array::index(file_, "$PhysicalNames");
                Indx offset = 1;
                StrgStream line_stream(file_.get(line_start + offset));
                Indx num_physical_names;
                line_stream >> num_physical_names;
                offset += 1;
                for (Indx i = 0; i < num_physical_names; ++i) {
                    line_stream = StrgStream(file_.get(line_start + offset));
                    Indx dim;
                    Indx tag;
                    Strg name;
                    line_stream >> dim >> tag >> name;
                    file::removeCharacter(name, '"');
                    physical_entities.data.push_back(PhysicalEntity{tag, dim, name});
                    offset += 1;
                }
                //return physical_entities;
            }

            Array<UnorderedSet<Indx>, 4>
            getSubGeometricalEntities(
                    Indx const & d,
                    Indx const & t,
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

            using Base = MeshBase<MeshFormatType::Gmsh, D, Med>;

            void
            makeNodes()
            {
                auto & nodes = this->elements_.template get<0>().template get<0>();
                Indx line_start = array::index(this->file_, "$Nodes");
                Indx offset = 1;
                StrgStream line_stream(this->file_.get(line_start + offset));
                Indx num_entity_block;
                Indx num_nodes;
                Indx min_node_tag;
                Indx max_node_tag;
                line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                offset += 1;
                for (Indx i = 0; i < num_entity_block; ++i) {
                    line_stream = StrgStream(this->file_.get(line_start + offset));
                    Indx entity_dim;
                    Indx entity_tag;
                    Indx parametric;
                    Indx num_nodes_in_block;
                    line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                    offset += 1;
                    for (Indx j = 0; j < num_nodes_in_block; ++j) {
                        line_stream = StrgStream(this->file_.get(line_start + offset));
                        Indx node_tag;
                        line_stream >> node_tag;
                        line_stream = StrgStream(this->file_.get(line_start + offset + num_nodes_in_block));
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
            makeNodeSet(
                    Strg
                    element_set_name
            )
            {
                auto const & physical_groups = this->module_.physical_groups;
                Indx physical_index;
                for (int i = 0; i < physical_groups.size(); ++i) {
                    if (physical_groups.get(i).name == element_set_name) {
                        physical_index = i;
                    }
                }
                auto & element_set_hashes = this->element_sets_.get(element_set_name).template get<0>().template get<0>();
                auto & nodes = this->elements_.template get<0>().template get<0>();
                for (Indx i = 0; i < 4; ++i) {
                    for (Indx j = 0; j < physical_groups.get(physical_index).geometrical_entities_tags.get(i).size(); ++j) {
                        Indx tag = physical_groups.get(physical_index).geometrical_entities_tags.get(i).get(j);
    //                    Indx line_start = mesh_file.getLineIndex("$Nodes");
                        Indx line_start = array::index(this->file_, "$Nodes");
                        Indx offset = 1;
                        StrgStream line_stream(this->file_.get(line_start + offset));
                        Indx num_entity_block;
                        Indx num_nodes;
                        Indx min_node_tag;
                        Indx max_node_tag;
                        line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                        offset += 1;
                        for (Indx k = 0; k < num_entity_block; ++k) {
                            line_stream = StrgStream(this->file_.get(line_start + offset));
                            Indx entity_dim;
                            Indx entity_tag;
                            Indx parametric;
                            Indx num_nodes_in_block;
                            line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                            offset += 1;
                            if (entity_dim == i && entity_tag == tag) {
                                for (Indx l = 0; l < num_nodes_in_block; ++l) {
                                    line_stream = StrgStream(this->file_.get(line_start + offset));
                                    Strg node_hash = line_stream.str();
                                    element_set_hashes.data.insert({node_hash, nodes.get(node_hash)});
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

            template<Element EE>
            void
            makeCells()
            {
                Indx line_start = array::index(this->file_, "$Elements");
                Indx offset = 1;
                StrgStream line_stream(this->file_.get(line_start + offset));
                Indx num_entity_blocks;
                Indx num_elements;
                Indx min_element_tag;
                Indx max_element_tag;
                line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
                offset += 1;
                for (Indx i = 0; i < num_entity_blocks; ++i) {
                    line_stream = StrgStream(this->file_.get(line_start + offset));
                    Indx entity_dim = 0;
                    Indx entity_tag = 0;
                    Indx element_type_tag = 0;
                    Indx num_elements_in_block = 0;
                    line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                    offset += 1;
                    Element const element_type_arg = getElementType(element_type_tag);
                    if (entity_dim == D.dim and getElementType(element_type_tag) == EE) {
                        for (Indx j = 0; j < num_elements_in_block; ++j) {
                            line_stream = StrgStream(this->file_.get(line_start + offset));
                            Indx tag;
                            line_stream >> tag;
                            Array<Indx, EE.num_nodes> node_tags;
                            for (Indx k = 0; k < EE.num_nodes; ++k) {
                                line_stream >> node_tags.get(k);
                            }
                            this->template makeCell<EE>(node_tags);
                            offset += 1;
                        }
                    } else {
                        offset += num_elements_in_block;
                    }
                }
            }

        };

    };

    template<typename M>
    static inline
    void
    printMesh(
            M const &
            mesh_arg
    )
    {
        auto elsets = mesh_arg.element_sets_;
        auto elems = mesh_arg.elements_;
        print("START ELEMS");
        print("NODES");
        for (auto const & [key, val] : elems.template get<0>().template get<0>().data) {
            print(key, val.get().tag);
            for (auto const & ghj: val.get().neighbours.template get<0>().template get<0>().data) {
                print(" --> ", ghj.ptr.get().hash());
            }
            for (auto const & ghj: val.get().neighbours.template get<1>().template get<0>().data) {
                print(" ----> ", ghj.ptr.get().hash());
            }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
        }
        print("NODES DONE");

        print("SEGS");
        for (auto const & [key, val] : elems.template get<1>().template get<0>().data) {
            print(key, val.get().tag);
            for (int i = 0; i < val.get().components.template get<0>().template get<0>().size(); ++i) {
                print("--> ", val.get().components.template get<0>().template get<0>().get(i).ptr.get().hash(), val.get().components.template get<0>().template get<0>().get(i).orientation);
            }
            for (auto const & ghj: val.get().neighbours.template get<1>().template get<0>().data) {
                print(" ----> ", ghj.ptr.get().hash());
            }
            for (auto const & ghj: val.get().neighbours.template get<0>().template get<0>().data) {
                print(" <-> ", ghj.ptr.get().hash());
            }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
        }
        print("SEGS DONE");

        print("QUAS");
        for (auto const & [key, val] : elems.template get<2>().template get<1>().data) {
            print(key, val.get().tag);
            for (int i = 0; i < val.get().components.template get<0>().template get<0>().size(); ++i) {
                print("--> ", val.get().components.template get<0>().template get<0>().get(i).ptr.get().hash(), val.get().components.template get<0>().template get<0>().get(i).orientation);
            }
            for (int i = 0; i < val.get().components.template get<1>().template get<0>().size(); ++i) {
                print("----> ", val.get().components.template get<1>().template get<0>().get(i).ptr.get().hash(), val.get().components.template get<1>().template get<0>().get(i).orientation);
            }
            for (auto const & ghj: val.get().neighbours.template get<0>().template get<0>().data) {
                print(" <-> ", ghj.ptr.get().hash());
            }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
        }
        print("QUAS DONE");

        print("TRIS");
        for (auto const & [key, val] : elems.template get<2>().template get<0>().data) {
            print(key, val.get().tag);
            for (int i = 0; i < val.get().components.template get<0>().template get<0>().size(); ++i) {
                print("--> ", val.get().components.template get<0>().template get<0>().get(i).ptr.get().hash(), val.get().components.template get<0>().template get<0>().get(i).orientation);
            }
            for (int i = 0; i < val.get().components.template get<1>().template get<0>().size(); ++i) {
                print("----> ", val.get().components.template get<1>().template get<0>().get(i).ptr.get().hash(), val.get().components.template get<1>().template get<0>().get(i).orientation);
            }
            for (auto const & ghj: val.get().neighbours.template get<0>().template get<0>().data) {
                print(" <-> ", ghj.ptr.get().hash());
            }
//                for (auto const & dom_label: val.get().domain_names.data) {
//                    print(" [-] ", dom_label);
//                }
        }
        print("TRIS DONE");

        print("END ELEMS");

        print("START ELSETS");
        for (auto const & [key, val] : elsets.data) {
            print(key);
            for (auto const & [elhash, elptr] : val.template get<0>().template get<0>().data) {
                print(elhash);
            }
        }
        for (auto const & [key, val] : elsets.data) {
            print(key);
            for (auto const & [elhash, elptr] : val.template get<1>().template get<0>().data) {
                print(elhash);
            }
        }
        for (auto const & [key, val] : elsets.data) {
            print(key);
            for (auto const & [elhash, elptr] : val.template get<2>().template get<0>().data) {
                print(elhash);
            }
        }
        print("END ELSETS");
    }

}

#endif //LOLITA_LOLITA_CORE_MESH_HXX
