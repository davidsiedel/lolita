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

    namespace elt = lolita::core::element;

    template<MeshFormatType, Domain D, auto Med>
    struct MeshModule;

    template<MeshFormatType Mft, Domain D, auto Med>
    struct MeshBase
    {

    protected:

        template<elt::Element E>
        using Element = elt::FiniteElement<E, D, Med>;

        template<elt::Element E, Indx I, Indx J>
        using ComponentPointer = typename Element<E>::Components::template Type<I>::template Type<J>::Type;

        template<elt::Element E, Indx I, Indx J>
        using Component = typename ComponentPointer<E, I, J>::Type;

        template<elt::Element E>
        using ElementPointerMap = UnorderedMap<Strg, SharedPointer<Element<E>>>;

        using Module = typename MeshModule<Mft, D, Med>::Module;

        using Implementation = typename MeshModule<Mft, D, Med>::Implementation;

    public:

        using ElementPointers = elt::Elements<D.dim, ElementPointerMap>;

        using ElementPointersSets = UnorderedMap<Strg, ElementPointers>;

        using DomainPointers = Array<SharedPointer<MeshInteriorDomain<Med>>>;

        ElementPointers element_collection;

        ElementPointersSets element_sets;

        DomainPointers mesh_domains;

        MeshBase(
                file::File const &
                mesh_file,
                Array<MeshInteriorDomain<Med>> &&
                domains
        )
        {
            Module const auxiliary_data(mesh_file);
            this->setDomainPointers(domains);
            this->setNodes(mesh_file);
            this->initializeElementSets(mesh_file, auxiliary_data);
            for (auto const & [key, val] : this->element_sets.data) {
                this->setNodesSet(mesh_file, key, auxiliary_data);
            }
            this->setCells(mesh_file);
            this->setNeighbours();
        }

//        Bool
//        operator==(
//                ParserBase const &
//                other
//        )
//        const = default;
//
//        Bool
//        operator!=(
//                ParserBase const &
//                other
//        )
//        const = default;

        void
        initializeElementSets(
                file::File const & mesh_file,
                Module const & helper
        )
        {
            static_cast<Implementation *>(this)->initializeElementSets(mesh_file, helper);
        }

        void
        setNodesSet(
                file::File const & mesh_file,
                Strg element_set_name,
                Module const & helper
        )
        {
            static_cast<Implementation *>(this)->setNodesSet(mesh_file, element_set_name, helper);
        }

        void
        setNodes(
                file::File const & mesh_file
        )
        {
            static_cast<Implementation *>(this)->setNodes(mesh_file);
        }

        template<Indx I = 0>
        void
        setCells(
                file::File const & mesh_file
        )
        {
            static_cast<Implementation *>(this)->template setCells<I>(mesh_file);
            auto const constexpr ddd = ElementPointers::template Type<D.dim>::size();
            if constexpr (I < ddd - 1) setCells<I + 1>(mesh_file);
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
                mesh_domains.data.push_back(SharedPointer<MeshInteriorDomain<Med>>(domain));
            }
        }

        template<elt::Element E>
        static
        auto
        getElementHash(
//                Array<Indx, E.num_nodes>
                auto
                node_tags_arg
        )
        {
            StrgStream element_hash;
            if (E.num_nodes > 1) {
                std::sort(std::execution::par_unseq, node_tags_arg.data.begin(), node_tags_arg.data.end());
            }
            for (Indx i = 0; i < E.num_nodes; ++i) {
                element_hash << node_tags_arg.get(i);
            }
            return element_hash.str();
        }

        template<elt::Element E, Indx I, Indx J>
        static
        auto
        getComponentNodeTags(
//                Array<Indx, E.num_nodes> const &
                auto const &
                node_tags_arg,
                Indx
                component_index_arg
        )
        {
//            Array<Indx, Component<E, I, J>::num_nodes> bound_node_tags;
            auto bound_node_tags = Array<Indx, Component<E, I, J>::element.num_nodes>();
            for (Indx j = 0; j < Component<E, I, J>::element.num_nodes; ++j) {
//                Indx k = Element<E>::node_connectivity.template get<I>().template get<J>()(component_index_arg)(j);
                auto const & ndscon = Element<E>::node_connectivity.template get<I>().template get<J>();
                auto const k = ndscon.get(component_index_arg).get(j);
//                Indx k = Element<E>::node_connectivity.template get<I>().template get<J>().get(component_index_arg, j);
                bound_node_tags.get(j) = node_tags_arg.get(k);
            }
            return bound_node_tags;
        }

        template<elt::Element E>
        auto
        getElementDomains(
                Array<Indx, E.num_nodes> const &
                node_tags_arg
        )
        {
            Array<Strg> domain_names;
            for (auto const & [set_name, element_tags] : element_sets.data) {
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
        }

    protected:

        void
        setElement(
                Indx const &
                node_tag_arg,
                Vector<Real, D.dim> const &
                coordinates_arg
        )
        {
            auto & nodes = element_collection.template get<0>().template get<0>();
            auto hash = getElementHash<elt::pnt_00>(Array<Indx, 1>{node_tag_arg});
            auto ptr_node = SharedPointer<Element<elt::pnt_00>>(Element<elt::pnt_00>());
            //ptr_node.get().coordinates = coordinates_arg;
            ptr_node.get().coordinates = SharedPointer<Vector<Real, D.dim>>(coordinates_arg);
            ptr_node.get().tag = nodes.size();
            nodes.data.insert({hash, ptr_node});
        }

        template<elt::Element E, Indx I = 0, Indx J = 0>
        void
        setElement(
                Array<Indx, E.num_nodes> const &
                element_node_tags_arg,
                typename Element<E>::Components &
                element_components_arg,
                SharedPointer<Element<E>> &
                ptr_element_arg
        )
        {
            /*
             * defining useful constants, to fetch both the map of all current structure elements, and that of all
             * is components.
             */
            auto const constexpr cmp = Component<E, I, J>::element;
            auto const constexpr bd = cmp.dim;
            auto const constexpr bt = elt::elementIndex<cmp>();
            auto const constexpr ed = E.dim;
            auto const constexpr et = elt::elementIndex<E>();
            /*
             * fetching all current structure elements, all current component elements, and all nodes in the mesh.
             * fetching the collection of components elements defining the Ith layer and Jth type of the current
             * structure element composition.
             * component_elements is the
             */
            auto & components = element_collection.template get<bd>().template get<bt>();
            auto & elements = element_collection.template get<ed>().template get<et>();
            auto & nodes = element_collection.template get<0>().template get<0>();
            auto & element_component_array = element_components_arg.template get<I>().template get<J>();
            auto element_hash = getElementHash<E>(element_node_tags_arg);
            for (Indx i = 0; i < element_component_array.size(); ++i) {
                auto component_node_tags = getComponentNodeTags<E, I, J>(element_node_tags_arg, i);
                auto component_hash = getElementHash<Component<E, I, J>::element>(component_node_tags);
                if constexpr (Component<E, I, J>::element == elt::pnt_00) {
                    /*
                     * if the component element is a point, it is already defined, since all nodes in the mesh are
                     * built prior to building structure elements.
                     */
                    element_component_array.get(i) = ComponentPointer<E, I, J>(nodes.get(component_hash), 1);
                }
                else {
                    /*
                     * else, one needs to check either the component element has already been defined.
                     */
                    if (components.data.contains(component_hash)) {
                        /*
                         * if the component element has already been defined, one just fetches it from the existing
                         * ones.
                         */
                        element_component_array.get(i) = ComponentPointer<E, I, J>(components.get(component_hash), -1);
                    }
                    else {
                        /*
                         * otherwise, one creates it by recursively calling the present function.
                         * A pointer to the component element is created, as well as a component element collection.
                         */
                        typename Element<cmp>::Components component_components;
                        auto ptr_component = SharedPointer<Element<cmp>>(Element<cmp>());
                        setElement<Component<E, I, J>::element>(
                                component_node_tags,
                                component_components,
                                ptr_component
                        );
                        /*
                         * Once the component element and its components defined, one adds it to existing ones, and
                         * assigns it to the current structure element
                         */
                        element_component_array.get(i) = ComponentPointer<E, I, J>(components.get(component_hash), 1);
                    }
                }
                auto & component_neighbours = components.get(component_hash).get().neighbours;
                if constexpr (Component<E, I, J>::element == elt::pnt_00) {
                    /*
                     * if the current component is a point
                     */
                    Indx const constexpr nd = ed - 1;
                    using ComponentNeighbours = typename Element<cmp>::Neighbours;
                    using ComponentNeighbour = typename ComponentNeighbours::template Type<nd>::template Type<et>::Type;
                    component_neighbours.template get<nd>().template get<et>().data.push_back(
                            ComponentNeighbour{ptr_element_arg}
                    );
                } else {
                    /*
                     * if the current component is not a point
                     */
                    Indx const constexpr nd = ed - bd;
                    using ComponentNeighbours = typename Element<cmp>::Neighbours;
                    using ComponentNeighbour = typename ComponentNeighbours::template Type<nd>::template Type<et>::Type;
                    component_neighbours.template get<nd>().template get<et>().data.push_back(
                            ComponentNeighbour{ptr_element_arg}
                    );
                }
            }
            if constexpr (J < Element<E>::Components::template Type<I>::size() - 1) {
                setElement<E, I, J + 1>(
                        element_node_tags_arg,
                        element_components_arg,
                        ptr_element_arg
                );
            }
            if constexpr(Intg(E.dim) - (1 + I) > 0) {
                setElement<E, I + 1, 0>(
                        element_node_tags_arg,
                        element_components_arg,
                        ptr_element_arg
                );
            }
            if constexpr((I == 0) && (J == 0)) {
                ptr_element_arg.get().tag = elements.size();
                ptr_element_arg.get().components = element_components_arg;
                //ptr_element_arg.get().initialize();
                elements.data.insert({element_hash, ptr_element_arg});
            }
        }

        template<Indx I = 1, Indx J = 0, Indx K = 0>
        void
        setNeighbours()
        {
            using Elementts = typename elt::Elements<D.dim, Element>;
            using Elementt = typename Elementts::template Type<I>::template Type<J>;
            using ElementtNeighbour = typename Elementt::Neighbours::template Type<0>::template Type<K>::Type;
            /*
             * for each element pointer in the mesh
             */
            auto & elements_selection = element_collection.template get<I>().template get<J>();
            for (auto & [element_label, ptr_element]: elements_selection.data) {
                auto const & element_nodes = ptr_element.get().components.template get<I - 1>().template get<0>();
                for (Indx i = 0; i < element_nodes.size(); ++i) {
                    auto const & ptr_node = element_nodes.get(i);
                    auto const & ptr_neighbours = ptr_node.get().neighbours.template get<I - 1>().template get<K>();
                    for (auto const & ptr_neighbour: ptr_neighbours.data) {
                        if (((ptr_neighbour.get().tag != ptr_element.get().tag) && K == J) || (K != J)) {
                            auto & neighbour_map = ptr_element.get().neighbours.template get<0>().template get<K>();
                            Bool found = false;
                            for (auto const & item: neighbour_map.data) {
                                if (item == ptr_neighbour) {
                                    found = true;
                                    break;
                                }
                            }
                            if (!found) {
                                neighbour_map.data.push_back(ElementtNeighbour{ptr_neighbour});
                            }
                        }
                    }
                }
            }
            if constexpr (K < ElementPointers::template Type<I>::size() - 1) {
                setNeighbours<I, J, K + 1>();
            }
            if constexpr (J < ElementPointers::template Type<I>::size() - 1) {
                setNeighbours<I, J + 1, 0>();
            }
            if constexpr (I < D.dim) {
                setNeighbours<I + 1, 0, 0>();
            }
        }

        template<Indx I = 0>
        Bool
        isNodeIsolated(
                Strg const &
                node_hash_arg,
                Bool
                is_isolated_arg = true
        )
        {
            auto & ptr_node = element_collection.template operator ()<0>().template operator ()<0>()(node_hash_arg);
            auto & node_neighbours = ptr_node().neighbours.template get<0>().template get<I>();
            if (node_neighbours.size() > 0) {
                is_isolated_arg = false;
            }
            if constexpr(I < ElementPointers::template Type<1>::size() - 1) {
                is_isolated_arg = isNodeIsolated<I + 1>(node_hash_arg, is_isolated_arg);
            }
            return is_isolated_arg;
        }

    };

    template<Domain D, auto Med>
    struct MeshModule<MeshFormatType::Gmsh, D, Med>
    {

        auto const static constexpr element_table = Array<elt::Element, 10>{
            elt::seg_02,
            elt::tri_03,
            elt::qua_04,
            elt::pnt_00,
        };

        static
        elt::Element
        getElementType(
                Indx
                tag_arg
        )
        {
            if (tag_arg == 15) {
                return elt::pnt_00;
            }
            else if (tag_arg == 1) {
                return elt::seg_02;
            }
            else if (tag_arg == 2) {
                return elt::tri_03;
            }
            else if (tag_arg == 3) {
                return elt::qua_04;
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
                    geometrical_entities(getGeometricalEntities(mesh_file)),
                    physical_entities(getPhysicalEntities(mesh_file)),
                    physical_groups(getPhysicalGroups(mesh_file))
            {}

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

            static
            GeometricalEntities
            getGeometricalEntities(
                    file::File const & mesh_file
            )
            {
                GeometricalEntities geometrical_entities;
    //            Indx line_start = mesh_file.getLineIndex("$Entities");
                Indx line_start = array::index(mesh_file, "$Entities");
                Indx offset = 1;
                StrgStream line_stream(mesh_file.get(line_start + offset));
                Indx num_points;
                Indx num_curves;
                Indx num_surfaces;
                Indx num_volumes;
                line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
                Array<Indx, 4> num_domains = {
                        num_points,
                        num_curves,
                        num_surfaces,
                        num_volumes
                };
                offset += 1;
                for (Indx i = 0; i < 4; ++i) {
                    for (Indx j = 0; j < num_domains.get(i); ++j) {
                        line_stream = StrgStream(mesh_file.get(line_start + offset));
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
                return geometrical_entities;
            }

            static
            PhysicalEntities
            getPhysicalEntities(
                    file::File const & mesh_file
            )
            {
                PhysicalEntities physical_entities;
    //            Indx line_start = mesh_file.getLineIndex("$PhysicalNames");
                Indx line_start = array::index(mesh_file, "$PhysicalNames");
                Indx offset = 1;
                StrgStream line_stream(mesh_file.get(line_start + offset));
                Indx num_physical_names;
                line_stream >> num_physical_names;
                offset += 1;
                for (Indx i = 0; i < num_physical_names; ++i) {
                    line_stream = StrgStream(mesh_file.get(line_start + offset));
                    Indx dim;
                    Indx tag;
                    Strg name;
                    line_stream >> dim >> tag >> name;
                    file::removeCharacter(name, '"');
                    physical_entities.data.push_back(PhysicalEntity{tag, dim, name});
                    offset += 1;
                }
                return physical_entities;
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

            PhysicalGroups
            getPhysicalGroups(
                    file::File const & mesh_file
            )
            {
                PhysicalGroups groups;
                for (Indx i = 0; i < physical_entities.size(); ++i) {
                    Array<UnorderedSet<Indx>, 4> tags;
                    for (Indx j = 0; j < 4; ++j) {
                        for (Indx k = 0; k < geometrical_entities.get(j).size(); ++k) {
                            for (Indx l = 0; l < geometrical_entities.get(j).get(k).physical_entities_tags.size(); ++l) {
                                Indx const & t0 = geometrical_entities.get(j).get(k).physical_entities_tags.get(l);
                                Indx const & t1 = physical_entities.get(i).tag;
                                if (t0 == t1) {
                                    tags = getSubGeometricalEntities(j, k + 1, tags);
                                }
                            }
                        }
                    }
                    Array<Array<Indx>, 4> group_tags;
                    for (Indx j = 0; j < 4; ++j) {
                        group_tags.get(j).data.assign(tags.get(j).data.begin(), tags.get(j).data.end());
                    }
                    groups.data.push_back(PhysicalGroup{physical_entities.get(i).name, group_tags});
                }
                return groups;
            }

        public:

            GeometricalEntities geometrical_entities;
            PhysicalEntities physical_entities;
            PhysicalGroups physical_groups;

        };

        struct Implementation : public MeshBase<MeshFormatType::Gmsh, D, Med>
        {

            using Base = MeshBase<MeshFormatType::Gmsh, D, Med>;

//            using Base::Base;

            template<elt::Element F>
            using EEl = typename Base::template Element<F>;

            void
            setNodes(
                    file::File const &
                    mesh_file
            )
            {
                auto & nodes = this->element_collection.template get<0>().template get<0>();
                Indx line_start = array::index(mesh_file, "$Nodes");
                Indx offset = 1;
                StrgStream line_stream(mesh_file.get(line_start + offset));
                Indx num_entity_block;
                Indx num_nodes;
                Indx min_node_tag;
                Indx max_node_tag;
                line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                offset += 1;
                for (Indx i = 0; i < num_entity_block; ++i) {
                    line_stream = StrgStream(mesh_file.get(line_start + offset));
                    Indx entity_dim;
                    Indx entity_tag;
                    Indx parametric;
                    Indx num_nodes_in_block;
                    line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                    offset += 1;
                    for (Indx j = 0; j < num_nodes_in_block; ++j) {
                        line_stream = StrgStream(mesh_file.get(line_start + offset));
                        Indx node_tag;
                        line_stream >> node_tag;
                        line_stream = StrgStream(mesh_file.get(line_start + offset + num_nodes_in_block));
                        Vector<Real, D.dim> coordinates;
                        for (int k = 0; k < D.dim; ++k) {
                            line_stream >> coordinates(k);
                        }
                        this->setElement(node_tag, coordinates);
                        offset += 1;
                    }
                    offset += num_nodes_in_block;
                }
            }

            void
            initializeElementSets(
                    file::File const &
                    mesh_file,
//                    ParserHelper1<MeshFormatType::Gmsh, D, Med> const &
                    Module const &
                    helper
            )
            {
                auto const & physical_groups = helper.physical_groups;
                for (int i = 0; i < physical_groups.size(); ++i) {
                    this->element_sets.get(physical_groups.get(i).name);
                }
            }

            void
            setNodesSet(
                    file::File const &
                    mesh_file,
                    Strg
                    element_set_name,
//                    ParserHelper1<MeshFormatType::Gmsh, D, Med> const &
                    Module const &
                    helper
            )
            {
                auto const & physical_groups = helper.physical_groups;
                Indx physical_index;
                for (int i = 0; i < physical_groups.size(); ++i) {
                    if (physical_groups.get(i).name == element_set_name) {
                        physical_index = i;
                    }
                }
                auto & element_set_hashes = this->element_sets.get(element_set_name).template get<0>().template get<0>();
                auto & nodes = this->element_collection.template get<0>().template get<0>();
                for (Indx i = 0; i < 4; ++i) {
                    for (Indx j = 0; j < physical_groups.get(physical_index).geometrical_entities_tags.get(i).size(); ++j) {
                        Indx tag = physical_groups.get(physical_index).geometrical_entities_tags.get(i).get(j);
    //                    Indx line_start = mesh_file.getLineIndex("$Nodes");
                        Indx line_start = array::index(mesh_file, "$Nodes");
                        Indx offset = 1;
                        StrgStream line_stream(mesh_file.get(line_start + offset));
                        Indx num_entity_block;
                        Indx num_nodes;
                        Indx min_node_tag;
                        Indx max_node_tag;
                        line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                        offset += 1;
                        for (Indx k = 0; k < num_entity_block; ++k) {
                            line_stream = StrgStream(mesh_file.get(line_start + offset));
                            Indx entity_dim;
                            Indx entity_tag;
                            Indx parametric;
                            Indx num_nodes_in_block;
                            line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                            offset += 1;
                            if (entity_dim == i && entity_tag == tag) {
                                for (Indx l = 0; l < num_nodes_in_block; ++l) {
                                    line_stream = StrgStream(mesh_file.get(line_start + offset));
    //                                hashes.data.insert(line_stream.str());
                                    Strg node_hash = line_stream.str();
                                    element_set_hashes.data.insert({node_hash, nodes.get(node_hash)});
    //                                nodes(node_hash)().domain_names.data.push_back(element_set_name);
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

            template<Indx I = 0>
            void
            setCells(
                    file::File const & mesh_file
            )
            {
                auto const constexpr dim_euclidean = D.dim;
                using Elements = typename elt::Elements<dim_euclidean, EEl>;
                using ElementTraits = typename Elements::template Type<dim_euclidean>::template Type<I>;
                auto const constexpr elt = ElementTraits::element;
    //            Indx line_start = mesh_file.getLineIndex("$Elements");
                Indx line_start = array::index(mesh_file, "$Elements");
                Indx offset = 1;
                StrgStream line_stream(mesh_file.get(line_start + offset));
                Indx num_entity_blocks;
                Indx num_elements;
                Indx min_element_tag;
                Indx max_element_tag;
                line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
                offset += 1;
                for (Indx i = 0; i < num_entity_blocks; ++i) {
                    line_stream = StrgStream(mesh_file.get(line_start + offset));
                    Indx entity_dim = 0;
                    Indx entity_tag = 0;
                    Indx element_type_tag = 0;
                    Indx num_elements_in_block = 0;
                    line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                    offset += 1;
                    elt::Element const constexpr element_type = ElementTraits::element;
                    elt::Element const element_type_arg = getElementType(element_type_tag);
                    if (entity_dim == dim_euclidean and getElementType(element_type_tag) == element_type) {
                        for (Indx j = 0; j < num_elements_in_block; ++j) {
                            line_stream = StrgStream(mesh_file.get(line_start + offset));
                            Indx tag;
                            line_stream >> tag;
                            Array<Indx, element_type.num_nodes> node_tags;
                            for (Indx k = 0; k < element_type.num_nodes; ++k) {
                                line_stream >> node_tags.get(k);
                            }
                            typename EEl<elt>::Components cell_boundary;
                            SharedPointer<EEl<elt>> ptr_cell = SharedPointer<EEl<elt>>(EEl<elt>());
                            this->template setElement<elt>(node_tags, cell_boundary, ptr_cell);
    //                        this->template setElement<CurrentElement::type>(node_tags, cell_boundary);
                            offset += 1;
                        }
                    } else {
                        offset += num_elements_in_block;
                    }
                }
            }

        };

//        explicit
//        ParserHelper1(
//                file::File const &
//                mesh_file
//        )
//                :
//                geometrical_entities(getGeometricalEntities(mesh_file)),
//                physical_entities(getPhysicalEntities(mesh_file)),
//                physical_groups(getPhysicalGroups(mesh_file))
//        {}
//
//        Bool
//        operator==(
//                ParserHelper1 const &
//                other
//        )
//        const = default;
//
//        Bool
//        operator!=(
//                ParserHelper1 const &
//                other
//        )
//        const = default;
//
//        struct PhysicalEntity
//        {
//
//            Bool
//            operator==(
//                    PhysicalEntity const &
//                    other
//            )
//            const = default;
//
//            Bool
//            operator!=(
//                    PhysicalEntity const &
//                    other
//            )
//            const = default;
//
//            Indx tag;
//
//            Indx dim;
//
//            Strg name;
//
//        };
//
//        struct GeometricalEntity
//        {
//
//            Bool
//            operator==(
//                    GeometricalEntity const &
//                    other
//            )
//            const = default;
//
//            Bool
//            operator!=(
//                    GeometricalEntity const &
//                    other
//            )
//            const = default;
//
//            Indx tag;
//
//            Indx dim;
//
//            Array<Indx> physical_entities_tags;
//
//            Array<Indx> bounding_entities_tags;
//
//        };
//
//        struct PhysicalGroup
//        {
//
//            Bool
//            operator==(
//                    PhysicalGroup const &
//                    other
//            )
//            const = default;
//
//            Bool
//            operator!=(
//                    PhysicalGroup const &
//                    other
//            )
//            const = default;
//
//            Strg name;
//
//            Array<Array<Indx>, 4> geometrical_entities_tags;
//
//        };
//
//        using PhysicalGroups = Array<PhysicalGroup>;
//        using GeometricalEntities = Array<Array<GeometricalEntity>, 4>;
//        using PhysicalEntities = Array<PhysicalEntity>;
//
//    private:
//
//        static
//        GeometricalEntities
//        getGeometricalEntities(
//                file::File const & mesh_file
//        )
//        {
//            GeometricalEntities geometrical_entities;
////            Indx line_start = mesh_file.getLineIndex("$Entities");
//            Indx line_start = array::index(mesh_file, "$Entities");
//            Indx offset = 1;
//            StrgStream line_stream(mesh_file.get(line_start + offset));
//            Indx num_points;
//            Indx num_curves;
//            Indx num_surfaces;
//            Indx num_volumes;
//            line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
//            Array<Indx, 4> num_domains = {
//                    num_points,
//                    num_curves,
//                    num_surfaces,
//                    num_volumes
//            };
//            offset += 1;
//            for (Indx i = 0; i < 4; ++i) {
//                for (Indx j = 0; j < num_domains.get(i); ++j) {
//                    line_stream = StrgStream(mesh_file.get(line_start + offset));
//                    Indx tag;
//                    line_stream >> tag;
//                    if (i == 0) {
//                        for (Indx k = 0; k < 3; ++k) {
//                            Real a;
//                            line_stream >> a;
//                        }
//                    } else {
//                        for (Indx k = 0; k < 6; ++k) {
//                            Real a;
//                            line_stream >> a;
//                        }
//                    }
//                    Indx num_physical_entities;
//                    line_stream >> num_physical_entities;
//                    Array<Indx> physical_entities_tags;
//                    for (Indx k = 0; k < num_physical_entities; ++k) {
//                        Intg physical_entity_tag;
//                        line_stream >> physical_entity_tag;
//                        physical_entities_tags.data.push_back(physical_entity_tag);
//                    }
//                    Array<Indx> bounding_entities_tags = {};
//                    if (i > 0) {
//                        Indx num_bounding_entities;
//                        line_stream >> num_bounding_entities;
//                        for (Indx k = 0; k < num_bounding_entities; ++k) {
//                            Intg bounding_entity_tag;
//                            line_stream >> bounding_entity_tag;
//                            bounding_entities_tags.data.push_back(std::abs(bounding_entity_tag));
//                        }
//                    }
//                    geometrical_entities.get(i).data.push_back(GeometricalEntity{
//                            tag,
//                            i,
//                            physical_entities_tags,
//                            bounding_entities_tags
//                    });
//                    offset += 1;
//                }
//            }
//            return geometrical_entities;
//        }
//
//        static
//        PhysicalEntities
//        getPhysicalEntities(
//                file::File const & mesh_file
//        )
//        {
//            PhysicalEntities physical_entities;
////            Indx line_start = mesh_file.getLineIndex("$PhysicalNames");
//            Indx line_start = array::index(mesh_file, "$PhysicalNames");
//            Indx offset = 1;
//            StrgStream line_stream(mesh_file.get(line_start + offset));
//            Indx num_physical_names;
//            line_stream >> num_physical_names;
//            offset += 1;
//            for (Indx i = 0; i < num_physical_names; ++i) {
//                line_stream = StrgStream(mesh_file.get(line_start + offset));
//                Indx dim;
//                Indx tag;
//                Strg name;
//                line_stream >> dim >> tag >> name;
//                file::removeCharacter(name, '"');
//                physical_entities.data.push_back(PhysicalEntity{tag, dim, name});
//                offset += 1;
//            }
//            return physical_entities;
//        }
//
//        Array<UnorderedSet<Indx>, 4>
//        getSubGeometricalEntities(
//                Indx const & d,
//                Indx const & t,
//                Array<UnorderedSet<Indx>, 4> & a
//        )
//        {
//            a.get(d).data.insert(t);
//            for (Indx i = 0; i < geometrical_entities.get(d).get(t - 1).bounding_entities_tags.size(); ++i) {
//                Indx const d2 = d - 1;
//                Indx const t2 = geometrical_entities.get(d).get(t - 1).bounding_entities_tags.get(i);
//                a = getSubGeometricalEntities(d2, t2, a);
//            }
//            return a;
//        }
//
//        PhysicalGroups
//        getPhysicalGroups(
//                file::File const & mesh_file
//        )
//        {
//            PhysicalGroups groups;
//            for (Indx i = 0; i < physical_entities.size(); ++i) {
//                Array<UnorderedSet<Indx>, 4> tags;
//                for (Indx j = 0; j < 4; ++j) {
//                    for (Indx k = 0; k < geometrical_entities.get(j).size(); ++k) {
//                        for (Indx l = 0; l < geometrical_entities.get(j).get(k).physical_entities_tags.size(); ++l) {
//                            Indx const & t0 = geometrical_entities.get(j).get(k).physical_entities_tags.get(l);
//                            Indx const & t1 = physical_entities.get(i).tag;
//                            if (t0 == t1) {
//                                tags = getSubGeometricalEntities(j, k + 1, tags);
//                            }
//                        }
//                    }
//                }
//                Array<Array<Indx>, 4> group_tags;
//                for (Indx j = 0; j < 4; ++j) {
//                    group_tags.get(j).data.assign(tags.get(j).data.begin(), tags.get(j).data.end());
//                }
//                groups.data.push_back(PhysicalGroup{physical_entities.get(i).name, group_tags});
//            }
//            return groups;
//        }
//
//    public:
//
//        GeometricalEntities geometrical_entities;
//        PhysicalEntities physical_entities;
//        PhysicalGroups physical_groups;

    };

    template<typename M>
    static inline
    void
    printMesh(
            M const &
            mesh_arg
    )
    {
        auto elsets = mesh_arg.element_sets;
        auto elems = mesh_arg.element_collection;
        print("START ELEMS");
        print("NODES");
        for (auto const & [key, val] : elems.template get<0>().template get<0>().data) {
            print(key, val.get().tag);
            for (auto const & ghj: val.get().neighbours.template get<0>().template get<0>().data) {
                print(" --> ", ghj.get().getHash());
            }
            for (auto const & ghj: val.get().neighbours.template get<1>().template get<0>().data) {
                print(" ----> ", ghj.get().getHash());
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
                print("--> ", val.get().components.template get<0>().template get<0>().get(i).get().getHash(), val.get().components.template get<0>().template get<0>().get(i).orientation);
            }
            for (auto const & ghj: val.get().neighbours.template get<1>().template get<0>().data) {
                print(" ----> ", ghj.get().getHash());
            }
            for (auto const & ghj: val.get().neighbours.template get<0>().template get<0>().data) {
                print(" <-> ", ghj.get().getHash());
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
                print("--> ", val.get().components.template get<0>().template get<0>().get(i).get().getHash(), val.get().components.template get<0>().template get<0>().get(i).orientation);
            }
            for (int i = 0; i < val.get().components.template get<1>().template get<0>().size(); ++i) {
                print("----> ", val.get().components.template get<1>().template get<0>().get(i).get().getHash(), val.get().components.template get<1>().template get<0>().get(i).orientation);
            }
            for (auto const & ghj: val.get().neighbours.template get<0>().template get<0>().data) {
                print(" <-> ", ghj.get().getHash());
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
                print("--> ", val.get().components.template get<0>().template get<0>().get(i).get().getHash(), val.get().components.template get<0>().template get<0>().get(i).orientation);
            }
            for (int i = 0; i < val.get().components.template get<1>().template get<0>().size(); ++i) {
                print("----> ", val.get().components.template get<1>().template get<0>().get(i).get().getHash(), val.get().components.template get<1>().template get<0>().get(i).orientation);
            }
            for (auto const & ghj: val.get().neighbours.template get<0>().template get<0>().data) {
                print(" <-> ", ghj.get().getHash());
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
