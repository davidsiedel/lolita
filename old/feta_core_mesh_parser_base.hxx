//
// Created by dsiedel on 15/03/2022.
//

#ifndef FETA_FETA_CORE_MESH_PARSER_BASE_HXX
#define FETA_FETA_CORE_MESH_PARSER_BASE_HXX

#include <execution>

//#include "feta/feta/static_array.hxx"
//#include "feta/feta/dynamic_array.hxx"
//#include "feta/feta/unordered_map.hxx"
//#include "feta/feta/unordered_set.hxx"
//#include "feta/feta/collection.hxx"
//#include "feta/feta/aliases.hxx"
//#include "feta/feta/matrix.hxx"
//#include "feta/feta/files.hxx"
////#include "feta/feta/functions.hxx"
//#include "feta/core/finite_element_connectivity.hxx"

#include "new/feta_core_element_element.hxx"
#include "new/_feta_unordered_map.hxx"

namespace feta::core
{

    template<auto Med>
    struct ParserBase
    {

    protected:

        template<elt::ElementDescription E>
        using Element = elt::Element<E, Med>;

        template<elt::ElementDescription E, Indx I, Indx J>
        using ComponentPointer = typename Element<E>::Components::template Type<I>::template Type<J>::Type;

        template<elt::ElementDescription E, Indx I, Indx J>
        using Component = typename ComponentPointer<E, I, J>::Type;

        template<elt::ElementDescription E>
        using ElementPointerMap = UnorderedMap<Strg, SharedPointer<Element<E>>>;

    public:

        using ElementPointers = core::Elements<Med.dim_euclidean, ElementPointerMap>;

        using ElementPointersSets = UnorderedMap<Strg, ElementPointers>;

        using DomainPointers = Array<SharedPointer<MeshInteriorDomain<Med>>>;

        ElementPointers element_collection;

        ElementPointersSets element_sets;

        DomainPointers mesh_domains;

        Bool
        operator==(
                ParserBase const &
                other
        )
        const
        {
            auto eq_0 = element_collection == other.element_collection;
            auto eq_1 = element_sets == other.element_sets;
            auto eq_2 = mesh_domains == other.mesh_domains;
            return eq_0 && eq_1 && eq_2;
        }

        Bool
        operator!=(
                ParserBase const &
                other
        )
        const
        {
            return !(other == * this);
        }

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

        template<elt::ElementDescription E>
        static
        auto
        getElementHash(
                Array<Indx, Element<E>::num_nodes>
                node_tags_arg
        )
        {
            StrgStream element_hash;
            if (Element<E>::num_nodes > 1) {
                std::sort(std::execution::par_unseq, node_tags_arg.data.begin(), node_tags_arg.data.end());
            }
            for (Indx i = 0; i < Element<E>::num_nodes; ++i) {
                element_hash << node_tags_arg(i);
            }
            return element_hash.str();
        }

        template<elt::ElementDescription E, Indx I, Indx J>
        static
        auto
        getComponentNodeTags(
                Array<Indx, Element<E>::num_nodes> const &
                node_tags_arg,
                Indx
                component_index_arg
        )
        {
            Array<Indx, Component<E, I, J>::num_nodes> bound_node_tags;
            for (Indx j = 0; j < Component<E, I, J>::num_nodes; ++j) {
//                Indx k = Element<E>::node_connectivity.template get<I>().template get<J>()(component_index_arg)(j);
                auto const & ndscon = Element<E>::node_connectivity.template get<I>().template get<J>();
                auto const k = ndscon.get(component_index_arg, j);
//                Indx k = Element<E>::node_connectivity.template get<I>().template get<J>().get(component_index_arg, j);
                bound_node_tags(j) = node_tags_arg(k);
            }
            return bound_node_tags;
        }

        template<elt::ElementDescription E>
        auto
        getElementDomains(
                Array<Indx, Element<E>::num_nodes> const &
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
                Matrix<Real, Med.dim_euclidean> const &
                coordinates_arg
        )
        {
            auto & nodes = element_collection.template get<0>().template get<0>();
            auto hash = getElementHash<elt::pnt_0>({node_tag_arg});
            auto ptr_node = SharedPointer<Element<elt::pnt_0>>(Element<elt::pnt_0>());
            ptr_node.get().coordinates = coordinates_arg;
            ptr_node.get().tag = nodes.getSize();
            nodes.data.insert({hash, ptr_node});
        }

        template<elt::ElementDescription E, Indx I = 0, Indx J = 0>
        void
        setElement(
                Array<Indx, Element<E>::num_nodes> const &
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
            auto const constexpr cmp = Component<E, I, J>::element_description;
            auto const constexpr bd = Component<E, I, J>::element_description.shape_description.ord_shape;
            auto const constexpr bt = Component<E, I, J>::reference_tag;
            auto const constexpr ed = Element<E>::element_description.shape_description.ord_shape;
            auto const constexpr et = Element<E>::reference_tag;
//            using ElementComponentArray = typename Element<E>::Components::template Type<I>::template Type<J>;
//            using ElementComponent = typename Element<E>::template Component<component_description>;
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
            for (Indx i = 0; i < element_component_array.getSize(); ++i) {
                auto component_node_tags = getComponentNodeTags<E, I, J>(element_node_tags_arg, i);
                auto component_hash = getElementHash<Component<E, I, J>::element_description>(component_node_tags);
                if constexpr (Component<E, I, J>::element_description == elt::pnt_0) {
                    /*
                     * if the component element is a point, it is already defined, since all nodes in the mesh are
                     * built prior to building structure elements.
                     */
                    element_component_array.get(i) = ComponentPointer<E, I, J>(nodes.get(component_hash), 1);
                } else {
                    /*
                     * else, one needs to check either the component element has already been defined.
                     */
                    if (components.data.contains(component_hash)) {
                        /*
                         * if the component element has already been defined, one just fetches it from the existing
                         * ones.
                         */
                        element_component_array.get(i) = ComponentPointer<E, I, J>(components.get(component_hash), -1);
                    } else {
                        /*
                         * otherwise, one creates it by recursively calling the present function.
                         * A pointer to the component element is created, as well as a component element collection.
                         */
                        typename Element<cmp>::Components component_components;
                        auto ptr_component = SharedPointer<Element<cmp>>(Element<cmp>());
                        setElement<Component<E, I, J>::element_description>(
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
                if constexpr (Component<E, I, J>::element_description == elt::pnt_0) {
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
            if constexpr (J < Element<E>::Components::template Type<I>::getSize() - 1) {
                setElement<E, I, J + 1>(
                        element_node_tags_arg,
                        element_components_arg,
                        ptr_element_arg
                );
            }
            if constexpr(Intg(Element<E>::element_description.shape_description.ord_shape) - (1 + I) > 0) {
                setElement<E, I + 1, 0>(
                        element_node_tags_arg,
                        element_components_arg,
                        ptr_element_arg
                );
            }
            if constexpr((I == 0) && (J == 0)) {
                ptr_element_arg.get().tag = elements.getSize();
                ptr_element_arg.get().components = element_components_arg;
                elements.data.insert({element_hash, ptr_element_arg});
            }
        }

        template<Indx I = 1, Indx J = 0, Indx K = 0>
        void
        setNeighbours()
        {
            using Elementts = typename core::Elements<Med.dim_euclidean, Element>;
            using Elementt = typename Elementts::template Type<I>::template Type<J>;
            using ElementtNeighbour = typename Elementt::Neighbours::template Type<0>::template Type<K>::Type;
            /*
             * for each element pointer in the mesh
             */
            auto & elements_selection = element_collection.template get<I>().template get<J>();
            for (auto & [element_label, ptr_element]: elements_selection.data) {
                auto const & element_nodes = ptr_element.get().components.template get<I - 1>().template get<0>();
                for (Indx i = 0; i < element_nodes.getSize(); ++i) {
                    auto const & ptr_node = element_nodes(i);
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
            if constexpr (K < ElementPointers::template Type<I>::getSize() - 1) {
                setNeighbours<I, J, K + 1>();
            }
            if constexpr (J < ElementPointers::template Type<I>::getSize() - 1) {
                setNeighbours<I, J + 1, 0>();
            }
            if constexpr (I < Med.dim_euclidean) {
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
            if (node_neighbours.getSize() > 0) {
                is_isolated_arg = false;
            }
            if constexpr(I < ElementPointers::template Type<1>::getSize() - 1) {
                is_isolated_arg = isNodeIsolated<I + 1>(node_hash_arg, is_isolated_arg);
            }
            return is_isolated_arg;
        }

    };

}

#endif //FETA_FETA_CORE_MESH_PARSER_BASE_HXX
