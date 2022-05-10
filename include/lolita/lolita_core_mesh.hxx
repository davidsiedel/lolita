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

    template<Element E, lolita::geometry::Domain D, auto Med>
    using FiniteElementE = core::element::FiniteElementF<E, D, Med>;

    template<lolita::mesh::MeshFormatType, lolita::geometry::Domain D, auto Med>
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

    template<lolita::mesh::MeshFormatType Mft, lolita::geometry::Domain D, auto Med>
    struct MeshBase
    {

    private:

        template<Element __element, lolita::geometry::Domain __domain, auto __arg>
        using ElementPointerMap = std::unordered_map<std::basic_string<lolita::character>, std::shared_ptr<FiniteElementE<__element, __domain, __arg>>>;

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

        using DegreeOfFreedomIndices = std::array<DegreeOfFreedomIndex, Med.size()>;

        using Loads = std::vector<lolita::finite_element::Load>;

        using Behaviours = std::vector<lolita::behaviour::MgisBehaviour>;

        std::shared_ptr<lolita::finite_element::LoadComponent> const null_load_;

    public:

        using ElementPointers = lolita::core::element::ElementCollection<ElementPointerMap, D, Med>;

        using ElementPointersSets = std::unordered_map<std::basic_string<lolita::character>, ElementPointers>;

        DegreeOfFreedomIndices unknown_indices;

        Module module_;

        ElementPointers elements_;

        ElementPointersSets element_sets_;

        Loads loads_;

        Behaviours behaviours_;

        MeshBase(
                lolita::utility::File const &
                mesh_file,
                Loads const &
                loads_arg,
                Behaviours const &
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
//            makeNeighbourhood();
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

//        void
//        makeNode(
//                lolita::numerics::RealConcept auto
//                node_tag_arg,
//                auto const &
//                node_coordinates_arg
//        )
//        {
//            auto hash = std::to_string(node_tag_arg);
////            auto ptr_node = std::shared_ptr<FiniteElementE<element::pnt_00, D, Med>>(FiniteElementE<element::pnt_00, D, Med>());
//            auto ptr_node = std::make_shared<FiniteElementE<element::pnt_00, D, Med>>(FiniteElementE<element::pnt_00, D, Med>());
////            auto ptr_coordinates = std::shared_ptr<lolita::geometry::Point>(node_coordinates_arg);
//            auto ptr_coordinates = std::make_shared<lolita::geometry::Point>(node_coordinates_arg);
//            /*
//             *
//             */
//            auto make_node = [&] <auto K = 0> (auto & self) {
////                pointer::make(ptr_node.get().template get<K>());
////                std::get<K>(* ptr_node) = std::make_shared<>();
//                ptr_node.get().template get<K>().get().coordinates = ptr_coordinates;
////                if constexpr (K < FiniteElementE<element::pnt_00, D, Med>::size() - 1) {
////                    self.template operator()<K + 1>(self);
////                }
//            };
//            /*
//             *
//             */
////            make_node(make_node);
////            auto & nodes = elements_.template get<0>().template get<0>();
//            auto & nodes = std::get<0>(std::get<0>(elements_));
//            ptr_node->tag = nodes.size();
//            ptr_node->coordinates_ = ptr_coordinates;
//            nodes.insert({hash, ptr_node});
//        }

        void
        makeNode2(
                lolita::numerics::RealConcept auto
                node_tag_arg,
                auto const &
                node_coordinates_arg
        )
        {
            auto hash = std::to_string(node_tag_arg);
            auto ptr_node = std::make_shared<FiniteElementE<element::pnt_00, D, Med>>(FiniteElementE<element::pnt_00, D, Med>());
            auto ptr_coordinates = std::make_shared<lolita::geometry::Point>(node_coordinates_arg);
//            auto & nodes = std::get<0>(std::get<0>(elements_));
            auto & nodes = elements_.template getElements<0, 0>();
            ptr_node->tag = nodes.size();
            ptr_node->coordinates_ = ptr_coordinates;
            nodes.insert({hash, ptr_node});
        }

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

        template<Element EE>
        void
        makeElement2(
                std::ranges::range auto const &
                cell_node_tags_arg
        )
        {
            /*
             *
             */
            auto get_element_domains = [&] <Element __element> (
                    std::array<lolita::index, __element.num_nodes> const &
                    node_tags
            )
            {
                std::vector<std::basic_string<lolita::character>> domain_names;
                for (auto const & [set_name, element_tags] : element_sets_) {
//                    auto const & node_hashes = std::get<0>(std::get<0>(element_tags));
                    auto const & node_hashes = element_tags.template getElements<0, 0>();
                    lolita::boolean belongs = true;
                    for (auto node_tag : node_tags) {
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
            };
            /*
             *
             */
            auto get_component_node_tags = [] <Element __element, auto _i, auto _j> (
                    std::array<lolita::index, __element.num_nodes> const &
                    node_tags,
                    lolita::index
                    index
            )
            {
                auto component_node_tags = std::array<lolita::index, element::neighbour<__element, 0, _i, _j>().num_nodes>();
                for (lolita::index j = 0; j < element::neighbour<__element, 0, _i, _j>().num_nodes; ++j) {
//                    auto const & element_node_connectivity = FiniteElementE<__element, D, Med>::node_connectivity;
//                    auto const & bound_node_connectivity = std::get<_j>(std::get<_i>(element_node_connectivity));
//                    auto const k = bound_node_connectivity[index][j];
                    auto const k = FiniteElementE<__element, D, Med>::template getNeighboursNodeConnectivity<_i, _j>(index, j);
                    component_node_tags[j] = node_tags[k];
                }
                return component_node_tags;
            };
            /*
             *
             */
            auto get_element_hash = [] <Element __element> (
                    std::array<lolita::index, __element.num_nodes>
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
            auto set_element_n = [&] <Element __element, lolita::index _i = 0, lolita::index _j = 0> (
                    std::array<lolita::index, __element.num_nodes> const &
                    element_node_tags_arg,
                    std::shared_ptr<FiniteElementE<__element, D, Med>> &
                    ptr_element_arg,
                    auto &
                    self
            )
            mutable
            {
                if constexpr (!element::PointConcept<__element>) {
                    auto const constexpr _neighbour = element::neighbour<__element, 0, _i, _j>();
                    auto const constexpr _pos_element = element::elementPosition<__element>();
                    auto const constexpr _pos_neighbour = element::elementPosition<_neighbour>();
//                    auto & components = std::get<_pos_neighbour[1]>(std::get<_pos_neighbour[0]>(elements_));
//                    auto & elements = std::get<_pos_element[1]>(std::get<_pos_element[0]>(elements_));
                    auto & components = elements_.template getElements<_pos_neighbour[0], _pos_neighbour[1]>();
                    auto & elements = elements_.template getElements<_pos_element[0], _pos_element[1]>();
                    auto & element_component_array = ptr_element_arg->template getNeighbours<0, _i, _j>();
//                    auto & element_component_array = std::get<_j>(std::get<_i>(std::get<0>(ptr_element_arg->neighbours)));
                    auto element_hash = get_element_hash.template operator ()<__element>(element_node_tags_arg);
                    for (auto i = 0; i < element::numNeighbours<__element, 0, _i, _j>(); ++i) {
                        auto component_node_tags = get_component_node_tags.template operator ()<__element, _i, _j>(element_node_tags_arg, i);
                        auto component_hash = get_element_hash.template operator ()<_neighbour>(component_node_tags);
                        if (!components.contains(component_hash)) {
                            auto ptr_component = std::make_shared<FiniteElementE<_neighbour, D, Med>>(FiniteElementE<_neighbour, D, Med>());
                            self.template operator ()<_neighbour>(component_node_tags, ptr_component, self);
                        }
                        element_component_array[i] = components[component_hash];
                        auto const constexpr poss = element::neighbourPosition<_neighbour, __element>();
                        components[component_hash]->template getNeighbours<poss[0], poss[1], poss[2]>().push_back(ptr_element_arg);
//                        std::get<poss[2]>(std::get<poss[1]>(std::get<poss[0]>(components[component_hash]->neighbours))).push_back(ptr_element_arg);
//                        components[component_hash]->template getNeighbours<poss[0], poss[1], poss[2]>().push_back(ptr_element_arg);
                    }
                    if constexpr (_j < element::numNeighbours<__element, 0, _i>() - 1) {
                        self.template operator()<__element, _i, _j + 1>(element_node_tags_arg, ptr_element_arg, self);
                    }
                    else if constexpr (_i < element::numNeighbours<__element, 0>() - 1) {
                        self.template operator()<__element, _i + 1, 0>(element_node_tags_arg, ptr_element_arg, self);
                    }
                    if constexpr (_i == 0 && _j == 0) {
                        ptr_element_arg->tag = elements.size();
                        auto element_domains = get_element_domains.template operator()<__element>(element_node_tags_arg);
                        elements.insert({element_hash, ptr_element_arg});
                    }
                }
            };
            /*
             *
             */
            auto ptr_cell = std::make_shared<FiniteElementE<EE, D, Med>>(FiniteElementE<EE, D, Med>());
            set_element_n.template operator ()<EE>(cell_node_tags_arg, ptr_cell, set_element_n);
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
                static_cast<Implementation *>(this)->template makeCells<element::element<D.dim_, I>()>();
                if constexpr (I < element::numElements<D.dim_>() - 1) {
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
            auto make_element_neighbourhood = [&] <Element E, auto K = 0> (auto & self)
            mutable
            {
//                auto const constexpr ed = E.dim;
//                auto const constexpr et = element::elementIndex<E>();
                auto const constexpr poss = element::elementPosition<E>();
//                auto & element_map = elements_.template get<ed>().template get<et>();
                auto & element_map = std::get<poss[1]>(std::get<poss[0]>(elements_));//.template get<ed>().template get<et>();
                for (auto & element_map_item: element_map) {
                    auto & element = element_map_item.second;
//                    auto const & element_nds = std::get<0>(std::get<ed - 1>(element->neighbours));//.template get<ed - 1>().template get<0>();
                    auto const & element_nds = element->template getNeighbours<0, poss[0] - 1, 0>();//.template get<ed - 1>().template get<0>();
//                    auto const & element_nds = std::get<0>(std::get<poss[0] - 1>(std::get<0>(element->neighbours)));//.template get<ed - 1>().template get<0>();
                    for (int i = 0; i < element_nds.size(); ++i) {
                        auto const & nde = element_nds[i];
//                        auto const & ngs = nde.get().neighbours.template get<ed - 1>().template get<K>();
                        auto const & ngs = nde->template getNeighbours<1, poss[0] - 1, K>();
//                        auto const & ngs = std::get<K>(std::get<poss[0] - 1>(std::get<1>(nde->neighbours)));
                        for (auto const & neighbour: ngs) {
                            if (((neighbour->tag != element->tag) && K == poss[1]) || (K != poss[1])) {
//                                auto & element_ngs = element.get().neighbours.template get<0>().template get<K>();
                                auto & element_ngs = element->template getNeighbours<1, 0, K>();
//                                auto & element_ngs = std::get<K>(std::get<0>(std::get<1>(element->neighbours)));
                                lolita::boolean found = false;
                                for (auto const & ngb: element_ngs) {
                                    if (ngb->tag == neighbour->tag) {
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
                if constexpr (K < element::numElements<poss[0]>() - 1) {
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
                else if constexpr (I < D.dim_) {
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
//                    for (auto const & c_ : elem_arg.get().components.template get<I>().template get<J>().data) {
                    for (auto const & c_ : elem_arg->template getNeighbours<0, I, J>()) {
                        os << "layer : " << I << " type : " << J << " <-- " << c_->hash() << std::endl;
                    }
                    if constexpr (J < element::numNeighbours<EE, 0, I>() - 1) {
                        self.template operator()<EE, I, J + 1>(elem_arg, self);
                    }
                    else if constexpr (I < element::numNeighbours<EE, 0>() - 1) {
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
//                for (auto const & c_ : elem_arg.get().neighbours.template get<I>().template get<J>().data) {
                for (auto const & c_ : elem_arg->template getNeighbours<1, I, J>()) {
                    if constexpr (!element::PointConcept<EE> && I == 0) {
                        os << "layer : " << I << " type : " << J << " <-> " << c_->hash() << std::endl;
                    }
                    else {
                        os << "layer : " << I << " type : " << J << " --> " << c_->hash() << std::endl;
                    }
                }
                if constexpr (J < element::numNeighbours<EE, 0, D.dim_, I>() - 1) {
                    self.template operator()<EE, I, J + 1>(elem_arg, self);
                }
                else if constexpr (I < element::numNeighbours<EE, 0, D.dim_>() - 1) {
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
//                for (auto const & item : set_arg.template get<I>().template get<J>().data) {
                for (auto const & item : std::get<J>(std::get<I>(set_arg))) {
                    os << "layer : " << I << " type : " << J << " <-- " << item.second->hash() << std::endl;
                }
                if constexpr (J < element::numElements<I>() - 1) {
                    self.template operator()<I, J + 1>(set_arg, self);
                }
                else if constexpr (I < D.dim_) {
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
            auto print_elements = [&] <auto I = 0, auto J = 0> (
                    auto & self
            )
            mutable
            {
                auto constexpr elt = element::element<I, J>();
                if constexpr (I == 0 && J == 0) {
                    os << "*** Elements : " << std::endl;
                }
//                for (auto const & element : base.elements_.template get<I>().template get<J>().data) {
                for (auto const & element : std::get<J>(std::get<I>(base.elements_))) {
                    os << "* Element : " << element.second->hash() << std::endl;
                    print_element_components.template operator()<elt>(element.second, print_element_components);
                    print_element_neighbours.template operator()<elt>(element.second, print_element_neighbours);
                }
                if constexpr (J < element::numElements<I>() - 1) {
                    self.template operator()<I, J + 1>(self);
                }
                else if constexpr (I < D.dim_) {
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

    template<lolita::geometry::Domain D, auto Med>
    struct MeshModule<lolita::mesh::MeshFormatType::Gmsh, D, Med>
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

                lolita::index tag;

                lolita::index dim;

                std::basic_string<lolita::character> name;

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

                lolita::index tag;

                lolita::index dim;

                std::vector<lolita::index> physical_entities_tags;

                std::vector<lolita::index> bounding_entities_tags;

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

                std::basic_string<lolita::character> name;

                std::array<std::vector<lolita::index>, 4> geometrical_entities_tags;

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
                for (lolita::index i = 0; i < geometrical_entities[d][t - 1].bounding_entities_tags.size(); ++i) {
                    lolita::index const d2 = d - 1;
                    lolita::index const t2 = geometrical_entities[d][t - 1].bounding_entities_tags[i];
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
                            for (lolita::index l = 0; l < geometrical_entities[j][k].physical_entities_tags.size(); ++l) {
                                lolita::index const & t0 = geometrical_entities[j][k].physical_entities_tags[l];
                                lolita::index const & t1 = physical_entities[i].tag;
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
                    physical_groups.push_back(PhysicalGroup{physical_entities[i].name, group_tags});
                }
                //return groups;
            }

        public:

            GeometricalEntities geometrical_entities;
            PhysicalEntities physical_entities;
            PhysicalGroups physical_groups;

        };

        struct Implementation : public MeshBase<lolita::mesh::MeshFormatType::Gmsh, D, Med>
        {

            //using Base = MeshBase<lolita::mesh::MeshFormatType::Gmsh, D, Med>;

            void
            makeNodes()
            {
//                auto & nodes = this->elements_.template get<0>().template get<0>();
//                auto & nodes = std::get<0>(std::get<0>(this->elements_));
                auto & nodes = this->elements_.template getElements<0, 0>();
                lolita::index line_start = std::distance(this->module_.file_.lines_.begin(), std::find(this->module_.file_.lines_.begin(), this->module_.file_.lines_.end(), "$Nodes"));
//                lolita::index line_start = array::index(this->module_.file_, "$Nodes");
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(this->module_.file_.lines_[line_start + offset]);
                lolita::index num_entity_block;
                lolita::index num_nodes;
                lolita::index min_node_tag;
                lolita::index max_node_tag;
                line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                offset += 1;
                for (lolita::index i = 0; i < num_entity_block; ++i) {
//                    line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                    line_stream = std::basic_stringstream<lolita::character>(this->module_.file_.lines_[line_start + offset]);
                    lolita::index entity_dim;
                    lolita::index entity_tag;
                    lolita::index parametric;
                    lolita::index num_nodes_in_block;
                    line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                    offset += 1;
                    for (lolita::index j = 0; j < num_nodes_in_block; ++j) {
//                        line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                        line_stream = std::basic_stringstream<lolita::character>(this->module_.file_.lines_[line_start + offset]);
                        lolita::natural node_tag;
                        line_stream >> node_tag;
//                        line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset + num_nodes_in_block));
                        line_stream = std::basic_stringstream<lolita::character>(this->module_.file_.lines_[line_start + offset + num_nodes_in_block]);
                        lolita::geometry::Point coordinates = lolita::geometry::Point::Zero();
                        for (int k = 0; k < D.dim_; ++k) {
                            line_stream >> coordinates(k);
                        }
//                        this->setElement(node_tag, coordinates);
//                        this->makeNode(node_tag, coordinates);
                        this->makeNode2(node_tag, coordinates);
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
                    this->element_sets_[physical_groups[i].name];
                }
            }

            void
            makeNodeSets()
            {
                auto const & physical_groups = this->module_.physical_groups;
                for (int m = 0; m < physical_groups.size(); ++m) {
                    auto & node_set_name = physical_groups[m].name;
//                    auto & node_set = std::get<0>(std::get<0>(this->element_sets_[node_set_name]));
                    auto & node_set = this->element_sets_[node_set_name].template getElements<0, 0>();
//                    auto & nodes = std::get<0>(std::get<0>(this->elements_));
                    auto & nodes = this->elements_.template getElements<0, 0>();
                    for (lolita::index i = 0; i < 4; ++i) {
                        for (lolita::index j = 0; j < physical_groups[m].geometrical_entities_tags[i].size(); ++j) {
                            lolita::index tag = physical_groups[m].geometrical_entities_tags[i][j];
//                            lolita::index line_start = array::index(this->module_.file_, "$Nodes");
                            lolita::index line_start = std::distance(this->module_.file_.lines_.begin(), std::find(this->module_.file_.lines_.begin(), this->module_.file_.lines_.end(), "$Nodes"));
                            lolita::index offset = 1;
//                            std::basic_string<lolita::character>Stream line_stream(this->module_.file_.get(line_start + offset));
                            std::basic_stringstream<lolita::character> line_stream(this->module_.file_.lines_[line_start + offset]);
                            lolita::index num_entity_block;
                            lolita::index num_nodes;
                            lolita::index min_node_tag;
                            lolita::index max_node_tag;
                            line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                            offset += 1;
                            for (lolita::index k = 0; k < num_entity_block; ++k) {
//                                line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                                line_stream = std::basic_stringstream<lolita::character>(this->module_.file_.lines_[line_start + offset]);
                                lolita::index entity_dim;
                                lolita::index entity_tag;
                                lolita::index parametric;
                                lolita::index num_nodes_in_block;
                                line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                                offset += 1;
                                if (entity_dim == i && entity_tag == tag) {
                                    for (lolita::index l = 0; l < num_nodes_in_block; ++l) {
//                                        line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                                        line_stream = std::basic_stringstream<lolita::character>(this->module_.file_.lines_[line_start + offset]);
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

            template<Element EE>
            void
            makeCells()
            {
//                lolita::index line_start = array::index(this->module_.file_, "$Elements");
                lolita::index line_start = std::distance(this->module_.file_.lines_.begin(), std::find(this->module_.file_.lines_.begin(), this->module_.file_.lines_.end(), "$Elements"));
                lolita::index offset = 1;
//                std::basic_string<lolita::character>Stream line_stream(this->module_.file_.get(line_start + offset));
                std::basic_stringstream<lolita::character> line_stream(this->module_.file_.lines_[line_start + offset]);
                lolita::index num_entity_blocks;
                lolita::index num_elements;
                lolita::index min_element_tag;
                lolita::index max_element_tag;
                line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
                offset += 1;
                for (lolita::index i = 0; i < num_entity_blocks; ++i) {
//                    line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                    line_stream = std::basic_stringstream<lolita::character>(this->module_.file_.lines_[line_start + offset]);
                    lolita::index entity_dim = 0;
                    lolita::index entity_tag = 0;
                    lolita::index element_type_tag = 0;
                    lolita::index num_elements_in_block = 0;
                    line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                    offset += 1;
                    Element const element_type_arg = getElementType(element_type_tag);
                    if (entity_dim == D.dim_ and getElementType(element_type_tag) == EE) {
                        for (lolita::index j = 0; j < num_elements_in_block; ++j) {
//                            line_stream = std::basic_string<lolita::character>Stream(this->module_.file_.get(line_start + offset));
                            line_stream = std::basic_stringstream<lolita::character>(this->module_.file_.lines_[line_start + offset]);
                            lolita::index tag;
                            line_stream >> tag;
                            std::array<lolita::index, EE.num_nodes> node_tags;
                            for (lolita::index k = 0; k < EE.num_nodes; ++k) {
                                line_stream >> node_tags[k];
                            }
//                            this->template makeElement<EE>(node_tags);
                            this->template makeElement2<EE>(node_tags);
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
