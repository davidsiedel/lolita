#ifndef BC3CC377_0788_4FE2_B68B_684298E3EEEA
#define BC3CC377_0788_4FE2_B68B_684298E3EEEA

#include "lolita/lolita_core_n_11.hxx"

namespace lolita2::geometry
{
    
    template<Domain t_domain>
    struct Mesh
    {

    private:
    
        template<Element t_element, Domain t__domain>
        using _ElementPointerMap = std::unordered_map<std::basic_string<lolita::character>, std::shared_ptr<FiniteElement<t_element, t__domain>>>;

    public:
    
        friend
        std::ostream &
        operator<<(
                std::ostream & os,
                Mesh const & mesh
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
            
            auto print_elements = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
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
                    os << "* Element : " << t_element << " " << element.second->hash() << " is in domains : ";
                    for (auto const & element_domain : element.second->domains_) {
                        os << * element_domain << " ";
                    }
                    os << std::endl;
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
            print_elements(print_elements);
            return os;
        }
        
        ElementCollection<_ElementPointerMap, t_domain> elements_;
        
        std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;

    };
    
    template<MeshData _mesh, Domain t_domain>
    struct MeshParserModule;
    
    template<MeshData _mesh, Domain t_domain>
    struct MeshParser
    {

    private:

        using t_Module = typename MeshParserModule<_mesh, t_domain>::Module;

        using t_Implementation = typename MeshParserModule<_mesh, t_domain>::Implementation;

    public:
    
        MeshParser(
            std::basic_string_view<lolita::character> mesh_file_path
        )
        {
            
            auto const module = t_Module(mesh_file_path);
            makeDomains(module);
            auto set_elements = [&] <lolita::index t_i = 0u, lolita::index t_j = 0u> (
                auto & set_elements_imp
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                makeElement<t_element>(module);
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1) {
                    set_elements_imp.template operator()<t_i, t_j + 1u>(set_elements_imp);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1) {
                    set_elements_imp.template operator()<t_i + 1u, 0u>(set_elements_imp);
                }
            };
            set_elements(set_elements);
            auto make_elements = [&] <lolita::index t_i = 0u, lolita::index t_j = 0u> (
                auto & self
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                buildElementN<t_element>();
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1) {
                    self.template operator()<t_i, t_j + 1u>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1) {
                    self.template operator()<t_i + 1u, 0u>(self);
                }
            };
            make_elements(make_elements);
        }

        static
        Mesh<t_domain>
        makeMesh(
            std::basic_string_view<lolita::character> mesh_file_path
        )
        {
            return MeshParser(mesh_file_path).mesh_data_;
        }

    protected:

        template<Element t_element>
        static
        std::basic_string<lolita::character>
        getElementHash(
            std::array<lolita::index, t_element.num_nodes_> node_tags
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

        template<Element t_element>
        void
        makeElement(
            std::shared_ptr<FiniteElement<t_element, t_domain>> & ptr_element,
            std::array<lolita::index, t_element.num_nodes_> node_tags
        )
        requires(!t_element.isNode())
        {
            using t_ElementDescription = ElementTraits<t_element, t_domain>;
            using t_MeshDescription = DomainTraits<t_domain>;
            auto const constexpr t_element_coordinates = t_MeshDescription::template getElementCoordinates<t_element>();
            auto const constexpr t_node_coordinates = t_ElementDescription::template getComponentCoordinates<Element::node()>();
            auto & elements = mesh_data_.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
            auto const & nodes = ptr_element->template getComponents<t_node_coordinates.dim_, t_node_coordinates.tag_>();
            auto domains = std::unordered_set<std::shared_ptr<std::basic_string<lolita::character>>>();
            for (auto const & domain : nodes[0]->domains_)
            {
                auto has_domain = true;
                for (lolita::index j = 1; j < t_element.num_nodes_; ++j)
                {
                    has_domain = std::find(nodes[j]->domains_.begin(), nodes[j]->domains_.end(), domain) != nodes[j]->domains_.end();
                }
                if (has_domain)
                {
                    domains.insert(domain);
                }
            }
            ptr_element->tag_ = elements.size();
            ptr_element->domains_.assign(domains.begin(), domains.end());
            elements.insert({getElementHash<t_element>(node_tags), ptr_element});
        }

        template<Element t_element, lolita::integer t_i = 0, lolita::integer t_j = 0>
        void
        makeElement2(
            std::shared_ptr<FiniteElement<t_element, t_domain>> & ptr_element,
            std::array<lolita::index, t_element.num_nodes_> node_tags
        )
        requires(!t_element.isNode())
        {
            using t_Element = FiniteElement<t_element, t_domain>;
            using t_ElementDescription = ElementTraits<t_element, t_domain>;
            using t_ComponentDescription = ElementTraits<t_ElementDescription::template getComponent<t_i, t_j>(), t_domain>;
            using t_MeshDescription = DomainTraits<t_domain>;
            auto const constexpr _is_initialized = t_i == 0 && t_j == 0;
            auto const constexpr _component = t_ElementDescription::template getComponent<t_i, t_j>();
            auto const constexpr _component_coordinates = t_MeshDescription::template getElementCoordinates<_component>();
            auto const constexpr _neighbour_coordinates = t_ComponentDescription::template getNeighbourCoordinates<t_element>();
            using t_Component = FiniteElement<_component, t_domain>;
            auto & components = mesh_data_.elements_.template getElements<_component_coordinates.dim_, _component_coordinates.tag_>();
            auto & element_component_array = ptr_element->template getComponents<t_i, t_j>();
            for (auto i = 0; i < element_component_array.size(); ++i)
            {
                auto component_hash = std::basic_string<lolita::character>();
                if constexpr(!_component.isNode())
                {
                    auto component_node_tags = std::array<lolita::index, _component.num_nodes_>();
                    for (int j = 0; j < _component.num_nodes_; ++j)
                    {
                        auto const k = t_Element::template getComponentNodeConnection<t_i, t_j>(i, j);
                        component_node_tags[j] = node_tags[k];
                    }
                    component_hash = getElementHash<_component>(component_node_tags);
                    if (!components.contains(component_hash))
                    {
                        auto ptr_component = std::make_shared<t_Component>(t_Component());
                        makeElement2<_component, 0, 0>(ptr_component, component_node_tags);
                    }
                }
                else
                {
                    component_hash = std::to_string(node_tags[t_Element::template getComponentNodeConnection<t_i, t_j>(i, 0)]);
                }
                element_component_array[i] = components[component_hash];
                components[component_hash]->template getNeighbours<_neighbour_coordinates.dim_, _neighbour_coordinates.tag_>().push_back(ptr_element);
            }
            if constexpr (t_j < t_ElementDescription::template getNumComponents<t_i>() - 1)
            {
                makeElement2<t_element, t_i, t_j + 1u>(ptr_element, node_tags);
            }
            else if constexpr (t_i < t_ElementDescription::getNumComponents() - 1)
            {
                makeElement2<t_element, t_i + 1u, 0u>(ptr_element, node_tags);
            }
            if constexpr (_is_initialized)
            {
                makeElement<t_element>(ptr_element, node_tags);
            }
        }

        template<Element t_element>
        void
        makeElement2(
            std::array<lolita::index, t_element.num_nodes_> node_tags
        )
        requires(!t_element.isNode())
        {
            using t_Element = FiniteElement<t_element, t_domain>;
            auto ptr_element = std::make_shared<t_Element>(t_Element());
            makeElement2<t_element, 0, 0>(ptr_element, node_tags);
        }

        template<Element t_element>
        void
        makeElement2(
                lolita::natural const & tag,
                lolita::domain::Point const & coordinates,
                std::vector<std::shared_ptr<std::basic_string<lolita::character>>> const & domains
        )
        requires(t_element.isNode())
        {
            using t_Element = FiniteElement<t_element, t_domain>;
            auto ptr_element = std::make_shared<t_Element>(t_Element());
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            ptr_element->tag_ = tag;
            ptr_element->domains_ = domains;
            ptr_element->coordinates_ = std::make_shared<lolita::domain::Point>(coordinates);
            auto elem_hash = std::to_string(tag);
            mesh_data_.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>().insert({elem_hash, ptr_element});
        }
    
        void
        makeDomains(
            t_Module const & module
        )
        {
            static_cast<t_Implementation *>(this)->makeDomains(module);
        }
        
        template<Element t_element>
        void
        makeElement(
            t_Module const & module
        )
        {
            static_cast<t_Implementation *>(this)->template makeElement<t_element>(module);
        }

        template<Element t_element, lolita::integer _k = 0>
        void
        buildElementN()
        {
            auto const constexpr t_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto & element_map = mesh_data_.elements_.template getElements<t_coordinates.dim_, t_coordinates.tag_>();
            for (auto & element_map_item : element_map)
            {
                auto & element = element_map_item.second;
                if constexpr (!t_element.isNode())
                {
                    auto const & element_nds = element->template getComponents<t_coordinates.dim_ - 1, 0>();
                    for (int i = 0; i < element_nds.size(); ++i)
                    {
                        auto const & nde = element_nds[i];
                        auto const & ngs = nde->template getNeighbours<t_coordinates.dim_ - 1, _k>();
                        for (auto const & neighbour : ngs)
                        {
                            if (((neighbour->tag_ != element->tag_) && _k == t_coordinates.tag_) || (_k != t_coordinates.tag_))
                            {
                                auto & element_ngs = element->template getNeighbours<0, _k>();
                                auto found = false;
                                for (auto const & ngb: element_ngs)
                                {
                                    if (ngb->tag_ == neighbour->tag_)
                                    {
                                        found = true;
                                        break;
                                    }
                                }
                                if (!found)
                                {
                                    element_ngs.push_back(neighbour);
                                }
                            }
                        }
                    }
                    if constexpr (_k < DomainTraits<t_domain>::template getNumElements<t_coordinates.dim_>() - 1)
                    {
                        buildElementN<t_element, _k + 1>();
                    }
                }
            }
        }

    public:
    
        Mesh<t_domain> mesh_data_;

    };

    template<MeshData _mesh, Domain t_domain>
    requires(_mesh.isGmsh())
    struct MeshParserModule<_mesh, t_domain>
    {

        lolita::integer static constexpr node_index_offset = 1;

        static constexpr
        Element
        getElemType(
                lolita::index tag
        )
        {
            if (tag == 15) return Element::node();
            else if (tag == 01) return Element::segment(1);
            else if (tag == 02) return Element::triangle(1);
            else if (tag == 03) return Element::quadrangle(1);
            else return Element::node();
        }

        struct Module
        {

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

            explicit
            Module(
                std::basic_string_view<lolita::character> mesh_file_path
            )
            :
            file_(mesh_file_path)
            {
                setGeometricalEntities();
                setPhysicalEntities();
                setPhysicalGroups();
            }

            std::array<std::vector<GeometricalEntity>, 4>
            makeGeometricalEntities()
            {
                auto geometrical_entities = std::array<std::vector<GeometricalEntity>, 4>();
                auto const & file_lines = file_.lines_;
                lolita::index line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Entities"));
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::index num_points;
                lolita::index num_curves;
                lolita::index num_surfaces;
                lolita::index num_volumes;
                line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
                auto num_domains = std::array<lolita::index, 4>{num_points, num_curves, num_surfaces, num_volumes};
                offset += 1;
                for (lolita::index i = 0; i < 4; ++i) {
                    for (lolita::index j = 0; j < num_domains[i]; ++j) {
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
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
                return geometrical_entities;
            }

            std::vector<PhysicalEntity>
            makePhysicalEntities()
            {
                auto physical_entities = std::vector<PhysicalEntity>();
                auto const & file_lines = this->file_.lines_;
                lolita::index line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::index num_physical_names;
                line_stream >> num_physical_names;
                offset += 1;
                for (lolita::index i = 0; i < num_physical_names; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    lolita::index dim;
                    lolita::index tag;
                    std::basic_string<lolita::character> name;
                    line_stream >> dim >> tag >> name;
                    lolita::utility::removeCharacter(name, '"');
                    physical_entities.push_back(PhysicalEntity{tag, dim, name});
                    offset += 1;
                }
                return physical_entities;
            }

            std::array<std::unordered_set<lolita::index>, 4>
            getSubGeometricalEntities(
                    lolita::index d,
                    lolita::index t,
                    std::array<std::unordered_set<lolita::index>, 4> & a
            )
            {
                a[d].insert(t);
                for (lolita::index i = 0; i < geometrical_entities_[d][t - 1].bounding_entities_tags_.size(); ++i) {
                    lolita::index const d2 = d - 1;
                    lolita::index const t2 = geometrical_entities_[d][t - 1].bounding_entities_tags_[i];
                    a = getSubGeometricalEntities(d2, t2, a);
                }
                return a;
            }

            std::vector<PhysicalGroup>
            makePhysicalGroups()
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
                    for (lolita::index i = 0; i < geometrical_entities_[d][t - 1].bounding_entities_tags.size(); ++i) {
                        a = self(d - 1, geometrical_entities_[d][t - 1].bounding_entities_tags[i], a, self);
                    }
                    return a;
                };
                /*
                 *
                 */
                auto physical_groups = std::vector<PhysicalGroup>();
                for (lolita::index i = 0; i < physical_entities_.size(); ++i) {
                    std::array<std::unordered_set<lolita::index>, 4> tags;
                    for (lolita::index j = 0; j < 4; ++j) {
                        for (lolita::index k = 0; k < geometrical_entities_[j].size(); ++k) {
                            for (lolita::index l = 0; l < geometrical_entities_[j][k].physical_entities_tags_.size(); ++l) {
                                lolita::index const & t0 = geometrical_entities_[j][k].physical_entities_tags_[l];
                                lolita::index const & t1 = physical_entities_[i].tag_;
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
                    physical_groups.push_back(PhysicalGroup{physical_entities_[i].name_, group_tags});
                }
                return physical_groups;
            }

            void
            setGeometricalEntities()
            {
                auto const & file_lines = file_.lines_;
                lolita::index line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Entities"));
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::index num_points;
                lolita::index num_curves;
                lolita::index num_surfaces;
                lolita::index num_volumes;
                line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
                auto num_domains = std::array<lolita::index, 4>{num_points, num_curves, num_surfaces, num_volumes};
                offset += 1;
                for (lolita::index i = 0; i < 4; ++i) {
                    for (lolita::index j = 0; j < num_domains[i]; ++j) {
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
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
                        geometrical_entities_[i].push_back(GeometricalEntity{
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
                auto const & file_lines = file_.lines_;
                lolita::index line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::index num_physical_names;
                line_stream >> num_physical_names;
                offset += 1;
                for (lolita::index i = 0; i < num_physical_names; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    lolita::index dim;
                    lolita::index tag;
                    std::basic_string<lolita::character> name;
                    line_stream >> dim >> tag >> name;
                    lolita::utility::removeCharacter(name, '"');
                    physical_entities_.push_back(PhysicalEntity{tag, dim, name});
                    offset += 1;
                }
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
                    for (lolita::index i = 0; i < geometrical_entities_[d][t - 1].bounding_entities_tags.size(); ++i) {
                        a = self(d - 1, geometrical_entities_[d][t - 1].bounding_entities_tags[i], a, self);
                    }
                    return a;
                };
                /*
                 *
                 */
                for (lolita::index i = 0; i < physical_entities_.size(); ++i) {
                    std::array<std::unordered_set<lolita::index>, 4> tags;
                    for (lolita::index j = 0; j < 4; ++j) {
                        for (lolita::index k = 0; k < geometrical_entities_[j].size(); ++k) {
                            for (lolita::index l = 0; l < geometrical_entities_[j][k].physical_entities_tags_.size(); ++l) {
                                lolita::index const & t0 = geometrical_entities_[j][k].physical_entities_tags_[l];
                                lolita::index const & t1 = physical_entities_[i].tag_;
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
                    physical_groups_.push_back(PhysicalGroup{physical_entities_[i].name_, group_tags});
                }
            }

            lolita::utility::File file_;

            std::array<std::vector<GeometricalEntity>, 4> geometrical_entities_;

            std::vector<PhysicalEntity> physical_entities_;

            std::vector<PhysicalGroup> physical_groups_;

        };

        struct Implementation : public MeshParser<_mesh, t_domain>
        {

        private:

            using t_MeshParser = MeshParser<_mesh, t_domain>;

        public:

            void
            makeDomains(
                    Module const & module
            )
            {
                auto const & physical_groups = module.physical_groups_;
                for (int k = 0; k < physical_groups.size(); ++k) {
                    this->mesh_data_.domains_.push_back(std::make_shared<std::basic_string<lolita::character>>(physical_groups[k].name_));
                }
            }

            template<Element t_element>
            void
            makeElement(
                    Module const & module
            )
            requires(t_element.isNode())
            {
                auto const & file_lines = module.file_.lines_;
                auto const & physical_groups = module.physical_groups_;
                auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Nodes"));
                auto offset = 1;
                auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto num_entity_block = lolita::index();
                auto num_nodes = lolita::index();
                auto min_node_tag = lolita::index();
                auto max_node_tag = lolita::index();
                line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                offset += 1;
                for (int i = 0; i < num_entity_block; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    auto entity_dim = lolita::index();
                    auto entity_tag = lolita::index();
                    auto parametric = lolita::index();
                    auto num_nodes_in_block = lolita::index();
                    line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                    offset += 1;
                    for (int j = 0; j < num_nodes_in_block; ++j) {
                        auto tag_ = lolita::natural();
                        auto coordinates_ = lolita::domain::Point();
                        auto domains_ = std::vector<std::shared_ptr<std::basic_string<lolita::character>>>();
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                        line_stream >> tag_;
                        tag_ --;
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset + num_nodes_in_block]);
                        coordinates_.setZero();
                        for (int k = 0; k < 3; ++k) {
                            line_stream >> coordinates_(k);
                        }
                        for (int k = 0; k < physical_groups.size(); ++k) {
                            auto const & node_set_name = physical_groups[k].name_;
                            for (int l = 0; l < physical_groups[k].geometrical_entities_tags_[entity_dim].size(); ++l) {
                                lolita::index tag = physical_groups[k].geometrical_entities_tags_[entity_dim][l];
                                if (entity_tag == tag) {
                                    auto is_equal = [&] (std::shared_ptr<std::basic_string<lolita::character>> const & ptr_domain) {
                                        return * ptr_domain == node_set_name;
                                    };
                                    auto const & domains = this->mesh_data_.domains_;
                                    auto domain_index = std::distance(domains.begin(), std::find_if(domains.begin(), domains.end(), is_equal));
                                    domains_.push_back(domains[domain_index]);
                                }
                            }
                        }
                        t_MeshParser::template makeElement2<t_element>(tag_, coordinates_, domains_);
                        offset += 1;
                    }
                    offset += num_nodes_in_block;
                }
            }

            template<Element t_element>
            void
            makeElement(
                    Module const &
                    module
            )
            requires(!t_element.isNode())
            {
                auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
                auto const & file_lines = module.file_.lines_;
                auto & elements = this->mesh_data_.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
                auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Elements"));
                auto offset = 1;
                auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto num_entity_blocks = 0;
                auto num_elements = 0;
                auto min_element_tag = 0;
                auto max_element_tag = 0;
                line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
                offset += 1;
                for (int i = 0; i < num_entity_blocks; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    auto entity_dim = 0;
                    auto entity_tag = 0;
                    auto element_type_tag = 0;
                    auto num_elements_in_block = 0;
                    line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                    offset += 1;
                    auto const element = getElemType(element_type_tag);
                    if (entity_dim == t_domain.dim_ && element == t_element) {
                        for (int j = 0; j < num_elements_in_block; ++j) {
                            auto node_tags = std::array<lolita::index, t_element.num_nodes_>();
                            line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                            auto tag = lolita::index();
                            line_stream >> tag;
                            for (lolita::index k = 0; k < element.num_nodes_; ++k) {
                                auto node_tag = lolita::integer();
                                line_stream >> node_tag;
                                node_tags[k] = node_tag - 1;
                            }
                            t_MeshParser::template makeElement2<t_element>(node_tags);
                            offset += 1;
                        }
                    }
                    else {
                        offset += num_elements_in_block;
                    }
                }
            }

        };

    };

}

#endif /* BC3CC377_0788_4FE2_B68B_684298E3EEEA */
