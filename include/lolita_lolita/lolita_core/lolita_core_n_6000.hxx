#ifndef B2C662AD_4CF7_4877_B631_1E1921DEC692
#define B2C662AD_4CF7_4877_B631_1E1921DEC692

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_2000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_3000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4001.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4002.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4003.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4004.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4005.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_5000.hxx"

namespace lolita
{

    template<Domain t_domain>
    struct FiniteElementMap : ElementMap<FiniteElementHolder, t_domain>
    {

        std::unique_ptr<FiniteElementSet<t_domain>>
        makeFiniteElementSet()
        const
        {
            auto finite_element_set = std::make_unique<FiniteElementSet<t_domain>>();
            auto make_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    finite_element_set->template getElements<t_i, t_j>().push_back(element.second);
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            }; 
            make_elements(make_elements);
            return finite_element_set;
        }

    };
    
    template<Element t_element, Domain t_domain>
    struct MeshElement;
    
    template<Element t_element, Domain t_domain>
    requires(!t_element.isNode())
    struct MeshElement<t_element, t_domain>
    {
        
        template<Integer t_i, Integer t_j>
        static constexpr
        Integer
        getInnerNeighborNodeConnection(
            Integer i,
            Integer j
        )
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element, t_domain>::node_connectivity_))[i][j];
        }

        static
        std::basic_string<Character>
        getHash(
            std::array<Natural, t_element.num_nodes_> node_tags
        )
        {
            auto hash = std::basic_stringstream<Character>();
            std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
            for (auto node_tag : node_tags)
            {
                hash << node_tag;
            }
            return hash.str();
        }

        MeshElement()
        :
        node_tags_()
        {}

        template<Integer t_i, Integer t_j>
        void
        setElement(
            FiniteElementMap<t_domain> & element_map,
            std::shared_ptr<FiniteElementHolder<t_element, t_domain>> & ptr_element
        )
        const
        {
            auto const constexpr t_component = ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>();
            auto const constexpr t_is_initialized = t_i == 0 && t_j == 0;
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto const constexpr t_node_coordinates = ElementTraits<t_element, t_domain>::template getInnerNeighborCoordinates<Element::node()>();
            auto const constexpr t_component_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_component>();
            auto const constexpr t_neighbour_coordinates = ElementTraits<t_component, t_domain>::template getOuterNeighborCoordinates<t_element>();
            auto & components = element_map.template getElements<t_component_coordinates.dim_, t_component_coordinates.tag_>();
            for (Integer i = 0; i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<t_i, t_j>(); ++i)
            {
                auto component_hash = std::basic_string<Character>();
                if constexpr(!t_component.isNode())
                {
                    auto cpt = MeshElement<t_component, t_domain>();
                    for (Integer j = 0; j < t_component.num_nodes_; ++j)
                    {
                        cpt.node_tags_[j] = node_tags_[getInnerNeighborNodeConnection<t_i, t_j>(i, j)];
                    }
                    component_hash = MeshElement<t_component, t_domain>::getHash(cpt.node_tags_);
                    if (!components.contains(component_hash))
                    {
                        auto ptr_component = std::make_shared<FiniteElementHolder<t_component, t_domain>>();
                        cpt.template setElement<0, 0>(element_map, ptr_component);
                    }
                }
                else
                {
                    component_hash = MeshElement<t_component, t_domain>::getHash(node_tags_[getInnerNeighborNodeConnection<t_i, t_j>(i, 0)]);
                }
                ptr_element->template getInnerNeighbors<t_i, t_j>()[i] = components[component_hash];
                components[component_hash]->template getOuterNeighbors<t_neighbour_coordinates.dim_, t_neighbour_coordinates.tag_>().push_back(ptr_element);
            }
            if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<t_i>() - 1)
            {
                setElement<t_i, t_j + 1>(element_map, ptr_element);
            }
            else if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<>() - 1)
            {
                setElement<t_i + 1, 0>(element_map, ptr_element);
            }
            if constexpr (t_is_initialized)
            {
                auto const & nodes = ptr_element->template getInnerNeighbors<t_node_coordinates.dim_, t_node_coordinates.tag_>();
                auto domains = std::set<std::shared_ptr<MeshDomain>>();
                for (auto const & domain : nodes[0]->domains_)
                {
                    auto has_domain = true;
                    for (Integer j = 1; j < t_element.num_nodes_; ++j)
                    {
                        has_domain = std::find(nodes[j]->domains_.begin(), nodes[j]->domains_.end(), domain) != nodes[j]->domains_.end();
                    }
                    if (has_domain)
                    {
                        domains.insert(domain);
                    }
                }
                ptr_element->tag_ = element_map.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>().size();
                ptr_element->domains_.assign(domains.begin(), domains.end());
                ptr_element->coordinates_ = std::make_shared<Point>(ptr_element->getCurrentCentroid());
                element_map.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>()[getHash(node_tags_)] = ptr_element;
            }
        }
        
        void
        makeElement(
            FiniteElementMap<t_domain> & element_map
        )
        const
        {
            auto ptr_element = std::make_shared<FiniteElementHolder<t_element, t_domain>>(FiniteElementHolder<t_element, t_domain>());
            setElement<0, 0>(element_map, ptr_element);
        }

        template<Integer t_i>
        static
        void
        setElementOuterNeighborhood(
            std::shared_ptr<FiniteElementHolder<t_element, t_domain>> & ptr_element
        )
        {
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto const constexpr t_node_coordinates = ElementTraits<t_element, t_domain>::template getInnerNeighborCoordinates<Element::node()>();
            auto const & element_nds = ptr_element->template getInnerNeighbors<t_node_coordinates.dim_, t_node_coordinates.tag_>();
            for (auto i = 0; i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<t_node_coordinates.dim_, t_node_coordinates.tag_>(); ++i)
            {
                auto const & nde = element_nds[i];
                auto const & ngs = nde->template getOuterNeighbors<t_element_coordinates.dim_ - 1, t_i>();
                for (auto const & neighbour : ngs)
                {
                    if (((neighbour->tag_ != ptr_element->tag_) && t_i == t_element_coordinates.tag_) || (t_i != t_element_coordinates.tag_))
                    {
                        auto & element_ngs = ptr_element->template getOuterNeighbors<0, t_i>();
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
            if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<t_element_coordinates.dim_>() - 1)
            {
                setElementOuterNeighborhood<t_i + 1>(ptr_element);
            }
        }
        
        static inline
        void
        initialize(
            std::shared_ptr<FiniteElementHolder<t_element, t_domain>> & ptr_element
        )
        {
            setElementOuterNeighborhood<0>(ptr_element);
        }

        std::array<Natural, t_element.num_nodes_> node_tags_;

    };

    template<Element t_element, Domain t_domain>
    requires(t_element.isNode())
    struct MeshElement<t_element, t_domain>
    {

        static
        std::basic_string<Character>
        getHash(
            Natural tag
        )
        {
            return std::to_string(tag);
        }

        MeshElement()
        :
        tag_(),
        coordinates_(std::make_shared<Point>(Point::Zero())),
        domains_()
        {}
        
        void
        setElement(
            FiniteElementMap<t_domain> & element_set,
            std::shared_ptr<FiniteElementHolder<t_element, t_domain>> & ptr_element
        )
        const
        {
            auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            ptr_element->tag_ = tag_;
            ptr_element->domains_ = domains_;
            ptr_element->coordinates_ = coordinates_;
            element_set.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>()[getHash(tag_)] = ptr_element;
        }
        
        void
        makeElement(
            FiniteElementMap<t_domain> & element_set
        )
        const
        {
            auto ptr_element = std::make_shared<FiniteElementHolder<t_element, t_domain>>();
            setElement(element_set, ptr_element);
        }
        
        static inline
        void
        initialize(
            std::shared_ptr<FiniteElementHolder<t_element, t_domain>> & ptr_element
        )
        {}
        
        Natural tag_;
        
        std::shared_ptr<Point> coordinates_;

        std::vector<std::shared_ptr<MeshDomain>> domains_;

    };
        
    template<Domain t_domain>
    struct MeshElementSet : ElementSet<MeshElement, t_domain>
    {

        MeshElementSet()
        {}
        
        std::unique_ptr<FiniteElementSet<t_domain>>
        makeFiniteElementSet()
        const
        {
            auto element_map = std::make_unique<FiniteElementMap<t_domain>>();
            auto make_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    element->template makeElement(* element_map);
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1u>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1)
                {
                    self.template operator()<t_i + 1u, 0u>(self);
                }
            };
            make_elements(make_elements);
            auto make_elements_outer_neighborhood = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto & element : element_map->template getElements<t_i, t_j>())
                {
                    MeshElement<t_element, t_domain>::initialize(element.second);
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
            make_elements_outer_neighborhood(make_elements_outer_neighborhood);
            return element_map->makeFiniteElementSet();
        }

    };    
    
    struct MeshFileParserBase
    {

        MeshFileParserBase(
            std::basic_string_view<Character> str
        )
        :
        file_(str)
        {}

        lolita::utility::File file_;
        
    };

    struct GmshFileParser : MeshFileParserBase
    {

    private:

        static constexpr
        Element
        getElemType(
            Integer tag
        )
        {
            if (tag == 15) return Element::node();
            else if (tag == 01) return Element::segment(1);
            else if (tag == 02) return Element::triangle(1);
            else if (tag == 03) return Element::quadrangle(1);
            else return Element::node();
        }

        static constexpr
        Integer
        getTagFromElement(
            Element element
        )
        {
            if (element.isNode()) return 15;
            else if (element.isSegment(1)) return 1;
            else if (element.isTriangle(1)) return 2;
            else if (element.isQuadrangle(1)) return 3;
            else return -1;
        }

        struct PhysicalEntity
        {

            static
            std::basic_string<Character>
            getHash(
                Integer tag
            )
            {
                return std::to_string(tag);
            }

            PhysicalEntity(
                std::basic_string<Character> const & name,
                Integer dim
            )
            :
            mesh_domain_(std::make_shared<MeshDomain>(dim, name))
            {}

            std::shared_ptr<MeshDomain> mesh_domain_;

        };

        struct GeometricalEntity
        {

            static
            std::basic_string<Character>
            getHash(
                Integer dim,
                Integer tag
            )
            {
                auto hash = std::basic_stringstream<Character>();
                hash << dim;
                hash << tag;
                return hash.str();
            }

            GeometricalEntity()
            :
            physical_entities_(),
            bounding_entities_()
            {}

            void
            setPhysicalEntity(
                std::shared_ptr<PhysicalEntity> const & physical_entity
            )
            {
                physical_entities_.insert(physical_entity);
                for (auto & bounding_entity : bounding_entities_)
                {
                    bounding_entity->setPhysicalEntity(physical_entity);
                }
            }

            void
            setBoundingEntity(
                std::shared_ptr<GeometricalEntity> & bounding_entity
            )
            {
                bounding_entities_.insert(bounding_entity);
                for (auto const & physical_entity : physical_entities_)
                {
                    bounding_entity->setPhysicalEntity(physical_entity);
                }
            }

            std::set<std::shared_ptr<PhysicalEntity>> physical_entities_;

            std::set<std::shared_ptr<GeometricalEntity>> bounding_entities_;

        };
        
        void
        setPhysicalEntities()
        {
            auto const & file_lines = this->file_.lines_;
            auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
            auto offset = 1;
            auto line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
            auto num_physical_names = Integer();
            line_stream >> num_physical_names;
            offset += 1;
            for (Integer i = 0; i < num_physical_names; ++i)
            {
                line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                auto dim = Integer();
                auto name = std::basic_string<Character>();
                auto physical_entity_tag = Integer();
                line_stream >> dim >> physical_entity_tag >> name;
                name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
                auto hash = PhysicalEntity::getHash(physical_entity_tag);
                physical_entities_[hash] = std::make_shared<PhysicalEntity>(PhysicalEntity(name, dim));
                offset += 1;
            }
        }

        void
        setGeometricalEntities()
        {
            auto const & file_lines = this->file_.lines_;
            auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Entities"));
            auto offset = 1;
            auto line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
            auto num_points = Integer();
            auto num_curves = Integer();
            auto num_surfaces = Integer();
            auto num_volumes = Integer();
            line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
            auto num_domains = std::array<Integer, 4>{num_points, num_curves, num_surfaces, num_volumes};
            offset += 1;
            for (Integer i = 0; i < 4; ++i)
            {
                for (Integer j = 0; j < num_domains[i]; ++j)
                {
                    line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                    auto tag = Integer();
                    line_stream >> tag;
                    auto geometrical_entity_hash = GeometricalEntity::getHash(i, tag);
                    geometrical_entities_[geometrical_entity_hash] = std::make_shared<GeometricalEntity>(GeometricalEntity());
                    if (i == 0)
                    {
                        for (Integer k = 0; k < 3; ++k)
                        {
                            auto a = Real();
                            line_stream >> a;
                        }
                    }
                    else
                    {
                        for (Integer k = 0; k < 6; ++k)
                        {
                            auto a = Real();
                            line_stream >> a;
                        }
                    }
                    auto num_physical_entities = Integer();
                    line_stream >> num_physical_entities;
                    for (Integer k = 0; k < num_physical_entities; ++k)
                    {
                        auto physical_entity_tag = Integer();
                        line_stream >> physical_entity_tag;
                        auto physical_entity_hash = PhysicalEntity::getHash(physical_entity_tag);
                        geometrical_entities_[geometrical_entity_hash]->setPhysicalEntity(physical_entities_[physical_entity_hash]);
                    }
                    if (i > 0)
                    {
                        auto num_bounding_entities = Integer();
                        line_stream >> num_bounding_entities;
                        for (Integer k = 0; k < num_bounding_entities; ++k)
                        {
                            auto bounding_entity_tag = Integer();
                            line_stream >> bounding_entity_tag;
                            auto bounding_entity_hash = GeometricalEntity::getHash(i - 1, std::abs(bounding_entity_tag));
                            geometrical_entities_[geometrical_entity_hash]->setBoundingEntity(geometrical_entities_[bounding_entity_hash]);
                        }
                    }
                    offset += 1;
                }
            }
        }

    public:

        GmshFileParser(
            std::basic_string_view<Character> str
        )
        :
        MeshFileParserBase(str)
        {
            setPhysicalEntities();
            setGeometricalEntities();
        }

        template<Domain t_domain, Element t_element>
        void
        setMeshElement(
            MeshElementSet<t_domain> & mesh_element_set
        )
        const
        requires(t_element.isNode())
        {
            auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto & elements = mesh_element_set.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
            auto const & file_lines = this->file_.lines_;
            auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Nodes"));
            auto offset = 1;
            auto line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
            auto num_entity_block = Integer();
            auto num_nodes = Integer();
            auto min_node_tag = Integer();
            auto max_node_tag = Integer();
            line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
            offset += 1;
            for (Integer i = 0; i < num_entity_block; ++i)
            {
                line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                auto entity_dim = Integer();
                auto entity_tag = Integer();
                auto parametric = Integer();
                auto num_nodes_in_block = Integer();
                line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                auto const entity_hash = GeometricalEntity::getHash(entity_dim, entity_tag);
                offset += 1;
                for (Integer j = 0; j < num_nodes_in_block; ++j)
                {
                    auto mesh_element = std::make_shared<MeshElement<t_element, t_domain>>();
                    auto tag = Natural();
                    auto coordinates_ = Point();
                    auto domains_ = std::vector<std::shared_ptr<std::basic_string<Character>>>();
                    line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                    line_stream >> tag;
                    mesh_element->tag_ = tag - 1;
                    line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset + num_nodes_in_block]);
                    for (Integer k = 0; k < 3; ++k)
                    {
                        line_stream >> mesh_element->coordinates_->operator()(k);
                    }
                    for (auto const & physical_entity : geometrical_entities_.at(entity_hash)->physical_entities_)
                    {
                        mesh_element->domains_.push_back(physical_entity->mesh_domain_);
                    }
                    // auto const element_hash = MeshElement<t_element, t_domain>::getHash(mesh_element->tag_);
                    // elements[element_hash] = mesh_element;
                    elements.push_back(mesh_element);
                    offset += 1;
                }
                offset += num_nodes_in_block;
            }
        }
        
        template<Domain t_domain, Element t_element>
        void
        setMeshElement(
            MeshElementSet<t_domain> & mesh_element_set
        )
        const
        requires(!t_element.isNode())
        {
            auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto & elements = mesh_element_set.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
            auto const & file_lines = this->file_.lines_;
            auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Elements"));
            auto offset = 1;
            auto line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
            auto num_entity_blocks = 0;
            auto num_elements = 0;
            auto min_element_tag = 0;
            auto max_element_tag = 0;
            line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
            offset += 1;
            for (Integer i = 0; i < num_entity_blocks; ++i)
            {
                line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                auto entity_dim = 0;
                auto entity_tag = 0;
                auto element_type_tag = 0;
                auto num_elements_in_block = 0;
                line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                offset += 1;
                auto const element = getElemType(element_type_tag);
                if (element == t_element)
                {
                    for (Integer j = 0; j < num_elements_in_block; ++j)
                    {
                        auto mesh_element = std::make_shared<MeshElement<t_element, t_domain>>();
                        line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                        auto tag = Integer();
                        line_stream >> tag;
                        for (Integer k = 0; k < t_element.num_nodes_; ++k)
                        {
                            auto node_tag = Integer();
                            line_stream >> node_tag;
                            mesh_element->node_tags_[k] = node_tag - 1;
                        }
                        // auto const element_hash = MeshElement<t_element, t_domain>::getHash(mesh_element->node_tags_);
                        // elements[element_hash] = mesh_element;
                        elements.push_back(mesh_element);
                        offset += 1;
                    }
                }
                else
                {
                    offset += num_elements_in_block;
                }
            }
        }

        template<Domain t_domain>
        static
        void
        setOutput(
            std::basic_string<Character> const & file_path,
            std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
            auto... behavior_label
        )
        {
            auto labels = std::array<std::basic_string<Character>, sizeof...(behavior_label)>{behavior_label...};
            auto outfile = std::ofstream();
            outfile.open(file_path);
            outfile << "$MeshFormat\n";
            outfile << "2.2 0 8\n";
            outfile << "$EndMeshFormat\n";
            outfile << "$Nodes\n";
            outfile << element_set->template getNumElements<0>() + numerics::sum(element_set->template getNumIntegrationPoints<>(behavior_label)...) << "\n";
            auto c_node = Natural(1);
            for (auto const & node : element_set->template getElements<0, 0>())
            {
                auto const & coordinates = node->getCurrentCoordinates();
                outfile << c_node << " " << coordinates(0) << " " << coordinates(1) << " " << coordinates(2) << "\n";
                c_node ++;
            }
            auto set_integration_nodes = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : element_set->template getElements<t_i, t_j>())
                {
                    for (auto const & label : labels)
                    {
                        if (element->quadrature_.contains(label))
                        {
                            for (auto const & integration_point : element->quadrature_.at(label).ips_)
                            {
                                auto const & coordinates = integration_point.getCurrentCoordinates();
                                outfile << c_node << " " << coordinates(0) << " " << coordinates(1) << " " << coordinates(2) << "\n";
                                c_node ++;
                            }
                        }
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            }; 
            set_integration_nodes(set_integration_nodes);
            outfile << "$EndNodes\n";
            outfile << "$Elements\n";
            outfile << element_set->template getNumElements<>() + numerics::sum(element_set->template getNumIntegrationPoints<>(behavior_label)...) << "\n";
            auto c_element = Natural(1);
            for (auto const & node : element_set->template getElements<0, 0>())
            {
                outfile << c_element << " 15 2 0 0 " << c_element << "\n";
                c_element ++;
            }
            auto set_integration_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : element_set->template getElements<t_i, t_j>())
                {
                    for (auto const & label : labels)
                    {
                        if (element->quadrature_.contains(label))
                        {
                            for (auto const & integration_point : element->quadrature_.at(label).ips_)
                            {
                                outfile << c_element << " 15 2 1 1 " << c_element << "\n";
                                c_element ++;
                            }
                        }
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            set_integration_elements(set_integration_elements);
            auto set_elements = [&] <Integer t_i = 1, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : element_set->template getElements<t_i, t_j>())
                {
                    outfile << c_element << " " << getTagFromElement(DomainTraits<t_domain>::template getElement<t_i, t_j>()) << " 2 0 0";
                    for (auto const & node : element->template getInnerNeighbors<t_i - 1, 0>())
                    {
                        outfile << " " << node->getTag() + 1;
                    }
                    outfile << "\n";
                    c_element ++;
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            set_elements(set_elements);
            outfile << "$EndElements\n";
        }

        template<Integer t_coordinate, Domain t_domain>
        static
        void
        addQuadratureStrainOutput(
            std::basic_string<Character> const & file_path,
            std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
            Integer time_step_index,
            Real time_step_value,
            std::basic_string<Character> behavior_label,
            Integer row
        )
        {
            if (!std::filesystem::exists(file_path))
            {
                throw std::runtime_error("File does not exist");
            }
            auto outfile = std::ofstream();
            auto c_element = Natural();
            outfile.open(file_path, std::ios_base::app);
            outfile << "$NodeData\n";
            outfile << "1\n";
            outfile << "\"" << behavior_label << " " << row << " Strain\"\n";
            outfile << "1\n";
            outfile << time_step_value << "\n";
            outfile << "3\n";
            outfile << time_step_index << "\n";
            outfile << 1 << "\n"; // size of gradients
            outfile << element_set->template getNumIntegrationPoints<>(behavior_label) << "\n"; // number of quad pts
            c_element = element_set->template getNumElements<0>() + 1;
            auto set_strain = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : element_set->template getElements<t_coordinate, t_j>())
                {
                    if (element->quadrature_.contains(behavior_label))
                    {
                        for (auto const & integration_point : element->quadrature_.at(behavior_label).ips_)
                        {
                            outfile << c_element << " " << integration_point.behavior_data_->s1.gradients[row] << "\n";
                            c_element ++;
                        }
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_coordinate>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            };
            set_strain(set_strain);
            outfile << "$EndNodeData\n";
        }

        template<Integer t_coordinate, Domain t_domain>
        static
        void
        addQuadratureStressOutput(
            std::basic_string<Character> const & file_path,
            std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
            Integer time_step_index,
            Real time_step_value,
            std::basic_string<Character> behavior_label,
            Integer row
        )
        {
            if (!std::filesystem::exists(file_path))
            {
                throw std::runtime_error("File does not exist");
            }
            auto outfile = std::ofstream();
            auto c_element = Natural();
            outfile.open(file_path, std::ios_base::app);
            outfile << "$NodeData\n";
            outfile << "1\n";
            outfile << "\"" << behavior_label << " " << row << " Stress\"\n";
            outfile << "1\n";
            outfile << time_step_value << "\n";
            outfile << "3\n";
            outfile << time_step_index << "\n";
            outfile << 1 << "\n"; // size of gradients
            outfile << element_set->template getNumIntegrationPoints<>(behavior_label) << "\n"; // number of quad pts
            c_element = element_set->template getNumElements<0>() + 1;
            auto set_strain = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : element_set->template getElements<t_coordinate, t_j>())
                {
                    if (element->quadrature_.contains(behavior_label))
                    {
                        for (auto const & integration_point : element->quadrature_.at(behavior_label).ips_)
                        {
                            outfile << c_element << " " << integration_point.behavior_data_->s1.thermodynamic_forces[row] << "\n";
                            c_element ++;
                        }
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_coordinate>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            };
            set_strain(set_strain);
            outfile << "$EndNodeData\n";
        }
        
        template<Integer t_coordinate, Domain t_domain, auto... t_args>
        static
        void
        addNodalDofOutput(
            std::basic_string<Character> const & file_path,
            std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
            Integer time_step_index,
            Real time_step_value,
            std::basic_string<Character> unknown_label,
            Integer row,
            Integer col
        )
        {
            if (!std::filesystem::exists(file_path))
            {
                throw std::runtime_error("File does not exist");
            }
            auto outfile = std::ofstream();
            auto c_element = Natural();
            outfile.open(file_path, std::ios_base::app);
            // writing strain
            outfile << "$NodeData\n";
            outfile << "1\n";
            // outfile << "\"" << label << "NodalValues\"\n";
            outfile << "\"" << unknown_label << " " << row << " " << col << " NodalValues\"\n";
            outfile << "1\n";
            outfile << time_step_value << "\n";
            outfile << "3\n";
            outfile << time_step_index << "\n";
            outfile << 1 << "\n"; // size of gradients
            outfile << element_set->template getNumElements<0>() << "\n"; // number of nodes
            c_element = 1;
            auto nodal_values = element_set->template getNodalValues<t_coordinate, t_args...>(unknown_label, row, col);
            for (auto nodal_value : nodal_values)
            {
                outfile << c_element << " " << nodal_value << "\n";
                c_element ++;
            }
            outfile << "$EndNodeData\n";
        }
        
        template<Integer t_coordinate, Domain t_domain, auto... t_args>
        static
        void
        addQuadratureDofOutput(
            std::basic_string<Character> const & file_path,
            std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
            Integer time_step_index,
            Real time_step_value,
            std::basic_string<Character> unknown_label,
            std::basic_string<Character> quadrature_label,
            Integer row,
            Integer col
        )
        {
            if (!std::filesystem::exists(file_path))
            {
                throw std::runtime_error("File does not exist");
            }
            auto outfile = std::ofstream();
            auto c_element = Natural();
            outfile.open(file_path, std::ios_base::app);
            // writing strain
            outfile << "$NodeData\n";
            outfile << "1\n";
            outfile << "\"" << unknown_label << " " << row << " " << col << " QuadratureValues\"\n";
            outfile << "1\n";
            outfile << time_step_value << "\n";
            outfile << "3\n";
            outfile << time_step_index << "\n";
            outfile << 1 << "\n"; // size of gradients
            std::cout << "!!!!\n";
            std::cout << element_set->template getNumIntegrationPoints<>(quadrature_label);
            std::cout << "\n!!!!\n";
            outfile << element_set->template getNumIntegrationPoints<>(quadrature_label) << "\n"; // number of nodes
            c_element = element_set->template getNumElements<0>() + 1;
            auto quadrature_values = element_set->template getQuadratureValues<t_coordinate, t_args...>(unknown_label, quadrature_label, row, col);
            for (auto quadrature_value : quadrature_values)
            {
                outfile << c_element << " " << quadrature_value << "\n";
                c_element ++;
            }
            outfile << "$EndNodeData\n";
        }
        
    private:
        
        std::map<std::basic_string<Character>, std::shared_ptr<PhysicalEntity>> physical_entities_;

        std::map<std::basic_string<Character>, std::shared_ptr<GeometricalEntity>> geometrical_entities_;

    };
    
    struct MeshFileParser
    {

        MeshFileParser(
            std::basic_string_view<Character> str
        )
        :
        file_path_(str)
        {}

        template<Domain t_domain>
        std::unique_ptr<FiniteElementSet<t_domain>>
        makeFiniteElementSet()
        const
        {
            auto gmsh_file_parser = GmshFileParser(file_path_);
            auto mesh_element_set = MeshElementSet<t_domain>();
            auto make_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                gmsh_file_parser.template setMeshElement<t_domain, t_element>(mesh_element_set);
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            make_elements(make_elements);
            return mesh_element_set.makeFiniteElementSet();
        }

        std::basic_string_view<Character> file_path_;

    };

} // namespace lolita


#endif /* B2C662AD_4CF7_4877_B631_1E1921DEC692 */
