#ifndef CB51704E_040A_4B3D_AE3C_46C0DE3B8543
#define CB51704E_040A_4B3D_AE3C_46C0DE3B8543

#include "lolita.hxx"
#include "core/000_physics_traits.hxx"
#include "core/100_geometry.hxx"
#include "core/linear_system.hxx"
#include "core/200_quadrature.hxx"
#include "core/201_finite_element_dof.hxx"
#include "core/202_finite_element_frm.hxx"
#include "core/300_finite_element.hxx"
#include "core/400_finite_element_basis.hxx"
#include "core/500_finite_element_hdg_discretization.hxx"
#include "core/600_finite_element_set.hxx"

namespace lolita
{

    template<Domain t_domain>
    struct FiniteElementMap : ElementMap<FiniteElement, t_domain>, DomainMap<FiniteDomain, t_domain>
    {

        std::unique_ptr<FiniteElementSet<t_domain>>
        makeFiniteElementSet()
        const
        {
            auto finite_element_set = std::make_unique<FiniteElementSet<t_domain>>();
            auto make_sets = [&] <Integer t_i = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & dm : this->template getDomains<t_i>())
                {
                    finite_element_set->template getDomains<t_i>().push_back(dm.second);
                }
                if constexpr (t_i < t_domain.getDim())
                {
                    self.template operator()<t_i + 1>(self);
                }
            };
            make_sets(make_sets);
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

    template<Integer t_dim, Domain t_domain>
    struct MeshDomain
    {
        
        static
        std::basic_string<Character>
        getHash(
            std::basic_string<Character> && label
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << t_dim;
            hash << std::forward<std::basic_string<Character>>(label);
            return hash.str();
        }
        
        static
        std::basic_string<Character>
        getHash(
            Integer tag
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << t_dim;
            hash << std::setfill('0') << std::setw(10) << tag;
            return hash.str();
        }

        MeshDomain()
        :
        label_()
        {}

        MeshDomain(
            std::basic_string<Character> const & label
        )
        :
        label_(label)
        {}

        MeshDomain(
            std::basic_string<Character> && label
        )
        :
        label_(std::move(label))
        {}

        std::basic_string<Character> const &
        getLabel()
        const
        {
            return label_;
        }

        std::shared_ptr<FiniteDomain<t_dim, t_domain>>
        makeDom()
        const
        {
            return std::make_shared<FiniteDomain<t_dim, t_domain>>(label_);
        }

        std::basic_string<Character> label_;
    
    };
    
    /**
     * @brief A helper class to create a Finite Element object for a node.
     * An object of this class is intended to be built by any mesh format parser.
     * 
     * @tparam t_element The element object defining the element geometry
     * @tparam t_domain The domain object defining the global domain geometry
     */
    template<Element t_element, Domain t_domain>
    struct MeshElement;
    
    /**
     * @brief A helper class to create a Finite Element object for a node.
     * An object of this class is intended to be built by any mesh format parser.
     * For a generic element that is not a node, it consists in defining
     * the list of all node tags the element is connected to, and the list of all domains it belongs to.
     * 
     * @tparam t_element The element object defining the element geometry
     * @tparam t_domain The domain object defining the global domain geometry
     */
    template<Element t_element, Domain t_domain>
    requires(!t_element.isNode())
    struct MeshElement<t_element, t_domain>
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element>;

        using t_Domains = typename t_ElementTraits::template Domains<MeshDomain, t_domain>;

    public:
        
        /**
         * @brief Get index of the jth node tag of the ith inner neighbor in the element node tag array.
         * 
         * @tparam t_i The dimension of the inner neighbor.
         * @tparam t_j The tag of the inner neighbor.
         * @param i The index of the inner neighbor.
         * @param j The index of the inner neighbor node.
         * @return The index of the jth node tag of the ith inner neighbor in the element node tag array.
         */
        template<Integer t_i, Integer t_j>
        static constexpr
        Integer
        getInnerNeighborNodeConnection(
            Integer i,
            Integer j
        )
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element>::node_connectivity_))[i][j];
        }

        /**
         * @brief Get an unique identifier for the element, based on its node tags. Node tags are sorted, and concatenated into a string object.
         * 
         * @param node_tags An array containing all node tags of the element 
         * @return A string object the uniquely identifies the element.
         */
        static
        std::basic_string<Character>
        getHash(
            std::array<Natural, t_element.getNumNodes()> node_tags
        )
        {
            auto hash = std::basic_stringstream<Character>();
            std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
            for (auto node_tag : node_tags)
            {
                hash << std::setfill('0') << std::setw(10) << node_tag;
            }
            return hash.str();
        }

        /**
         * @brief Construct a new Mesh Element object
         * 
         */
        MeshElement()
        :
        node_tags_(),
        domains_(),
        finite_element_()
        {}

        /**
         * @brief Add a domain pointer that the element is connected to. If the domain pointer is already present in the list, nothing happens.
         * 
         * @param domain The domain pointer to be added.
         */
        void
        addDomain(
            std::shared_ptr<MeshDomain<t_element.getDim(), t_domain>> const & domain
        )
        {
            for (auto const & dom : domains_)
            {
                if (dom == domain)
                {
                    return;
                }
            }
            domains_.push_back(domain);
        }

        /**
         * @brief Set the ith node tag of the element
         * 
         * @param index The ith node tag index of the element
         * @param tag The ith node tag value of the element
         */
        void
        setNodeTag(
            Integer index,
            Natural const & tag
        )
        {
            node_tags_[index] = tag;
        }

        /**
         * @brief Get the ith node tag of the element
         * 
         * @param index The node tag index of the element
         * @return A reference to the node tag of the element
         */
        Natural const &
        getNodeTag(
            Integer index
        )
        const
        {
            return node_tags_[index];
        }

        /**
         * @brief Get the array containing all node tags of the element
         * 
         * @return A reference to the array containing all node tags of the element
         */
        std::array<Natural, t_element.getNumNodes()> const &
        getNodeTags()
        const
        {
            return node_tags_;
        }

        /**
         * @brief Get an unique identifier for the element, based on its node tags. Node tags are sorted, and concatenated into a string object.
         * 
         * @return A string object the uniquely identifies the element.
         */
        std::basic_string<Character>
        getHash()
        const
        {
            auto node_tags = node_tags_;
            auto hash = std::basic_stringstream<Character>();
            std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
            for (auto node_tag : node_tags)
            {
                hash << std::setfill('0') << std::setw(10) << node_tag;
            }
            return hash.str();
        }

        /**
         * @brief Creates a finite element object from the data members, that are the node tags and the list of domains containing the element. A finite
         * element is created, and all its inner neighbors are sought in the finite element set given in argument. If an inner neighbor is not found, i.e. if
         * it does not exists, it is created by calling this function recursively. Once all inner neighbors are either found or created, they are linked to
         * the finite element, which is then initialized by linking it to all the domains it belongs to, and by setting its coordinates, which is the finite
         * element barycenter if the element is not a node.
         * 
         * @param element_map The finite element list of the mesh.
         */
        void
        makeElement(
            FiniteElementMap<t_domain> & element_map
        )
        {
            auto make_element = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto const constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<t_i, t_j>();
                auto const constexpr t_is_initialized = t_i == 0 && t_j == 0;
                auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
                auto const constexpr t_node_coordinates = ElementTraits<t_element>::template getInnerNeighborCoordinates<Element::node()>();
                auto const constexpr t_inner_neighbor_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_inner_neighbor>();
                auto const constexpr t_neighbour_coordinates = ElementTraits<t_inner_neighbor>::template getOuterNeighborCoordinates<t_domain, t_element>();
                auto & inner_neighbors = element_map.template getElements<t_inner_neighbor_coordinates.getDim(), t_inner_neighbor_coordinates.getTag()>();
                if constexpr (t_is_initialized)
                {
                    auto tag = element_map.template getElements<t_element_coordinates.getDim(), t_element_coordinates.getTag()>().size();
                    finite_element_ = std::make_shared<FiniteElement<t_element, t_domain>>(FiniteElement<t_element, t_domain>(tag));
                }
                for (auto i = 0; i < ElementTraits<t_element>::template getNumInnerNeighbors<t_i, t_j>(); ++i)
                {
                    auto inner_neighbor_hash = std::basic_string<Character>();
                    if constexpr(!t_inner_neighbor.isNode())
                    {
                        auto mesh_inner_neighbor = MeshElement<t_inner_neighbor, t_domain>();
                        for (auto j = 0; j < t_inner_neighbor.num_nodes_; ++j)
                        {
                            mesh_inner_neighbor.setNodeTag(j, getNodeTag(getInnerNeighborNodeConnection<t_i, t_j>(i, j)));
                        }
                        inner_neighbor_hash = mesh_inner_neighbor.getHash();
                        if (!inner_neighbors.contains(inner_neighbor_hash))
                        {
                            mesh_inner_neighbor.makeElement(element_map);
                        }
                    }
                    else
                    {
                        inner_neighbor_hash = MeshElement<t_inner_neighbor, t_domain>::getHash(getNodeTag(getInnerNeighborNodeConnection<t_i, t_j>(i, 0)));
                    }
                    auto & inner_neighbor = inner_neighbors.at(inner_neighbor_hash);
                    finite_element_->template getInnerNeighbors<t_i, t_j>()[i] = inner_neighbor;
                    inner_neighbor->template getOuterNeighbors<t_neighbour_coordinates.dim_, t_neighbour_coordinates.tag_>().push_back(finite_element_);
                }
                if constexpr (t_j < ElementTraits<t_element>::template getNumInnerNeighbors<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<>() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
                if constexpr (t_is_initialized)
                {
                    for (auto const & domain : domains_)
                    {
                        finite_element_->addDomain(element_map.template getDomains<t_element.getDim()>().at(domain->getLabel()));
                    }
                    finite_element_->setCoordinates(finite_element_->getCurrentCentroid());
                    element_map.template getElements<t_element_coordinates.getDim(), t_element_coordinates.getTag()>()[getHash(node_tags_)] = finite_element_;
                }
            };
            make_element(make_element);
        }

    private:

        /**
         * @brief The tag of each node the element is connected to.
         * 
         */
        ElementNodeConnectivity<t_element> node_tags_;

        /**
         * @brief A vector containing each MeshDomain the element is connected to.
         * 
         */
        t_Domains domains_;

        /**
         * @brief A pointer to the FiniteElement to be build.
         * 
         */
        std::shared_ptr<FiniteElement<t_element, t_domain>> finite_element_;

    };

    /**
     * @brief A helper class to create a Finite Element object for a node.
     * An object of this class is intended to be built by any mesh format parser.
     * For a node, it consists in defining a tag that unambiguously identifies it,
     * its coordinates, and the list of all domains it belongs to.
     * 
     * @tparam t_element 
     * @tparam t_domain 
     */
    template<Element t_element, Domain t_domain>
    requires(t_element.isNode())
    struct MeshElement<t_element, t_domain>
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element>;

        using t_Domains = typename t_ElementTraits::template Domains<MeshDomain, t_domain>;

    public:

        static
        std::basic_string<Character>
        getHash(
            Natural tag
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << tag;
            return hash.str();
        }

        MeshElement()
        :
        tag_(),
        coordinates_(Point::Zero()),
        domains_(),
        finite_element_()
        {}

        void
        addDomain(
            std::shared_ptr<MeshDomain<t_element.getDim(), t_domain>> const & domain
        )
        {
            for (auto const & dom : domains_)
            {
                if (dom == domain)
                {
                    return;
                }
            }
            domains_.push_back(domain);
        }

        void
        setTag(
            Natural tag
        )
        {
            tag_ = tag;
        }

        Natural const &
        getNodeTag()
        const
        {
            return tag_;
        }

        void
        setCoordinate(
            Integer i,
            Real const & coordinate
        )
        {
            coordinates_(i) = coordinate;
        }
        
        std::basic_string<Character>
        getHash()
        const
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << tag_;
            return hash.str();
        }
        
        void
        makeElement(
            FiniteElementMap<t_domain> & element_set
        )
        {
            auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            finite_element_ = std::make_shared<FiniteElement<t_element, t_domain>>(FiniteElement<t_element, t_domain>(tag_));
            finite_element_->setCoordinates(coordinates_);
            for (auto const & domain : domains_)
            {
                finite_element_->addDomain(element_set.template getDomains<t_element.getDim()>().at(domain->getLabel()));
            }
            element_set.template getElements<t_element_coordinates.getDim(), t_element_coordinates.getTag()>()[getHash(tag_)] = finite_element_;
        }

    private:
        
        /**
         * @brief A tag that unambiguously identifies the node.
         * 
         */
        Natural tag_;
        
        /**
         * @brief The coordinates of the node in the mesh.
         * 
         */
        Point coordinates_;

        /**
         * @brief A vector containing each MeshDomain the element is connected to.
         * 
         */
        t_Domains domains_;

        /**
         * @brief A pointer to the FiniteElement to be build.
         * 
         */
        std::shared_ptr<FiniteElement<t_element, t_domain>> finite_element_;

    };
        
    template<Domain t_domain>
    struct MeshElementSet : ElementMap<MeshElement, t_domain>, DomainMap<MeshDomain, t_domain>
    {

    private:

        template<Element t_element>
        static
        void
        setElementOuterNeighborhood(
            std::shared_ptr<FiniteElement<t_element, t_domain>> const & finite_element
        )
        requires(t_element.isNode())
        {}

        template<Element t_element, Integer t_i = 0>
        static
        void
        setElementOuterNeighborhood(
            std::shared_ptr<FiniteElement<t_element, t_domain>> const & finite_element
        )
        requires(!t_element.isNode())
        {
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto const constexpr t_node_coordinates = ElementTraits<t_element>::template getInnerNeighborCoordinates<Element::node()>();
            for (auto const & nde : finite_element->template getInnerNeighbors<t_node_coordinates.getDim(), t_node_coordinates.getTag()>())
            {
                auto const & ngs = nde->template getOuterNeighbors<t_element_coordinates.getDim() - 1, t_i>();
                for (auto const & neighbour : ngs)
                {
                    if (((neighbour->getTag() != finite_element->getTag()) && t_i == t_element_coordinates.getTag()) || (t_i != t_element_coordinates.getTag()))
                    {
                        auto & element_ngs = finite_element->template getOuterNeighbors<0, t_i>();
                        auto found = false;
                        for (auto const & ngb: element_ngs)
                        {
                            if (ngb->getTag() == neighbour->getTag())
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
            if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<t_element_coordinates.getDim()>() - 1)
            {
                setElementOuterNeighborhood<t_element, t_i + 1>(finite_element);
            }
        }

    public:

        MeshElementSet()
        {}
        
        std::unique_ptr<FiniteElementSet<t_domain>>
        makeFiniteElementSet()
        const
        {
            auto element_map = std::make_unique<FiniteElementMap<t_domain>>();
            auto make_sets = [&] <Integer t_i = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & dm : this->template getDomains<t_i>())
                {
                    element_map->template getDomains<t_i>()[dm.second->getLabel()] = dm.second->makeDom();
                }
                if constexpr (t_i < t_domain.getDim())
                {
                    self.template operator()<t_i + 1>(self);
                }
            };
            make_sets(make_sets);
            auto make_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    element.second->template makeElement(* element_map);
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
            make_elements(make_elements);
            auto make_elements_outer_neighborhood = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto & element : element_map->template getElements<t_i, t_j>())
                {
                    setElementOuterNeighborhood<t_element>(element.second);
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

    template<Domain t_domain>
    struct MshFileParser : MeshFileParserBase, MeshElementSet<t_domain>
    {

        template<Integer t_dim, auto...>
        struct GeometricEntity
        {
        
            static
            std::basic_string<Character>
            getHash(
                Integer tag
            )
            {
                auto hash = std::basic_stringstream<Character>();
                hash << std::setfill('0') << std::setw(10) << t_dim;
                hash << std::setfill('0') << std::setw(10) << tag;
                return hash.str();
            }

            GeometricEntity()
            :
            tag_(),
            domains_()
            {}

            GeometricEntity(
                Integer tag
            )
            :
            tag_(tag),
            domains_()
            {}

            void
            addDomain(
                std::shared_ptr<MeshDomain<t_dim, t_domain>> const & d
            )
            {
                domains_.push_back(d);
            }

            std::vector<std::shared_ptr<MeshDomain<t_dim, t_domain>>> const &
            getDomains()
            const
            {
                return domains_;
            }

        private:

            Integer tag_;

            std::vector<std::shared_ptr<MeshDomain<t_dim, t_domain>>> domains_;

        };

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

        static
        std::basic_string<Character>
        getElemHash(
            Natural const & tag
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(20) << tag;
            return hash.str();
        }

        MshFileParser(
            std::basic_string_view<Character> str
        )
        :
        MeshFileParserBase(str)
        {
            setMeshDomain();
            setGeometricalEntities();
            setMeshElements();
        }
        
        void
        setMeshDomain()
        {
            auto const & file_lines = this->file_.lines_;
            auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
            auto offset = 1;
            auto line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
            auto num_physical_names = Integer();
            line_stream >> num_physical_names;
            offset += 1;
            for (auto i = 0; i < num_physical_names; ++i)
            {
                line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                auto dim = Integer();
                auto name = std::basic_string<Character>();
                auto physical_entity_tag = Integer();
                line_stream >> dim >> physical_entity_tag >> name;
                name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
                auto set_domain = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    if (t_i == dim)
                    {
                        auto physical_entity_hash = MeshDomain<t_i, t_domain>::getHash(physical_entity_tag);
                        this->template getDomains<t_i>()[physical_entity_hash] = std::make_shared<MeshDomain<t_i, t_domain>>(name);
                    }
                    if constexpr (t_i < t_domain.getDim())
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_domain(set_domain);
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
            for (auto i = 0; i < 4; ++i)
            {
                for (auto j = 0; j < num_domains[i]; ++j)
                {
                    auto set_geometric_entity = [&] <Integer t_i = 0> (
                        auto & self
                    )
                    constexpr mutable
                    {
                        if (t_i == i)
                        {
                            line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                            auto tag = Integer();
                            line_stream >> tag;
                            auto geometric_entity = std::make_shared<GeometricEntity<t_i, t_domain>>(tag);
                            if (i == 0)
                            {
                                for (auto k = 0; k < 3; ++k)
                                {
                                    auto a = Real();
                                    line_stream >> a;
                                }
                            }
                            else
                            {
                                for (auto k = 0; k < 6; ++k)
                                {
                                    auto a = Real();
                                    line_stream >> a;
                                }
                            }
                            auto num_physical_entities = Integer();
                            line_stream >> num_physical_entities;
                            for (auto k = 0; k < num_physical_entities; ++k)
                            {
                                auto physical_entity_tag = Integer();
                                line_stream >> physical_entity_tag;
                                auto physical_entity_hash = MeshDomain<t_i, t_domain>::getHash(physical_entity_tag);
                                geometric_entity->addDomain(this->template getDomains<t_i>().at(physical_entity_hash));
                            }
                            geometric_entities_.template getDomains<t_i>()[GeometricEntity<t_i, t_domain>::getHash(tag)] = geometric_entity;
                        }
                        if constexpr (t_i < t_domain.getDim())
                        {
                            self.template operator ()<t_i + 1>(self);
                        }
                    };
                    set_geometric_entity(set_geometric_entity);
                    offset += 1;
                }
            }
        }

        template<Element t_element>
        void
        setMeshElement()
        requires(t_element.isNode())
        {
            auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto & elements = this->template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
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
            for (auto i = 0; i < num_entity_block; ++i)
            {
                line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                auto entity_dim = Integer();
                auto entity_tag = Integer();
                auto parametric = Integer();
                auto num_nodes_in_block = Integer();
                line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                offset += 1;
                for (auto j = 0; j < num_nodes_in_block; ++j)
                {
                    auto mesh_element = std::make_shared<MeshElement<t_element, t_domain>>();
                    auto tag = Natural();
                    auto domains_ = std::vector<std::shared_ptr<std::basic_string<Character>>>();
                    line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                    line_stream >> tag;
                    mesh_element->setTag(tag - 1);
                    line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset + num_nodes_in_block]);
                    auto coordinate = Real();
                    for (auto k = 0; k < 3; ++k)
                    {
                        line_stream >> coordinate;
                        mesh_element->setCoordinate(k, coordinate);
                    }
                    if (entity_dim == t_element.getDim())
                    {
                        auto const entity_hash = GeometricEntity<t_element.getDim(), t_domain>::getHash(entity_tag);
                        for (auto const & physical_entity : geometric_entities_.template getDomains<t_element.getDim()>().at(entity_hash)->getDomains())
                        {
                            mesh_element->addDomain(physical_entity);
                        }
                    }
                    elements[getElemHash(tag)] = mesh_element;
                    offset += 1;
                }
                offset += num_nodes_in_block;
            }
        }
        
        template<Element t_element>
        void
        setMeshElement()
        requires(!t_element.isNode())
        {
            auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto & elements = this->template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
            auto const & file_lines = this->file_.lines_;
            auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Elements"));
            auto offset = 1;
            auto line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
            auto num_entity_blocks = Integer();
            auto num_elements = Integer();
            auto min_element_tag = Integer();
            auto max_element_tag = Integer();
            line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
            offset += 1;
            for (auto i = 0; i < num_entity_blocks; ++i)
            {
                line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                auto entity_dim = Integer();
                auto entity_tag = Integer();
                auto element_type_tag = Integer();
                auto num_elements_in_block = Integer();
                line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                offset += 1;
                auto const element = getElemType(element_type_tag);
                if (element == t_element)
                {
                    auto const entity_hash = GeometricEntity<t_element.getDim(), t_domain>::getHash(entity_tag);
                    for (auto j = 0; j < num_elements_in_block; ++j)
                    {
                        auto mesh_element = std::make_shared<MeshElement<t_element, t_domain>>();
                        line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                        auto tag = Natural();
                        line_stream >> tag;
                        auto node_tag = Integer();
                        for (auto k = 0; k < t_element.getNumNodes(); ++k)
                        {
                            line_stream >> node_tag;
                            mesh_element->setNodeTag(k, node_tag - 1);
                        }
                        for (auto const & physical_entity : geometric_entities_.template getDomains<t_element.getDim()>().at(entity_hash)->getDomains())
                        {
                            mesh_element->addDomain(physical_entity);
                        }
                        elements[getElemHash(tag)] = mesh_element;
                        offset += 1;
                    }
                }
                else
                {
                    offset += num_elements_in_block;
                }
            }
        }

        void
        setMeshElements()
        {

            auto set_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                setMeshElement<t_element>();
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
        }

    private:

        DomainMap<GeometricEntity, t_domain> geometric_entities_;

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
            return MshFileParser<t_domain>(file_path_).makeFiniteElementSet();
        }

        std::basic_string_view<Character> file_path_;

    };

    // template<Domain t_domain>
    // struct GmshFileParser : MeshFileParserBase<t_domain>
    // {

    // private:

    //     static constexpr
    //     Element
    //     getElemType(
    //         Integer tag
    //     )
    //     {
    //         if (tag == 15) return Element::node();
    //         else if (tag == 01) return Element::segment(1);
    //         else if (tag == 02) return Element::triangle(1);
    //         else if (tag == 03) return Element::quadrangle(1);
    //         else return Element::node();
    //     }

    //     static constexpr
    //     Integer
    //     getTagFromElement(
    //         Element element
    //     )
    //     {
    //         if (element.isNode()) return 15;
    //         else if (element.isSegment(1)) return 1;
    //         else if (element.isTriangle(1)) return 2;
    //         else if (element.isQuadrangle(1)) return 3;
    //         else return -1;
    //     }

    //     struct GeomEnt : DomainSet<MeshDomain, t_domain>
    //     {
            
    //         static
    //         std::basic_string<Character>
    //         getHash(
    //             Integer dim,
    //             Integer tag
    //         )
    //         {
    //             auto hash = std::basic_stringstream<Character>();
    //             hash << std::setfill('0') << std::setw(10) << dim;
    //             hash << std::setfill('0') << std::setw(10) << tag;
    //             return hash.str();
    //         }

    //         GeometricalEntity()
    //         {}

    //         template<Integer t_i>
    //         void
    //         add(
    //             std::shared_ptr<MeshDomain<t_i, t_domain>> const & physical_entity
    //         )
    //         {
    //             for (auto const & m : this->template getDomains<t_i>())
    //             {
    //                 if (m == physical_entity)
    //                 {
    //                     return;
    //                 }                    
    //             }
    //             this->template getDomains<t_i>().push_back(physical_entity);
    //         }

    //     };

    //     struct PhysicalEntity
    //     {

    //         static
    //         std::basic_string<Character>
    //         getHash(
    //             Integer tag
    //         )
    //         {
    //             auto hash = std::basic_stringstream<Character>();
    //             hash << std::setfill('0') << std::setw(10) << tag;
    //             return hash.str();
    //         }

    //         PhysicalEntity(
    //             Integer dim,
    //             std::basic_string<Character> const & label
    //         )
    //         :
    //         mesh_domain_(std::make_shared<MeshDomain>(dim, label))
    //         {}

    //         std::shared_ptr<MeshDomain> mesh_domain_;

    //     };

    //     struct GeometricalEntity
    //     {

    //         static
    //         std::basic_string<Character>
    //         getHash(
    //             Integer dim,
    //             Integer tag
    //         )
    //         {
    //             auto hash = std::basic_stringstream<Character>();
    //             hash << std::setfill('0') << std::setw(10) << dim;
    //             hash << std::setfill('0') << std::setw(10) << tag;
    //             return hash.str();
    //         }

    //         GeometricalEntity()
    //         :
    //         physical_entities_(),
    //         bounding_entities_()
    //         {}

    //         void
    //         addPhysicalEntity(
    //             std::shared_ptr<PhysicalEntity> const & physical_entity
    //         )
    //         {
    //             physical_entities_.push_back(physical_entity);
    //         }

    //         std::vector<std::shared_ptr<PhysicalEntity>> const &
    //         getPhysicalEntities()
    //         const
    //         {
    //             return physical_entities_;
    //         }

    //         // void
    //         // addPhysicalEntity(
    //         //     std::shared_ptr<PhysicalEntity> const & physical_entity
    //         // )
    //         // {
    //         //     physical_entities_.insert(physical_entity);
    //         //     for (auto & bounding_entity : bounding_entities_)
    //         //     {
    //         //         bounding_entity->addPhysicalEntity(physical_entity);
    //         //     }
    //         // }

    //         void
    //         addBoundingEntity(
    //             std::shared_ptr<GeometricalEntity> & bounding_entity
    //         )
    //         {
    //             bounding_entities_.push_back(bounding_entity);
    //         }

    //         std::vector<std::shared_ptr<GeometricalEntity>> const &
    //         getBoundingEntities()
    //         const
    //         {
    //             return bounding_entities_;
    //         }

    //         // void
    //         // addBoundingEntity(
    //         //     std::shared_ptr<GeometricalEntity> & bounding_entity
    //         // )
    //         // {
    //         //     bounding_entities_.insert(bounding_entity);
    //         //     for (auto const & physical_entity : physical_entities_)
    //         //     {
    //         //         bounding_entity->addPhysicalEntity(physical_entity);
    //         //     }
    //         // }

    //         std::vector<std::shared_ptr<PhysicalEntity>> physical_entities_;

    //         std::vector<std::shared_ptr<GeometricalEntity>> bounding_entities_;

    //         // std::set<std::shared_ptr<PhysicalEntity>> physical_entities_;

    //         // std::set<std::shared_ptr<GeometricalEntity>> bounding_entities_;

    //     };

    //     inline
    //     std::basic_string<Character>
    //     getMeshDomainHash(
    //         Integer tag
    //     )
    //     {
    //         auto hash = std::basic_stringstream<Character>();
    //         hash << std::setfill('0') << std::setw(10) << tag;
    //         return hash.str();
    //     }
        
    //     template<Integer t_i = 0>
    //     void
    //     setMeshDomain()
    //     {
    //         auto const & file_lines = this->file_.lines_;
    //         auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
    //         auto offset = 1;
    //         auto line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //         auto num_physical_names = Integer();
    //         line_stream >> num_physical_names;
    //         offset += 1;
    //         for (auto i = 0; i < num_physical_names; ++i)
    //         {
    //             line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //             auto dim = Integer();
    //             auto name = std::basic_string<Character>();
    //             auto physical_entity_tag = Integer();
    //             line_stream >> dim >> physical_entity_tag >> name;
    //             name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
    //             auto hash = getMeshDomainHash(physical_entity_tag);
    //             if (t_i == dim)
    //             {
    //                 this->template getDomains<t_i>()[hash] = std::make_shared<MeshDomain<t_i, t_domain>>(name);
    //             }
    //             offset += 1;
    //         }
    //         if constexpr (t_i < t_domain.getDim())
    //         {
    //             setMeshDomain<t_i + 1>();
    //         }
    //     }
        
    //     void
    //     setPhysicalEntities()
    //     {
    //         auto const & file_lines = this->file_.lines_;
    //         auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
    //         auto offset = 1;
    //         auto line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //         auto num_physical_names = Integer();
    //         line_stream >> num_physical_names;
    //         offset += 1;
    //         for (auto i = 0; i < num_physical_names; ++i)
    //         {
    //             line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //             auto dim = Integer();
    //             auto name = std::basic_string<Character>();
    //             auto physical_entity_tag = Integer();
    //             line_stream >> dim >> physical_entity_tag >> name;
    //             name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
    //             auto hash = PhysicalEntity::getHash(physical_entity_tag);
    //             physical_entities_[hash] = std::make_shared<PhysicalEntity>(dim, name);
    //             offset += 1;
    //         }
    //     }

    //     void
    //     setGeometricalEntities()
    //     {
    //         auto const & file_lines = this->file_.lines_;
    //         auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Entities"));
    //         auto offset = 1;
    //         auto line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //         auto num_points = Integer();
    //         auto num_curves = Integer();
    //         auto num_surfaces = Integer();
    //         auto num_volumes = Integer();
    //         line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
    //         auto num_domains = std::array<Integer, 4>{num_points, num_curves, num_surfaces, num_volumes};
    //         offset += 1;
    //         for (auto i = 0; i < 4; ++i)
    //         {
    //             for (auto j = 0; j < num_domains[i]; ++j)
    //             {
    //                 line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //                 auto tag = Integer();
    //                 line_stream >> tag;
    //                 auto geometrical_entity_hash = GeometricalEntity::getHash(i, tag);
    //                 geometrical_entities_[geometrical_entity_hash] = std::make_shared<GeometricalEntity>();
    //                 if (i == 0)
    //                 {
    //                     for (auto k = 0; k < 3; ++k)
    //                     {
    //                         auto a = Real();
    //                         line_stream >> a;
    //                     }
    //                 }
    //                 else
    //                 {
    //                     for (auto k = 0; k < 6; ++k)
    //                     {
    //                         auto a = Real();
    //                         line_stream >> a;
    //                     }
    //                 }
    //                 auto num_physical_entities = Integer();
    //                 line_stream >> num_physical_entities;
    //                 for (auto k = 0; k < num_physical_entities; ++k)
    //                 {
    //                     auto physical_entity_tag = Integer();
    //                     line_stream >> physical_entity_tag;
    //                     auto physical_entity_hash = PhysicalEntity::getHash(physical_entity_tag);
    //                     geometrical_entities_.at(geometrical_entity_hash)->addPhysicalEntity(physical_entities_.at(physical_entity_hash));
    //                 }
    //                 if (i > 0)
    //                 {
    //                     auto num_bounding_entities = Integer();
    //                     line_stream >> num_bounding_entities;
    //                     for (auto k = 0; k < num_bounding_entities; ++k)
    //                     {
    //                         auto bounding_entity_tag = Integer();
    //                         line_stream >> bounding_entity_tag;
    //                         auto bounding_entity_hash = GeometricalEntity::getHash(i - 1, std::abs(bounding_entity_tag));
    //                         geometrical_entities_.at(geometrical_entity_hash)->addBoundingEntity(geometrical_entities_.at(bounding_entity_hash));
    //                     }
    //                 }
    //                 offset += 1;
    //             }
    //         }
    //     }

    // public:

    //     GmshFileParser(
    //         std::basic_string_view<Character> str
    //     )
    //     :
    //     MeshFileParserBase<t_domain>(str)
    //     {
    //         setMeshDomain();
    //         setPhysicalEntities();
    //         setGeometricalEntities();
    //     }

    //     template<Element t_element>
    //     void
    //     setMeshElement(
    //         MeshElementSet<t_domain> & mesh_element_set
    //     )
    //     const
    //     requires(t_element.isNode())
    //     {
    //         auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
    //         auto & elements = mesh_element_set.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
    //         auto const & file_lines = this->file_.lines_;
    //         auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Nodes"));
    //         auto offset = 1;
    //         auto line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //         auto num_entity_block = Integer();
    //         auto num_nodes = Integer();
    //         auto min_node_tag = Integer();
    //         auto max_node_tag = Integer();
    //         line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
    //         offset += 1;
    //         for (auto i = 0; i < num_entity_block; ++i)
    //         {
    //             line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //             auto entity_dim = Integer();
    //             auto entity_tag = Integer();
    //             auto parametric = Integer();
    //             auto num_nodes_in_block = Integer();
    //             line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
    //             auto const entity_hash = GeometricalEntity::getHash(entity_dim, entity_tag);
    //             offset += 1;
    //             for (auto j = 0; j < num_nodes_in_block; ++j)
    //             {
    //                 auto mesh_element = std::make_shared<MeshElement<t_element, t_domain>>();
    //                 auto tag = Natural();
    //                 auto coordinates_ = Point();
    //                 auto domains_ = std::vector<std::shared_ptr<std::basic_string<Character>>>();
    //                 line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //                 line_stream >> tag;
    //                 mesh_element->tag_ = tag - 1;
    //                 line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset + num_nodes_in_block]);
    //                 for (auto k = 0; k < 3; ++k)
    //                 {
    //                     line_stream >> mesh_element->coordinates_->operator()(k);
    //                 }
    //                 for (auto const & physical_entity : geometrical_entities_.at(entity_hash)->getPhysicalEntities())
    //                 {
    //                     if (physical_entity->mesh_domain_->dim_ == t_element.getDim())
    //                     {
    //                         mesh_element->domains_.push_back(physical_entity->mesh_domain_);
    //                     }
    //                 }
    //                 elements.push_back(mesh_element);
    //                 // elements[mesh_element->getHash()] = mesh_element;
    //                 offset += 1;
    //             }
    //             offset += num_nodes_in_block;
    //         }
    //     }
        
    //     template<Element t_element>
    //     void
    //     setMeshElement(
    //         MeshElementSet<t_domain> & mesh_element_set
    //     )
    //     const
    //     requires(!t_element.isNode())
    //     {
    //         auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
    //         auto & elements = mesh_element_set.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
    //         auto const & file_lines = this->file_.lines_;
    //         auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Elements"));
    //         auto offset = 1;
    //         auto line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //         auto num_entity_blocks = Integer();
    //         auto num_elements = Integer();
    //         auto min_element_tag = Integer();
    //         auto max_element_tag = Integer();
    //         line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
    //         offset += 1;
    //         for (auto i = 0; i < num_entity_blocks; ++i)
    //         {
    //             line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //             auto entity_dim = Integer();
    //             auto entity_tag = Integer();
    //             auto element_type_tag = Integer();
    //             auto num_elements_in_block = Integer();
    //             line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
    //             auto const entity_hash = GeometricalEntity::getHash(entity_dim, entity_tag);
    //             offset += 1;
    //             auto const element = getElemType(element_type_tag);
    //             if (element == t_element)
    //             {
    //                 for (auto j = 0; j < num_elements_in_block; ++j)
    //                 {
    //                     auto mesh_element = std::make_shared<MeshElement<t_element, t_domain>>();
    //                     line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
    //                     auto tag = Integer();
    //                     line_stream >> tag;
    //                     for (auto k = 0; k < t_element.num_nodes_; ++k)
    //                     {
    //                         auto node_tag = Integer();
    //                         line_stream >> node_tag;
    //                         mesh_element->node_tags_[k] = node_tag - 1;
    //                     }
    //                     for (auto const & physical_entity : geometrical_entities_.at(entity_hash)->getPhysicalEntities())
    //                     {
    //                         if (physical_entity->mesh_domain_->dim_ == t_element.getDim())
    //                         {
    //                             mesh_element->domains_.push_back(physical_entity->mesh_domain_);
    //                         }
    //                     }
    //                     elements.push_back(mesh_element);
    //                     // elements[mesh_element->getHash()] = mesh_element;
    //                     offset += 1;
    //                 }
    //             }
    //             else
    //             {
    //                 offset += num_elements_in_block;
    //             }
    //         }
    //     }

    //     static
    //     void
    //     setOutput(
    //         std::basic_string<Character> const & file_path,
    //         std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
    //         auto... behavior_label
    //     )
    //     {
    //         auto labels = std::array<std::basic_string<Character>, sizeof...(behavior_label)>{behavior_label...};
    //         auto outfile = std::basic_ofstream<Character>();
    //         outfile.open(file_path);
    //         outfile << std::fixed << std::setprecision(17);
    //         outfile << "$MeshFormat\n";
    //         outfile << "2.2 0 8\n";
    //         outfile << "$EndMeshFormat\n";
    //         outfile << "$Nodes\n";
    //         outfile << element_set->template getNumElements<0>() + numerics::sum(element_set->template getNumIntegrationPoints<>(behavior_label)...) << "\n";
    //         auto c_node = Natural(1);
    //         for (auto const & node : element_set->template getElements<0, 0>())
    //         {
    //             auto const & coordinates = node->getCurrentCoordinates();
    //             outfile << c_node << " " << coordinates(0) << " " << coordinates(1) << " " << coordinates(2) << "\n";
    //             c_node ++;
    //         }
    //         auto set_integration_nodes = [&] <Integer t_i = 0, Integer t_j = 0> (
    //             auto & self
    //         )
    //         mutable
    //         {
    //             for (auto const & element : element_set->template getElements<t_i, t_j>())
    //             {
    //                 for (auto const & label : labels)
    //                 {
    //                     if (element->quadrature_.contains(label))
    //                     {
    //                         for (auto const & integration_point : element->quadrature_.at(label).ips_)
    //                         {
    //                             auto const & coordinates = integration_point.getCurrentCoordinates();
    //                             outfile << c_node << " " << coordinates(0) << " " << coordinates(1) << " " << coordinates(2) << "\n";
    //                             c_node ++;
    //                         }
    //                     }
    //                 }
    //             }
    //             if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
    //             {
    //                 self.template operator()<t_i, t_j + 1>(self);
    //             }
    //             else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
    //             {
    //                 self.template operator()<t_i + 1, 0>(self);
    //             }
    //         }; 
    //         set_integration_nodes(set_integration_nodes);
    //         outfile << "$EndNodes\n";
    //         outfile << "$Elements\n";
    //         outfile << element_set->template getNumElements<0>() + numerics::sum(element_set->template getNumIntegrationPoints<>(behavior_label)...) << "\n";
    //         auto c_element = Natural(1);
    //         for (auto const & node : element_set->template getElements<0, 0>())
    //         {
    //             outfile << c_element << " 15 2 0 0 " << c_element << "\n";
    //             c_element ++;
    //         }
    //         auto set_integration_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
    //             auto & self
    //         )
    //         mutable
    //         {
    //             for (auto const & element : element_set->template getElements<t_i, t_j>())
    //             {
    //                 for (auto const & label : labels)
    //                 {
    //                     if (element->quadrature_.contains(label))
    //                     {
    //                         for (auto const & integration_point : element->quadrature_.at(label).ips_)
    //                         {
    //                             outfile << c_element << " 15 2 1 1 " << c_element << "\n";
    //                             c_element ++;
    //                         }
    //                     }
    //                 }
    //             }
    //             if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
    //             {
    //                 self.template operator()<t_i, t_j + 1>(self);
    //             }
    //             else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
    //             {
    //                 self.template operator()<t_i + 1, 0>(self);
    //             }
    //         };
    //         set_integration_elements(set_integration_elements);
    //         auto set_elements = [&] <Integer t_i = 1, Integer t_j = 0> (
    //             auto & self
    //         )
    //         mutable
    //         {
    //             for (auto const & element : element_set->template getElements<t_i, t_j>())
    //             {
    //                 outfile << c_element << " " << getTagFromElement(DomainTraits<t_domain>::template getElement<t_i, t_j>()) << " 2 0 0";
    //                 for (auto const & node : element->template getInnerNeighbors<t_i - 1, 0>())
    //                 {
    //                     outfile << " " << node->getTag() + 1;
    //                 }
    //                 outfile << "\n";
    //                 c_element ++;
    //             }
    //             if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
    //             {
    //                 self.template operator()<t_i, t_j + 1>(self);
    //             }
    //             else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
    //             {
    //                 self.template operator()<t_i + 1, 0>(self);
    //             }
    //         };
    //         // set_elements(set_elements);
    //         outfile << "$EndElements\n";
    //     }

    //     template<Integer t_coordinate>
    //     static
    //     void
    //     addQuadratureStrainOutput(
    //         std::basic_string<Character> const & file_path,
    //         std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
    //         Integer time_step_index,
    //         Real time_step_value,
    //         std::basic_string<Character> behavior_label,
    //         Integer row
    //     )
    //     {
    //         if (!std::filesystem::exists(file_path))
    //         {
    //             throw std::runtime_error("File does not exist");
    //         }
    //         auto outfile = std::ofstream();
    //         auto c_element = Natural();
    //         outfile.open(file_path, std::ios_base::app);
    //         outfile << std::fixed << std::setprecision(17);
    //         outfile << "$NodeData\n";
    //         outfile << "1\n";
    //         outfile << "\"" << behavior_label << " " << row << " Strain\"\n";
    //         outfile << "1\n";
    //         outfile << time_step_value << "\n";
    //         outfile << "3\n";
    //         outfile << time_step_index << "\n";
    //         outfile << 1 << "\n"; // size of gradients
    //         outfile << element_set->template getNumIntegrationPoints<>(behavior_label) << "\n"; // number of quad pts
    //         c_element = element_set->template getNumElements<0>() + 1;
    //         auto set_strain = [&] <Integer t_j = 0> (
    //             auto & self
    //         )
    //         mutable
    //         {
    //             for (auto const & element : element_set->template getElements<t_coordinate, t_j>())
    //             {
    //                 if (element->quadrature_.contains(behavior_label))
    //                 {
    //                     for (auto const & integration_point : element->quadrature_.at(behavior_label).ips_)
    //                     {
    //                         outfile << c_element << " " << integration_point.behavior_data_->s1.gradients[row] << "\n";
    //                         c_element ++;
    //                     }
    //                 }
    //             }
    //             if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_coordinate>() - 1)
    //             {
    //                 self.template operator()<t_j + 1>(self);
    //             }
    //         };
    //         set_strain(set_strain);
    //         outfile << "$EndNodeData\n";
    //     }

    //     template<Integer t_coordinate>
    //     static
    //     void
    //     addQuadratureStressOutput(
    //         std::basic_string<Character> const & file_path,
    //         std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
    //         Integer time_step_index,
    //         Real time_step_value,
    //         std::basic_string<Character> behavior_label,
    //         Integer row
    //     )
    //     {
    //         if (!std::filesystem::exists(file_path))
    //         {
    //             throw std::runtime_error("File does not exist");
    //         }
    //         auto outfile = std::ofstream();
    //         auto c_element = Natural();
    //         outfile.open(file_path, std::ios_base::app);
    //         outfile << std::fixed << std::setprecision(17);
    //         outfile << "$NodeData\n";
    //         outfile << "1\n";
    //         outfile << "\"" << behavior_label << " " << row << " Stress\"\n";
    //         outfile << "1\n";
    //         outfile << time_step_value << "\n";
    //         outfile << "3\n";
    //         outfile << time_step_index << "\n";
    //         outfile << 1 << "\n"; // size of gradients
    //         outfile << element_set->template getNumIntegrationPoints<>(behavior_label) << "\n"; // number of quad pts
    //         c_element = element_set->template getNumElements<0>() + 1;
    //         auto set_strain = [&] <Integer t_j = 0> (
    //             auto & self
    //         )
    //         mutable
    //         {
    //             for (auto const & element : element_set->template getElements<t_coordinate, t_j>())
    //             {
    //                 if (element->quadrature_.contains(behavior_label))
    //                 {
    //                     for (auto const & integration_point : element->quadrature_.at(behavior_label).ips_)
    //                     {
    //                         outfile << c_element << " " << integration_point.behavior_data_->s1.thermodynamic_forces[row] << "\n";
    //                         c_element ++;
    //                     }
    //                 }
    //             }
    //             if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_coordinate>() - 1)
    //             {
    //                 self.template operator()<t_j + 1>(self);
    //             }
    //         };
    //         set_strain(set_strain);
    //         outfile << "$EndNodeData\n";
    //     }

    //     template<Integer t_coordinate>
    //     static
    //     void
    //     addQuadratureInternalVariableOutput(
    //         std::basic_string<Character> const & file_path,
    //         std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
    //         Integer time_step_index,
    //         Real time_step_value,
    //         std::basic_string<Character> behavior_label,
    //         Integer row
    //     )
    //     {
    //         if (!std::filesystem::exists(file_path))
    //         {
    //             throw std::runtime_error("File does not exist");
    //         }
    //         auto outfile = std::ofstream();
    //         auto c_element = Natural();
    //         outfile.open(file_path, std::ios_base::app);
    //         outfile << std::fixed << std::setprecision(17);
    //         outfile << "$NodeData\n";
    //         outfile << "1\n";
    //         outfile << "\"" << behavior_label << " " << row << " InternalVariable\"\n";
    //         outfile << "1\n";
    //         outfile << time_step_value << "\n";
    //         outfile << "3\n";
    //         outfile << time_step_index << "\n";
    //         outfile << 1 << "\n"; // size of gradients
    //         outfile << element_set->template getNumIntegrationPoints<>(behavior_label) << "\n"; // number of quad pts
    //         c_element = element_set->template getNumElements<0>() + 1;
    //         auto set_strain = [&] <Integer t_j = 0> (
    //             auto & self
    //         )
    //         mutable
    //         {
    //             for (auto const & element : element_set->template getElements<t_coordinate, t_j>())
    //             {
    //                 if (element->quadrature_.contains(behavior_label))
    //                 {
    //                     for (auto const & integration_point : element->quadrature_.at(behavior_label).ips_)
    //                     {
    //                         outfile << c_element << " " << integration_point.behavior_data_->s1.internal_state_variables[row] << "\n";
    //                         c_element ++;
    //                     }
    //                 }
    //             }
    //             if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_coordinate>() - 1)
    //             {
    //                 self.template operator()<t_j + 1>(self);
    //             }
    //         };
    //         set_strain(set_strain);
    //         outfile << "$EndNodeData\n";
    //     }

    //     template<Integer t_coordinate>
    //     static
    //     void
    //     addQuadratureExternalVariableOutput(
    //         std::basic_string<Character> const & file_path,
    //         std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
    //         Integer time_step_index,
    //         Real time_step_value,
    //         std::basic_string<Character> behavior_label,
    //         Integer row
    //     )
    //     {
    //         if (!std::filesystem::exists(file_path))
    //         {
    //             throw std::runtime_error("File does not exist");
    //         }
    //         auto outfile = std::ofstream();
    //         auto c_element = Natural();
    //         outfile.open(file_path, std::ios_base::app);
    //         outfile << std::fixed << std::setprecision(17);
    //         outfile << "$NodeData\n";
    //         outfile << "1\n";
    //         outfile << "\"" << behavior_label << " " << row << " InternalVariable\"\n";
    //         outfile << "1\n";
    //         outfile << time_step_value << "\n";
    //         outfile << "3\n";
    //         outfile << time_step_index << "\n";
    //         outfile << 1 << "\n"; // size of gradients
    //         outfile << element_set->template getNumIntegrationPoints<>(behavior_label) << "\n"; // number of quad pts
    //         c_element = element_set->template getNumElements<0>() + 1;
    //         auto set_strain = [&] <Integer t_j = 0> (
    //             auto & self
    //         )
    //         mutable
    //         {
    //             for (auto const & element : element_set->template getElements<t_coordinate, t_j>())
    //             {
    //                 if (element->quadrature_.contains(behavior_label))
    //                 {
    //                     for (auto const & integration_point : element->quadrature_.at(behavior_label).ips_)
    //                     {
    //                         outfile << c_element << " " << integration_point.behavior_data_->s1.external_state_variables[row] << "\n";
    //                         c_element ++;
    //                     }
    //                 }
    //             }
    //             if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_coordinate>() - 1)
    //             {
    //                 self.template operator()<t_j + 1>(self);
    //             }
    //         };
    //         set_strain(set_strain);
    //         outfile << "$EndNodeData\n";
    //     }

    //     template<Integer t_coordinate>
    //     static
    //     void
    //     addQuadratureDissipatedEnergyOutput(
    //         std::basic_string<Character> const & file_path,
    //         std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
    //         Integer time_step_index,
    //         Real time_step_value,
    //         std::basic_string<Character> behavior_label
    //     )
    //     {
    //         if (!std::filesystem::exists(file_path))
    //         {
    //             throw std::runtime_error("File does not exist");
    //         }
    //         auto outfile = std::ofstream();
    //         auto c_element = Natural();
    //         outfile.open(file_path, std::ios_base::app);
    //         outfile << std::fixed << std::setprecision(17);
    //         outfile << "$NodeData\n";
    //         outfile << "1\n";
    //         outfile << "\"" << behavior_label << " DissipatedEnergy\"\n";
    //         outfile << "1\n";
    //         outfile << time_step_value << "\n";
    //         outfile << "3\n";
    //         outfile << time_step_index << "\n";
    //         outfile << 1 << "\n"; // size of gradients
    //         outfile << element_set->template getNumIntegrationPoints<>(behavior_label) << "\n"; // number of quad pts
    //         c_element = element_set->template getNumElements<0>() + 1;
    //         auto set_strain = [&] <Integer t_j = 0> (
    //             auto & self
    //         )
    //         mutable
    //         {
    //             for (auto const & element : element_set->template getElements<t_coordinate, t_j>())
    //             {
    //                 if (element->quadrature_.contains(behavior_label))
    //                 {
    //                     for (auto const & integration_point : element->quadrature_.at(behavior_label).ips_)
    //                     {
    //                         outfile << c_element << " " << integration_point.behavior_data_->s1.dissipated_energy << "\n";
    //                         c_element ++;
    //                     }
    //                 }
    //             }
    //             if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_coordinate>() - 1)
    //             {
    //                 self.template operator()<t_j + 1>(self);
    //             }
    //         };
    //         set_strain(set_strain);
    //         outfile << "$EndNodeData\n";
    //     }

    //     template<Integer t_coordinate>
    //     static
    //     void
    //     addQuadratureStoredEnergyOutput(
    //         std::basic_string<Character> const & file_path,
    //         std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
    //         Integer time_step_index,
    //         Real time_step_value,
    //         std::basic_string<Character> behavior_label
    //     )
    //     {
    //         if (!std::filesystem::exists(file_path))
    //         {
    //             throw std::runtime_error("File does not exist");
    //         }
    //         auto outfile = std::ofstream();
    //         auto c_element = Natural();
    //         outfile.open(file_path, std::ios_base::app);
    //         outfile << std::fixed << std::setprecision(17);
    //         outfile << "$NodeData\n";
    //         outfile << "1\n";
    //         outfile << "\"" << behavior_label << " StoredEnergy\"\n";
    //         outfile << "1\n";
    //         outfile << time_step_value << "\n";
    //         outfile << "3\n";
    //         outfile << time_step_index << "\n";
    //         outfile << 1 << "\n"; // size of gradients
    //         outfile << element_set->template getNumIntegrationPoints<>(behavior_label) << "\n"; // number of quad pts
    //         c_element = element_set->template getNumElements<0>() + 1;
    //         auto set_strain = [&] <Integer t_j = 0> (
    //             auto & self
    //         )
    //         mutable
    //         {
    //             for (auto const & element : element_set->template getElements<t_coordinate, t_j>())
    //             {
    //                 if (element->quadrature_.contains(behavior_label))
    //                 {
    //                     for (auto const & integration_point : element->quadrature_.at(behavior_label).ips_)
    //                     {
    //                         outfile << c_element << " " << integration_point.behavior_data_->s1.stored_energy << "\n";
    //                         c_element ++;
    //                     }
    //                 }
    //             }
    //             if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_coordinate>() - 1)
    //             {
    //                 self.template operator()<t_j + 1>(self);
    //             }
    //         };
    //         set_strain(set_strain);
    //         outfile << "$EndNodeData\n";
    //     }
        
    //     template<Integer t_coordinate, auto... t_args>
    //     static
    //     void
    //     addNodalDofOutput(
    //         std::basic_string<Character> const & file_path,
    //         std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
    //         Integer time_step_index,
    //         Real time_step_value,
    //         std::basic_string<Character> unknown_label,
    //         Integer row,
    //         Integer col
    //     )
    //     {
    //         if (!std::filesystem::exists(file_path))
    //         {
    //             throw std::runtime_error("File does not exist");
    //         }
    //         auto outfile = std::ofstream();
    //         auto c_element = Natural();
    //         outfile.open(file_path, std::ios_base::app);
    //         outfile << std::fixed << std::setprecision(17);
    //         // writing strain
    //         outfile << "$NodeData\n";
    //         outfile << "1\n";
    //         outfile << "\"" << unknown_label << " " << row << " " << col << " NodalValues\"\n";
    //         outfile << "1\n";
    //         outfile << time_step_value << "\n";
    //         outfile << "3\n";
    //         outfile << time_step_index << "\n";
    //         outfile << 1 << "\n"; // size of gradients
    //         outfile << element_set->template getNumElements<0>() << "\n"; // number of nodes
    //         c_element = 1;
    //         auto nodal_values = element_set->template getNodalValues<t_coordinate, t_args...>(unknown_label, row, col);
    //         for (auto nodal_value : nodal_values)
    //         {
    //             outfile << c_element << " " << nodal_value << "\n";
    //             c_element ++;
    //         }
    //         outfile << "$EndNodeData\n";
    //     }
        
    //     template<Integer t_coordinate, auto... t_args>
    //     static
    //     void
    //     addQuadratureDofOutput(
    //         std::basic_string<Character> const & file_path,
    //         std::unique_ptr<FiniteElementSet<t_domain>> const & element_set,
    //         Integer time_step_index,
    //         Real time_step_value,
    //         std::basic_string<Character> unknown_label,
    //         std::basic_string<Character> quadrature_label,
    //         Integer row,
    //         Integer col
    //     )
    //     {
    //         if (!std::filesystem::exists(file_path))
    //         {
    //             throw std::runtime_error("File does not exist");
    //         }
    //         auto outfile = std::ofstream();
    //         auto c_element = Natural();
    //         outfile.open(file_path, std::ios_base::app);
    //         outfile << std::fixed << std::setprecision(17);
    //         // writing strain
    //         outfile << "$NodeData\n";
    //         outfile << "1\n";
    //         outfile << "\"" << unknown_label << " " << row << " " << col << " QuadratureValues\"\n";
    //         outfile << "1\n";
    //         outfile << time_step_value << "\n";
    //         outfile << "3\n";
    //         outfile << time_step_index << "\n";
    //         outfile << 1 << "\n"; // size of gradients
    //         outfile << element_set->template getNumIntegrationPoints<t_coordinate>(quadrature_label) << "\n"; // number of nodes
    //         c_element = element_set->template getNumElements<0>() + 1;
    //         auto quadrature_values = element_set->template getQuadratureValues<t_coordinate, t_args...>(unknown_label, quadrature_label, row, col);
    //         for (auto quadrature_value : quadrature_values)
    //         {
    //             outfile << c_element << " " << quadrature_value << "\n";
    //             c_element ++;
    //         }
    //         outfile << "$EndNodeData\n";
    //     }
        
    // private:

    //     std::vector<std::shared_ptr<GeomEnt>> geom_entities_;
        
    //     std::map<std::basic_string<Character>, std::shared_ptr<PhysicalEntity>> physical_entities_;

    //     std::map<std::basic_string<Character>, std::shared_ptr<GeometricalEntity>> geometrical_entities_;

    // };
    
} // namespace lolita



#endif /* CB51704E_040A_4B3D_AE3C_46C0DE3B8543 */
