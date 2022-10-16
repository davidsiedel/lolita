#ifndef E47C7236_24B0_44FD_8B1A_C03859454D62
#define E47C7236_24B0_44FD_8B1A_C03859454D62

#include "2/core/_include.hxx"
#include "2/core/region.hxx"
#include "2/core/element.hxx"
#include "2/core/mesh.hxx"

namespace lolita::core
{

    /**
     * @brief 
     * 
     * @tparam t_domain 
     */
    template<MeshConcept auto t_domain>
    struct FiniteElementMap : ElementMap<FiniteElement, t_domain>, DomainMap<FiniteDomain, t_domain>
    {

    private:

        /**
         * @brief 
         * 
         */
        using FiniteElementSet_ = FiniteElementSet<t_domain>;

    public:
    
        /**
         * @brief 
         * 
         * @return std::unique_ptr<FiniteElementSet_> 
         */
        std::unique_ptr<FiniteElementSet_>
        makeFiniteElementSet()
        const
        {
            auto finite_element_set = std::make_unique<FiniteElementSet_>();
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
                if constexpr (t_j < MeshTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < MeshTraits<t_domain>::template getNumElements<>() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            }; 
            make_elements(make_elements);
            return finite_element_set;
        }

    };

    /**
     * @brief A helper class to create a Domain object.
     * An object of this class is intended to be built by any mesh format parser.
     * 
     * @tparam t_dim The dimension of the domain to be built
     * @tparam t_domain The domain object defining the global domain geometry
     */
    template<DomainConcept auto t_dim, MeshConcept auto t_domain>
    struct MeshDomain
    {

    private:

        /**
         * @brief Alias for FiniteDomain 
         * 
         */
        using Domain_ = FiniteDomain<t_dim, t_domain>;

    public:

        /**
         * @brief Get the Hash object
         * 
         * @param label 
         * @return std::basic_string<Character> 
         */
        static
        std::basic_string<Character>
        getHash(
            std::basic_string<Character> && label
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << t_dim.getDim();
            hash << std::forward<std::basic_string<Character>>(label);
            return hash.str();
        }
        
        /**
         * @brief Get the Hash object
         * 
         * @param tag 
         * @return std::basic_string<Character> 
         */
        static
        std::basic_string<Character>
        getHash(
            Integer tag
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << t_dim.getDim();
            hash << std::setfill('0') << std::setw(10) << tag;
            return hash.str();
        }
        
        /**
         * @brief Construct a new Mesh Domain object
         * 
         */
        MeshDomain()
        :
        label_()
        {}

        /**
         * @brief Construct a new Mesh Domain object
         * 
         * @param label 
         */
        MeshDomain(
            std::basic_string<Character> const & label
        )
        :
        label_(label)
        {}

        /**
         * @brief Construct a new Mesh Domain object
         * 
         * @param label 
         */
        MeshDomain(
            std::basic_string<Character> && label
        )
        :
        label_(std::move(label))
        {}

        /**
         * @brief Get the Label object
         * 
         * @return std::basic_string<Character> const& 
         */
        std::basic_string<Character> const &
        getLabel()
        const
        {
            return label_;
        }

        /**
         * @brief 
         * 
         * @return std::shared_ptr<FiniteDomain<t_dim, t_domain>> 
         */
        std::shared_ptr<Domain_>
        letDomain()
        const
        {
            return std::make_shared<Domain_>(label_);
        }

        /**
         * @brief The name of the domain to be built
         * 
         */
        std::basic_string<Character> label_;
    
    };
    
    /**
     * @brief A helper class to create a Finite Element object for an element.
     * An object of this class is intended to be built by any mesh format parser.
     * 
     * @tparam t_element The element object defining the element geometry
     * @tparam t_domain The domain object defining the global domain geometry
     */
    template<ShapeConcept auto t_element, MeshConcept auto t_domain>
    struct MeshElement;
    
    /**
     * @brief A helper class to create a Finite Element object for an element.
     * An object of this class is intended to be built by any mesh format parser.
     * For a generic element that is not a node, it consists in defining
     * the list of all node tags the element is connected to, and the list of all domains it belongs to.
     * 
     * @tparam t_element The element object defining the element geometry
     * @tparam t_domain The domain object defining the global domain geometry
     */
    template<ShapeConcept auto t_element, MeshConcept auto t_domain>
    requires(t_element != Node())
    struct MeshElement<t_element, t_domain>
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element>;

        using FiniteElement_ = FiniteElement<t_element, t_domain>;

        using NodeConnectivity_ = ElementNodeConnectivity<t_element>;
        
        template<DomainConcept auto t_dim, MeshConcept auto t__domain>
        using t_Dom = MeshDomain<t_dim, t__domain>;

        /**
         * @brief The MeshDomain the element is connected to.
         * 
         */
        using MeshDomain_ = typename t_ElementTraits::template DomainConnectivity<t_Dom, t_domain>;

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
        finite_element_()
        {}

        /**
         * @brief Set the domain pointer that the element is connected to.
         * 
         * @param domain The domain pointer to be set.
         */
        void
        setDomain(
            std::shared_ptr<MeshDomain_> const & domain
        )
        {
            domain_ = domain;
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
                auto const constexpr t_element_coordinates = MeshTraits<t_domain>::template getElementCoordinates<t_element>();
                auto const constexpr t_node_coordinates = ElementTraits<t_element>::template getInnerNeighborCoordinates<Node{}>();
                auto const constexpr t_inner_neighbor_coordinates = MeshTraits<t_domain>::template getElementCoordinates<t_inner_neighbor>();
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
                    if constexpr(t_inner_neighbor != Node())
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
                    // for (auto const & domain : domains_)
                    // {
                    //     finite_element_->addDomain(element_map.template getDomains<t_element.getDim()>().at(domain->getLabel()));
                    // }
                    if (domain_ != nullptr)
                    {
                        finite_element_->setDomain(element_map.template getDomains<t_element.getDim()>().at(domain_->getLabel()));
                    }                    
                    // finite_element_->setDomain(element_map.template getDomains<t_element.getDim()>().at(domain_->getLabel()));
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
        NodeConnectivity_ node_tags_;

        /**
         * @brief A pointer to the FiniteElement to be build.
         * 
         */
        std::shared_ptr<FiniteElement_> finite_element_;

        /**
         * @brief The MeshDomain the element is connected to.
         * 
         */
        std::shared_ptr<MeshDomain_> domain_;

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
    template<ShapeConcept auto t_element, MeshConcept auto t_domain>
    requires(t_element == Node())
    struct MeshElement<t_element, t_domain>
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element>;
        
        template<DomainConcept auto t_dim, MeshConcept auto t__domain>
        using t_Dom = MeshDomain<t_dim, t__domain>;

        /**
         * @brief The MeshDomain the element is connected to.
         * 
         */
        using MeshDomain_ = typename t_ElementTraits::template DomainConnectivity<t_Dom, t_domain>;

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
        finite_element_()
        {}
        
        /**
         * @brief Set the domain pointer that the element is connected to.
         * 
         * @param domain The domain pointer to be set.
         */
        void
        setDomain(
            std::shared_ptr<MeshDomain_> const & domain
        )
        {
            domain_ = domain;
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
            auto constexpr t_element_coordinates = MeshTraits<t_domain>::template getElementCoordinates<t_element>();
            finite_element_ = std::make_shared<FiniteElement<t_element, t_domain>>(FiniteElement<t_element, t_domain>(tag_));
            finite_element_->setCoordinates(coordinates_);
            if (domain_ != nullptr)
            {
                finite_element_->setDomain(element_set.template getDomains<t_element.getDim()>().at(domain_->getLabel()));
            }
            // finite_element_->setDomain(element_set.template getDomains<t_element.getDim()>().at(domain_->getLabel()));
            // for (auto const & domain : domains_)
            // {
            //     finite_element_->addDomain(element_set.template getDomains<t_element.getDim()>().at(domain->getLabel()));
            // }
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
         * @brief A pointer to the FiniteElement to be build.
         * 
         */
        std::shared_ptr<FiniteElement<t_element, t_domain>> finite_element_;

        /**
         * @brief The MeshDomain the element is connected to.
         * 
         */
        std::shared_ptr<MeshDomain_> domain_;

    };
        
    template<MeshConcept auto t_domain>
    struct MeshElementSet : ElementMap<MeshElement, t_domain>, DomainMap<MeshDomain, t_domain>
    {

    private:

        template<ShapeConcept auto t_element>
        static
        void
        setElementOuterNeighborhood(
            std::shared_ptr<FiniteElement<t_element, t_domain>> const & finite_element
        )
        requires(t_element == Node())
        {}

        template<ShapeConcept auto t_element, Integer t_i = 0>
        static
        void
        setElementOuterNeighborhood(
            std::shared_ptr<FiniteElement<t_element, t_domain>> const & finite_element
        )
        requires(t_element != Node())
        {
            auto const constexpr t_element_coordinates = MeshTraits<t_domain>::template getElementCoordinates<t_element>();
            auto const constexpr t_node_coordinates = ElementTraits<t_element>::template getInnerNeighborCoordinates<Node{}>();
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
            if constexpr (t_i < MeshTraits<t_domain>::template getNumElements<t_element_coordinates.getDim()>() - 1)
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
                    element_map->template getDomains<t_i>()[dm.second->getLabel()] = dm.second->letDomain();
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
                auto constexpr t_element = MeshTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    element.second->template makeElement(* element_map);
                }
                if constexpr (t_j < MeshTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < MeshTraits<t_domain>::getNumElements() - 1)
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
                auto const constexpr t_element = MeshTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto & element : element_map->template getElements<t_i, t_j>())
                {
                    setElementOuterNeighborhood<t_element>(element.second);
                }
                if constexpr (t_j < MeshTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < MeshTraits<t_domain>::getNumElements() - 1)
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

    template<MeshConcept auto t_domain>
    struct MshFileParser : MeshFileParserBase, MeshElementSet<t_domain>
    {

        template<DomainConcept auto t_dim, auto...>
        struct GeometricEntity
        {

            using MeshDomain_ = MeshDomain<t_dim, t_domain>;
        
            static
            std::basic_string<Character>
            getHash(
                Integer tag
            )
            {
                auto hash = std::basic_stringstream<Character>();
                hash << std::setfill('0') << std::setw(10) << t_dim.getDim();
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
                std::shared_ptr<MeshDomain_> const & d
            )
            {
                domains_.push_back(d);
            }

            std::vector<std::shared_ptr<MeshDomain_>> const &
            getDomains()
            const
            {
                return domains_;
            }

            void
            setDomain(
                std::shared_ptr<MeshDomain_> const & d
            )
            {
                domain_ = d;
            }

            std::shared_ptr<MeshDomain_> const &
            getDomain()
            const
            {
                return domain_;
            }

        private:

            Integer tag_;

            std::vector<std::shared_ptr<MeshDomain_>> domains_;

            std::shared_ptr<MeshDomain_> domain_;

        };

        // static constexpr
        // ShapeConcept auto
        // getElemType(
        //     Integer tag
        // )
        // {
        //     if (tag == 15) return Node();
        //     else if (tag == 01) return LinearSegment();
        //     else if (tag == 02) return LinearTriangle();
        //     else if (tag == 03) return LinearQuadrangle();
        //     else return Node();
        // }

        // static constexpr
        // Integer
        // getTagFromElement(
        //     ShapeConcept auto const & element
        // )
        // {
        //     if (element == Node()) return 15;
        //     else if (element == LinearSegment()) return 1;
        //     else if (element == LinearTriangle()) return 2;
        //     else if (element == LinearQuadrangle()) return 3;
        //     else return -1;
        // }

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
                        auto constexpr dom_n = MeshTraits<t_domain>::template getDomain<t_i>();
                        auto physical_entity_hash = MeshDomain<dom_n, t_domain>::getHash(physical_entity_tag);
                        this->template getDomains<t_i>()[physical_entity_hash] = std::make_shared<MeshDomain<dom_n, t_domain>>(name);
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
                            auto constexpr dom_n = MeshTraits<t_domain>::template getDomain<t_i>();
                            line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                            auto tag = Integer();
                            line_stream >> tag;
                            auto geometric_entity = std::make_shared<GeometricEntity<dom_n, t_domain>>(tag);
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
                                if (k == 0)
                                {
                                    auto physical_entity_tag = Integer();
                                    line_stream >> physical_entity_tag;
                                    auto physical_entity_hash = MeshDomain<dom_n, t_domain>::getHash(physical_entity_tag);
                                    geometric_entity->setDomain(this->template getDomains<t_i>().at(physical_entity_hash));
                                }
                                // auto physical_entity_tag = Integer();
                                // line_stream >> physical_entity_tag;
                                // auto physical_entity_hash = MeshDomain<t_i, t_domain>::getHash(physical_entity_tag);
                                // geometric_entity->addDomain(this->template getDomains<t_i>().at(physical_entity_hash));
                            }
                            geometric_entities_.template getDomains<t_i>()[GeometricEntity<dom_n, t_domain>::getHash(tag)] = geometric_entity;
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

        template<ShapeConcept auto t_element>
        void
        setMeshElement()
        requires(t_element == Node())
        {
            auto constexpr t_element_coordinates = MeshTraits<t_domain>::template getElementCoordinates<t_element>();
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
                        auto constexpr dom_n = MeshTraits<t_domain>::template getDomain<t_element.getDim()>();
                        auto const entity_hash = GeometricEntity<dom_n, t_domain>::getHash(entity_tag);
                        auto const & domain = geometric_entities_.template getDomains<t_element.getDim()>().at(entity_hash)->getDomain();
                        mesh_element->setDomain(domain);
                        // for (auto const & physical_entity : geometric_entities_.template getDomains<t_element.getDim()>().at(entity_hash)->getDomains())
                        // {
                        //     mesh_element->addDomain(physical_entity);
                        // }
                    }
                    elements[getElemHash(tag)] = mesh_element;
                    offset += 1;
                }
                offset += num_nodes_in_block;
            }
        }
        
        template<ShapeConcept auto t_element>
        void
        setMeshElement()
        requires(t_element != Node())
        {
            auto constexpr t_element_coordinates = MeshTraits<t_domain>::template getElementCoordinates<t_element>();
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
                // auto const element = getElemType(element_type_tag);
                if (element_type_tag == t_element.getMshTag())
                {
                    auto constexpr dom_n = MeshTraits<t_domain>::template getDomain<t_element.getDim()>();
                    auto const entity_hash = GeometricEntity<dom_n, t_domain>::getHash(entity_tag);
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
                        auto const & domain = geometric_entities_.template getDomains<t_element.getDim()>().at(entity_hash)->getDomain();
                        mesh_element->setDomain(domain);
                        // for (auto const & physical_entity : geometric_entities_.template getDomains<t_element.getDim()>().at(entity_hash)->getDomains())
                        // {
                        //     mesh_element->addDomain(physical_entity);
                        // }
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
                auto const constexpr t_element = MeshTraits<t_domain>::template getElement<t_i, t_j>();
                setMeshElement<t_element>();
                if constexpr (t_j < MeshTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < MeshTraits<t_domain>::template getNumElements<>() - 1)
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

        template<MeshConcept auto t_domain>
        std::unique_ptr<FiniteElementSet<t_domain>>
        makeFiniteElementSet()
        const
        {
            return MshFileParser<t_domain>(file_path_).makeFiniteElementSet();
        }

        std::basic_string_view<Character> file_path_;

    };

} // namespace lolita::core

#endif /* E47C7236_24B0_44FD_8B1A_C03859454D62 */
