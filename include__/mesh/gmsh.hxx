#ifndef C5421678_DA6F_4EA0_AFC9_73889C52D424
#define C5421678_DA6F_4EA0_AFC9_73889C52D424

#include "geometry/domain.hxx"
#include "geometry/point.hxx"
#include "geometry/shape.hxx"
#include "geometry/frame.hxx"
#include "mesh/element.hxx"
#include "mesh/region.hxx"
#include "mesh/factory.hxx"

namespace lolita::mesh
{

    struct Msh
    {

        static constexpr
        std::basic_string_view<Character>
        getFileNameExtension()
        {
            return "msh";
        }

    };

    template<>
    struct Table<geometry::Node, Msh>
    {

        static constexpr
        Integer tag_ = 15;

    };

    template<>
    struct Table<geometry::Segment, Msh>
    {

        static constexpr
        Integer tag_ = 1;

        static constexpr
        ShapeInnerNeighborhoodNodeConnectivity<geometry::Segment> node_connectivity_ = {
            {
                {
                    0,
                    1,
                }
            }
        };

    };

    template<>
    struct Table<geometry::Triangle, Msh>
    {

        static constexpr
        Integer tag_ = 2;

        static constexpr
        ShapeInnerNeighborhoodNodeConnectivity<geometry::Triangle> node_connectivity_ = {
            {
                {
                    0, 1,
                    1, 2,
                    2, 0,
                }
            },
            {
                {
                    0,
                    1,
                    2,
                }
            }
        };

    };

    template<>
    struct Table<geometry::Quadrangle, Msh>
    {

        static constexpr
        Integer tag_ = 3;

        static constexpr
        ShapeInnerNeighborhoodNodeConnectivity<geometry::Quadrangle> node_connectivity_ = {
            {
                {
                    0, 1,
                    1, 2,
                    2, 3,
                    3, 0,
                }
            },
            {
                {
                    0,
                    1,
                    2,
                    3,
                }
            }
        };

    };

    template<typename T_>
    concept ReaderConcept = requires(
        std::remove_reference_t<T_> const & t
    )
    {
        // TODO : express the fact that is has these member functions...
        // { t.template setMeshDomain<Domain_>(m) } -> std::same_as<void>;
        // { t.template setMeshElement<Shape_>(m) } -> std::same_as<void>;
        // { t.file_ } -> std::same_as<File>;
        // requires(std::same_as<decltype(t.file_), File>);
        requires(true);
    };

    template<geometry::DomainConcept Domain_, geometry::FrameConcept Frame_>
    struct GeometricEntity
    {
    
        static
        std::basic_string<Character>
        makeHash(
            Integer tag
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << Domain_::getDimDomain();
            hash << std::setfill('0') << std::setw(10) << tag;
            return hash.str();
        }

        GeometricEntity()
        :
        physical_entity_tags_()
        {}

        void
        addPhysicalEntityTag(
            Integer tag
        )
        {
            physical_entity_tags_.push_back(tag);
        }

        std::vector<Integer> const &
        getPhysicalEntityTags()
        const
        {
            return physical_entity_tags_;
        }

    private:

        std::vector<Integer> physical_entity_tags_;

    };

    template<geometry::FrameConcept Frame_>
    struct MshOne<Frame_, Msh>
    {

        template<typename... T_>
        using GeometricEntities1 = std::unordered_map<std::basic_string<Character>, GeometricEntity<T_...>>;

        using GeometricEntities_ = geometry::DomainCollection<Frame_, GeometricEntities1, Frame_>;

        explicit
        MshOne(
            std::basic_string<Character> const & file_path
        )
        :
        file_(file_path)
        {
            setGeometricalEntities();
        }

        explicit
        MshOne(
            std::basic_string<Character> && file_path
        )
        :
        file_(std::move(file_path))
        {
            setGeometricalEntities();
        }

        template<geometry::DomainConcept Domain_>
        void
        addGeometricEntity(
            Integer tag,
            GeometricEntity<Domain_, Frame_> const & geometric_entity
        )
        {
            geometric_entities_.template getComponent<Domain_>()[GeometricEntity<Domain_, Frame_>::makeHash(tag)] = geometric_entity;
        }

        template<geometry::DomainConcept Domain_>
        GeometricEntity<Domain_, Frame_> const &
        getGeometricEntity(
            Integer tag
        )
        const
        {
            return geometric_entities_.template getComponent<Domain_>().at(GeometricEntity<Domain_, Frame_>::makeHash(tag));
        }
        
        void
        setGeometricalEntities()
        {
            auto const & file_lines = this->file_.getContent();
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
                    auto set_geometric_entity = [&] <Integer i_ = 0> (
                        auto & set_geometric_entity_
                    )
                    constexpr mutable
                    {
                        if (i_ == i)
                        {
                            using Domain_ = geometry::Domain<i_>;
                            line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                            auto geometric_entity_tag = Integer();
                            line_stream >> geometric_entity_tag;
                            auto geometric_entity = GeometricEntity<Domain_, Frame_>();
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
                                geometric_entity.addPhysicalEntityTag(physical_entity_tag);
                            }
                            addGeometricEntity(geometric_entity_tag, geometric_entity);
                        }
                        if constexpr (i_ < Frame_::getDimEuclidean())
                        {
                            set_geometric_entity_.template operator ()<i_ + 1>(set_geometric_entity_);
                        }
                    };
                    set_geometric_entity(set_geometric_entity);
                    offset += 1;
                }
            }
        }

        template<geometry::DomainConcept Domain_>
        void
        setMeshDomain(
            ElementFactoryMap<Frame_> & my_map_
        )
        const
        {
            auto const & file_lines = this->file_.getContent();
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
                if (dim == Domain_::getDimDomain())
                {
                    my_map_.template addRegion<Domain_>(physical_entity_tag, name);
                }
                offset += 1;
            }
        }

        template<geometry::ShapeConcept Shape_>
        void
        setMeshElement(
            ElementFactoryMap<Frame_> & my_map_
        )
        const
        requires(std::same_as<Shape_, geometry::Node>)
        {
            auto const & file_lines = this->file_.getContent();
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
                    auto node_tag = Natural();
                    auto domains_ = std::vector<std::shared_ptr<std::basic_string<Character>>>();
                    line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                    line_stream >> node_tag;
                    line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset + num_nodes_in_block]);
                    auto node_coordinates = geometry::Point<Frame_::getDimEuclidean()>();
                    for (auto k = 0; k < Frame_::getDimEuclidean(); ++k)
                    {
                        auto coordinate = Real();
                        line_stream >> coordinate;
                        node_coordinates.setCoordinate(k, coordinate);
                    }
                    my_map_.template addElement<Shape_>(node_tag - 1, node_coordinates);
                    if (entity_dim == Shape_::getDimShape())
                    {
                        using Domain_ = geometry::Domain<Shape_::getDimShape()>;
                        auto const & geometric_entity = getGeometricEntity<Domain_>(entity_tag);
                        for (auto const & physical_entity_tag : geometric_entity.getPhysicalEntityTags())
                        {
                            my_map_.template addRegionToElement<Shape_>(physical_entity_tag, node_tag - 1);
                        }
                    }
                    offset += 1;
                }
                offset += num_nodes_in_block;
            }
        }

        template<geometry::ShapeConcept Shape_>
        void
        setMeshElement(
            ElementFactoryMap<Frame_> & my_map_
        )
        const
        {
            auto const & file_lines = this->file_.getContent();
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
                if (element_type_tag == Table<Shape_, Msh>::tag_)
                {
                    using Domain_ = geometry::Domain<Shape_::getDimShape()>;
                    for (auto j = 0; j < num_elements_in_block; ++j)
                    {
                        line_stream = std::basic_stringstream<Character>(file_lines[line_start + offset]);
                        auto element_tag = Natural();
                        line_stream >> element_tag;
                        auto node_tags = ShapeNodeConnectivity<Shape_>();
                        for (auto k = 0; k < Shape_::getNumNodes(); ++k)
                        {
                            auto node_tag = Integer();
                            line_stream >> node_tag;
                            node_tags[k] = node_tag - 1;
                        }
                        my_map_.template addElement<Shape_>(element_tag, node_tags);
                        auto const & geometric_entity = getGeometricEntity<Domain_>(entity_tag);
                        for (auto const & physical_entity_tag : geometric_entity.getPhysicalEntityTags())
                        {
                            my_map_.template addRegionToElement<Shape_>(physical_entity_tag, element_tag);
                        }
                        offset += 1;
                    }
                }
                else
                {
                    offset += num_elements_in_block;
                }
            }
        }

    private:

        File file_;

        GeometricEntities_ geometric_entities_;

    };

    // template<geometry::FrameConcept Frame_>
    // struct MyElementMap
    // {

    // private:

    //     template<typename... T_>
    //     using ElementPointer_ = std::shared_ptr<Element<T_...>>;

    //     template<typename... T_>
    //     using RegionPointer_ = std::shared_ptr<Region<T_...>>;

    //     using Elements_ = ElementMap<ElementPointer_, Frame_>;

    //     using Regions_ = RegionMap<RegionPointer_, Frame_>;
    
    //     using MyElementSet_ = MyElementSet<Frame_>;

    // public:

    //     MyElementMap()
    //     :
    //     my_map_(),
    //     elements_(),
    //     regions_()
    //     {}
        
    //     std::unique_ptr<MyElementSet_>
    //     makeElementSet()
    //     const
    //     {
    //         auto finite_element_set = std::make_unique<MyElementSet_>();
    //         auto make_sets = [&] <Integer i_ = 0> (
    //             auto & make_sets_
    //         )
    //         mutable
    //         {
    //             for (auto const & dm : regions_->template getDomains<i_>())
    //             {
    //                 finite_element_set->regions_.template getDomains<i_>().push_back(dm.second);
    //             }
    //             if constexpr (i_ < Frame_::getDimEuclidean())
    //             {
    //                 make_sets_.template operator()<i_ + 1>(make_sets_);
    //             }
    //         };
    //         make_sets(make_sets);
    //         auto make_elements = [&] <Integer i_ = 0, Integer j_ = 0> (
    //             auto & make_elements_
    //         )
    //         mutable
    //         {
    //             for (auto const & element : elements_->template getElements<i_, j_>())
    //             {
    //                 finite_element_set->elements_.template getElements<i_, j_>().push_back(element.second);
    //             }
    //             if constexpr (j_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<i_>() - 1)
    //             {
    //                 make_elements_.template operator()<i_, j_ + 1>(make_elements_);
    //             }
    //             else if constexpr (i_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<>() - 1)
    //             {
    //                 make_elements_.template operator()<i_ + 1, 0>(make_elements_);
    //             }
    //         }; 
    //         make_elements(make_elements);
    //         return finite_element_set;
    //     }

    //     template<typename Mesh>
    //     void
    //     setIt(
    //         std::basic_string<Character> && file_path
    //     )
    //     {
    //         auto my_one = MshOne<Frame_, Mesh>(std::forward<std::basic_string<Character>>(file_path));
    //         auto set_regions = [&] <Integer i_ = 0> (
    //             auto & set_regions_
    //         )
    //         mutable
    //         {
    //             using Domain_ = typename geometry::DomainLibrary<Frame_>::template Domain<i_>;
    //             my_one.template setMeshDomain<Domain_>(my_map_);
    //             if constexpr (i_ < geometry::DomainLibrary<Frame_>::getNumDomains() - 1)
    //             {
    //                 set_regions_.template operator()<i_ + 1>(set_regions_);
    //             }
    //         };
    //         set_regions(set_regions);
    //         auto set_elements = [&] <Integer i_ = 0, Integer j_ = 0> (
    //             auto & set_elements_
    //         )
    //         mutable
    //         {
    //             using Shape_ = geometry::ShapeLibraryTraits<Frame_>::template Shape<i_, j_>;
    //             my_one.template setMeshElement<Shape_>(my_map_);
    //             if constexpr (j_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<i_>() - 1)
    //             {
    //                 set_elements_.template operator()<i_, j_ + 1>(set_elements_);
    //             }
    //             else if constexpr (i_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<>() - 1)
    //             {
    //                 set_elements_.template operator()<i_ + 1, 0>(set_elements_);
    //             }
    //         };
    //         set_elements(set_elements);
    //     }
        
    //     // void
    //     // makeElement(
    //     //     FiniteElementMap<t_domain> & element_map
    //     // )
    //     // {
    //     //     auto make_element = [&] <Integer t_i = 0, Integer t_j = 0> (
    //     //         auto & self
    //     //     )
    //     //     mutable
    //     //     {
    //     //         auto const constexpr t_inner_neighbor = ShapeTraits<t_element>::template getInnerNeighbor<t_i, t_j>();
    //     //         auto const constexpr t_is_initialized = t_i == 0 && t_j == 0;
    //     //         auto const constexpr t_element_coordinates = MeshTraits<t_domain>::template getElementCoordinates<t_element>();
    //     //         auto const constexpr t_node_coordinates = ShapeTraits<t_element>::template getInnerNeighborCoordinates<Node{}>();
    //     //         auto const constexpr t_inner_neighbor_coordinates = MeshTraits<t_domain>::template getElementCoordinates<t_inner_neighbor>();
    //     //         auto const constexpr t_neighbour_coordinates = ShapeTraits<t_inner_neighbor>::template getOuterNeighborCoordinates<t_domain, t_element>();
    //     //         auto & inner_neighbors = element_map.template getElements<t_inner_neighbor_coordinates.getDim(), t_inner_neighbor_coordinates.getTag()>();
    //     //         if constexpr (t_is_initialized)
    //     //         {
    //     //             auto tag = element_map.template getElements<t_element_coordinates.getDim(), t_element_coordinates.getTag()>().size();
    //     //             finite_element_ = std::make_shared<FiniteElement<t_element, t_domain>>(FiniteElement<t_element, t_domain>(tag));
    //     //         }
    //     //         for (auto i = 0; i < ShapeTraits<t_element>::template getNumInnerNeighbors<t_i, t_j>(); ++i)
    //     //         {
    //     //             auto inner_neighbor_hash = std::basic_string<Character>();
    //     //             if constexpr(t_inner_neighbor != Node())
    //     //             {
    //     //                 auto mesh_inner_neighbor = MeshElement<t_inner_neighbor, t_domain>();
    //     //                 for (auto j = 0; j < t_inner_neighbor.num_nodes_; ++j)
    //     //                 {
    //     //                     mesh_inner_neighbor.setNodeTag(j, getNodeTag(getInnerNeighborNodeConnection<t_i, t_j>(i, j)));
    //     //                 }
    //     //                 inner_neighbor_hash = mesh_inner_neighbor.getHash();
    //     //                 if (!inner_neighbors.contains(inner_neighbor_hash))
    //     //                 {
    //     //                     mesh_inner_neighbor.makeElement(element_map);
    //     //                 }
    //     //             }
    //     //             else
    //     //             {
    //     //                 inner_neighbor_hash = MeshElement<t_inner_neighbor, t_domain>::getHash(getNodeTag(getInnerNeighborNodeConnection<t_i, t_j>(i, 0)));
    //     //             }
    //     //             auto & inner_neighbor = inner_neighbors.at(inner_neighbor_hash);
    //     //             finite_element_->template getInnerNeighbors<t_i, t_j>()[i] = inner_neighbor;
    //     //             inner_neighbor->template getOuterNeighbors<t_neighbour_coordinates.dim_, t_neighbour_coordinates.tag_>().push_back(finite_element_);
    //     //         }
    //     //         if constexpr (t_j < ShapeTraits<t_element>::template getNumInnerNeighbors<t_i>() - 1)
    //     //         {
    //     //             self.template operator()<t_i, t_j + 1>(self);
    //     //         }
    //     //         else if constexpr (t_i < ShapeTraits<t_element>::template getNumInnerNeighbors<>() - 1)
    //     //         {
    //     //             self.template operator()<t_i + 1, 0>(self);
    //     //         }
    //     //         if constexpr (t_is_initialized)
    //     //         {
    //     //             // for (auto const & domain : domains_)
    //     //             // {
    //     //             //     finite_element_->addDomain(element_map.template getDomains<t_element.getDim()>().at(domain->getLabel()));
    //     //             // }
    //     //             if (domain_ != nullptr)
    //     //             {
    //     //                 finite_element_->setDomain(element_map.template getDomains<t_element.getDim()>().at(domain_->getLabel()));
    //     //             }                    
    //     //             // finite_element_->setDomain(element_map.template getDomains<t_element.getDim()>().at(domain_->getLabel()));
    //     //             finite_element_->setCoordinates(finite_element_->getCurrentCentroid());
    //     //             element_map.template getElements<t_element_coordinates.getDim(), t_element_coordinates.getTag()>()[getHash(node_tags_)] = finite_element_;
    //     //         }
    //     //     };
    //     //     make_element(make_element);
    //     // }

    // private:

    //     ElementFactoryMap<Frame_> my_map_;

    //     Elements_ elements_;

    //     Regions_ regions_;

    // };
 
} // namespace lolita::mesh


#endif /* C5421678_DA6F_4EA0_AFC9_73889C52D424 */
