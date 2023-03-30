#ifndef C5421678_DA6F_4EA0_AFC9_73889C52D424
#define C5421678_DA6F_4EA0_AFC9_73889C52D424

#include "geometry/domain.hxx"
#include "geometry/point.hxx"
#include "geometry/shape.hxx"
#include "geometry/frame.hxx"
#include "mesh/element.hxx"
#include "mesh/element_set.hxx"
#include "mesh/region.hxx"
#include "mesh/region_set.hxx"
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
    struct Table<Msh, geometry::Node>
    {

        static constexpr
        Integer tag_ = 15;

    };

    template<>
    struct Table<Msh, geometry::Segment>
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
    struct Table<Msh, geometry::Triangle>
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
    struct Table<Msh, geometry::Quadrangle>
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
    struct MshOne
    {

        MshOne(
            std::basic_string<Character> && file_path
        )
        :
        file_(std::forward<std::basic_string<Character>>(file_path))
        {
            ElementFactoryMap<Frame_> my_map_;
            setGeometricalEntities();
            setMeshDomain(my_map_);
            setMeshElement<geometry::Node>(my_map_);
            setMeshElement<geometry::Triangle>(my_map_);
        }

        template<geometry::DomainConcept Domain_>
        void
        addGeometricEntity(
            Integer tag,
            GeometricEntity<Domain_, Frame_> const & geometric_entity
        )
        {
            geometric_entities_.template getDomains<Domain_::getDimDomain()>()[GeometricEntity<Domain_, Frame_>::makeHash(tag)] = geometric_entity;
        }

        template<geometry::DomainConcept Domain_>
        GeometricEntity<Domain_, Frame_> const &
        getGeometricEntity(
            Integer tag
        )
        const
        {
            return geometric_entities_.template getDomains<Domain_::getDimDomain()>().at(GeometricEntity<Domain_, Frame_>::makeHash(tag));
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
        
        void
        setMeshDomain(
            ElementFactoryMap<Frame_> & my_map_
        )
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
                auto set_domain = [&] <Integer i_ = 0> (
                    auto & set_domain_
                )
                constexpr mutable
                {
                    if (i_ == dim)
                    {
                        using Domain_ = geometry::Domain<i_>;
                        my_map_.template addRegion<Domain_>(physical_entity_tag, name);
                    }
                    if constexpr (i_ < Frame_::getDimEuclidean())
                    {
                        set_domain_.template operator ()<i_ + 1>(set_domain_);
                    }
                };
                set_domain(set_domain);
                offset += 1;
            }
        }

        template<geometry::ShapeConcept Shape_>
        void
        setMeshElement(
            ElementFactoryMap<Frame_> & my_map_
        )
        requires(std::same_as<Shape_, geometry::Node>)
        {
            auto constexpr t_element_coordinates = geometry::ShapeLibraryTraits<Frame_>::template getShapeCoordinates<Shape_>();
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
                    my_map_.template addElement<Shape_>(node_tag, node_coordinates);
                    if (entity_dim == Shape_::getDimShape())
                    {
                        using Domain_ = geometry::Domain<Shape_::getDimShape()>;
                        auto const & geometric_entity = getGeometricEntity<Domain_>(entity_tag);
                        for (auto const & physical_entity_tag : geometric_entity.getPhysicalEntityTags())
                        {
                            my_map_.template addElementInDomain<Shape_>(node_tag, physical_entity_tag);
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
        requires(!std::same_as<Shape_, geometry::Node>)
        {
            auto constexpr t_element_coordinates = geometry::ShapeLibraryTraits<Frame_>::template getShapeCoordinates<Shape_>();
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
                if (element_type_tag == Table<Msh, Shape_>::tag_)
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
                            my_map_.template addElementInDomain<Shape_>(element_tag, physical_entity_tag);
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

        RegionMap<GeometricEntity, Frame_> geometric_entities_;

        File file_;

    };
 
} // namespace lolita::mesh


#endif /* C5421678_DA6F_4EA0_AFC9_73889C52D424 */
