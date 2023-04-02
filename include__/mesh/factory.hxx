/**
 * @file factory.hxx
 * @author David Siedel (davidsiedel@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2023-03-29
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef BA740F29_E221_4994_BF1D_7892518F367D
#define BA740F29_E221_4994_BF1D_7892518F367D

#include "geometry/frame.hxx"
#include "geometry/domain.hxx"
#include "geometry/point.hxx"
#include "geometry/shape.hxx"
#include "mesh/region.hxx"
#include "mesh/element.hxx"
#include "mesh/gmsh.hxx"

namespace lolita::mesh
{
    
    template<geometry::ShapeConcept Shape_>
    using ShapeNodeConnectivity = std::array<Natural, Shape_::getNumNodes()>;

    template<geometry::ShapeConcept Shape_>
    using ShapeInnerNeighborhoodNodeConnectivity = typename geometry::ShapeInnerNeighborhood<Shape_, ShapeNodeConnectivity>::Components;

    template<typename T_>
    concept MeshFormatConcept = requires
    {

        typename T_::NodeTag;

        { T_::getFileNameExtension() } -> std::same_as<std::basic_string_view<Character>>;

    };

    template<geometry::ShapeConcept Shape_, MeshFormatConcept>
    struct Table;

    template<geometry::FrameConcept Frame_, MeshFormatConcept MeshFormat_>
    struct MshOne;
    
    template<geometry::FrameConcept Frame_, typename Mesh_>
    struct ElementFactoryMap
    {

    private:

        template<typename... T_>
        using ElementPointer2_ = std::unordered_map<std::basic_string<Character>, std::shared_ptr<Element<T_...>>>;

        template<typename... T_>
        using RegionPointer2_ = std::unordered_map<std::basic_string<Character>, std::shared_ptr<Region<T_...>>>;

        using Elements2_ = geometry::ShapeCollection<Frame_, ElementPointer2_, Frame_>;

        using Regions2_ = geometry::DomainCollection<Frame_, RegionPointer2_, Frame_>;

        // using NodeTag = typename Mesh_::NodeTag;
        
        template<geometry::DomainConcept Domain_>
        static
        std::basic_string<Character>
        makeRegionHash(
            std::basic_string<Character> && label
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << Domain_::getDimDomain();
            hash << std::forward<std::basic_string<Character>>(label);
            return hash.str();
        }
        
        template<geometry::DomainConcept Domain_>
        static
        std::basic_string<Character>
        makeRegionHash(
            Integer tag
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << Domain_::getDimDomain();
            hash << std::setfill('0') << std::setw(10) << tag;
            return hash.str();
        }
        
        template<geometry::PointConcept Shape_>
        static
        std::basic_string<Character>
        makeElementHash(
            auto const & tag
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << tag;
            return hash.str();
        }
        
        template<geometry::ShapeConcept Shape_>
        static
        std::basic_string<Character>
        makeElementHash(
            ShapeNodeConnectivity<Shape_> node_tags
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

        template<geometry::ShapeConcept Shape_, geometry::ShapeConcept Neighbor_>
        static constexpr
        Integer
        getInnerNeighborNode(
            Integer i,
            Integer j
        )
        {
            auto constexpr cc = geometry::ShapeInnerNeighborhood<Shape_>::template ShapeTraits<Neighbor_>::coordinates_;
            return std::get<cc[1]>(std::get<cc[0]>(Table<Shape_, Mesh_>::node_connectivity_))[i][j];
        }

    public:

        ElementFactoryMap()
        {}

        template<geometry::DomainConcept Domain_>
        void
        addRegion(
            auto const & domain_tag,
            std::basic_string<Character> const & domain_label
        )
        {
            auto & regions = regions2_.template getComponent<Domain_>();
            auto region_hash = makeRegionHash<Domain_>(domain_tag);
            regions[region_hash] = std::make_shared<Region<Domain_, Frame_>>(domain_label);
            std::cout << "adding region : " << domain_label << std::endl;
        }

        template<geometry::DomainConcept Domain_>
        std::shared_ptr<Region<Domain_, Frame_>> const &
        getRegion(
            auto const & domain_tag
        )
        const
        {
            auto region_hash = makeRegionHash<Domain_>(domain_tag);
            return regions2_.template getComponent<Domain_>().at(region_hash);
        }

        template<geometry::ShapeConcept Shape_>
        void
        addElement(
            ShapeNodeConnectivity<Shape_> const & node_tags
        )
        {
            auto & elements = elements2_.template getComponent<Shape_>();
            auto element = std::make_shared<Element<Shape_, Frame_>>(elements.size());
            elements[makeElementHash<Shape_>(node_tags)] = element;
            std::cout << "making element : " << element->getTag() << " : " << makeElementHash<Shape_>(node_tags) << std::endl;
            auto make_neighborhood = [&] <geometry::ShapeConcept Neighbor_> ()
            {
                auto & inner_neighbors = elements2_.template getComponent<Neighbor_>();
                for (auto i = 0; i < geometry::ShapeInnerNeighborhood<Shape_>::template getNumComponents<Neighbor_>(); ++i)
                {
                    auto inner_neighbor_hash = std::basic_string<Character>();
                    if constexpr (!std::same_as<Neighbor_, geometry::Node>)
                    {
                        auto inner_neighbor_node_tags = ShapeNodeConnectivity<Neighbor_>();
                        for (auto j = 0; j < Neighbor_::getNumNodes(); ++j)
                        {
                            inner_neighbor_node_tags[j] = node_tags[getInnerNeighborNode<Shape_, Neighbor_>(i, j)];
                        }
                        inner_neighbor_hash = makeElementHash<Neighbor_>(inner_neighbor_node_tags);
                        if (!inner_neighbors.contains(inner_neighbor_hash))
                        {
                            addElement<Neighbor_>(inner_neighbor_node_tags);
                        }
                    }
                    else
                    {
                        auto node_tag = node_tags[getInnerNeighborNode<Shape_, Neighbor_>(i, 0)];
                        inner_neighbor_hash = makeElementHash<Neighbor_>(node_tag);
                    }
                    auto & inner_neighbor = inner_neighbors.at(inner_neighbor_hash);
                    element->getInnerNeighborhood().template getComponent<Neighbor_>()[i] = inner_neighbor;
                    if constexpr (std::same_as<Neighbor_, geometry::Node>)
                    {
                        inner_neighbor->template getOuterNeighborhood().template getComponent<Shape_>().push_back(element);
                    }
                }
            };
            geometry::ShapeInnerNeighborhood<Shape_>::make(make_neighborhood);
        }

        template<geometry::ShapeConcept Shape_>
        void
        addElement(
            auto const & element_tag,
            geometry::Point<Frame_::getDimEuclidean()> const & node_coordinates
        )
        requires(std::same_as<Shape_, geometry::Node>)
        {
            auto & elements = elements2_.template getComponent<Shape_>();
            auto element = std::make_shared<Element<Shape_, Frame_>>(element_tag);
            elements[makeElementHash<Shape_>(element_tag)] = element;
            element->setCoordinates(node_coordinates);
            std::cout << "making node : " << element->getTag() << " : " << makeElementHash<Shape_>(element_tag) << std::endl;
        }

        template<geometry::ShapeConcept Shape_>
        void
        addRegionToElement(
            auto const & region_tag,
            ShapeNodeConnectivity<Shape_> const & node_tags
        )
        {
            auto & elements = elements2_.template getComponent<Shape_>();
            auto element_hash = makeElementHash<Shape_>(node_tags);
            auto const & region = this->template getRegion<typename Shape_::Domain>(region_tag);
            elements.at(element_hash)->addRegion(region);
            std::cout << "adding region : " << region->getLabel() << " to element : " << makeElementHash<Shape_>(node_tags) << std::endl;
        }

        template<geometry::ShapeConcept Shape_>
        void
        addRegionToElement(
            auto const & region_tag,
            auto const & node_tag
        )
        requires(std::same_as<Shape_, geometry::Node>)
        {
            auto & elements = elements2_.template getComponent<Shape_>();
            auto element_hash = makeElementHash<Shape_>(node_tag);
            auto const & region = this->template getRegion<typename Shape_::Domain>(region_tag);
            elements.at(element_hash)->addRegion(region);
            std::cout << "adding region : " << region->getLabel() << " to element : " << makeElementHash<Shape_>(node_tag) << std::endl;
        }
        
        void
        setIt(
            std::basic_string<Character> && file_path
        )
        {
            auto my_one = MshOne<Frame_, Mesh_>(std::forward<std::basic_string<Character>>(file_path));
            auto set_regions = [&] <geometry::DomainConcept Domain_> ()
            {
                my_one.template setMeshDomain<Domain_>(* this);
            };
            geometry::DomainCollection<Frame_>::make(set_regions);
            auto set_elements = [&] <geometry::ShapeConcept Shape_> ()
            {
                my_one.template setMeshElement<Shape_>(* this);
            };
            geometry::ShapeCollection<Frame_>::make(set_elements);
        }

    private:

        Elements2_ elements2_;

        Regions2_ regions2_;

    };

} // namespace lolita::mesh


#endif /* BA740F29_E221_4994_BF1D_7892518F367D */
