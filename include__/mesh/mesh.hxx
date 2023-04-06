/**
 * @file mesh.hxx
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

namespace lolita::mesh
{
    
    template<geometry::ShapeConcept Shape_>
    using ShapeNodeConnectivity = std::array<Natural, Shape_::getNumNodes()>;

    template<geometry::ShapeConcept Shape_>
    using ShapeInnerNeighborhoodNodeConnectivity = typename geometry::ShapeInnerNeighborhoodTraits<Shape_>::template InnerNeighborhood<ShapeNodeConnectivity>;

    template<typename T_>
    concept MeshFormatConcept = requires
    {

        typename T_::NodeTag;

        { T_::getFileNameExtension() } -> std::same_as<std::basic_string_view<Character>>;

    };

    struct Msh
    {

        using NodeTag = Integer;

        static constexpr
        std::basic_string_view<Character>
        getFileNameExtension()
        {
            return "msh";
        }

    };

    using MeshFormats = std::tuple<
        Msh
    >;

    struct MeshFormatLibrary
    {
        
        template<Integer i_>
        using MeshFormat = std::tuple_element_t<i_, MeshFormats>;

        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<MeshFormats>;
        }

        static
        void
        apply(
            auto const & fun
        )
        {
            auto apply = [&] <Integer i_ = 0> (
                auto & apply_
            )
            mutable
            {
                fun.template operator()<MeshFormat<i_>>();
                if constexpr (i_ < getNumComponents() - 1)
                {
                    apply_.template operator()<i_ + 1>(apply_);
                }
            };
            apply(apply);
        }

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

    template<geometry::ShapeConcept Shape_, MeshFormatConcept>
    struct Table;

    template<geometry::FrameConcept Frame_, MeshFormatConcept MeshFormat_>
    struct MeshFileParser;
    
    template<geometry::FrameConcept Frame_, MeshFormatConcept MeshFormat_>
    struct MeshFactory;

    struct MeshFile
    {

    public:

        explicit
        MeshFile(
            std::basic_string<Character> && file_path
        )
        :
        file_path_(std::move(file_path)),
        lines_(readLines())
        {}

        explicit
        MeshFile(
            std::basic_string<Character> const & file_path
        )
        :
        file_path_(file_path),
        lines_(readLines())
        {}

        inline
        std::vector<std::basic_string<Character>> const &
        getContent()
        const
        {
            return lines_;
        }

    private:

        inline
        std::vector<std::basic_string<Character>>
        readLines()
        {
            auto lines = std::vector<std::basic_string<Character>>();
            auto file = std::basic_ifstream<Character>(std::forward<std::basic_string<Character>>(file_path_));
            if (file)
            {
                for (std::basic_string<Character> line; std::getline(file, line); )
                {
                    lines.push_back(line);
                }
                return lines;
            }
            throw std::runtime_error("Could not open file");
        }

        std::basic_string<Character> file_path_;

        std::vector<std::basic_string<Character>> const lines_;

    };

    template<geometry::FrameConcept Frame_, typename... T_>
    struct RegionCollection
    {

        template<typename... U_>
        using Regions = std::vector<std::shared_ptr<Region<U_...>>>;

        using Collection = typename geometry::DomainCollectionTraits<Frame_>::template Collection<Regions, Frame_>;

        RegionCollection()
        {}
    
        template<geometry::DomainConcept Domain_>
        Regions<Domain_, Frame_> const &
        getRegions()
        const
        {
            auto constexpr i_ = geometry::DomainCollectionTraits<Frame_>::template getDomainIndex<Domain_>();
            return std::get<i_>(collection_);
        }
    
        template<geometry::DomainConcept Domain_>
        Regions<Domain_, Frame_> &
        getRegions()
        {
            auto constexpr i_ = geometry::DomainCollectionTraits<Frame_>::template getDomainIndex<Domain_>();
            return std::get<i_>(collection_);
        }

    private:

        Collection collection_;

    };

    template<geometry::FrameConcept Frame_, typename... T_>
    struct ElementCollection
    {

        template<typename... U_>
        using Elements = std::vector<std::shared_ptr<Element<U_...>>>;

        using Collection = typename geometry::ShapeCollectionTraits<Frame_>::template Collection<Elements, Frame_>;

        ElementCollection()
        {}
    
        template<geometry::ShapeConcept Shape_>
        Elements<Shape_, Frame_> const &
        getElements()
        const
        {
            auto constexpr i_ = geometry::ShapeCollectionTraits<Frame_>::template getShapeIndex<Shape_>(0);
            auto constexpr j_ = geometry::ShapeCollectionTraits<Frame_>::template getShapeIndex<Shape_>(1);
            return std::get<j_>(std::get<i_>(collection_));
        }
    
        template<geometry::ShapeConcept Shape_>
        Elements<Shape_, Frame_> &
        getElements()
        {
            auto constexpr i_ = geometry::ShapeCollectionTraits<Frame_>::template getShapeIndex<Shape_>(0);
            auto constexpr j_ = geometry::ShapeCollectionTraits<Frame_>::template getShapeIndex<Shape_>(1);
            return std::get<j_>(std::get<i_>(collection_));
        }

    private:

        Collection collection_;

    };
    
    template<geometry::FrameConcept Frame_>
    struct Mesh
    {

        using RegionCollection_ = RegionCollection<Frame_>;

        using ElementCollection_ = ElementCollection<Frame_>;

        static
        std::unique_ptr<Mesh>
        makeFromMeshFile(
            std::basic_string<Character> && file_path
        )
        {
            auto mesh = std::unique_ptr<Mesh>();
            auto set_mesh = [&] <mesh::ReaderConcept MeshReader_> ()
            {
                auto mesh_format_tag = "." + std::basic_string<Character>(MeshReader_::getFileNameExtension());
                if (std::forward<std::basic_string<Character>>(file_path).find(mesh_format_tag) != std::string::npos)
                {
                    mesh = std::make_unique<Mesh>(MeshFactory<Frame_, MeshReader_>(std::forward<std::basic_string<Character>>(file_path)).makeMesh());
                }
            };
            MeshFormatLibrary::apply(set_mesh);
            if (mesh == nullptr)
            {
                throw std::runtime_error("Mesh format not implemented");
            }
            return mesh;
        }

        Mesh(
            RegionCollection_ && region_collection,
            ElementCollection_ && element_collection
        )
        :
        region_collection_(std::move(region_collection)),
        element_collection_(std::move(element_collection))
        {}

        void
        working()
        {
            std::cout << "first node tag : " << element_collection_.template getElements<geometry::Node>()[0]->getTag() << std::endl;
        }

    private:

        RegionCollection_ region_collection_;

        ElementCollection_ element_collection_;

    };
    
    template<geometry::FrameConcept Frame_, MeshFormatConcept MeshFormat_>
    struct MeshFactory
    {

    private:

        template<typename... T_>
        using ElementTags_ = std::unordered_map<std::basic_string<Character>, Natural>;

        template<typename... T_>
        using RegionTags_ = std::unordered_map<std::basic_string<Character>, Natural>;

        using ElementTagsCollection_ = typename geometry::ShapeCollectionTraits<Frame_>::template Collection<ElementTags_, Frame_>;

        using RegionTagsCollection_ = typename geometry::DomainCollectionTraits<Frame_>::template Collection<RegionTags_, Frame_>;

        using ElementCollection_ = ElementCollection<Frame_>;

        using RegionCollection_ = RegionCollection<Frame_>;

        using MeshFileParser_ = MeshFileParser<Frame_, MeshFormat_>;
        
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
            Natural const & tag
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
            auto constexpr i_ = geometry::ShapeInnerNeighborhoodTraits<Shape_>::template getShapeIndex<Neighbor_>(0);
            auto constexpr j_ = geometry::ShapeInnerNeighborhoodTraits<Shape_>::template getShapeIndex<Neighbor_>(1);
            return std::get<j_>(std::get<i_>(Table<Shape_, MeshFormat_>::node_connectivity_))[i][j];
        }
    
        template<geometry::DomainConcept Domain_>
        RegionTags_<Domain_, Frame_> const &
        getRegionTags()
        const
        {
            auto constexpr i_ = geometry::DomainCollectionTraits<Frame_>::template getDomainIndex<Domain_>();
            return std::get<i_>(region_tags_collection_);
        }
    
        template<geometry::DomainConcept Domain_>
        RegionTags_<Domain_, Frame_> &
        getRegionTags()
        {
            auto constexpr i_ = geometry::DomainCollectionTraits<Frame_>::template getDomainIndex<Domain_>();
            return std::get<i_>(region_tags_collection_);
        }
    
        template<geometry::ShapeConcept Shape_>
        ElementTags_<Shape_, Frame_> const &
        getElementTags()
        const
        {
            auto constexpr i_ = geometry::ShapeCollectionTraits<Frame_>::template getShapeIndex<Shape_>(0);
            auto constexpr j_ = geometry::ShapeCollectionTraits<Frame_>::template getShapeIndex<Shape_>(1);
            return std::get<j_>(std::get<i_>(element_tags_collection_));
        }
    
        template<geometry::ShapeConcept Shape_>
        ElementTags_<Shape_, Frame_> &
        getElementTags()
        {
            auto constexpr i_ = geometry::ShapeCollectionTraits<Frame_>::template getShapeIndex<Shape_>(0);
            auto constexpr j_ = geometry::ShapeCollectionTraits<Frame_>::template getShapeIndex<Shape_>(1);
            return std::get<j_>(std::get<i_>(element_tags_collection_));
        }

    public:

        explicit
        MeshFactory(
            std::basic_string<Character> && file_path
        )
        requires(ReaderConcept<MeshFileParser_>)
        :
        mesh_file_parser_(MeshFileParser_(std::forward<std::basic_string<Character>>(file_path))),
        region_tags_collection_(),
        element_tags_collection_(),
        region_collection_(),
        element_collection_()
        {}

        template<geometry::DomainConcept Domain_>
        void
        addRegion(
            auto const & domain_tag,
            std::basic_string<Character> const & domain_label
        )
        {
            getRegionTags<Domain_>()[makeRegionHash<Domain_>(domain_tag)] = region_collection_.template getRegions<Domain_>().size();
            region_collection_.template getRegions<Domain_>().push_back(std::make_shared<Region<Domain_, Frame_>>(domain_label));
            std::cout << "adding region : " << domain_label << std::endl;
        }

        template<geometry::ShapeConcept Shape_>
        void
        addElement(
            ShapeNodeConnectivity<Shape_> const & node_tags
        )
        {
            getElementTags<Shape_>()[makeElementHash<Shape_>(node_tags)] = element_collection_.template getElements<Shape_>().size();
            auto element = std::make_shared<Element<Shape_, Frame_>>(element_collection_.template getElements<Shape_>().size());
            element_collection_.template getElements<Shape_>().push_back(element);
            std::cout << "making element : " << element->getTag() << " : " << makeElementHash<Shape_>(node_tags) << std::endl;
            auto add_neighbors = [&] <geometry::ShapeConcept Neighbor_> ()
            {
                for (auto i = 0; i < geometry::ShapeInnerNeighborhoodTraits<Shape_>::template getNumComponents<Neighbor_>(); ++i)
                {
                    auto inner_neighbor_hash = std::basic_string<Character>();
                    if constexpr (geometry::PointConcept<Neighbor_>)
                    {
                        auto node_tag = node_tags[getInnerNeighborNode<Shape_, Neighbor_>(i, 0)];
                        inner_neighbor_hash = makeElementHash<Neighbor_>(node_tag);
                    }
                    else
                    {
                        auto inner_neighbor_node_tags = ShapeNodeConnectivity<Neighbor_>();
                        for (auto j = 0; j < Neighbor_::getNumNodes(); ++j)
                        {
                            inner_neighbor_node_tags[j] = node_tags[getInnerNeighborNode<Shape_, Neighbor_>(i, j)];
                        }
                        inner_neighbor_hash = makeElementHash<Neighbor_>(inner_neighbor_node_tags);
                        if (!getElementTags<Neighbor_>().contains(inner_neighbor_hash))
                        {
                            addElement<Neighbor_>(inner_neighbor_node_tags);
                        }
                    }
                    auto & inner_neighbor = element_collection_.template getElements<Neighbor_>()[getElementTags<Neighbor_>().at(inner_neighbor_hash)];
                    element->getConnectivity().template getInnerNeighbors<Neighbor_>()[i] = inner_neighbor;
                    inner_neighbor->getConnectivity().template getOuterNeighbors<Shape_>().push_back(element);
                }
            };
            geometry::ShapeInnerNeighborhoodTraits<Shape_>::apply(add_neighbors);
        }

        template<geometry::ShapeConcept Shape_>
        void
        addElement(
            auto const & element_tag,
            geometry::Point<Frame_::getDimEuclidean()> const & node_coordinates
        )
        requires(geometry::PointConcept<Shape_>)
        {
            getElementTags<Shape_>()[makeElementHash<Shape_>(element_tag)] = element_collection_.template getElements<Shape_>().size();
            auto element = std::make_shared<Element<Shape_, Frame_>>(element_collection_.template getElements<Shape_>().size());
            element_collection_.template getElements<Shape_>().push_back(element);
            element->getConnectivity().setCoordinates(node_coordinates);
            std::cout << "making node : " << element->getTag() << " : " << makeElementHash<Shape_>(element_tag) << std::endl;
        }

        template<geometry::ShapeConcept Shape_>
        void
        addRegionToElement(
            auto const & region_tag,
            ShapeNodeConnectivity<Shape_> const & node_tags
        )
        {
            auto element_index = getElementTags<Shape_>().at(makeElementHash<Shape_>(node_tags));
            auto region_index = getRegionTags<typename Shape_::Domain>().at(makeRegionHash<typename Shape_::Domain>(region_tag));
            auto const & element = element_collection_.template getElements<Shape_>()[element_index];
            auto const & region = region_collection_.template getRegions<typename Shape_::Domain>()[region_index];
            element->addRegion(region);
            std::cout << "adding region : " << region->getLabel() << " to element : " << makeElementHash<Shape_>(node_tags) << std::endl;
        }

        template<geometry::ShapeConcept Shape_>
        void
        addRegionToElement(
            auto const & region_tag,
            auto const & node_tag
        )
        requires(geometry::PointConcept<Shape_>)
        {
            auto element_index = getElementTags<Shape_>().at(makeElementHash<Shape_>(node_tag));
            auto region_index = getRegionTags<typename Shape_::Domain>().at(makeRegionHash<typename Shape_::Domain>(region_tag));
            auto const & element = element_collection_.template getElements<Shape_>()[element_index];
            auto const & region = region_collection_.template getRegions<typename Shape_::Domain>()[region_index];
            element->addRegion(region);
            std::cout << "adding region : " << region->getLabel() << " to element : " << makeElementHash<Shape_>(node_tag) << std::endl;
        }
        
        Mesh<Frame_>
        makeMesh()
        requires(ReaderConcept<MeshFileParser_>)
        {
            auto set_regions = [&] <geometry::DomainConcept Domain_> ()
            {
                mesh_file_parser_.template setMeshDomain<Domain_>(* this);
            };
            geometry::DomainCollectionTraits<Frame_>::apply(set_regions);
            auto set_elements = [&] <geometry::ShapeConcept Shape_> ()
            {
                mesh_file_parser_.template setMeshElement<Shape_>(* this);
            };
            geometry::ShapeCollectionTraits<Frame_>::apply(set_elements);
            return Mesh<Frame_>(std::move(region_collection_), std::move(element_collection_));
        }

    private:

        MeshFileParser_ mesh_file_parser_;

        RegionCollection_ region_collection_;

        ElementCollection_ element_collection_;

        RegionTagsCollection_ region_tags_collection_;

        ElementTagsCollection_ element_tags_collection_;

    };

} // namespace lolita::mesh


#endif /* BA740F29_E221_4994_BF1D_7892518F367D */