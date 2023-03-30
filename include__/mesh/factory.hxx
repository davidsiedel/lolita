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
#include "mesh/region_set.hxx"
#include "mesh/element.hxx"
#include "mesh/element_set.hxx"
#include "mesh/gmsh.hxx"

namespace lolita::mesh
{
    
    template<geometry::ShapeConcept Shape_>
    using ShapeNodeConnectivity = std::array<Natural, Shape_::getNumNodes()>;

    template<geometry::ShapeConcept Shape_>
    using ShapeInnerNeighborhoodNodeConnectivity = typename geometry::ShapeInnerNeighborhoodTraits<Shape_>::template InnerNeighborhood<ShapeNodeConnectivity>;

    template<typename T_>
    concept MeshFormatConcept = requires
    {

        { T_::getFileNameExtension() } -> std::same_as<std::basic_string_view<Character>>;

    };

    template<MeshFormatConcept, geometry::ShapeConcept Shape_>
    struct Table;

    template<geometry::FrameConcept Frame_>
    struct MyElementSet
    {

        template<typename... T_>
        using ElementPointer_ = std::shared_ptr<Element<T_...>>;

        template<typename... T_>
        using RegionPointer_ = std::shared_ptr<Region<T_...>>;

        using Elements_ = ElementSet<ElementPointer_, Frame_>;

        using Regions_ = RegionSet<RegionPointer_, Frame_>;

        Elements_ elements_;

        Regions_ regions_;
    
    };

    template<geometry::FrameConcept Frame_>
    struct MyElementMap
    {

    private:

        template<typename... T_>
        using ElementPointer_ = std::shared_ptr<Element<T_...>>;

        template<typename... T_>
        using RegionPointer_ = std::shared_ptr<Region<T_...>>;

        using Elements_ = ElementMap<ElementPointer_, Frame_>;

        using Regions_ = RegionMap<RegionPointer_, Frame_>;
    
        using MyElementSet_ = MyElementSet<Frame_>;

    public:

        MyElementMap()
        :
        elements_(),
        regions_()
        {}
        
        std::unique_ptr<MyElementSet_>
        makeElementSet()
        const
        {
            auto finite_element_set = std::make_unique<MyElementSet_>();
            auto make_sets = [&] <Integer i_ = 0> (
                auto & make_sets_
            )
            mutable
            {
                for (auto const & dm : regions_->template getDomains<i_>())
                {
                    finite_element_set->regions_.template getDomains<i_>().push_back(dm.second);
                }
                if constexpr (i_ < Frame_::getDimEuclidean())
                {
                    make_sets_.template operator()<i_ + 1>(make_sets_);
                }
            };
            make_sets(make_sets);
            auto make_elements = [&] <Integer i_ = 0, Integer j_ = 0> (
                auto & make_elements_
            )
            mutable
            {
                for (auto const & element : elements_->template getElements<i_, j_>())
                {
                    finite_element_set->elements_.template getElements<i_, j_>().push_back(element.second);
                }
                if constexpr (j_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<i_>() - 1)
                {
                    make_elements_.template operator()<i_, j_ + 1>(make_elements_);
                }
                else if constexpr (i_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<>() - 1)
                {
                    make_elements_.template operator()<i_ + 1, 0>(make_elements_);
                }
            }; 
            make_elements(make_elements);
            return finite_element_set;
        }

    private:

        Elements_ elements_;

        Regions_ regions_;

    };

    /**
     * @brief A helper class to create a Domain object.
     * An object of this class is intended to be built by any mesh format parser.
     * 
     * @tparam Domain_ 
     * @tparam Frame_ 
     */
    template<geometry::DomainConcept Domain_, geometry::FrameConcept Frame_>
    struct RegionFactory
    {

    private:

        /**
         * @brief Alias for FiniteDomain 
         * 
         */
        using Region_ = Region<Domain_, Frame_>;

    public:
        
        static
        std::basic_string<Character>
        makeHash(
            std::basic_string<Character> && label
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << Frame_::getDimEuclidean();
            hash << std::forward<std::basic_string<Character>>(label);
            return hash.str();
        }
        
        static
        std::basic_string<Character>
        makeHash(
            Integer tag
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << Frame_::getDimEuclidean();
            hash << std::setfill('0') << std::setw(10) << tag;
            return hash.str();
        }
        
        // /**
        //  * @brief Construct a new Mesh Domain object
        //  * 
        //  */
        // RegionFactory()
        // :
        // label_()
        // {}

        /**
         * @brief Construct a new Mesh Domain object
         * 
         * @param label 
         */
        RegionFactory(
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
        RegionFactory(
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
         * @return std::shared_ptr<Region_> 
         */
        std::shared_ptr<Region_>
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

    template<geometry::ShapeConcept Shape_, geometry::FrameConcept Frame_>
    struct ElementFactory
    {

    private:
        
        using ElementPointer_ = std::shared_ptr<Element<Shape_, Frame_>>;

        using Element_ = Element<Shape_, Frame_>;

        using NodeConnectivity_ = ShapeNodeConnectivity<Shape_>;
        
        using Region_ = std::shared_ptr<RegionFactory<geometry::Domain<Shape_::getDimShape()>, Frame_>>;
        
        using Regions_ = std::vector<std::shared_ptr<RegionFactory<geometry::Domain<Shape_::getDimShape()>, Frame_>>>;

    public:

        template<typename Mesh_, Integer i_, Integer j_>
        static constexpr
        Integer
        getInnerNeighborNodeConnection(
            Integer i,
            Integer j
        )
        requires(!geometry::PointConcept<Shape_>)
        {
            return std::get<j_>(std::get<i_>(Table<Mesh_, Shape_>::node_connectivity_))[i][j];
        }
        
        static
        std::basic_string<Character>
        makeHash(
            auto const & tag
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << tag;
            return hash.str();
        }
        
        static
        std::basic_string<Character>
        makeHash(
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
        
        // ElementFactory()
        // :
        // node_tags_(),
        // finite_element_()
        // {}
        
        explicit
        ElementFactory(
            ShapeNodeConnectivity<Shape_> const & node_tags
        )
        :
        node_tags_(node_tags),
        finite_element_()
        {}
        
        void
        setDomain(
            Region_ const & domain
        )
        {
            domain_ = domain;
        }

        void
        addRegion(
            Region_ const & domain
        )
        {
            regions_.push_back(domain);
        }
        
        // void
        // setNodeTag(
        //     Integer index,
        //     Natural const & tag
        // )
        // {
        //     node_tags_[index] = tag;
        // }
        
        // Natural const &
        // getNodeTag(
        //     Integer index
        // )
        // const
        // {
        //     return node_tags_[index];
        // }
        
        ShapeNodeConnectivity<Shape_> const &
        getNodeTags()
        const
        {
            return node_tags_;
        }
        
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

        Region_ domain_;

        Regions_ regions_;

        ShapeNodeConnectivity<Shape_> node_tags_;

        ElementPointer_ finite_element_;
    
    };

    template<geometry::PointConcept Shape_, geometry::FrameConcept Frame_>
    struct ElementFactory<Shape_, Frame_>
    {

    private:
        
        using ElementPointer_ = std::shared_ptr<Element<Shape_, Frame_>>;

        using Element_ = Element<Shape_, Frame_>;

        // using NodeConnectivity_ = ShapeNodeConnectivity<Shape_>;
        
        using Region_ = std::shared_ptr<RegionFactory<geometry::Domain<Shape_::getDimShape()>, Frame_>>;
        
        using Regions_ = std::vector<std::shared_ptr<RegionFactory<geometry::Domain<Shape_::getDimShape()>, Frame_>>>;

    public:

    // private:
    
    //     using t_ElementTraits = ShapeTraits<t_element>;
        
    //     template<DomainConcept auto t_dim, MeshConcept auto t__domain>
    //     using t_Dom = MeshDomain<t_dim, t__domain>;

    //     /**
    //      * @brief The MeshDomain the element is connected to.
    //      * 
    //      */
    //     using MeshDomain_ = typename t_ElementTraits::template DomainConnectivity<t_Dom, t_domain>;

    // public:

        static
        std::basic_string<Character>
        makeHash(
            Natural tag
        )
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << tag;
            return hash.str();
        }

        // MeshElement()
        // :
        // tag_(),
        // coordinates_(),
        // finite_element_()
        // {}

        explicit
        ElementFactory(
            geometry::Point<Frame_::getDimEuclidean()> coordinates
        )
        :
        tag_(),
        coordinates_(coordinates),
        finite_element_()
        {}

        void
        addRegion(
            Region_ const & domain
        )
        {
            regions_.push_back(domain);
        }
        
        // /**
        //  * @brief Set the domain pointer that the element is connected to.
        //  * 
        //  * @param domain The domain pointer to be set.
        //  */
        // void
        // setDomain(
        //     std::shared_ptr<MeshDomain_> const & domain
        // )
        // {
        //     domain_ = domain;
        // }

        // void
        // setTag(
        //     Natural tag
        // )
        // {
        //     tag_ = tag;
        // }

        // Natural const &
        // getNodeTag()
        // const
        // {
        //     return tag_;
        // }

        // void
        // setCoordinate(
        //     Integer i,
        //     Real const & coordinate
        // )
        // {
        //     coordinates_(i) = coordinate;
        // }
        
        std::basic_string<Character>
        getHash()
        const
        {
            auto hash = std::basic_stringstream<Character>();
            hash << std::setfill('0') << std::setw(10) << tag_;
            return hash.str();
        }

    private:

        ElementPointer_ finite_element_;

        Regions_ regions_;
        
        /**
         * @brief A tag that unambiguously identifies the node.
         * 
         */
        Natural tag_;
        
        /**
         * @brief The coordinates of the node in the mesh.
         * 
         */
        geometry::Point<Frame_::getDimEuclidean()> coordinates_;

        /**
         * @brief A pointer to the FiniteElement to be build.
         * 
         */
        // std::shared_ptr<FiniteElement<t_element, t_domain>> finite_element_;

        /**
         * @brief The MeshDomain the element is connected to.
         * 
         */
        // std::shared_ptr<MeshDomain_> domain_;
        
    };
    
    template<geometry::FrameConcept Frame_>
    struct ElementFactoryMap
    {

    private:

        template<typename... T_>
        using ElementPointer_ = std::shared_ptr<ElementFactory<T_...>>;

        template<typename... T_>
        using RegionPointer_ = std::shared_ptr<RegionFactory<T_...>>;

        using Elements_ = ElementMap<ElementPointer_, Frame_>;

        using Regions_ = RegionMap<RegionPointer_, Frame_>;

        template<geometry::ShapeConcept Shape_>
        static
        void
        setElementOuterNeighborhood(
            std::shared_ptr<Element<Shape_, Frame_>> const & finite_element
        )
        requires(std::same_as<Shape_, geometry::Node>)
        {}

        template<geometry::ShapeConcept Shape_, Integer i_ = 0>
        static
        void
        setElementOuterNeighborhood(
            std::shared_ptr<Element<Shape_, Frame_>> const & finite_element
        )
        requires(!std::same_as<Shape_, geometry::Node>)
        {
            auto constexpr shape_coordinates = geometry::ShapeLibraryTraits<Frame_>::template getShapeCoordinates<Shape_>();
            auto constexpr node_coordinates = geometry::ShapeInnerNeighborhoodTraits<Shape_>::template getShapeCoordinates<geometry::Node>();
            for (auto const & nde : finite_element->template getInnerNeighbors<node_coordinates.i_, node_coordinates.j_>())
            {
                auto const & ngs = nde->template getOuterNeighbors<shape_coordinates.i_ - 1, i_>();
                for (auto const & neighbour : ngs)
                {
                    if (((neighbour->getTag() != finite_element->getTag()) && i_ == shape_coordinates.j_) || (i_ != shape_coordinates.j_))
                    {
                        auto & element_ngs = finite_element->template getOuterNeighbors<0, i_>();
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
            if constexpr (i_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<shape_coordinates.i_>() - 1)
            {
                setElementOuterNeighborhood<Shape_, i_ + 1>(finite_element);
            }
        }

    public:

        ElementFactoryMap()
        {}
        
        std::unique_ptr<MyElementSet<Frame_>>
        makeElementSet()
        const
        {
            auto element_map = std::make_unique<MyElementMap<Frame_>>();
            auto make_sets = [&] <Integer i_ = 0> (
                auto & make_sets_
            )
            mutable
            {
                for (auto const & dm : regions_->template getDomains<i_>())
                {
                    element_map->template getDomains<i_>()[dm.second->getLabel()] = dm.second->letDomain();
                }
                if constexpr (i_ < Frame_::getDimEuclidean())
                {
                    make_sets_.template operator()<i_ + 1>(make_sets_);
                }
            };
            make_sets(make_sets);
            auto make_elements = [&] <Integer i_ = 0, Integer j_ = 0> (
                auto & make_elements_
            )
            mutable
            {
                auto constexpr t_element = geometry::ShapeLibraryTraits<Frame_>::template getElement<i_, j_>();
                for (auto const & element : elements_->template getElements<i_, j_>())
                {
                    element.second->template makeElement(* element_map);
                }
                if constexpr (j_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<i_>() - 1)
                {
                    make_elements_.template operator()<i_, j_ + 1>(make_elements_);
                }
                else if constexpr (i_ < geometry::ShapeLibraryTraits<Frame_>::getNumElements() - 1)
                {
                    make_elements_.template operator()<i_ + 1, 0>(make_elements_);
                }
            };
            make_elements(make_elements);
            auto make_elements_outer_neighborhood = [&] <Integer i_ = 0, Integer j_ = 0> (
                auto & make_elements_outer_neighborhood_
            )
            mutable
            {
                auto const constexpr t_element = geometry::ShapeLibraryTraits<Frame_>::template getElement<i_, j_>();
                for (auto & element : element_map->template getElements<i_, j_>())
                {
                    setElementOuterNeighborhood<t_element>(element.second);
                }
                if constexpr (j_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<i_>() - 1)
                {
                    make_elements_outer_neighborhood_.template operator()<i_, j_ + 1>(make_elements_outer_neighborhood_);
                }
                else if constexpr (i_ < geometry::ShapeLibraryTraits<Frame_>::getNumElements() - 1)
                {
                    make_elements_outer_neighborhood_.template operator()<i_ + 1, 0>(make_elements_outer_neighborhood_);
                }
            };
            make_elements_outer_neighborhood(make_elements_outer_neighborhood);
            return element_map->makeElementSet();
        }

        template<geometry::DomainConcept Domain_>
        void
        addRegion(
            auto const & domain_tag,
            std::basic_string<Character> const & domain_label
        )
        {
            auto & regions = regions_.template getDomains<Domain_::getDimDomain()>();
            auto region_hash = RegionFactory<Domain_, Frame_>::makeHash(domain_tag);
            regions[region_hash] = std::make_shared<RegionFactory<Domain_, Frame_>>(domain_label);
        }

        template<geometry::DomainConcept Domain_>
        std::shared_ptr<RegionFactory<Domain_, Frame_>> const &
        getRegion(
            auto const & domain_tag
        )
        const
        {
            auto region_hash = RegionFactory<Domain_, Frame_>::makeHash(domain_tag);
            return regions_.template getDomains<Domain_::getDimDomain()>().at(region_hash);
        }

        template<geometry::ShapeConcept Shape_>
        void
        addElement(
            auto const & element_tag,
            ShapeNodeConnectivity<Shape_> const & node_tags
        )
        {
            auto constexpr shape_coordinates = geometry::ShapeLibraryTraits<Frame_>::template getShapeCoordinates<Shape_>();
            auto & elements = elements_.template getElements<shape_coordinates.i_, shape_coordinates.j_>();
            auto element_hash = ElementFactory<Shape_, Frame_>::makeHash(element_tag);
            // auto const & region = this->template getRegion<geometry::Domain<Shape_::getDimShape()>>(domain_tag);
            elements[element_hash] = std::make_shared<ElementFactory<Shape_, Frame_>>(node_tags);
        }

        template<geometry::ShapeConcept Shape_>
        void
        addElement(
            auto const & element_tag,
            geometry::Point<Frame_::getDimEuclidean()> const & node_coordinates
        )
        requires(std::same_as<Shape_, geometry::Node>)
        {
            auto constexpr shape_coordinates = geometry::ShapeLibraryTraits<Frame_>::template getShapeCoordinates<Shape_>();
            auto & elements = elements_.template getElements<shape_coordinates.i_, shape_coordinates.j_>();
            auto element_hash = ElementFactory<Shape_, Frame_>::makeHash(element_tag);
            // auto const & region = this->template getRegion<geometry::Domain<Shape_::getDimShape()>>(domain_tag);
            elements[element_hash] = std::make_shared<ElementFactory<Shape_, Frame_>>(node_coordinates);
        }

        template<geometry::ShapeConcept Shape_>
        void
        addElementInDomain(
            auto const & element_tag,
            auto const & domain_tag
        )
        {
            auto constexpr shape_coordinates = geometry::ShapeLibraryTraits<Frame_>::template getShapeCoordinates<Shape_>();
            auto & elements = elements_.template getElements<shape_coordinates.i_, shape_coordinates.j_>();
            auto element_hash = ElementFactory<Shape_, Frame_>::makeHash(element_tag);
            auto const & region = this->template getRegion<geometry::Domain<Shape_::getDimShape()>>(domain_tag);
            elements.at(element_hash)->addRegion(region);
        }

    private:

        Elements_ elements_;

        Regions_ regions_;

    };

} // namespace lolita::mesh


#endif /* BA740F29_E221_4994_BF1D_7892518F367D */
