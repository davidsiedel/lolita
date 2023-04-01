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

        { T_::getFileNameExtension() } -> std::same_as<std::basic_string_view<Character>>;

    };

    template<geometry::ShapeConcept Shape_, MeshFormatConcept>
    struct Table;

    template<geometry::FrameConcept Frame_, MeshFormatConcept MeshFormat_>
    struct MshOne;

    template<geometry::FrameConcept Frame_>
    struct ElementFactoryMap;

    template<geometry::DomainConcept Domain_, geometry::FrameConcept Frame_>
    struct RegionFactory
    {

    private:
    
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
        
        explicit
        RegionFactory(
            ElementFactoryMap<Frame_> & element_map,
            std::basic_string<Character> const & label
        )
        :
        label_(label),
        element_map_(element_map),
        region_(std::make_shared<Region_>(label_))
        {}
        
        std::basic_string<Character> const &
        getLabel()
        const
        {
            return label_;
        }
        
        void
        makeDomain()
        const
        {
            element_map_.regions2_.template getComponent<Domain_>()[label_] = region_;
        }

    // private:
        
        std::basic_string<Character> label_;

        std::shared_ptr<Region_> region_;

        ElementFactoryMap<Frame_> & element_map_;
    
    };

    template<geometry::ShapeConcept Shape_, geometry::FrameConcept Frame_>
    struct ElementFactory
    {

    private:

        using Element_ = Element<Shape_, Frame_>;

        using NodeConnectivity_ = ShapeNodeConnectivity<Shape_>;

        using Region_ = RegionFactory<typename Shape_::Domain, Frame_>;

    public:

        template<Integer i_, Integer j_, typename Mesh_>
        static constexpr
        Integer
        getInnerNeighborNodeConnection(
            Integer i,
            Integer j
        )
        {
            return std::get<j_>(std::get<i_>(Table<Shape_, Mesh_>::node_connectivity_))[i][j];
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
        
        ElementFactory(
            ElementFactoryMap<Frame_> & element_map,
            ShapeNodeConnectivity<Shape_> const & node_tags
        )
        :
        node_tags_(node_tags),
        regions_(),
        finite_element_(std::make_shared<Element_>()),
        element_map_(element_map)
        {}

        void
        addRegion(
            std::shared_ptr<Region_> const & region
        )
        {
            regions_.push_back(region);
        }

        template<typename Mesh_>
        void
        makeElement()
        {
            auto make_element = [&] <Integer i_ = 0, Integer j_ = 0> (
                auto & make_element_
            )
            mutable
            {
                using t_inner_neighbor = typename geometry::ShapeInnerNeighborhood<Shape_>::template Component<i_, j_>;
                auto & inner_neighbors = element_map_.elements2_.template getComponent<t_inner_neighbor>();
                for (auto i = 0; i < geometry::ShapeInnerNeighborhood<Shape_>::getNumComponents(i_, j_); ++i)
                {
                    auto inner_neighbor_hash = std::basic_string<Character>();
                    if constexpr(!std::same_as<t_inner_neighbor, geometry::Node>)
                    {
                        auto inner_neighbor_node_tags = ShapeNodeConnectivity<t_inner_neighbor>();
                        for (auto j = 0; j < t_inner_neighbor::getNumNodes(); ++j)
                        {
                            inner_neighbor_node_tags[j] = node_tags_[getInnerNeighborNodeConnection<i_, j_, Mesh_>(i, j)];
                        }
                        inner_neighbor_hash = ElementFactory<t_inner_neighbor, Frame_>::makeHash(inner_neighbor_node_tags);
                        if (!inner_neighbors.contains(inner_neighbor_hash))
                        {
                            ElementFactory<t_inner_neighbor, Frame_>(element_map_, inner_neighbor_node_tags).template makeElement<Mesh_>();
                        }
                    }
                    else
                    {
                        auto node_tag = node_tags_[getInnerNeighborNodeConnection<i_, j_, Mesh_>(i, 0)];
                        inner_neighbor_hash = ElementFactory<t_inner_neighbor, Frame_>::makeHash(node_tag);
                    }
                    auto & inner_neighbor = inner_neighbors.at(inner_neighbor_hash);
                    finite_element_->getInnerNeighborhood().template getComponent<t_inner_neighbor>()[i] = inner_neighbor;
                    inner_neighbor->template getOuterNeighborhood().template getComponent<Shape_>().push_back(finite_element_);
                }
                if constexpr (j_ < geometry::ShapeInnerNeighborhood<Shape_>::getNumComponents(i_) - 1)
                {
                    make_element_.template operator()<i_, j_ + 1>(make_element_);
                }
                else if constexpr (i_ < geometry::ShapeInnerNeighborhood<Shape_>::getNumComponents() - 1)
                {
                    make_element_.template operator()<i_ + 1, 0>(make_element_);
                }
                if constexpr (i_ == 0 && j_ == 0)
                {
                    for (auto const & region : regions_)
                    {
                        finite_element_->addRegion(element_map_.regions2_.template getComponent<typename Shape_::Domain>().at(region->getLabel()));
                    }
                    finite_element_->setTag(element_map_.elements2_.template getComponent<Shape_>().size());
                    element_map_.elements2_.template getComponent<Shape_>()[makeHash(node_tags_)] = finite_element_;
                }
            };
            make_element(make_element);
        }

        ShapeNodeConnectivity<Shape_> node_tags_;

        std::vector<std::shared_ptr<Region_>> regions_;

        std::shared_ptr<Element_> finite_element_;
        
        ElementFactoryMap<Frame_> & element_map_;
    
    };

    template<geometry::PointConcept Shape_, geometry::FrameConcept Frame_>
    struct ElementFactory<Shape_, Frame_>
    {

    private:

        using Element_ = Element<Shape_, Frame_>;

        using Region_ = RegionFactory<typename Shape_::Domain, Frame_>;

    public:

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
        
        ElementFactory(
            ElementFactoryMap<Frame_> & element_map,
            Natural const & tag,
            geometry::Point<Frame_::getDimEuclidean()> const & coordinates
        )
        :
        tag_(tag),
        coordinates_(coordinates),
        regions_(),
        finite_element_(std::make_shared<Element_>()),
        element_map_(element_map)
        {}

        void
        addRegion(
            std::shared_ptr<Region_> const & domain
        )
        {
            regions_.push_back(domain);
        }
        
        template<typename Mesh_>
        void
        makeElement()
        {
            finite_element_->setCoordinates(coordinates_);
            for (auto const & region : regions_)
            {
                finite_element_->addRegion(element_map_.regions2_.template getComponent<typename Shape_::Domain>().at(region->getLabel()));
            }
            finite_element_->setTag(tag_);
            element_map_.elements2_.template getComponent<Shape_>()[makeHash(tag_)] = finite_element_;
        }

    private:
        
        Natural tag_;
        
        geometry::Point<Frame_::getDimEuclidean()> coordinates_;

        std::vector<std::shared_ptr<Region_>> regions_;

        std::shared_ptr<Element_> finite_element_;
        
        ElementFactoryMap<Frame_> & element_map_;
        
    };
    
    template<geometry::FrameConcept Frame_>
    struct ElementFactoryMap
    {

    private:

        template<typename... T_>
        using ElementPointer_ = std::unordered_map<std::basic_string<Character>, std::shared_ptr<ElementFactory<T_...>>>;

        template<typename... T_>
        using RegionPointer_ = std::unordered_map<std::basic_string<Character>, std::shared_ptr<RegionFactory<T_...>>>;

        template<typename... T_>
        using ElementPointer2_ = std::unordered_map<std::basic_string<Character>, std::shared_ptr<Element<T_...>>>;

        template<typename... T_>
        using RegionPointer2_ = std::unordered_map<std::basic_string<Character>, std::shared_ptr<Region<T_...>>>;

        using Elements_ = geometry::ShapeCollection<Frame_, ElementPointer_, Frame_>;

        using Regions_ = geometry::DomainCollection<Frame_, RegionPointer_, Frame_>;

        using Elements2_ = geometry::ShapeCollection<Frame_, ElementPointer2_, Frame_>;

        using Regions2_ = geometry::DomainCollection<Frame_, RegionPointer2_, Frame_>;

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
            // auto constexpr shape_coordinates = geometry::ShapeLibraryTraits<Frame_>::template getShapeCoordinates<Shape_>();
            // auto constexpr node_coordinates = geometry::ShapeInnerNeighborhoodTraits<Shape_>::template getShapeCoordinates<geometry::Node>();
            // for (auto const & nde : finite_element->template getInnerNeighbors<node_coordinates.i_, node_coordinates.j_>())
            // {
            //     auto const & ngs = nde->template getOuterNeighbors<shape_coordinates.i_ - 1, i_>();
            //     for (auto const & neighbour : ngs)
            //     {
            //         if (((neighbour->getTag() != finite_element->getTag()) && i_ == shape_coordinates.j_) || (i_ != shape_coordinates.j_))
            //         {
            //             auto & element_ngs = finite_element->template getOuterNeighbors<0, i_>();
            //             auto found = false;
            //             for (auto const & ngb: element_ngs)
            //             {
            //                 if (ngb->getTag() == neighbour->getTag())
            //                 {
            //                     found = true;
            //                     break;
            //                 }
            //             }
            //             if (!found)
            //             {
            //                 element_ngs.push_back(neighbour);
            //             }
            //         }
            //     }
            // }
            // if constexpr (i_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<shape_coordinates.i_>() - 1)
            // {
            //     setElementOuterNeighborhood<Shape_, i_ + 1>(finite_element);
            // }
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
            auto & regions = regions_.template getComponent<Domain_>();
            auto region_hash = RegionFactory<Domain_, Frame_>::makeHash(domain_tag);
            regions[region_hash] = std::make_shared<RegionFactory<Domain_, Frame_>>(* this, domain_label);
        }

        template<geometry::DomainConcept Domain_>
        std::shared_ptr<RegionFactory<Domain_, Frame_>> const &
        getRegion(
            auto const & domain_tag
        )
        const
        {
            auto region_hash = RegionFactory<Domain_, Frame_>::makeHash(domain_tag);
            return regions_.template getComponent<Domain_>().at(region_hash);
        }

        template<geometry::ShapeConcept Shape_>
        void
        addElement(
            auto const & element_tag,
            ShapeNodeConnectivity<Shape_> const & node_tags
        )
        {
            auto & elements = elements_.template getComponent<Shape_>();
            auto element_hash = ElementFactory<Shape_, Frame_>::makeHash(element_tag);
            elements[element_hash] = std::make_shared<ElementFactory<Shape_, Frame_>>(* this, node_tags);
        }

        template<geometry::ShapeConcept Shape_>
        void
        addElement(
            auto const & element_tag,
            geometry::Point<Frame_::getDimEuclidean()> const & node_coordinates
        )
        requires(std::same_as<Shape_, geometry::Node>)
        {
            auto & elements = elements_.template getComponent<Shape_>();
            auto element_hash = ElementFactory<Shape_, Frame_>::makeHash(element_tag);
            elements[element_hash] = std::make_shared<ElementFactory<Shape_, Frame_>>(* this, element_tag, node_coordinates);
        }

        template<geometry::ShapeConcept Shape_>
        void
        addRegionToElement(
            auto const & region_tag,
            auto const & element_tag
        )
        {
            auto & elements = elements_.template getComponent<Shape_>();
            auto element_hash = ElementFactory<Shape_, Frame_>::makeHash(element_tag);
            auto const & region = this->template getRegion<typename Shape_::Domain>(region_tag);
            elements.at(element_hash)->addRegion(region);
        }
        
        // template<typename Domain_>
        // void
        // makeDomain(
        //     std::shared_ptr<RegionFactory<Domain_, Frame_>> const & element_fac
        // )
        // const
        // {
        //     regions2_.template getComponent<Domain_>()[element_fac->getLabel()] = element_fac->region_;
        // }

        // template<geometry::ShapeConcept Shape_, typename Mesh_>
        // void
        // makeElement(
        //     std::shared_ptr<ElementFactory<Shape_, Frame_>> const & element_fac
        // )
        // {
        //     auto make_element = [&] <Integer i_ = 0, Integer j_ = 0> (
        //         auto & make_element_
        //     )
        //     mutable
        //     {
        //         // auto const constexpr t_inner_neighbor = ShapeTraits<t_element>::template getInnerNeighbor<i_, j_>();
        //         // using InnerNeighborhoodTraits = typename geometry::ShapeInnerNeighborhood<Shape_>;
        //         // using ICI = geometry::ShapeLibraryTraits<Frame_>;
        //         // auto const constexpr t_element_coordinates = geometry::ShapeInnerNeighborhoodTraits<Frame_>::template getElementCoordinates<Shape_>();
        //         // auto const constexpr t_node_coordinates = ShapeTraits<t_element>::template getInnerNeighborCoordinates<Node{}>();
        //         // auto const constexpr t_inner_neighbor_coordinates = MeshTraits<Frame_>::template getElementCoordinates<t_inner_neighbor>();
        //         // auto const constexpr t_neighbour_coordinates = ShapeTraits<t_inner_neighbor>::template getOuterNeighborCoordinates<t_domain, Shape_>();
        //         // auto & inner_neighbors = element_map.template getElements<ICI2::template getCoordinate<t_inner_neighbor>(0), ICI2::template getCoordinate<t_inner_neighbor>(1)>();
        //         // ::template getCoordinate<t_inner_neighbor>(0);
        //         // mesh_inner_neighbor.setNodeTag(j, getNodeTag(getInnerNeighborNodeConnection<i_, j_>(i, j)));
        //                 // auto mesh_inner_neighbor = ElementFactory<t_inner_neighbor, Frame_>(inner_neighbor_node_tags);
        //                 // mesh_inner_neighbor.getHash();
        //         // if constexpr (t_is_initialized)
        //         // {
        //         //     auto tag = element_map_.elements2_.template getComponent<Shape_>().size();
        //         //     finite_element_ = std::make_shared<Element_>(Element_(tag));
        //         // }
        //         //
        //         //
        //         //
        //         using t_inner_neighbor = typename geometry::ShapeInnerNeighborhood<Shape_>::template Component<i_, j_>;
        //         // auto constexpr t_is_initialized = ;
        //         auto & inner_neighbors = elements2_.template getComponent<t_inner_neighbor>();
        //         for (auto i = 0; i < geometry::ShapeInnerNeighborhood<Shape_>::getNumComponents(i_, j_); ++i)
        //         {
        //             auto inner_neighbor_hash = std::basic_string<Character>();
        //             if constexpr(!std::same_as<t_inner_neighbor, geometry::Node>)
        //             {
        //                 auto inner_neighbor_node_tags = ShapeNodeConnectivity<t_inner_neighbor>();
        //                 for (auto j = 0; j < t_inner_neighbor::getNumNodes(); ++j)
        //                 {
        //                     auto ici = ElementFactory<Shape_, Frame_>::getInnerNeighborNodeConnection<t_inner_neighbor, Mesh_>(i, j);
        //                     inner_neighbor_node_tags[j] = node_tags_[ici];
        //                     // inner_neighbor_node_tags[j] = node_tags_[getInnerNeighborNodeConnection<Mesh_, i_, j_>(i, j)];
        //                 }
        //                 inner_neighbor_hash = ElementFactory<t_inner_neighbor, Frame_>::makeHash(inner_neighbor_node_tags);
        //                 if (!inner_neighbors.contains(inner_neighbor_hash))
        //                 {
        //                     makeElement<t_inner_neighbor, Mesh_>(ElementFactory<t_inner_neighbor, Frame_>(element_map_, inner_neighbor_node_tags));
        //                     // ElementFactory<t_inner_neighbor, Frame_>(element_map_, inner_neighbor_node_tags).template makeElement<Mesh_>();
        //                 }
        //             }
        //             else
        //             {
        //                 auto ici = ElementFactory<Shape_, Frame_>::getInnerNeighborNodeConnection<geometry::Node, Mesh_>(i, 0);
        //                 auto node_tag = node_tags_[ici];
        //                 // auto node_tag = node_tags_[getInnerNeighborNodeConnection<Mesh_, i_, j_>(i, 0)];
        //                 inner_neighbor_hash = ElementFactory<t_inner_neighbor, Frame_>::makeHash(node_tag);
        //             }
        //             // std::cout << "inner_neighbor_hash : " << inner_neighbor_hash << std::endl;
        //             auto & inner_neighbor = inner_neighbors.at(inner_neighbor_hash);
        //             element_fac->finite_element_->getInnerNeighborhood().template getComponent<t_inner_neighbor>()[i] = inner_neighbor;
        //             inner_neighbor->template getOuterNeighborhood().template getComponent<Shape_>().push_back(element_fac->finite_element_);
        //         }
        //         if constexpr (j_ < geometry::ShapeInnerNeighborhood<Shape_>::getNumComponents(i_) - 1)
        //         {
        //             make_element_.template operator()<i_, j_ + 1>(make_element_);
        //         }
        //         else if constexpr (i_ < geometry::ShapeInnerNeighborhood<Shape_>::getNumComponents() - 1)
        //         {
        //             make_element_.template operator()<i_ + 1, 0>(make_element_);
        //         }
        //         if constexpr (i_ == 0 && j_ == 0)
        //         {
        //             for (auto const & region : regions_)
        //             {
        //                 element_fac->finite_element_->addRegion(regions2_.template getComponent<typename Shape_::Domain>().at(region->getLabel()));
        //             }
        //             element_fac->finite_element_->setTag(elements2_.template getComponent<Shape_>().size());
        //             elements2_.template getComponent<Shape_>()[makeHash(node_tags_)] = finite_element_;
        //         }
        //     };
        //     make_element(make_element);
        // }

        // template<geometry::PointConcept Shape_, typename Mesh_>
        // void
        // makeElement(
        //     std::shared_ptr<ElementFactory<Shape_, Frame_>> const & element_fac
        // )
        // {
        //     // auto constexpr t_element_coordinates = MeshTraits<t_domain>::template getElementCoordinates<t_element>();
        //     // finite_element_ = std::make_shared<Element<Shape_, Frame_>>(Element<Shape_, Frame_>(tag_));
        //     // finite_element_->setTag(element_map_.elements2_.template getComponent<Shape_>().size());
        //     // if (domain_ != nullptr)
        //     // {
        //     //     finite_element_->setDomain(element_set.regions2_.template getDomains<Shape_::getDimShape()>().at(domain_->getLabel()));
        //     // }
        //     // finite_element_->setDomain(element_set.template getDomains<t_element.getDim()>().at(domain_->getLabel()));
        //     element_fac->finite_element_->setCoordinates(element_fac->coordinates_);
        //     for (auto const & region : element_fac->regions_)
        //     {
        //         element_fac->finite_element_->addRegion(regions2_.template getComponent<typename Shape_::Domain>().at(region->getLabel()));
        //     }
        //     element_fac->finite_element_->setTag(element_fac->tag_);
        //     elements2_.template getComponent<Shape_>()[makeHash(tag_)] = element_fac->finite_element_;
        // }

        template<typename Mesh_>
        void
        setIt(
            std::basic_string<Character> && file_path
        )
        {
            /**
             * @brief 
             * 
             */
            auto my_one = MshOne<Frame_, Mesh_>(std::forward<std::basic_string<Character>>(file_path));
            /**
             * @brief 
             * 
             */
            auto set_regions = [&] <Integer i_ = 0> (
                auto & set_regions_
            )
            mutable
            {
                using Domain_ = typename geometry::DomainCollection<Frame_>::template Component<i_>;
                my_one.template setMeshDomain<Domain_>(* this);
                if constexpr (i_ < geometry::DomainCollection<Frame_>::getNumComponents() - 1)
                {
                    set_regions_.template operator()<i_ + 1>(set_regions_);
                }
            };
            set_regions(set_regions);
            /**
             * @brief 
             * 
             */
            auto set_elements = [&] <Integer i_ = 0, Integer j_ = 0> (
                auto & set_elements_
            )
            mutable
            {
                using Shape_ = geometry::ShapeCollection<Frame_>::template Component<i_, j_>;
                my_one.template setMeshElement<Shape_>(* this);
                if constexpr (j_ < geometry::ShapeCollection<Frame_>::getNumComponents(i_) - 1)
                {
                    set_elements_.template operator()<i_, j_ + 1>(set_elements_);
                }
                else if constexpr (i_ < geometry::ShapeCollection<Frame_>::getNumComponents() - 1)
                {
                    set_elements_.template operator()<i_ + 1, 0>(set_elements_);
                }
            };
            set_elements(set_elements);
            /**
             * @brief 
             * 
             */
            auto make_sets = [&] <Integer i_ = 0> (
                auto & make_sets_
            )
            mutable
            {
                using Domain_ = geometry::DomainCollection<Frame_>::template Component<i_>;
                for (auto const & dm : regions_.template getComponent<Domain_>())
                {
                    dm.second->makeDomain();
                }
                if constexpr (i_ < geometry::DomainCollection<Frame_>::getNumComponents() - 1)
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
                // using t_element = typename geometry::ShapeLibraryTraits<Frame_>::template Shape<i_, j_>;
                using Shape_ = geometry::ShapeCollection<Frame_>::template Component<i_, j_>;
                for (auto const & element : elements_.template getComponent<Shape_>())
                {
                    element.second->template makeElement<Mesh_>();
                }
                if constexpr (j_ < geometry::ShapeCollection<Frame_>::getNumComponents(i_) - 1)
                {
                    make_elements_.template operator()<i_, j_ + 1>(make_elements_);
                }
                else if constexpr (i_ < geometry::ShapeCollection<Frame_>::getNumComponents() - 1)
                {
                    make_elements_.template operator()<i_ + 1, 0>(make_elements_);
                }
                // if constexpr (j_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<i_>() - 1)
                // {
                //     make_elements_.template operator()<i_, j_ + 1>(make_elements_);
                // }
                // else if constexpr (i_ < geometry::ShapeLibraryTraits<Frame_>::getNumElements() - 1)
                // {
                //     make_elements_.template operator()<i_ + 1, 0>(make_elements_);
                // }
            };
            make_elements(make_elements);
        }
        
        void
        setIt2()
        const
        {
            // auto element_map = std::make_unique<MyElementMap<Frame_>>();
            auto make_sets = [&] <Integer i_ = 0> (
                auto & make_sets_
            )
            mutable
            {
                for (auto const & dm : regions_->template getDomains<i_>())
                {
                    regions2_.template getDomains<i_>()[dm.second->getLabel()] = dm.second->letDomain();
                }
                if constexpr (i_ < Frame_::getDimEuclidean())
                {
                    make_sets_.template operator()<i_ + 1>(make_sets_);
                }
            };
            make_sets(make_sets);
            // auto make_elements = [&] <Integer i_ = 0, Integer j_ = 0> (
            //     auto & make_elements_
            // )
            // mutable
            // {
            //     auto constexpr t_element = geometry::ShapeLibraryTraits<Frame_>::template getElement<i_, j_>();
            //     for (auto const & element : elements_->template getElements<i_, j_>())
            //     {
            //         element.second->template makeElement(* element_map);
            //     }
            //     if constexpr (j_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<i_>() - 1)
            //     {
            //         make_elements_.template operator()<i_, j_ + 1>(make_elements_);
            //     }
            //     else if constexpr (i_ < geometry::ShapeLibraryTraits<Frame_>::getNumElements() - 1)
            //     {
            //         make_elements_.template operator()<i_ + 1, 0>(make_elements_);
            //     }
            // };
            // make_elements(make_elements);
            // auto make_elements_outer_neighborhood = [&] <Integer i_ = 0, Integer j_ = 0> (
            //     auto & make_elements_outer_neighborhood_
            // )
            // mutable
            // {
            //     auto const constexpr t_element = geometry::ShapeLibraryTraits<Frame_>::template getElement<i_, j_>();
            //     for (auto & element : element_map->template getElements<i_, j_>())
            //     {
            //         setElementOuterNeighborhood<t_element>(element.second);
            //     }
            //     if constexpr (j_ < geometry::ShapeLibraryTraits<Frame_>::template getNumElements<i_>() - 1)
            //     {
            //         make_elements_outer_neighborhood_.template operator()<i_, j_ + 1>(make_elements_outer_neighborhood_);
            //     }
            //     else if constexpr (i_ < geometry::ShapeLibraryTraits<Frame_>::getNumElements() - 1)
            //     {
            //         make_elements_outer_neighborhood_.template operator()<i_ + 1, 0>(make_elements_outer_neighborhood_);
            //     }
            // };
            // make_elements_outer_neighborhood(make_elements_outer_neighborhood);
            // return element_map->makeElementSet();
        }

    // private:

        Elements_ elements_;

        Regions_ regions_;

        Elements2_ elements2_;

        Regions2_ regions2_;

    };

} // namespace lolita::mesh


#endif /* BA740F29_E221_4994_BF1D_7892518F367D */
