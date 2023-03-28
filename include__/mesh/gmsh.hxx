#ifndef C5421678_DA6F_4EA0_AFC9_73889C52D424
#define C5421678_DA6F_4EA0_AFC9_73889C52D424

#include "geometry/shape_.hxx"
#include "geometry/frame.hxx"

namespace lolita::mesh
{

    namespace internal
    {
    
        template<geometry::ShapeConcept Shape_>
        using NodeConnectivity = std::array<Integer, Shape_::getNumNodes()>;
    
        template<geometry::ShapeConcept Shape_>
        using NodeConnectivity2 = typename Shape_::template InnerNeighborhood<NodeConnectivity>;
        
    } // namespace internal

    struct GMSH
    {};

    template<typename, geometry::ShapeConcept Shape_>
    struct Table;

    template<>
    struct Table<GMSH, geometry::Segment>
    {

        static constexpr
        internal::NodeConnectivity2<geometry::Segment> node_connectivity_ = {
            {
                {
                    0,
                    1,
                }
            }
        };

    };

    template<>
    struct Table<GMSH, geometry::Triangle>
    {

        static constexpr
        internal::NodeConnectivity2<geometry::Triangle> node_connectivity_ = {
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
    struct Table<GMSH, geometry::Quadrangle>
    {

        static constexpr
        internal::NodeConnectivity2<geometry::Quadrangle> node_connectivity_ = {
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

    template<Integer dim_>
    struct Domain;

    template<typename...>
    struct FE;

    template<typename...>
    struct FEMap;

    template<typename Mesh_, geometry::ShapeConcept Shape_, geometry::FrameConcept Frame_>
    struct MeshElement
    {

    private:
    
        // using t_ElementTraits = ShapeTraits<t_element>;

        using FiniteElement_ = FE<Shape_, Frame_>;

        using NodeConnectivity_ = internal::NodeConnectivity<Shape_>;
        
        // template<DomainConcept auto t_dim, MeshConcept auto t__domain>
        // using t_Dom = MeshDomain<t_dim, t__domain>;

        // /**
        //  * @brief The MeshDomain the element is connected to.
        //  * 
        //  */
        // // using MeshDomain_ = typename t_ElementTraits::template DomainConnectivity<t_Dom, t_domain>;
        using MeshDomain_ = Domain<Shape_::getDimShape()>;

    public:

        template<Integer i_, Integer j_>
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
        getHash(
            std::array<Natural, Shape_::getNumNodes()> node_tags
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
        
        MeshElement()
        :
        node_tags_()
        // finite_element_()
        {}
        
        void
        setDomain(
            std::shared_ptr<MeshDomain_> const & domain
        )
        {
            domain_ = domain;
        }
        
        void
        setNodeTag(
            Integer index,
            Natural const & tag
        )
        {
            node_tags_[index] = tag;
        }
        
        Natural const &
        getNodeTag(
            Integer index
        )
        const
        {
            return node_tags_[index];
        }
        
        std::array<Natural, Shape_::getNumNodes()> const &
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
            FEMap<Frame_> & element_map
        )
        {
            auto make_element = [&] <Integer i_ = 0, Integer j_ = 0> (
                auto & make_element_
            )
            mutable
            {
                // auto const constexpr i_nner_neighbor = ShapeTraits<t_element>::template getInnerNeighbor<i_, j_>();
                using InnerNeighbor = typename geometry::ShapeInnerNeighborhoodTraits<Shape_>::template ShapeInnerNeighbor<i_, j_>;
                // auto const constexpr i_s_initialized = i_ == 0 && j_ == 0;
                auto constexpr i_s_initialized = i_ == 0 && j_ == 0;
                // auto const constexpr t_element_coordinates = MeshTraits<t_domain>::template getElementCoordinates<t_element>();
                auto constexpr t_element_coordinates = geometry::ShapeLibraryTraits::getShapeCoordinates<Shape_>();
                // auto const constexpr i_nner_neighbor_coordinates = MeshTraits<t_domain>::template getElementCoordinates<i_nner_neighbor>();
                auto constexpr i_nner_neighbor_coordinates = geometry::ShapeLibraryTraits::getShapeCoordinates<InnerNeighbor>();
                // auto const constexpr t_neighbour_coordinates = ShapeTraits<i_nner_neighbor>::template getOuterNeighborCoordinates<t_domain, t_element>();
                auto constexpr t_neighbour_coordinates = geometry::ShapeOuterNeighborhoodTraits<InnerNeighbor, Frame_>::template getShapeCoordinates<Shape_>();
                auto & inner_neighbors = element_map.template getElements<i_nner_neighbor_coordinates.i_, i_nner_neighbor_coordinates.j_>();
                if constexpr (i_s_initialized)
                {
                    auto tag = element_map.template getElements<t_element_coordinates.i_, t_element_coordinates.j_>().size();
                    finite_element_ = std::make_shared<FE<Shape_, Frame_>>(FE<Shape_, Frame_>(tag));
                }
                // for (auto i = 0; i < ShapeTraits<t_element>::template getNumInnerNeighbors<i_, j_>(); ++i)
                for (auto i = 0; i < geometry::ShapeInnerNeighborhoodTraits<Shape_>::template getNumInnerNeighbors<i_, j_>(); ++i)
                {
                    auto inner_neighbor_hash = std::basic_string<Character>();
                    // if constexpr(i_nner_neighbor != Node())
                    if constexpr(!std::same_as<InnerNeighbor, geometry::Node>)
                    {
                        auto mesh_inner_neighbor = MeshElement<Mesh_, InnerNeighbor, Frame_>();
                        for (auto j = 0; j < InnerNeighbor::getNumNodes(); ++j)
                        {
                            mesh_inner_neighbor.setNodeTag(j, getNodeTag(getInnerNeighborNodeConnection<i_, j_>(i, j)));
                        }
                        inner_neighbor_hash = mesh_inner_neighbor.getHash();
                        if (!inner_neighbors.contains(inner_neighbor_hash))
                        {
                            mesh_inner_neighbor.makeElement(element_map);
                        }
                    }
                    else
                    {
                        inner_neighbor_hash = MeshElement<Mesh_, InnerNeighbor, Frame_>::getHash(getNodeTag(getInnerNeighborNodeConnection<i_, j_>(i, 0)));
                    }
                    auto & inner_neighbor = inner_neighbors.at(inner_neighbor_hash);
                    finite_element_->template getInnerNeighbors<i_, j_>()[i] = inner_neighbor;
                    inner_neighbor->template getOuterNeighbors<t_neighbour_coordinates.i_, t_neighbour_coordinates.j_>().push_back(finite_element_);
                }
                if constexpr (j_ < geometry::ShapeInnerNeighborhoodTraits<Shape_>::template getNumInnerNeighbors<i_>() - 1)
                {
                    make_element_.template operator()<i_, j_ + 1>(make_element_);
                }
                else if constexpr (i_ < geometry::ShapeInnerNeighborhoodTraits<Shape_>::template getNumInnerNeighbors<>() - 1)
                {
                    make_element_.template operator()<i_ + 1, 0>(make_element_);
                }
                if constexpr (i_s_initialized)
                {
                    // for (auto const & domain : domains_)
                    // {
                    //     finite_element_->addDomain(element_map.template getDomains<t_element.getDim()>().at(domain->getLabel()));
                    // }
                    if (domain_ != nullptr)
                    {
                        finite_element_->setDomain(element_map.template getDomains<Shape_::getDimShape()>().at(domain_->getLabel()));
                    }                    
                    // finite_element_->setDomain(element_map.template getDomains<t_element.getDim()>().at(domain_->getLabel()));
                    // finite_element_->setCoordinates(finite_element_->getCurrentCentroid());
                    element_map.template getElements<t_element_coordinates.i_, t_element_coordinates.j_>()[getHash(node_tags_)] = finite_element_;
                }
            };
            make_element(make_element);
        }

        // /**
        //  * @brief Creates a finite element object from the data members, that are the node tags and the list of domains containing the element. A finite
        //  * element is created, and all its inner neighbors are sought in the finite element set given in argument. If an inner neighbor is not found, i.e. if
        //  * it does not exists, it is created by calling this function recursively. Once all inner neighbors are either found or created, they are linked to
        //  * the finite element, which is then initialized by linking it to all the domains it belongs to, and by setting its coordinates, which is the finite
        //  * element barycenter if the element is not a node.
        //  * 
        //  * @param element_map The finite element list of the mesh.
        //  */
        // void
        // makeElement(
        //     FiniteElementMap<t_domain> & element_map
        // )
        // {
        //     auto make_element = [&] <Integer i_ = 0, Integer j_ = 0> (
        //         auto & self
        //     )
        //     mutable
        //     {
        //         auto const constexpr i_nner_neighbor = ShapeTraits<t_element>::template getInnerNeighbor<i_, j_>();
        //         auto const constexpr i_s_initialized = i_ == 0 && j_ == 0;
        //         auto const constexpr t_element_coordinates = MeshTraits<t_domain>::template getElementCoordinates<t_element>();
        //         auto const constexpr t_node_coordinates = ShapeTraits<t_element>::template getInnerNeighborCoordinates<Node{}>();
        //         auto const constexpr i_nner_neighbor_coordinates = MeshTraits<t_domain>::template getElementCoordinates<i_nner_neighbor>();
        //         auto const constexpr t_neighbour_coordinates = ShapeTraits<i_nner_neighbor>::template getOuterNeighborCoordinates<t_domain, t_element>();
        //         auto & inner_neighbors = element_map.template getElements<i_nner_neighbor_coordinates.getDim(), i_nner_neighbor_coordinates.getTag()>();
        //         if constexpr (i_s_initialized)
        //         {
        //             auto tag = element_map.template getElements<t_element_coordinates.getDim(), t_element_coordinates.getTag()>().size();
        //             finite_element_ = std::make_shared<FiniteElement<t_element, t_domain>>(FiniteElement<t_element, t_domain>(tag));
        //         }
        //         for (auto i = 0; i < ShapeTraits<t_element>::template getNumInnerNeighbors<i_, j_>(); ++i)
        //         {
        //             auto inner_neighbor_hash = std::basic_string<Character>();
        //             if constexpr(i_nner_neighbor != Node())
        //             {
        //                 auto mesh_inner_neighbor = MeshElement<i_nner_neighbor, t_domain>();
        //                 for (auto j = 0; j < i_nner_neighbor.num_nodes_; ++j)
        //                 {
        //                     mesh_inner_neighbor.setNodeTag(j, getNodeTag(getInnerNeighborNodeConnection<i_, j_>(i, j)));
        //                 }
        //                 inner_neighbor_hash = mesh_inner_neighbor.getHash();
        //                 if (!inner_neighbors.contains(inner_neighbor_hash))
        //                 {
        //                     mesh_inner_neighbor.makeElement(element_map);
        //                 }
        //             }
        //             else
        //             {
        //                 inner_neighbor_hash = MeshElement<i_nner_neighbor, t_domain>::getHash(getNodeTag(getInnerNeighborNodeConnection<i_, j_>(i, 0)));
        //             }
        //             auto & inner_neighbor = inner_neighbors.at(inner_neighbor_hash);
        //             finite_element_->template getInnerNeighbors<i_, j_>()[i] = inner_neighbor;
        //             inner_neighbor->template getOuterNeighbors<t_neighbour_coordinates.dim_, t_neighbour_coordinates.tag_>().push_back(finite_element_);
        //         }
        //         if constexpr (j_ < ShapeTraits<t_element>::template getNumInnerNeighbors<i_>() - 1)
        //         {
        //             self.template operator()<i_, j_ + 1>(self);
        //         }
        //         else if constexpr (i_ < ShapeTraits<t_element>::template getNumInnerNeighbors<>() - 1)
        //         {
        //             self.template operator()<i_ + 1, 0>(self);
        //         }
        //         if constexpr (i_s_initialized)
        //         {
        //             // for (auto const & domain : domains_)
        //             // {
        //             //     finite_element_->addDomain(element_map.template getDomains<t_element.getDim()>().at(domain->getLabel()));
        //             // }
        //             if (domain_ != nullptr)
        //             {
        //                 finite_element_->setDomain(element_map.template getDomains<t_element.getDim()>().at(domain_->getLabel()));
        //             }                    
        //             // finite_element_->setDomain(element_map.template getDomains<t_element.getDim()>().at(domain_->getLabel()));
        //             finite_element_->setCoordinates(finite_element_->getCurrentCentroid());
        //             element_map.template getElements<t_element_coordinates.getDim(), t_element_coordinates.getTag()>()[getHash(node_tags_)] = finite_element_;
        //         }
        //     };
        //     make_element(make_element);
        // }

    // private:

    //     /**
    //      * @brief The tag of each node the element is connected to.
    //      * 
    //      */
        NodeConnectivity_ node_tags_;

    //     /**
    //      * @brief A pointer to the FiniteElement to be build.
    //      * 
    //      */
        std::shared_ptr<FiniteElement_> finite_element_;

    //     /**
    //      * @brief The MeshDomain the element is connected to.
    //      * 
    //      */
        std::shared_ptr<MeshDomain_> domain_;

    };

    // struct Gmesh
    // {

    //     // template<geometry::ShapeConcept Shape_>
    //     // struct Table;

    //     // template<>
    //     // struct Table<geometry::Segment>
    //     // {

    //     //     Hello<geometry::Segment> node_connectivity_ = {
    //     //         {
    //     //             {
    //     //                 {
    //     //                     0,
    //     //                     1,
    //     //                 }
    //     //             }
    //     //         }
    //     //     };

    //     // };

    //     static constexpr
    //     Hello<geometry::Segment> segment_connectivity_ = {
    //         {
    //             {
    //                 {
    //                     0,
    //                     1,
    //                 }
    //             }
    //         }
    //     };

    //     static constexpr
    //     Hello<geometry::Triangle> triangle_connectivity_ = {
    //         {
    //             {
    //                 {
    //                     0, 1,
    //                     1, 2,
    //                     2, 0,
    //                 }
    //             },
    //             {
    //                 {
    //                     0,
    //                     1,
    //                     2,
    //                 }
    //             }
    //         }
    //     };

    //     static constexpr
    //     Hello<geometry::Quadrangle> quadrangle_connectivity_ = {
    //         {
    //             {
    //                 {
    //                     0, 1,
    //                     1, 2,
    //                     2, 3,
    //                     3, 0,
    //                 }
    //             },
    //             {
    //                 {
    //                     0,
    //                     1,
    //                     2,
    //                     3,
    //                 }
    //             }
    //         }
    //     };

    //     static constexpr
    //     Hello<geometry::Tetrahedron> tetrahedron_connectivity_ = {
    //         {
    //             {
    //                 {
    //                     0, 1, 3,
    //                     1, 2, 3,
    //                     2, 3, 3,
    //                     3, 0, 3,
    //                 }
    //             },
    //             {
    //                 {
    //                     0, 1,
    //                     1, 2,
    //                     2, 3,
    //                     3, 0,
    //                     3, 0,
    //                     3, 0,
    //                 }
    //             },
    //             {
    //                 {
    //                     0,
    //                     1,
    //                     2,
    //                     3,
    //                 }
    //             }
    //         }
    //     };

    // };



    // template<geometry::ShapeConcept Shape_>
    // struct Hello
    // {

    //     using Table = typename Shape_::template InnerNeighborhood<internal::NodeConnectivity>;

    //     constexpr
    //     Hello(
    //         Table const & point_connectivity
    //     )
    //     :
    //     point_connectivity_(point_connectivity)
    //     {}

    //     constexpr
    //     Hello(
    //         Table && point_connectivity
    //     )
    //     :
    //     point_connectivity_(std::move(point_connectivity))
    //     {}

    //     template<Integer i_, Integer j_>
    //     constexpr
    //     Integer
    //     getNodeConnectivity(
    //         Integer i,
    //         Integer j
    //     )
    //     const
    //     {
    //         return std::get<j_>(std::get<i_>(point_connectivity_))[i][j];
    //     }

    // private:

    //     Table const point_connectivity_;

    // };
    
} // namespace lolita::mesh


#endif /* C5421678_DA6F_4EA0_AFC9_73889C52D424 */
