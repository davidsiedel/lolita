#ifndef F548FA01_2D9F_4847_B854_824601F31371
#define F548FA01_2D9F_4847_B854_824601F31371

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_000.hxx"
#include "lolita/lolita_core_n_001.hxx"
#include "lolita/lolita_core_n_002.hxx"
#include "lolita/lolita_core_n_003.hxx"
#include "lolita/lolita_core_n_004.hxx"

namespace lolita2::geometry
{

    template<Domain t_domain>
    struct FiniteElementSet : ElementSet<FiniteElementHolder, t_domain>
    {

        template<ElementType t_ii>
        void
        activate(
            std::basic_string_view<lolita::character> label,
            std::basic_string_view<lolita::character> domain
        )
        {
            auto activate_elements = [&] <lolita::integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element.second->isIn(domain))
                    {
                        element.second->getFiniteElement(label)->activate();
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
        }

        template<ElementType t_ii, Field t_field, Basis t_basis>
        void
        setDegreeOfFreedom(
            std::basic_string_view<lolita::character> label,
            std::basic_string_view<lolita::character> domain,
            std::shared_ptr<DegreeOfFreedom> & degree_of_freedom
        )
        {
            auto activate_elements = [&] <lolita::integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element.second->isIn(domain))
                    {
                        element.second->getFiniteElement(label)->template setDegreeOfFreedom<t_field, t_basis>(degree_of_freedom);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
        }

        // template<FiniteElementMethodConcept auto t_arg, ElementType t_ii>
        // void
        // setLoad(
        //     std::basic_string_view<lolita::character> domain,
        //     lolita::integer row,
        //     lolita::integer col,
        //     Loading && load
        // )
        // {
        //     auto activate_elements = [&] <lolita::integer t_j = 0> (
        //         auto & self
        //     )
        //     mutable
        //     {
        //         auto constexpr t_i = t_ii.getDim();
        //         auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
        //         auto ptr_load = std::make_shared<Loading>(std::forward<Loading>(load));
        //         for (auto const & element : this->template getElements<t_i, t_j>())
        //         {
        //             if (element.second->isIn(domain))
        //             {
        //                 element.second->template getFiniteElement<t_arg>()->template setLoad(row, col, ptr_load);
        //             }
        //         }
        //         if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
        //         {
        //             self.template operator()<t_j + 1>(self);
        //         }
        //     }; 
        //     activate_elements(activate_elements);
        // }
        
        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            FiniteElementSet const & finite_element_set
        )
        {
            auto print_element_inner_neighbors = [&] <Element t_element, lolita::integer t_i = 0, lolita::integer t_j = 0> (
                    auto const & element,
                    auto & self
            )
            mutable
            {
                if constexpr (!t_element.isNode())
                {
                    auto const constexpr t_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>();
                    for (auto const & c_ : element->template getInnerNeighbors<t_i, t_j>())
                    {
                        os << "layer : " << t_i << " type : " << t_j << " <-- " << t_neighbor << " " << c_->getHash() << std::endl;
                    }
                    if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<t_i>() - 1)
                    {
                        self.template operator()<t_element, t_i, t_j + 1>(element, self);
                    }
                    else if constexpr (t_i < ElementTraits<t_element, t_domain>::getNumInnerNeighbors() - 1)
                    {
                        self.template operator()<t_element, t_i + 1, 0>(element, self);
                    }
                }
            };
            auto print_element_outer_neighbors = [&] <Element t_element, lolita::integer t_i = 0, lolita::integer t_j = 0> (
                    auto const & element,
                    auto & self
            )
            mutable
            {
                auto const constexpr t_neighbor = ElementTraits<t_element, t_domain>::template getOuterNeighbor<t_i, t_j>();
                for (auto const & c_ : element->template getOuterNeighbors<t_i, t_j>())
                {
                    if constexpr (!t_element.isNode() && t_i == 0)
                    {
                        os << "layer : " << t_i << " type : " << t_j << " <-> " << t_neighbor << " " << c_->getHash() << std::endl;
                    }
                    else
                    {
                        os << "layer : " << t_i << " type : " << t_j << " --> " << t_neighbor << " " << c_->getHash() << std::endl;
                    }
                }
                if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumOuterNeighbors<t_i>() - 1)
                {
                    self.template operator()<t_element, t_i, t_j + 1>(element, self);
                }
                else if constexpr (t_i < ElementTraits<t_element, t_domain>::getNumOuterNeighbors() - 1)
                {
                    self.template operator()<t_element, t_i + 1, 0>(element, self);
                }
            };
            auto print_elements = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                if constexpr (t_i == 0 && t_j == 0)
                {
                    os << "*** Elements : " << std::endl;
                }
                for (auto const & element : finite_element_set.template getElements<t_i, t_j>())
                {
                    os << "* Element : " << t_element << " " << element.second->getHash() << std::endl;
                    os << "* Domains : ";
                    for (auto const & domain : element.second->domains_)
                    {
                        os << domain->tag_ << " ";
                    }
                    os << std::endl;
                    print_element_inner_neighbors.template operator()<t_element>(element.second, print_element_inner_neighbors);
                    print_element_outer_neighbors.template operator()<t_element>(element.second, print_element_outer_neighbors);
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
            print_elements(print_elements);
            return os;
        }

    };

    template<Element t_element, Domain t_domain>
    struct MeshElement;
    
    template<Element t_element, Domain t_domain>
    requires(!t_element.isNode())
    struct MeshElement<t_element, t_domain>
    {
        
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita::integer
        getInnerNeighborNodeConnection(
            lolita::integer i,
            lolita::integer j
        )
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element, t_domain>::node_connectivity_))[i][j];
        }

        static
        std::basic_string<lolita::character>
        getHash(
            std::array<lolita::natural, t_element.num_nodes_> node_tags
        )
        {
            auto hash = std::basic_stringstream<lolita::character>();
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

        template<lolita::integer t_i, lolita::integer t_j, typename... t_Labels>
        void
        setElement(
            FiniteElementSet<t_domain> & element_set,
            std::shared_ptr<FiniteElementHolder<t_element, t_domain>> & ptr_element,
            t_Labels... labels
        )
        const
        {
            auto const constexpr t_component = ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>();
            auto const constexpr t_is_initialized = t_i == 0 && t_j == 0;
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto const constexpr t_node_coordinates = ElementTraits<t_element, t_domain>::template getInnerNeighborCoordinates<Element::node()>();
            auto const constexpr t_component_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_component>();
            auto const constexpr t_neighbour_coordinates = ElementTraits<t_component, t_domain>::template getOuterNeighborCoordinates<t_element>();
            auto & components = element_set.template getElements<t_component_coordinates.dim_, t_component_coordinates.tag_>();
            for (lolita::integer i = 0; i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<t_i, t_j>(); ++i)
            {
                auto component_hash = std::basic_string<lolita::character>();
                if constexpr(!t_component.isNode())
                {
                    auto cpt = MeshElement<t_component, t_domain>();
                    for (lolita::integer j = 0; j < t_component.num_nodes_; ++j)
                    {
                        cpt.node_tags_[j] = node_tags_[getInnerNeighborNodeConnection<t_i, t_j>(i, j)];
                    }
                    component_hash = MeshElement<t_component, t_domain>::getHash(cpt.node_tags_);
                    if (!components.contains(component_hash))
                    {
                        using t_Component = FiniteElementHolder<t_component, t_domain>;
                        auto ptr_component = std::make_shared<t_Component>(t_Component());
                        cpt.template setElement<0, 0>(element_set, ptr_component, labels...);
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
                setElement<t_i, t_j + 1>(element_set, ptr_element, labels...);
            }
            else if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<>() - 1)
            {
                setElement<t_i + 1, 0>(element_set, ptr_element, labels...);
            }
            if constexpr (t_is_initialized)
            {
                auto const & nodes = ptr_element->template getInnerNeighbors<t_node_coordinates.dim_, t_node_coordinates.tag_>();
                auto domains = std::set<std::shared_ptr<MeshDomain>>();
                for (auto const & domain : nodes[0]->domains_)
                {
                    auto has_domain = true;
                    for (lolita::integer j = 1; j < t_element.num_nodes_; ++j)
                    {
                        has_domain = std::find(nodes[j]->domains_.begin(), nodes[j]->domains_.end(), domain) != nodes[j]->domains_.end();
                    }
                    if (has_domain)
                    {
                        domains.insert(domain);
                    }
                }
                ptr_element->tag_ = element_set.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>().size();
                ptr_element->domains_.assign(domains.begin(), domains.end());
                auto labelss = std::array<std::basic_string_view<lolita::character>, sizeof...(t_Labels)>{labels...};
                for (auto i = 0; i < sizeof...(t_Labels); i++)
                {
                    auto fem = std::make_shared<FiniteElement<t_element, t_domain>>();
                    fem->label_ = labelss[i];
                    ptr_element->finite_elements_.push_back(fem);
                }
                element_set.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>()[getHash(node_tags_)] = ptr_element;
            }
        }

        template<typename... t_Labels>
        void
        makeElement(
            FiniteElementSet<t_domain> & element_set,
            t_Labels... labels
        )
        const
        {
            auto ptr_element = std::make_shared<FiniteElementHolder<t_element, t_domain>>(FiniteElementHolder<t_element, t_domain>());
            setElement<0, 0>(element_set, ptr_element, labels...);
        }

        template<lolita::integer t_i>
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
            for (auto a = 0; a < ptr_element->finite_elements_.size(); a++)
            {
                ptr_element->getFiniteElement(a)->tag_ = ptr_element->tag_;
                ptr_element->getFiniteElement(a)->domains_ = ptr_element->domains_;
                auto initialize_inner_neighbors = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                        auto & t_initialize_inner_neighbors
                )
                mutable
                {
                    for (lolita::integer i = 0; i < ptr_element->template getInnerNeighbors<t_i, t_j>().size(); ++i)
                    {
                        auto & rhs = ptr_element->template getInnerNeighbors<t_i, t_j>()[i]->getFiniteElement(a);
                        auto & lhs = ptr_element->getFiniteElement(a)->template getInnerNeighbors<t_i, t_j>()[i];
                        lhs = rhs;
                    }
                    if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<t_i>() - 1)
                    {
                        t_initialize_inner_neighbors.template operator()<t_i, t_j + 1>(t_initialize_inner_neighbors);
                    }
                    else if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<>() - 1)
                    {
                        t_initialize_inner_neighbors.template operator()<t_i + 1, 0>(t_initialize_inner_neighbors);
                    }
                };
                initialize_inner_neighbors(initialize_inner_neighbors);
                auto initialize_outer_neighbors = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                        auto & t_initialize_outer_neighbors
                )
                mutable
                {
                    for (lolita::integer i = 0; i < ptr_element->template getOuterNeighbors<t_i, t_j>().size(); ++i)
                    {
                        auto & rhs = ptr_element->template getOuterNeighbors<t_i, t_j>()[i]->getFiniteElement(a);
                        auto & lhs = ptr_element->getFiniteElement(a)->template getOuterNeighbors<t_i, t_j>();
                        lhs.push_back(rhs);
                    }
                    if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumOuterNeighbors<t_i>() - 1)
                    {
                        t_initialize_outer_neighbors.template operator()<t_i, t_j + 1>(t_initialize_outer_neighbors);
                    }
                    else if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumOuterNeighbors<>() - 1)
                    {
                        t_initialize_outer_neighbors.template operator()<t_i + 1, 0>(t_initialize_outer_neighbors);
                    }
                };
                initialize_outer_neighbors(initialize_outer_neighbors);
            }
        }

        std::array<lolita::natural, t_element.num_nodes_> node_tags_;

    };

    template<Element t_element, Domain t_domain>
    requires(t_element.isNode())
    struct MeshElement<t_element, t_domain>
    {

        static
        std::basic_string<lolita::character>
        getHash(
            lolita::natural tag
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

        template<typename... t_Labels>
        void
        setElement(
            FiniteElementSet<t_domain> & element_set,
            std::shared_ptr<FiniteElementHolder<t_element, t_domain>> & ptr_element,
            t_Labels... labels
        )
        const
        {
            auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            ptr_element->tag_ = tag_;
            ptr_element->domains_ = domains_;
            ptr_element->coordinates_ = coordinates_;
            auto labelss = std::array<std::basic_string_view<lolita::character>, sizeof...(t_Labels)>{labels...};
            for (auto i = 0; i < sizeof...(t_Labels); i++)
            {
                auto fem = std::make_shared<FiniteElement<t_element, t_domain>>();
                fem->label_ = labelss[i];
                ptr_element->finite_elements_.push_back(fem);
            }
            element_set.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>()[getHash(tag_)] = ptr_element;
        }

        template<typename... t_Labels>
        void
        makeElement(
            FiniteElementSet<t_domain> & element_set,
            t_Labels... labels
        )
        const
        {
            auto ptr_element = std::make_shared<FiniteElementHolder<t_element, t_domain>>();
            setElement(element_set, ptr_element, labels...);
        }
        
        static inline
        void
        initialize(
                std::shared_ptr<FiniteElementHolder<t_element, t_domain>> & ptr_element
        )
        {
            for (auto a = 0; a < ptr_element->finite_elements_.size(); ++a)
            {
                ptr_element->getFiniteElement(a)->inner_neighbors_ = ptr_element->inner_neighbors_;
                ptr_element->getFiniteElement(a)->coordinates_ = ptr_element->coordinates_;
                ptr_element->getFiniteElement(a)->tag_ = ptr_element->tag_;
                ptr_element->getFiniteElement(a)->domains_ = ptr_element->domains_;
                auto & rhs = ptr_element->template getOuterNeighbors<0, 0>()[0]->getFiniteElement(a);
                auto & lhs = ptr_element->getFiniteElement(a)->template getOuterNeighbors<0, 0>();
                auto initialize_outer_neighbors = [&] <lolita::integer t_i = 0u, lolita::integer t_j = 0u> (
                        auto & t_initialize_outer_neighbors
                )
                mutable
                {
                    for (lolita::integer i = 0; i < ptr_element->template getOuterNeighbors<t_i, t_j>().size(); ++i)
                    {
                        auto & rhs = ptr_element->template getOuterNeighbors<t_i, t_j>()[i]->getFiniteElement(a);
                        auto & lhs = ptr_element->getFiniteElement(a)->template getOuterNeighbors<t_i, t_j>();
                        lhs.push_back(rhs);
                    }
                    if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumOuterNeighbors<t_i>() - 1)
                    {
                        t_initialize_outer_neighbors.template operator()<t_i, t_j + 1>(t_initialize_outer_neighbors);
                    }
                    else if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumOuterNeighbors<>() - 1)
                    {
                        t_initialize_outer_neighbors.template operator()<t_i + 1, 0>(t_initialize_outer_neighbors);
                    }
                };
                initialize_outer_neighbors(initialize_outer_neighbors);
            }
            // auto initialize_element = [&] <lolita::integer t_k = 0u> (
            //         auto & t_initialize_element
            // )
            // constexpr mutable
            // {
            //     ptr_element->template getFiniteElement<t_k>()->inner_neighbors_ = ptr_element->inner_neighbors_;
            //     ptr_element->template getFiniteElement<t_k>()->coordinates_ = ptr_element->coordinates_;
            //     ptr_element->template getFiniteElement<t_k>()->tag_ = ptr_element->tag_;
            //     ptr_element->template getFiniteElement<t_k>()->domains_ = ptr_element->domains_;
            //     auto initialize_outer_neighbors = [&] <lolita::integer t_i = 0u, lolita::integer t_j = 0u> (
            //             auto & t_initialize_outer_neighbors
            //     )
            //     mutable
            //     {
            //         for (lolita::integer i = 0; i < ptr_element->template getOuterNeighbors<t_i, t_j>().size(); ++i)
            //         {
            //             auto & rhs = ptr_element->template getOuterNeighbors<t_i, t_j>()[i]->template getFiniteElement<t_k>();
            //             auto & lhs = ptr_element->template getFiniteElement<t_k>()->template getOuterNeighbors<t_i, t_j>();
            //             lhs.push_back(rhs);
            //         }
            //         if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumOuterNeighbors<t_i>() - 1)
            //         {
            //             t_initialize_outer_neighbors.template operator()<t_i, t_j + 1>(t_initialize_outer_neighbors);
            //         }
            //         else if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumOuterNeighbors<>() - 1)
            //         {
            //             t_initialize_outer_neighbors.template operator()<t_i + 1, 0>(t_initialize_outer_neighbors);
            //         }
            //     };
            //     initialize_outer_neighbors(initialize_outer_neighbors);
            //     if constexpr (t_k < std::tuple_size_v<typename FiniteElementHolder<t_element, t_domain, t_args...>::t_FiniteElements> - 1)
            //     {
            //         t_initialize_element.template operator()<t_k + 1>(t_initialize_element);
            //     }
            // };
            // initialize_element(initialize_element);
        }
        
        lolita::natural tag_;
        
        std::shared_ptr<Point> coordinates_;

        std::vector<std::shared_ptr<MeshDomain>> domains_;

    };
    
    template<Domain t_domain>
    struct MeshElementSet : ElementSet<MeshElement, t_domain>
    {

        MeshElementSet()
        {}
        
        template<typename... t_Labels>
        std::unique_ptr<FiniteElementSet<t_domain>>
        makeFiniteElementSet(
            t_Labels... labels
        )
        const
        {
            auto element_set = std::make_unique<FiniteElementSet<t_domain>>();
            auto make_elements = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    element.second->template makeElement(* element_set, labels...);
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
            auto make_elements_outer_neighborhood = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto & element : element_set->template getElements<t_i, t_j>())
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
            return element_set;
        }

    };    
    
    struct MeshFileParserBase
    {

        MeshFileParserBase(
            std::basic_string_view<lolita::character> str
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
                lolita::integer tag
        )
        {
            if (tag == 15) return Element::node();
            else if (tag == 01) return Element::segment(1);
            else if (tag == 02) return Element::triangle(1);
            else if (tag == 03) return Element::quadrangle(1);
            else return Element::node();
        }

        struct PhysicalEntity
        {

            static
            std::basic_string<lolita::character>
            getHash(
                lolita::integer tag
            )
            {
                return std::to_string(tag);
            }

            PhysicalEntity(
                std::basic_string<lolita::character> const & name,
                lolita::integer dim
            )
            :
            mesh_domain_(std::make_shared<MeshDomain>(dim, name))
            {}

            std::shared_ptr<MeshDomain> mesh_domain_;

        };

        struct GeometricalEntity
        {

            static
            std::basic_string<lolita::character>
            getHash(
                lolita::integer dim,
                lolita::integer tag
            )
            {
                auto hash = std::basic_stringstream<lolita::character>();
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
            addPhysicalEntity(
                std::shared_ptr<PhysicalEntity> const & physical_entity
            )
            {
                physical_entities_.insert(physical_entity);
                for (auto & bounding_entity : bounding_entities_)
                {
                    bounding_entity->addPhysicalEntity(physical_entity);
                }
            }

            void
            addBoundingEntity(
                std::shared_ptr<GeometricalEntity> & bounding_entity
            )
            {
                bounding_entities_.insert(bounding_entity);
                for (auto const & physical_entity : physical_entities_)
                {
                    bounding_entity->addPhysicalEntity(physical_entity);
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
            auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
            auto num_physical_names = lolita::integer();
            line_stream >> num_physical_names;
            offset += 1;
            for (lolita::integer i = 0; i < num_physical_names; ++i)
            {
                line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto dim = lolita::integer();
                auto name = std::basic_string<lolita::character>();
                auto physical_entity_tag = lolita::integer();
                line_stream >> dim >> physical_entity_tag >> name;
                lolita::utility::removeCharacter(name, '"');
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
            auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
            auto num_points = lolita::integer();
            auto num_curves = lolita::integer();
            auto num_surfaces = lolita::integer();
            auto num_volumes = lolita::integer();
            line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
            auto num_domains = std::array<lolita::integer, 4>{num_points, num_curves, num_surfaces, num_volumes};
            offset += 1;
            for (lolita::integer i = 0; i < 4; ++i)
            {
                for (lolita::integer j = 0; j < num_domains[i]; ++j)
                {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    auto tag = lolita::integer();
                    line_stream >> tag;
                    auto geometrical_entity_hash = GeometricalEntity::getHash(i, tag);
                    geometrical_entities_[geometrical_entity_hash] = std::make_shared<GeometricalEntity>(GeometricalEntity());
                    if (i == 0)
                    {
                        for (lolita::integer k = 0; k < 3; ++k)
                        {
                            auto a = lolita::real();
                            line_stream >> a;
                        }
                    }
                    else
                    {
                        for (lolita::integer k = 0; k < 6; ++k)
                        {
                            auto a = lolita::real();
                            line_stream >> a;
                        }
                    }
                    auto num_physical_entities = lolita::integer();
                    line_stream >> num_physical_entities;
                    for (lolita::integer k = 0; k < num_physical_entities; ++k)
                    {
                        auto physical_entity_tag = lolita::integer();
                        line_stream >> physical_entity_tag;
                        auto physical_entity_hash = PhysicalEntity::getHash(physical_entity_tag);
                        geometrical_entities_[geometrical_entity_hash]->addPhysicalEntity(physical_entities_[physical_entity_hash]);
                    }
                    if (i > 0)
                    {
                        auto num_bounding_entities = lolita::integer();
                        line_stream >> num_bounding_entities;
                        for (lolita::integer k = 0; k < num_bounding_entities; ++k)
                        {
                            auto bounding_entity_tag = lolita::integer();
                            line_stream >> bounding_entity_tag;
                            auto bounding_entity_hash = GeometricalEntity::getHash(i - 1, std::abs(bounding_entity_tag));
                            geometrical_entities_[geometrical_entity_hash]->addBoundingEntity(geometrical_entities_[bounding_entity_hash]);
                        }
                    }
                    offset += 1;
                }
            }
        }

    public:

        GmshFileParser(
            std::basic_string_view<lolita::character> str
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
            auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
            auto num_entity_block = lolita::integer();
            auto num_nodes = lolita::integer();
            auto min_node_tag = lolita::integer();
            auto max_node_tag = lolita::integer();
            line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
            offset += 1;
            for (lolita::integer i = 0; i < num_entity_block; ++i)
            {
                line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto entity_dim = lolita::integer();
                auto entity_tag = lolita::integer();
                auto parametric = lolita::integer();
                auto num_nodes_in_block = lolita::integer();
                line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                auto const entity_hash = GeometricalEntity::getHash(entity_dim, entity_tag);
                offset += 1;
                for (lolita::integer j = 0; j < num_nodes_in_block; ++j)
                {
                    auto mesh_element = std::make_shared<MeshElement<t_element, t_domain>>(MeshElement<t_element, t_domain>());
                    auto tag = lolita::natural();
                    auto coordinates_ = Point();
                    auto domains_ = std::vector<std::shared_ptr<std::basic_string<lolita::character>>>();
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    line_stream >> tag;
                    mesh_element->tag_ = tag - 1;
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset + num_nodes_in_block]);
                    for (lolita::integer k = 0; k < 3; ++k)
                    {
                        line_stream >> mesh_element->coordinates_->operator()(k);
                    }
                    for (auto const & physical_entity : geometrical_entities_.at(entity_hash)->physical_entities_)
                    {
                        mesh_element->domains_.push_back(physical_entity->mesh_domain_);
                    }
                    auto const element_hash = MeshElement<t_element, t_domain>::getHash(mesh_element->tag_);
                    elements[element_hash] = mesh_element;
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
            auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
            auto num_entity_blocks = 0;
            auto num_elements = 0;
            auto min_element_tag = 0;
            auto max_element_tag = 0;
            line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
            offset += 1;
            for (lolita::integer i = 0; i < num_entity_blocks; ++i)
            {
                line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto entity_dim = 0;
                auto entity_tag = 0;
                auto element_type_tag = 0;
                auto num_elements_in_block = 0;
                line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                offset += 1;
                auto const element = getElemType(element_type_tag);
                if (element == t_element)
                {
                    for (lolita::integer j = 0; j < num_elements_in_block; ++j)
                    {
                        auto mesh_element = std::make_shared<MeshElement<t_element, t_domain>>(MeshElement<t_element, t_domain>());
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                        auto tag = lolita::integer();
                        line_stream >> tag;
                        for (lolita::integer k = 0; k < t_element.num_nodes_; ++k)
                        {
                            auto node_tag = lolita::integer();
                            line_stream >> node_tag;
                            mesh_element->node_tags_[k] = node_tag - 1;
                        }
                        auto const element_hash = MeshElement<t_element, t_domain>::getHash(mesh_element->node_tags_);
                        elements[element_hash] = mesh_element;
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
        
        std::map<std::basic_string<lolita::character>, std::shared_ptr<PhysicalEntity>> physical_entities_;

        std::map<std::basic_string<lolita::character>, std::shared_ptr<GeometricalEntity>> geometrical_entities_;

    };
    
    struct MeshFileParser
    {

        MeshFileParser(
            std::basic_string_view<lolita::character> str
        )
        :
        file_path_(str)
        {}

        template<Domain t_domain, typename... t_Labels>
        std::unique_ptr<FiniteElementSet<t_domain>>
        makeFiniteElementSet(
            t_Labels... labels
        )
        const
        {
            auto gmsh_file_parser = GmshFileParser(file_path_);
            auto mesh_element_set = MeshElementSet<t_domain>();
            auto make_elements = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
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
            return mesh_element_set.makeFiniteElementSet(labels...);
        }

        std::basic_string_view<lolita::character> file_path_;

    };
    
}


#endif /* F548FA01_2D9F_4847_B854_824601F31371 */
