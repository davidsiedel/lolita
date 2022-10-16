#ifndef F30FA5DA_AED7_4C4E_A979_CB98F2F66287
#define F30FA5DA_AED7_4C4E_A979_CB98F2F66287

#include "2/core/_include.hxx"
#include "2/core/region.hxx"
#include "2/core/element.hxx"

namespace lolita::core
{
    
    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam t_domain 
     */
    template<template<DomainConcept auto, MeshConcept auto> typename T, MeshConcept auto t_domain>
    struct DomainMap
    {

    private:

        template<DomainConcept auto t__dim, MeshConcept auto t__domain>
        using t__RegionMap = std::unordered_map<std::basic_string<Character>, std::shared_ptr<T<t__dim, t__domain>>>;

        using t_RegionMap = lolita::utility::tuple_slice_t<typename DomainLibrary::Domains<t__RegionMap, t_domain>, 0, t_domain.getDim() + 1>;

    public:

        DomainMap()
        {}
    
        template<Integer t_i>
        std::tuple_element_t<t_i, t_RegionMap> const &
        getDomains()
        const
        {
            return std::get<t_i>(domains_);
        }
        
        template<Integer t_i>
        std::tuple_element_t<t_i, t_RegionMap> &
        getDomains()
        {
            return std::get<t_i>(domains_);
        }
        
        t_RegionMap domains_;
        
    };

    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam t_domain 
     */
    template<template<DomainConcept auto, MeshConcept auto> typename T, MeshConcept auto t_domain>
    struct DomainSet
    {

    private:

        template<DomainConcept auto t__dim, MeshConcept auto t__domain>
        using t__RegionSet = std::vector<std::shared_ptr<T<t__dim, t__domain>>>;

        using t_RegionSet = lolita::utility::tuple_slice_t<typename DomainLibrary::Domains<t__RegionSet, t_domain>, 0, t_domain.getDim() + 1>;

    public:

        DomainSet()
        {}
    
        template<Integer t_i>
        std::tuple_element_t<t_i, t_RegionSet> const &
        getDomains()
        const
        {
            return std::get<t_i>(domains_);
        }
        
        template<Integer t_i>
        std::tuple_element_t<t_i, t_RegionSet> &
        getDomains()
        {
            return std::get<t_i>(domains_);
        }
        
        t_RegionSet domains_;
        
    };
    
    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam t_domain 
     */
    template<template<ShapeConcept auto, MeshConcept auto> typename T, MeshConcept auto t_domain>
    struct ElementMap
    {

    private:

        template<ShapeConcept auto t_element, MeshConcept auto t__domain, auto... t__args>
        using t_ElementMap = std::unordered_map<std::basic_string<Character>, std::shared_ptr<T<t_element, t__domain, t__args...>>>;

        using t_Elements = lolita::utility::tuple_slice_t<typename ShapeLibrary::Elements<t_ElementMap, t_domain>, 0, t_domain.getDim() + 1>;

    public:

        ElementMap()
        {}
    
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> const &
        getElements()
        const
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> &
        getElements()
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        t_Elements elements_;

    };

    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam t_domain 
     */
    template<template<ShapeConcept auto, MeshConcept auto> typename T, MeshConcept auto t_domain>
    struct ElementSet
    {

    private:

        template<ShapeConcept auto t_element, MeshConcept auto t__domain>
        using t_ElementSet = std::vector<std::shared_ptr<T<t_element, t__domain>>>;

        using t_Elements = lolita::utility::tuple_slice_t<typename ShapeLibrary::Elements<t_ElementSet, t_domain>, 0, t_domain.getDim() + 1>;

    public:

        ElementSet()
        {}
    
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> const &
        getElements()
        const
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_Elements>> &
        getElements()
        {
            return std::get<t_j>(std::get<t_i>(elements_));
        }
        
        t_Elements elements_;

    };

    /**
     * @brief 
     * 
     * @tparam t_domain 
     */
    template<MeshConcept auto t_domain>
    struct FiniteElementSet : ElementSet<FiniteElement, t_domain>, DomainSet<FiniteDomain, t_domain>
    {

        using MeshTraits_ = MeshTraits<t_domain>;

        template<Integer t_i>
        void
        caller2(
            std::basic_string<Character> && domain,
            auto && fun,
            Integer num_threads = std::thread::hardware_concurrency()
        )
        const
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto evaluation_loop = [&] (
                    Integer start,
                    Integer stop
                )
                {
                    for (auto i = start; i < stop; ++i)
                    {
                        auto const & finite_element = this->template getElements<t_i, t_j>()[i];
                        if (finite_element->isIn(std::forward<std::basic_string<Character>>(domain)))
                        {
                            std::forward<decltype(fun)>(fun)(finite_element);
                        }
                    }
                };
                if (num_threads < 2)
                {
                    evaluation_loop(0, this->template getElements<t_i, t_j>().size());
                }
                else
                {
                    auto batch_size = this->template getElements<t_i, t_j>().size() / num_threads;
                    auto batch_remainder = this->template getElements<t_i, t_j>().size() % num_threads;
                    auto threads = std::vector<std::jthread>(num_threads);
                    for(auto i = 0; i < num_threads; ++i)
                    {
                        threads[i] = std::jthread(evaluation_loop, i * batch_size, (i + 1) * batch_size);
                    }
                    evaluation_loop(num_threads * batch_size, num_threads * batch_size + batch_remainder);
                    if constexpr (t_j < MeshTraits<t_domain>::template getNumElements<t_i>() - 1)
                    {
                        self.template operator()<t_j + 1>(self);
                    }
                }
            }; 
            activate_elements(activate_elements);
        }

        /**
         * Domain
         * *****************************************************************************************************************************************************
         */

        template<Integer t_i>
        void
        addDomain(
            std::basic_string<Character> && domain_label,
            std::basic_string<Character> && sub_domain_label,
            std::function<Boolean(Point const &)> && function
        )
        const
        {
            for (auto const & f : this->template getDomains<t_i>())
            {
                if (f->getLabel() == std::forward<std::basic_string<Character>>(sub_domain_label))
                {
                    return;
                }
            }
            auto constexpr t_dom = MeshTraits_::template getDomain<t_i>();
            auto sub_domain = std::make_shared<FiniteDomain<t_dom, t_domain>>(std::forward<std::basic_string<Character>>(sub_domain_label));
            auto fun = [&] (auto const & finite_element)
            {
                if (std::forward<std::function<Boolean(Point const &)>>(function)(finite_element->getCoordinates()))
                {
                    finite_element->addDomain(sub_domain);
                }
            };
            this->template getDomains<t_i>().push_back(sub_domain);
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i>
        FiniteDomain<MeshTraits_::template getDomain<t_i>(), t_domain> const &
        getFiniteDomain(
            std::basic_string<Character> && domain_label
        )
        const
        {
            for (auto const & f : this->template getDomains<t_i>())
            {
                if (f->getLabel() == std::forward<std::basic_string<Character>>(domain_label))
                {
                    return * f;
                }
            }
            throw std::runtime_error("No such domain in the mesh");
        }

        template<Integer t_i>
        FiniteDomain<MeshTraits_::template getDomain<t_i>(), t_domain> &
        getFiniteDomain(
            std::basic_string<Character> && domain_label
        )
        {
            for (auto const & f : this->template getDomains<t_i>())
            {
                if (f->getLabel() == std::forward<std::basic_string<Character>>(domain_label))
                {
                    return * f;
                }
            }
            throw std::runtime_error("No such domain in the mesh");
        }

        template<Integer t_i, FieldConcept auto t_field>
        void
        setDomainDiscreteField(
            std::basic_string<Character> && domain_label
        )
        {
            auto & finite_domain = this->template getFiniteDomain<t_i>(std::forward<std::basic_string<Character>>(domain_label));
            finite_domain.template getDiscreteField<t_field>().template setDiscreteField<t_field>();
        }

        template<Integer t_i, FieldConcept auto t_field>
        void
        addDomainDiscreteFieldLoad(
            std::basic_string<Character> && domain_label,
            Integer row,
            Integer col,
            std::function<Real(Point const &, Real const &)> && function
        )
        {
            auto & finite_domain = this->template getFiniteDomain<t_i>(std::forward<std::basic_string<Character>>(domain_label));
            finite_domain.template getDiscreteField<t_field>().addLoad(row, col, std::forward<std::function<Real(Point const &, Real const &)>>(function));
        }

        template<Integer t_i, LagrangianConcept auto t_lag>
        void
        setDomainLagrangian(
            std::basic_string<Character> && domain_label
        )
        {
            auto & finite_domain = this->template getFiniteDomain<t_i>(std::forward<std::basic_string<Character>>(domain_label));
            finite_domain.template setLagrangian<t_lag>();
        }

        template<Integer t_i, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
        void
        setDomainPotential(
            std::basic_string<Character> && domain_label,
            auto const &... args
        )
        {
            auto & finite_domain = this->template getFiniteDomain<t_i>(std::forward<std::basic_string<Character>>(domain_label));
            finite_domain.template getLagrangian<t_lag>().template setPotential<t_lag, t_potential>(args...);
        }

        /**
         * Element DOF
         * *****************************************************************************************************************************************************
         */

        template<Integer t_i, FieldConcept auto t_field>
        void
        setElementDiscreteField(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template setDiscreteField<t_field>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, FieldConcept auto t_field, auto a>
        void
        addElementDiscreteFieldToLinearSystem(
            std::basic_string<Character> && domain_label,
            LinearSystem<a> & linear_system
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template getDiscreteField<t_field>().template addToLinearSystem<t_field>(linear_system);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun, 0);
        }

        /**
         * Formulation
         * *****************************************************************************************************************************************************
         */
        
        template<Integer t_i, LagrangianConcept auto t_lag>
        void
        setLagrangian(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template setLagrangian<t_lag>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
        void
        setPotential(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template getLagrangian<t_lag>().template setPotential<t_lag, t_potential>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
        void
        setPotentialStrainOperators(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template getLagrangian<t_lag>().template getPotential<t_lag, t_potential>().setStrainOperators();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
        void
        setPotentialStrains(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template getLagrangian<t_lag>().template getPotential<t_lag, t_potential>().setStrains();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }
        
        template<Integer t_i, LagrangianConcept auto t_lag, PotentialConcept auto t_behavior, Label t_label>
        void
        setPotentialMaterialProperty(
            std::basic_string<Character> && domain_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto label = std::basic_string<Character>(t_label.view());
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template getLagrangian<t_lag>().template getPotential<t_lag, t_behavior>().setMaterialProperty(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }
        
        template<Integer t_i, LagrangianConcept auto t_lag, PotentialConcept auto t_behavior, Label t_label>
        void
        setPotentialExternalVariable(
            std::basic_string<Character> && domain_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto label = std::basic_string<Character>(t_label.view());
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template getLagrangian<t_lag>().template getPotential<t_lag, t_behavior>().setExternalVariable(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }
        
        template<Integer t_i, LagrangianConcept auto t_lag, PotentialConcept auto t_behavior>
        void
        integratePotentialConstitutiveEquation(
            std::basic_string<Character> && domain_label
        )
        {
            auto res = std::atomic<Boolean>(true);
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template getLagrangian<t_lag>().template getPotential<t_lag, t_behavior>().integrateConstitutiveEquation(res);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
            if (!res)
            {
                std::cout << "Integration Failure" << std::endl;
            }
        }

        template<Integer t_i, LagrangianConcept auto t_lag>
        void
        setLagrangianJacobianMatrix(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template getLagrangian<t_lag>().template setJacobianMatrix<t_lag>();
                finite_element->template getLagrangian<t_lag>().template setResidualVector<t_lag>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }
        
        template<Integer t_i, PotentialConcept auto t_behavior>
        void
        reservePotentialBehaviorData(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template reserveBehaviorData<t_behavior>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }
        
        template<Integer t_i, PotentialConcept auto t_behavior>
        void
        recoverPotentialBehaviorData(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template recoverBehaviorData<t_behavior>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }
        
        template<Integer t_i, auto... t_args>
        void
        CALL(
            std::basic_string<Character> && domain_label
        )
        const
        {
            auto fun = [&] (auto const & finite_element)
            {
                std::cout << "finite_element->template getBandWidth<t_args...>();" << std::endl;
                std::cout << finite_element->template getBandWidth<t_args...>() << std::endl;
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun, 0);
        }
        
        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            FiniteElementSet const & finite_element_set
        )
        {
            auto print_element_inner_neighbors = [&] <ShapeConcept auto t_element, Integer t_i = 0, Integer t_j = 0> (
                    auto const & element,
                    auto & self
            )
            mutable
            {
                if constexpr (t_element != Node())
                {
                    auto const constexpr t_neighbor = ElementTraits<t_element>::template getInnerNeighbor<t_i, t_j>();
                    for (auto const & c_ : element->template getInnerNeighbors<t_i, t_j>())
                    {
                        os << "layer : " << t_i << " type : " << t_j << " <-- " << t_neighbor << " " << c_->getHash() << std::endl;
                    }
                    if constexpr (t_j < ElementTraits<t_element>::template getNumInnerNeighbors<t_i>() - 1)
                    {
                        self.template operator()<t_element, t_i, t_j + 1>(element, self);
                    }
                    else if constexpr (t_i < ElementTraits<t_element>::getNumInnerNeighbors() - 1)
                    {
                        self.template operator()<t_element, t_i + 1, 0>(element, self);
                    }
                }
            };
            auto print_element_outer_neighbors = [&] <ShapeConcept auto t_element, Integer t_i = 0, Integer t_j = 0> (
                    auto const & element,
                    auto & self
            )
            mutable
            {
                auto const constexpr t_neighbor = ElementTraits<t_element>::template getOuterNeighbor<t_domain, t_i, t_j>();
                for (auto const & c_ : element->template getOuterNeighbors<t_i, t_j>())
                {
                    if constexpr (t_element != Node() && t_i == 0)
                    {
                        os << "layer : " << t_i << " type : " << t_j << " <-> " << t_neighbor << " " << c_->getHash() << std::endl;
                    }
                    else
                    {
                        os << "layer : " << t_i << " type : " << t_j << " --> " << t_neighbor << " " << c_->getHash() << std::endl;
                    }
                }
                if constexpr (t_j < ElementTraits<t_element>::template getNumOuterNeighbors<t_domain, t_i>() - 1)
                {
                    self.template operator()<t_element, t_i, t_j + 1>(element, self);
                }
                else if constexpr (t_i < ElementTraits<t_element>::template getNumOuterNeighbors<t_domain>() - 1)
                {
                    self.template operator()<t_element, t_i + 1, 0>(element, self);
                }
            };
            auto print_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_element = MeshTraits<t_domain>::template getElement<t_i, t_j>();
                if constexpr (t_i == 0 && t_j == 0)
                {
                    os << "*** Elements : " << std::endl;
                }
                for (auto const & element : finite_element_set.template getElements<t_i, t_j>())
                {
                    os << "* Element : " << t_element << " " << element->getHash() << std::endl;
                    if (element->hasDomain())
                    {
                        os << "* Domain : " << element->getDomain()->getLabel() << std::endl;
                    }
                    print_element_inner_neighbors.template operator()<t_element>(element, print_element_inner_neighbors);
                    print_element_outer_neighbors.template operator()<t_element>(element, print_element_outer_neighbors);
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
            print_elements(print_elements);
            return os;
        }

    };

} // namespace lolita::core


#endif /* F30FA5DA_AED7_4C4E_A979_CB98F2F66287 */
