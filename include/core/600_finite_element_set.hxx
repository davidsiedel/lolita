#ifndef AB07387F_F0B7_4F5A_875F_1CC5F0206148
#define AB07387F_F0B7_4F5A_875F_1CC5F0206148

#include "lolita.hxx"
#include "core/000_physics_traits.hxx"
#include "core/100_geometry.hxx"
#include "core/linear_system.hxx"
#include "core/200_quadrature.hxx"
#include "core/201_finite_element_dof.hxx"
#include "core/202_finite_element_frm.hxx"
#include "core/300_finite_element.hxx"
#include "core/400_finite_element_basis.hxx"
#include "core/500_finite_element_hdg_discretization.hxx"
#include "core/501_finite_element_hdg_discretization.hxx"
// #include "core/lolita_core_n_4001.hxx"

namespace lolita
{

    template<Domain t_domain>
    struct FiniteElementSet : ElementSet<FiniteElement, t_domain>, DomainSet<FiniteDomain, t_domain>
    {

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
                    evaluation_loop(0, this->template getElements<t_i, t_j>().size() - 1);
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
                    if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
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
            auto sub_domain = std::make_shared<FiniteDomain<t_i, t_domain>>(std::forward<std::basic_string<Character>>(sub_domain_label));
            auto fun = [&] (auto const & finite_element)
            {
                if (std::forward<std::function<Boolean(Point const &)>>(function)(finite_element->getCoordinates()))
                {
                    finite_element->addDomain(sub_domain);
                }
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i>
        FiniteDomain<t_i, t_domain> const &
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
        FiniteDomain<t_i, t_domain> &
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
        addDomainDiscreteField(
            std::basic_string<Character> && domain_label
        )
        {
            auto & finite_domain = this->template getFiniteDomain<t_i>(std::forward<std::basic_string<Character>>(domain_label));
            finite_domain.template addDiscreteField<t_field>();
        }

        template<Integer t_i, FieldConcept auto t_field, Integer t_size>
        void
        addDomainDiscreteFieldDegreeOfFreedom(
            std::basic_string<Character> && domain_label,
            auto const &... args
        )
        {
            auto & finite_domain = this->template getFiniteDomain<t_i>(std::forward<std::basic_string<Character>>(domain_label));
            finite_domain.template addDiscreteFieldDegreeOfFreedom<t_field, t_size>(args...);
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
            finite_domain.template addDiscreteFieldLoad<t_field>(row, col, std::forward<std::function<Real(Point const &, Real const &)>>(function));
        }

        /**
         * Element DOF
         * *****************************************************************************************************************************************************
         */

        template<Integer t_i, FieldConcept auto t_field>
        void
        addElementDiscreteField(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template addDiscreteField<t_field>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        // template<Integer t_i, FieldConcept auto t_field, auto t_arg>
        // void
        // addElementDiscreteFieldDegreeOfFreedom(
        //     std::basic_string<Character> && domain_label,
        //     auto const &... args
        // )
        // {
        //     auto fun = [&] (auto const & finite_element)
        //     {
        //         finite_element->template addDiscreteFieldDegreeOfFreedom<t_field, t_arg>(args...);
        //     };
        //     caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        // }

        template<Integer t_i, FieldConcept auto t_field>
        void
        addElementDiscreteFieldDegreeOfFreedom(
            std::basic_string<Character> && domain_label,
            auto const &... args
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template addDiscreteFieldDegreeOfFreedom<t_field>(args...);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        // template<Integer t_i, FieldConcept auto t_field, auto t_arg, Label t_label>
        // void
        // addElementDiscreteFieldOperator(
        //     std::basic_string<Character> && domain_label
        // )
        // {
        //     auto fun = [&] (auto const & finite_element)
        //     {
        //         finite_element->template addDiscreteFieldOperator<t_field, t_arg, t_label>();
        //     };
        //     caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        // }

        template<Integer t_i, FieldConcept auto t_field, Label t_label>
        void
        addElementDiscreteFieldOperator(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template addDiscreteFieldOperator<t_field, t_label>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        // template<Integer t_i, FieldConcept auto t_field, auto t_arg>
        // void
        // upgradeElementDiscreteFieldDegreeOfFreedom(
        //     std::basic_string<Character> && domain_label,
        //     auto const &... args
        // )
        // {
        //     auto fun = [&] (auto const & finite_element)
        //     {
        //         finite_element->template upgradeDiscreteFieldDegreeOfFreedom<t_field, t_arg>(args...);
        //     };
        //     caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        // }

        template<Integer t_i, FieldConcept auto t_field>
        void
        upgradeElementDiscreteFieldDegreeOfFreedom(
            std::basic_string<Character> && domain_label,
            auto const &... args
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template upgradeDiscreteFieldDegreeOfFreedom<t_field>(args...);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, FieldConcept auto t_field>
        void
        recoverElementDiscreteFieldDegreeOfFreedom(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template recoverDiscreteFieldDegreeOfFreedom<t_field>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, FieldConcept auto t_field>
        void
        reserveElementDiscreteFieldDegreeOfFreedom(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template reserveDiscreteFieldDegreeOfFreedom<t_field>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        /**
         * Formulation
         * *****************************************************************************************************************************************************
         */

        template<Integer t_i, PotentialConcept auto t_behavior>
        void
        addFormulation(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template addFormulation<t_behavior>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, PotentialConcept auto t_behavior>
        void
        addFormulation(
            std::basic_string<Character> && domain_label,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template addFormulation<t_behavior>(behavior);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, PotentialConcept auto t_behavior>
        void
        addFormulation(
            std::basic_string<Character> && domain_label,
            auto const &... args
        )
        {
            auto behavior = std::make_shared<mgis::behaviour::Behaviour>(mgis::behaviour::load(args...));
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template addFormulation<t_behavior>(behavior);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, PotentialConcept auto t_behavior, Quadrature t_quadrature>
        void
        addFormulation(
            std::basic_string<Character> && domain_label,
            auto const &... args
        )
        {
            auto behavior = std::make_shared<mgis::behaviour::Behaviour>(mgis::behaviour::load(args...));
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template addFormulation<t_behavior, t_quadrature>(behavior);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        // template<Integer t_i, PotentialConcept auto t_behavior, MappingConcept auto t_strain, DiscretizationConcept auto t_discretization>
        // void
        // addFormulationStrainOperator(
        //     std::basic_string<Character> && domain_label
        // )
        // {
        //     auto fun = [&] (auto const & finite_element)
        //     {
        //         finite_element->template addFormulationStrainOperator<t_behavior, t_strain, t_discretization>();
        //     };
        //     caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        // }

        template<Integer t_i, PotentialConcept auto t_behavior, MappingConcept auto t_strain>
        void
        addFormulationStrainOperator(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template addFormulationStrainOperator<t_behavior, t_strain>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, PotentialConcept auto... t_potential>
        void
        setStrainOperators(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template setElementStrainOperators<t_potential...>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, PotentialConcept auto... t_potential>
        void
        setStrains(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template setElementStrains<t_potential...>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }
        
        template<Integer t_i, PotentialConcept auto... t_potential>
        Real
        getInternalEnergy(
            std::basic_string<Character> && domain_label
        )
        const
        {
            auto energy = Real(0);
            auto fun = [&] (auto const & finite_element)
            {
                energy += finite_element->template getElementInternalEnergy<t_potential...>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun, 0);
            return energy;
        }
        
        template<Integer t_i, PotentialConcept auto... t_potential>
        Real
        getExternalEnergy(
            std::basic_string<Character> && domain_label,
            Real const & time
        )
        const
        {
            auto energy = Real(0);
            auto fun = [&] (auto const & finite_element)
            {
                energy += finite_element->template getElementExternalEnergy<t_potential...>(time);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun, 0);
            return energy;
        }
        
        template<Integer t_i, FieldConcept auto... t_fields>
        Real
        getExternalEnergy(
            std::basic_string<Character> && domain_label,
            Real const & time
        )
        const
        {
            auto energy = Real(0);
            auto fun = [&] (auto const & finite_element)
            {
                energy += finite_element->template getElementExternalEnergy<t_fields...>(time);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun, 0);
            return energy;
        }
        
        template<Integer t_i, PotentialConcept auto... t_potential>
        Real
        getResidualEnergy(
            std::basic_string<Character> && domain_label,
            Real const & time
        )
        const
        {
            auto energy = Real(0);
            auto fun = [&] (auto const & finite_element)
            {
                energy += finite_element->template getElementResidualEnergy<t_potential...>(time);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun, 0);
            return energy;
        }

        // template<Integer t_i, PotentialConcept auto t_behavior, MappingConcept auto t_strain, DiscretizationConcept auto t_discretization>
        // void
        // setFormulationStrain(
        //     std::basic_string<Character> && domain_label
        // )
        // {
        //     auto fun = [&] (auto const & finite_element)
        //     {
        //         finite_element->template setFormulationStrain<t_behavior, t_strain, t_discretization>();
        //     };
        //     caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        // }

        template<Integer t_i, PotentialConcept auto t_behavior, MappingConcept auto t_strain>
        void
        setFormulationStrain(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template setFormulationStrain<t_behavior, t_strain>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }
        
        template<Integer t_i, PotentialConcept auto t_behavior, Label t_label>
        void
        setFormulationMaterialProperty(
            std::basic_string<Character> && domain_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto label = std::basic_string<Character>(t_label.view());
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template setMaterialProperty<t_behavior>(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }
        
        template<Integer t_i, PotentialConcept auto t_behavior, Label t_label>
        void
        setFormulationExternalVariable(
            std::basic_string<Character> && domain_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto label = std::basic_string<Character>(t_label.view());
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template setExternalVariable<t_behavior>(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }
        
        template<Integer t_i, PotentialConcept auto t_behavior>
        void
        integrateFormulationConstitutiveEquation(
            std::basic_string<Character> && domain_label
        )
        {
            auto res = std::atomic<Boolean>(true);
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template integrateConstitutiveEquation<t_behavior>(res);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
            if (!res)
            {
                std::cout << "Integration Failure" << std::endl;
            }
        }
        
        template<Integer t_i, PotentialConcept auto t_behavior>
        void
        reserveFormulationBehaviorData(
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
        recoverFormulationBehaviorData(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template recoverBehaviorData<t_behavior>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }
        
        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            FiniteElementSet const & finite_element_set
        )
        {
            auto print_element_inner_neighbors = [&] <Element t_element, Integer t_i = 0, Integer t_j = 0> (
                    auto const & element,
                    auto & self
            )
            mutable
            {
                if constexpr (!t_element.isNode())
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
            auto print_element_outer_neighbors = [&] <Element t_element, Integer t_i = 0, Integer t_j = 0> (
                    auto const & element,
                    auto & self
            )
            mutable
            {
                auto const constexpr t_neighbor = ElementTraits<t_element>::template getOuterNeighbor<t_domain, t_i, t_j>();
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
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                if constexpr (t_i == 0 && t_j == 0)
                {
                    os << "*** Elements : " << std::endl;
                }
                for (auto const & element : finite_element_set.template getElements<t_i, t_j>())
                {
                    os << "* Element : " << t_element << " " << element->getHash() << std::endl;
                    os << "* Domains : ";
                    for (auto const & domain : element->getDomains())
                    {
                        os << domain->getLabel() << " ";
                    }
                    os << std::endl;
                    print_element_inner_neighbors.template operator()<t_element>(element, print_element_inner_neighbors);
                    print_element_outer_neighbors.template operator()<t_element>(element, print_element_outer_neighbors);
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
    
} // namespace lolita

#endif /* AB07387F_F0B7_4F5A_875F_1CC5F0206148 */
