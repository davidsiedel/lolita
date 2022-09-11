#ifndef AB07387F_F0B7_4F5A_875F_1CC5F0206148
#define AB07387F_F0B7_4F5A_875F_1CC5F0206148

#include "lolita.hxx"
#include "core/000_physics_traits.hxx"
#include "core/001_geometry.hxx"
#include "core/linear_system.hxx"
#include "core/003_quadrature.hxx"
#include "core/004_finite_element.hxx"
#include "core/005_finite_element_basis.hxx"
#include "core/006_finite_element_hdg_discretization.hxx"
// #include "core/lolita_core_n_4001.hxx"

namespace lolita
{

    template<Domain t_domain>
    struct FiniteElementSet : ElementSet<FiniteElement, t_domain>
    {

        template<Integer t_i>
        void
        caller(
            std::basic_string_view<Character> domain,
            auto & fun
        )
        const
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        fun(element);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
        }

        template<ElementType t_ii>
        void
        caller(
            std::basic_string_view<Character> domain,
            auto & fun
        )
        const
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        fun(element);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
        }

        template<Integer t_i>
        void
        caller2(
            std::basic_string_view<Character> domain,
            auto & fun,
            Integer num_threads = std::thread::hardware_concurrency()
        )
        const
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                if (num_threads == 0)
                {
                    caller<t_i>(domain, fun);
                }
                else
                {
                    auto batch_size = this->template getElements<t_i, t_j>().size() / num_threads;
                    auto batch_remainder = this->template getElements<t_i, t_j>().size() % num_threads;
                    auto threads = std::vector<std::jthread>(num_threads);
                    //
                    auto doit = [&] (
                        Integer start,
                        Integer stop
                    )
                    {
                        for (auto i = start; i < stop; ++i)
                        {
                            auto const & finite_element = this->template getElements<t_i, t_j>()[i];
                            if (finite_element->isIn(domain))
                            {
                                fun(finite_element);
                            }
                        }
                    };
                    //
                    for(auto i = 0; i < num_threads; ++i)
                    {
                        threads[i] = std::jthread(doit, i * batch_size, (i + 1) * batch_size);
                    }
                    doit(num_threads * batch_size, num_threads * batch_size + batch_remainder);
                    //
                    if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                    {
                        self.template operator()<t_j + 1>(self);
                    }
                }
            }; 
            activate_elements(activate_elements);
            // caller<t_ii>(domain, fun);
        }

        template<ElementType t_ii>
        void
        caller2(
            std::basic_string_view<Character> domain,
            auto & fun,
            Integer num_threads = std::thread::hardware_concurrency()
        )
        const
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                // auto num_threads = std::thread::hardware_concurrency() == 0 ? 1 : (std::thread::hardware_concurrency());
                // auto num_threads = std::thread::hardware_concurrency();
                // auto num_threads = 2;
                if (num_threads == 0)
                {
                    caller<t_ii>(domain, fun);
                }
                else
                {
                    auto batch_size = this->template getElements<t_i, t_j>().size() / num_threads;
                    auto batch_remainder = this->template getElements<t_i, t_j>().size() % num_threads;
                    auto threads = std::vector<std::jthread>(num_threads);
                    //
                    auto doit = [&] (
                        Integer start,
                        Integer stop
                    )
                    {
                        for (auto i = start; i < stop; ++i)
                        {
                            auto const & finite_element = this->template getElements<t_i, t_j>()[i];
                            if (finite_element->isIn(domain))
                            {
                                fun(finite_element);
                            }
                        }
                    };
                    //
                    for(auto i = 0; i < num_threads; ++i)
                    {
                        threads[i] = std::jthread(doit, i * batch_size, (i + 1) * batch_size);
                    }
                    doit(num_threads * batch_size, num_threads * batch_size + batch_remainder);
                    //
                    if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                    {
                        self.template operator()<t_j + 1>(self);
                    }
                }
            }; 
            activate_elements(activate_elements);
            // caller<t_ii>(domain, fun);
        }
        
        std::unique_ptr<FiniteElementSet>
        makeFiniteElementSubSet(
            std::basic_string_view<Character> domain
        )
        const
        {
            auto sub_set = std::make_unique<FiniteElementSet>();
            auto activate_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        sub_set->template getElements<t_i, t_j>().push_back(element);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            }; 
            activate_elements(activate_elements);
            return sub_set;
        }

        template<Integer t_i>
        void
        setDof(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template setDof();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, auto... t_args, Strategy t_s>
        void
        addDof(
            std::basic_string<Character> && domain_label,
            std::unique_ptr<LinearSystem<t_s>> const & linear_system
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template addDof<t_args...>(linear_system);
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<Integer t_i, auto... t_args>
        void
        addDof(
            std::basic_string<Character> && domain_label
        )
        {
            auto fun = [&] (auto const & finite_element)
            {
                finite_element->template addDof<t_args...>();
            };
            caller2<t_i>(std::forward<std::basic_string<Character>>(domain_label), fun);
        }

        template<ElementType t_ii, auto... t_args>
        void
        setElementarySystems(
            std::basic_string<Character> domain_label,
            std::basic_string<Character> behavior_label,
            std::basic_string<Character> unknown_label,
            std::unique_ptr<System> const & system
        )
        const
        {
            auto fun = [&] (auto const & element)
            {
                element->template setElementarySystem<t_args...>(behavior_label, unknown_label, system);
            };
            caller2<t_ii>(domain_label, fun);
        }

        template<ElementType t_ii>
        Natural
        getFormulationSize(
            std::basic_string<Character> domain_label,
            std::basic_string<Character> behavior_label
        )
        const
        {
            auto value = Natural(0);
            auto mutex = std::mutex();
            auto set_value = [&] (auto input_value)
            {
                auto lock = std::scoped_lock<std::mutex>(mutex);
                value += input_value;
            };
            auto fun = [&] (auto const & element)
            {
                set_value(std::pow(element->quadrature_.at(behavior_label).residual_vector_.size(), 2));
            };
            caller2<t_ii>(domain_label, fun);
            return value;
        }

        // template<ElementType t_ii>
        // void
        // addDomain(
        //     std::basic_string_view<Character> domain_label,
        //     std::basic_string<Character> && sub_domain_label,
        //     Integer sub_domain_dim,
        //     std::function<Boolean(Point const &)> && function
        // )
        // const
        // {
        //     auto sub_domain = std::make_shared<MeshDomain>(sub_domain_dim, std::forward<std::basic_string<Character>>(sub_domain_label));
        //     auto fun = [&] (auto const & element)
        //     {
        //         if (std::forward<std::function<Boolean(Point const &)>>(function)(* element->coordinates_))
        //         {
        //             element->domains_.push_back(sub_domain);
        //         }
        //     };
        //     caller2<t_ii>(domain_label, fun);
        // }

        template<ElementType t_ii>
        void
        addDomain(
            std::basic_string_view<Character> domain_label,
            std::basic_string<Character> && sub_domain_label,
            Integer sub_domain_dim,
            std::function<Boolean(Point const &)> const & function
        )
        const
        {
            auto sub_domain = std::make_shared<MeshDomain>(sub_domain_dim, std::forward<std::basic_string<Character>>(sub_domain_label));
            auto fun = [&] (auto const & element)
            {
                if (function(element->getCurrentCentroid()))
                {
                    element->domains_.push_back(sub_domain);
                }
            };
            caller<t_ii>(domain_label, fun);
        }

        template<ElementType t_ii>
        void
        setParameter(
            std::basic_string_view<Character> domain,
            std::basic_string<Character> && parameter_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto lab = std::basic_string<Character>(parameter_label);
            auto fun = [&] (auto const & element)
            {
                element->template setParameter(lab, std::forward<std::function<Real(Point const &)>>(function));
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, auto... t_args>
        std::shared_ptr<Vector<Real>>
        setDegreeOfFreedom(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label
        )
        {
            auto dof = std::make_shared<Vector<Real>>();
            auto lab = std::basic_string<Character>(label);
            auto fun = [&] (auto const & element)
            {
                element->template setDegreeOfFreedom<t_args...>(lab, dof);
            };
            caller<t_ii>(domain, fun);
            dof->setZero();
            return dof;
        }

        template<ElementType t_ii, auto... t_args>
        void
        updateUnknown(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label,
            std::unique_ptr<System> const & system
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template updateUnknown<t_args...>(label, system);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, auto... t_args>
        void
        updateBinding(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label,
            std::unique_ptr<System> const & system
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template updateBinding<t_args...>(label, system);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, auto... t_args>
        void
        setBandWidth(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label,
            std::unique_ptr<System> const & system
        )
        {
            auto fun = [&] (auto const & element)
            {
                system->setBandWidth(element->template getBandWidth<t_args...>(label));
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii>
        // std::shared_ptr<Loading>
        std::shared_ptr<Function>
        setLoad(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label,
            Loading const & loading,
            Integer row,
            Integer col
        )
        {
            // auto load = std::make_shared<Loading>(loading);
            auto load = std::make_shared<Function>(loading, row, col);
            auto lab = std::basic_string<Character>(label);
            auto fun = [&] (auto const & element)
            {
                element->setLoad(lab, load);
            };
            caller2<t_ii>(domain, fun);
            return load;
        }

        template<ElementType t_ii>
        // std::shared_ptr<Loading>
        std::shared_ptr<Function>
        setLoad(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label,
            Loading && loading,
            Integer row,
            Integer col
        )
        {
            // auto load = std::make_shared<Loading>(std::forward<Loading>(loading));
            auto load = std::make_shared<Function>(std::forward<Loading>(loading), row, col);
            auto lab = std::basic_string<Character>(label);
            auto fun = [&] (auto const & element)
            {
                element->setLoad(lab, load);
            };
            caller2<t_ii>(domain, fun);
            return load;
        }

        template<ElementType t_ii>
        // std::shared_ptr<Loading>
        std::shared_ptr<Function>
        setConstraint(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label,
            Loading const & loading,
            Integer row,
            Integer col
        )
        {
            // auto load = std::make_shared<Loading>(loading);
            auto load = std::make_shared<Function>(loading, row, col);
            auto lab = std::basic_string<Character>(label);
            auto fun = [&] (auto const & element)
            {
                element->setConstraint(lab, load);
            };
            caller2<t_ii>(domain, fun);
            return load;
        }

        template<ElementType t_ii>
        // std::shared_ptr<Loading>
        std::shared_ptr<Function>
        setConstraint(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label,
            Loading && loading,
            Integer row,
            Integer col
        )
        {
            // auto load = std::make_shared<Loading>(std::forward<Loading>(loading));
            auto load = std::make_shared<Function>(std::forward<Loading>(loading), row, col);
            auto lab = std::basic_string<Character>(label);
            auto fun = [&] (auto const & element)
            {
                element->setConstraint(lab, load);
            };
            caller2<t_ii>(domain, fun);
            return load;
        }

        template<ElementType t_ii, Quadrature t_quadrature>
        std::shared_ptr<mgis::behaviour::Behaviour>
        setBehavior(
            std::basic_string_view<Character> domain,
            auto const &... args
        )
        {
            auto behavior = std::make_shared<mgis::behaviour::Behaviour>(mgis::behaviour::load(args...));
            auto fun = [&] (auto const & element)
            {
                element->template setBehavior<t_quadrature>(behavior);
            };
            caller2<t_ii>(domain, fun);
            return behavior;
        }

        template<ElementType t_ii, Quadrature t_quadrature>
        void
        setBehavior(
            std::basic_string_view<Character> domain,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template setBehavior<t_quadrature>(behavior);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, auto t_arg, auto t_discretization>
        void
        setStrainOperators(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> finite_element_label,
            std::basic_string_view<Character> finite_element_label2
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template setStrainOperators<t_arg, t_discretization>(finite_element_label, finite_element_label2);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, auto t_arg, auto t_discretization>
        void
        setStrainValues(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> strain_label
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template setStrainValues<t_arg, t_discretization>(behavior_label, strain_label);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii>
        Output
        integrate(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label
        )
        {
            auto output_handler = OutputHandler();
            auto fun = [&] (auto const & element)
            {
                element->integrate(behavior_label, output_handler);
            };
            caller2<t_ii>(domain, fun);
            return output_handler.getOutput();
        }

        template<ElementType t_ii>
        Real
        getStoredEnergy(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label
        )
        const
        {
            auto value = Real(0);
            auto mutex = std::mutex();
            auto set_value = [&] (auto input_value)
            {
                auto lock = std::scoped_lock<std::mutex>(mutex);
                value += input_value;
            };
            auto fun = [&] (auto const & element)
            {
                auto val = Real(0);
                for (auto const & ip : element->quadrature_.at(std::string(behavior_label)).ips_)
                {
                    val += ip.weight_ * ip.behavior_data_->s1.stored_energy;
                }
                set_value(val);
            };
            caller2<t_ii>(domain, fun);
            return value;
        }

        template<ElementType t_ii>
        Real
        getDissipatedEnergy(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label
        )
        const
        {
            auto value = Real(0);
            auto mutex = std::mutex();
            auto set_value = [&] (auto input_value)
            {
                auto lock = std::scoped_lock<std::mutex>(mutex);
                value += input_value;
            };
            auto fun = [&] (auto const & element)
            {
                auto val = Real(0);
                for (auto const & ip : element->quadrature_.at(std::string(behavior_label)).ips_)
                {
                    val += ip.weight_ * ip.behavior_data_->s1.dissipated_energy;
                }
                set_value(val);
            };
            caller2<t_ii>(domain, fun);
            return value;
        }

        template<ElementType t_ii, auto... t_args>
        Real
        getBindingIntegral(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> binding_label,
            Integer row,
            Integer col
        )
        const
        {
            auto value = Real(0);
            auto mutex = std::mutex();
            auto set_value = [&] (auto input_value)
            {
                auto lock = std::scoped_lock<std::mutex>(mutex);
                value += input_value;
            };
            auto fun = [&] (auto const & element)
            {
                set_value(element->template getBindingIntegral<t_args...>(binding_label, row, col));
            };
            caller2<t_ii>(domain, fun);
            return value;
        }

        template<ElementType t_ii>
        void
        reserveBehaviorData(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->reserveBehaviorData(behavior_label);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii>
        void
        recoverBehaviorData(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->recoverBehaviorData(behavior_label);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, auto... t_args>
        void
        reserveUnknownCoefficients(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template reserveUnknownCoefficients<t_args...>(behavior_label);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, auto... t_args>
        void
        recoverUnknownCoefficients(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template recoverUnknownCoefficients<t_args...>(behavior_label);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, auto t_arg, auto t_discretization>
        void
        setElementOperators(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template setElementOperator<t_arg, t_discretization>(label);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii>
        void
        setMaterialProperty(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->setMaterialProperty(behavior_label, material_property_label, std::forward<std::function<Real(Point const &)>>(function));
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii>
        void
        setExternalVariable(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->setExternalVariable(behavior_label, material_property_label, std::forward<std::function<Real(Point const &)>>(function));
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        void
        assembleUnknownBlock(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> degree_of_freedom_label,
            std::unique_ptr<System> const & system
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template assembleUnknownBlock<t_finite_element_method, t_discretization>(behavior_label, degree_of_freedom_label, system);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        void
        assembleUnknownVector(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> degree_of_freedom_label,
            std::unique_ptr<System> const & system
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template assembleUnknownVector<t_finite_element_method, t_discretization>(behavior_label, degree_of_freedom_label, system);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        void
        assembleBindingBlock(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> binding_label,
            std::basic_string_view<Character> unknown_label,
            std::basic_string_view<Character> constraint_label,
            std::unique_ptr<System> const & system
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template assembleBindingBlock<t_finite_element_method, t_discretization>(binding_label, unknown_label, constraint_label, system);
            };
            caller2<t_ii>(domain, fun);
        }

        template<ElementType t_ii, FiniteElementMethodConcept auto t_finite_element_method, auto t_discretization>
        void
        assembleBindingVector(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> binding_label,
            std::basic_string_view<Character> unknown_label,
            std::basic_string_view<Character> constraint_label,
            std::unique_ptr<System> const & system,
            Real const & time
        )
        {
            auto fun = [&] (auto const & element)
            {
                element->template assembleBindingVector<t_finite_element_method, t_discretization>(binding_label, unknown_label, constraint_label, system, time);
            };
            caller2<t_ii>(domain, fun);
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
                    for (auto const & domain : element->domains_)
                    {
                        os << domain->tag_ << " ";
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

        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        // Num elements
        // -----------------------------------------------------------------------------------------------------------------------------------------------------

        template<Integer... t_coordinates>
        Natural
        getNumElements()
        const
        requires(sizeof...(t_coordinates) == 2)
        {
            auto constexpr coordinates = std::array<Integer, 2>{t_coordinates...};
            return this->template getElements<coordinates[0], coordinates[1]>().size();
        }

        template<Integer... t_coordinates>
        Natural
        getNumElements()
        const
        requires(sizeof...(t_coordinates) == 1)
        {
            auto constexpr coordinates = std::array<Integer, 1>{t_coordinates...};
            auto num_elements = Natural(0);
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                num_elements += this->template getElements<coordinates[0], t_j>().size();
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<coordinates[0]>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
            return num_elements;
        }

        template<Integer... t_coordinates>
        Natural
        getNumElements()
        const
        requires(sizeof...(t_coordinates) == 0)
        {
            auto num_elements = Natural(0);
            auto activate_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                num_elements += this->template getElements<t_i, t_j>().size();
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            }; 
            activate_elements(activate_elements);
            return num_elements;
        }

        template<Integer... t_coordinates>
        Natural
        getNumElements(
            std::basic_string_view<Character> domain
        )
        const
        requires(sizeof...(t_coordinates) == 2)
        {
            auto num_elements = Natural(0);
            auto constexpr coordinates = std::array<Integer, 2>{t_coordinates...};
            for (auto const & element : this->template getElements<coordinates[0], coordinates[1]>())
            {
                if (element->isIn(domain))
                {
                    num_elements ++;
                }
            }
            return num_elements;
        }

        template<Integer... t_coordinates>
        Natural
        getNumElements(
            std::basic_string_view<Character> domain
        )
        const
        requires(sizeof...(t_coordinates) == 1)
        {
            auto constexpr coordinates = std::array<Integer, 1>{t_coordinates...};
            auto num_elements = Natural(0);
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<coordinates[0], t_j>())
                {
                    if (element->isIn(domain))
                    {
                        num_elements ++;
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<coordinates[0]>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
            return num_elements;
        }

        template<Integer... t_coordinates>
        Natural
        getNumElements(
            std::basic_string_view<Character> domain
        )
        const
        requires(sizeof...(t_coordinates) == 0)
        {
            auto num_elements = Natural(0);
            auto activate_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        num_elements ++;
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            }; 
            activate_elements(activate_elements);
            return num_elements;
        }

        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        // Num integration points
        // -----------------------------------------------------------------------------------------------------------------------------------------------------

        template<Integer... t_coordinates>
        Natural
        getNumIntegrationPoints(
            std::basic_string_view<Character> behavior_label
        )
        const
        requires(sizeof...(t_coordinates) == 2)
        {
            auto constexpr coordinates = std::array<Integer, 2>{t_coordinates...};
            auto num_integration_points = Natural(0);
            for (auto const & element : this->template getElements<coordinates[0], coordinates[1]>())
            {
                num_integration_points += element->getNumIntegrationPoints(behavior_label);
            }
            return num_integration_points;
        }

        template<Integer... t_coordinates>
        Natural
        getNumIntegrationPoints(
            std::basic_string_view<Character> behavior_label
        )
        const
        requires(sizeof...(t_coordinates) == 1)
        {
            auto constexpr coordinates = std::array<Integer, 1>{t_coordinates...};
            auto num_integration_points = Natural(0);
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<coordinates[0], t_j>())
                {
                    num_integration_points += element->getNumIntegrationPoints(behavior_label);
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<coordinates[0]>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
            return num_integration_points;
        }

        template<Integer... t_coordinates>
        Natural
        getNumIntegrationPoints(
            std::basic_string_view<Character> behavior_label
        )
        const
        requires(sizeof...(t_coordinates) == 0)
        {
            auto num_integration_points = Natural(0);
            auto activate_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    num_integration_points += element->getNumIntegrationPoints(behavior_label);
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            }; 
            activate_elements(activate_elements);
            return num_integration_points;
        }

        template<Integer... t_coordinates>
        Natural
        getNumIntegrationPoints(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label
        )
        const
        requires(sizeof...(t_coordinates) == 2)
        {
            auto constexpr coordinates = std::array<Integer, 2>{t_coordinates...};
            auto num_integration_points = Natural(0);
            for (auto const & element : this->template getElements<coordinates[0], coordinates[1]>())
            {
                if (element->isIn(domain))
                {
                    num_integration_points += element->getNumIntegrationPoints(behavior_label);
                }
            }
            return num_integration_points;
        }

        template<Integer... t_coordinates>
        Natural
        getNumIntegrationPoints(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label
        )
        const
        requires(sizeof...(t_coordinates) == 1)
        {
            auto constexpr coordinates = std::array<Integer, 1>{t_coordinates...};
            auto num_integration_points = Natural(0);
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<coordinates[0], t_j>())
                {
                    if (element->isIn(domain))
                    {
                        num_integration_points += element->getNumIntegrationPoints(behavior_label);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<coordinates[0]>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
            return num_integration_points;
        }

        template<Integer... t_coordinates>
        Natural
        getNumIntegrationPoints(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label
        )
        const
        requires(sizeof...(t_coordinates) == 0)
        {
            auto num_integration_points = Natural(0);
            auto activate_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        num_integration_points += element->getNumIntegrationPoints(behavior_label);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            }; 
            activate_elements(activate_elements);
            return num_integration_points;
        }

        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        // Output
        // -----------------------------------------------------------------------------------------------------------------------------------------------------

        template<Integer t_coordinate, auto... t_args>
        std::vector<Real>
        getNodalValues(
            std::basic_string_view<Character> unknown_label,
            Integer row,
            Integer col
        )
        const
        {
            auto nodal_values = std::vector<Real>(getNumElements<0>(), 0);
            // //
            // for (auto const & node : this->template getElements<0, 0>())
            // {
            //     auto nodal_value = 0.0;
            //     auto set_nodal_values = [&] <Integer t_j = 0> (
            //         auto & self
            //     )
            //     mutable
            //     {
            //         auto c_element_node = 0;
            //         for (auto const & outer_neighbor : node->template getInnerNeighbors<t_coordinate - 1, t_j>())
            //         {
            //             auto point = outer_neighbor->getReferenceCoordinates().col(c_element_node);
            //             outer_neighbor.template getUnknownValue<t_args...>(unknown_label, point, row, col);
            //         }
            //     };
            //     // if (node->degrees_of_freedom_.contains(std::string(unknown_label)))
            //     // {
            //     //     auto point = node->getReferenceCoordinates();
            //     //     nodal_values[node->getTag()] += node->template getUnknownValue<t_args...>(unknown_label, point, row, col);
            //     //     nodal_indices[node->getTag()] += 1;
            //     // }
            // }
            // //
            // auto nodal_indices = std::vector<Integer>(getNumElements<0>(), 0);
            auto set_nodal_values = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<t_coordinate, t_j>())
                {
                    if (element->degrees_of_freedom_.contains(std::string(unknown_label)))
                    {
                        auto c_element_node = 0;
                        for (auto const & node : element->template getInnerNeighbors<t_coordinate - 1, 0>())
                        {
                            auto point = element->getReferenceCoordinates(c_element_node);
                            nodal_values[node->getTag()] += element->template getUnknownValue<t_args...>(unknown_label, point, row, col);
                            // nodal_indices[node->getTag()] += 1;
                            c_element_node ++;
                        }
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_coordinate>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            };
            set_nodal_values(set_nodal_values);
            // for (auto const & node : this->template getElements<0, 0>())
            // {
            //     if (node->degrees_of_freedom_.contains(std::string(unknown_label)))
            //     {
            //         auto point = node->getReferenceCoordinates();
            //         nodal_values[node->getTag()] += node->template getUnknownValue<t_args...>(unknown_label, point, row, col);
            //         // nodal_indices[node->getTag()] += 1;
            //     }
            // }
            auto c_nds = 0;
            for (auto const & node : this->template getElements<0, 0>())
            {
                // nodal_values[c_nds] *= (1.0 / node->template getInnerNeighbors<t_coordinate - 1, t_j>().size());
                auto sum = 0.0;
                auto set_nodal_values2 = [&] <Integer t_j = 0> (
                    auto & self
                )
                mutable
                {
                    sum += node->template getOuterNeighbors<t_coordinate - 1, t_j>().size();
                    if constexpr (t_j < ElementTraits<Element::node()>::template getNumOuterNeighbors<t_domain, t_coordinate - 1>() - 1)
                    {
                        self.template operator()<t_j + 1>(self);
                    }
                };
                set_nodal_values2(set_nodal_values2);
                nodal_values[c_nds] *= (1.0 / sum);
                c_nds ++;
            }
            // for (auto & val : nodal_values)
            // {
            //     val /= Real(nodal_indices[c_nds]);
            //     c_nds ++;
            // }
            return nodal_values;
        }

        template<Integer t_coordinate, auto... t_args>
        std::vector<Real>
        getQuadratureValues(
            std::basic_string_view<Character> unknown_label,
            std::basic_string_view<Character> quadrature_label,
            Integer row,
            Integer col
        )
        const
        {
            auto quadrature_values = std::vector<Real>(getNumIntegrationPoints<t_coordinate>(quadrature_label), 0);
            auto tag = 0;
            auto set_quadrature_values = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : this->template getElements<t_coordinate, t_j>())
                {
                    if (element->degrees_of_freedom_.contains(std::string(unknown_label)))
                    {
                        if (element->quadrature_.contains(std::string(quadrature_label)))
                        {
                            for (auto const & integration_point : element->quadrature_.at(std::string(quadrature_label)).ips_)
                            {
                                // auto const & point = integration_point.getCurrentCoordinates();
                                auto const & point = integration_point.getReferenceCoordinates();
                                quadrature_values[tag] += element->template getUnknownValue<t_args...>(unknown_label, point, row, col);
                                tag ++;
                            }
                        }
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_coordinate>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            };
            set_quadrature_values(set_quadrature_values);
            return quadrature_values;
        }

        // void
        // setOutput(
        //     std::basic_string<Character> && file_path,
        //     auto... behavior_label
        // )
        // {
        //     auto labels = std::array<std::basic_string<Character>, sizeof...(behavior_label)>{behavior_label...};
        //     auto outfile = std::ofstream();
        //     outfile.open(std::forward<std::basic_string<Character>>(file_path));
        //     outfile << "$MeshFormat\n";
        //     outfile << "2.2 0 8\n";
        //     outfile << "$EndMeshFormat\n";
        //     outfile << "$Nodes\n";
        //     outfile << getNumElements<0>() + numerics::sum(getNumIntegrationPoints<>(behavior_label)...) << "\n";
        //     auto c_node = Natural(1);
        //     for (auto const & node : this->template getElements<0, 0>())
        //     {
        //         auto const & coordinates = node->getCurrentCoordinates();
        //         outfile << c_node << " " << coordinates(0) << " " << coordinates(1) << " " << coordinates(2) << "\n";
        //         c_node ++;
        //     }
        //     auto set_integration_nodes = [&] <Integer t_i = 0, Integer t_j = 0> (
        //         auto & self
        //     )
        //     mutable
        //     {
        //         for (auto const & element : this->template getElements<t_i, t_j>())
        //         {
        //             for (auto const & label : labels)
        //             {
        //                 if (element->quadrature_.contains(label))
        //                 {
        //                     for (auto const & integration_point : element->quadrature_.at(label).ips_)
        //                     {
        //                         auto const & coordinates = integration_point.getCurrentCoordinates();
        //                         outfile << c_node << " " << coordinates(0) << " " << coordinates(1) << " " << coordinates(2) << "\n";
        //                         c_node ++;
        //                     }
        //                 }
        //             }
        //         }
        //         if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
        //         {
        //             self.template operator()<t_i, t_j + 1>(self);
        //         }
        //         else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
        //         {
        //             self.template operator()<t_i + 1, 0>(self);
        //         }
        //     };
        //     set_integration_nodes(set_integration_nodes);
        //     outfile << "$EndNodes\n";
        //     outfile << "$Elements\n";
        //     outfile << getNumElements<>() + numerics::sum(getNumIntegrationPoints<>(behavior_label)...) << "\n";
        //     auto c_element = Natural(1);
        //     for (auto const & node : this->template getElements<0, 0>())
        //     {
        //         outfile << c_element << " 15 2 0 0 " << c_element << "\n";
        //         c_element ++;
        //     }
        //     auto set_integration_elements = [&] <Integer t_i = 0, Integer t_j = 0> (
        //         auto & self
        //     )
        //     mutable
        //     {
        //         for (auto const & element : this->template getElements<t_i, t_j>())
        //         {
        //             for (auto const & label : labels)
        //             {
        //                 if (element->quadrature_.contains(label))
        //                 {
        //                     for (auto const & integration_point : element->quadrature_.at(label).ips_)
        //                     {
        //                         outfile << c_element << " 15 2 1 1 " << c_element << "\n";
        //                         c_element ++;
        //                     }
        //                 }
        //             }
        //         }
        //         if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
        //         {
        //             self.template operator()<t_i, t_j + 1>(self);
        //         }
        //         else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
        //         {
        //             self.template operator()<t_i + 1, 0>(self);
        //         }
        //     };
        //     set_integration_elements(set_integration_elements);
        //     auto set_elements = [&] <Integer t_i = 1, Integer t_j = 0> (
        //         auto & self
        //     )
        //     mutable
        //     {
        //         for (auto const & element : this->template getElements<t_i, t_j>())
        //         {
        //             outfile << c_element << " " << tag << " 2 0 0 " << c_element;
        //             for (auto const & node : element->template getInnerNeighbors<t_i - 1, 0>())
        //             {
        //                 outfile << " " << node->getTag() + 1;
        //             }
        //             outfile << "\n";
        //             c_element ++;
        //         }
        //         if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
        //         {
        //             self.template operator()<t_i, t_j + 1>(self);
        //         }
        //         else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
        //         {
        //             self.template operator()<t_i + 1, 0>(self);
        //         }
        //     };
        //     set_elements(set_elements);
        //     outfile << "$EndElements\n";
        // }

        // void
        // addOutput(
        //     std::basic_string<Character> && file_path,
        //     Integer time_step_index,
        //     Real time_step_value,
        //     auto... behavior_label
        // )
        // {
        //     // if (!std::filesystem::exists(std::forward<std::basic_string<Character>>(file_path)))
        //     // {
        //     //     throw std::runtime_error("File does not exist");
        //     // }
        //     auto labels = std::array<std::basic_string<Character>, sizeof...(behavior_label)>{behavior_label...};
        //     auto outfile = std::ofstream();
        //     auto c_element = Natural();
        //     outfile.open(std::forward<std::basic_string<Character>>(file_path), std::ios_base::app);
        //     // writing strain
        //     outfile << "$NodeData\n";
        //     outfile << "1\n";
        //     outfile << "\"" << labels[0] << "Strain\"\n";
        //     outfile << "1\n";
        //     outfile << time_step_value << "\n";
        //     outfile << "3\n";
        //     outfile << time_step_index << "\n";
        //     outfile << 4 << "\n"; // size of gradients
        //     outfile << getNumIntegrationPoints<>(labels[0]) << "\n"; // number of quad pts
        //     c_element = getNumElements<0>() + 1;
        //     auto set_strain = [&] <Integer t_i = 0, Integer t_j = 0> (
        //         auto & self
        //     )
        //     mutable
        //     {
        //         for (auto const & element : this->template getElements<t_i, t_j>())
        //         {
        //             if (element->quadrature_.contains(labels[0]))
        //             {
        //                 for (auto const & integration_point : element->quadrature_.at(labels[0]).ips_)
        //                 {
        //                     outfile << c_element;
        //                     for (auto val : integration_point.behavior_data_->s1.gradients)
        //                     {
        //                         outfile << " " << val;
        //                     }
        //                     outfile << "\n";
        //                     c_element ++;
        //                 }
        //             }
        //         }
        //         if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
        //         {
        //             self.template operator()<t_i, t_j + 1>(self);
        //         }
        //         else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
        //         {
        //             self.template operator()<t_i + 1, 0>(self);
        //         }
        //     };
        //     set_strain(set_strain);
        //     outfile << "$EndNodeData\n";
        //     // writing stress
        //     outfile << "$NodeData\n";
        //     outfile << "1\n";
        //     outfile << "\"" << labels[0] << "Stress\"";
        //     outfile << "1\n";
        //     outfile << time_step_value << "\n";
        //     outfile << "3\n";
        //     outfile << time_step_index << "\n";
        //     outfile << 4 << "\n"; // size of gradients
        //     outfile << getNumIntegrationPoints<>(labels[0]) << "\n"; // number of quad pts
        //     c_element = getNumElements<0>() + 1;
        //     auto set_stress = [&] <Integer t_i = 0, Integer t_j = 0> (
        //         auto & self
        //     )
        //     mutable
        //     {
        //         for (auto const & element : this->template getElements<t_i, t_j>())
        //         {
        //             if (element->quadrature_.contains(labels[0]))
        //             {
        //                 for (auto const & integration_point : element->quadrature_.at(labels[0]).ips_)
        //                 {
        //                     outfile << c_element;
        //                     for (auto val : integration_point.behavior_data_->s1.thermodynamic_forces)
        //                     {
        //                         outfile << " " << val;
        //                     }
        //                     outfile << "\n";
        //                     c_element ++;
        //                 }
        //             }
        //         }
        //         if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
        //         {
        //             self.template operator()<t_i, t_j + 1>(self);
        //         }
        //         else if constexpr (t_i < DomainTraits<t_domain>::template getNumElements<>() - 1)
        //         {
        //             self.template operator()<t_i + 1, 0>(self);
        //         }
        //     };
        //     set_stress(set_stress);
        //     outfile << "$EndNodeData\n";
        // }

    };
    
} // namespace lolita

#endif /* AB07387F_F0B7_4F5A_875F_1CC5F0206148 */
