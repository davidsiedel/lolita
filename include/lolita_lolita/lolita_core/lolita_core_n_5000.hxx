#ifndef B1612904_E704_486C_8F8D_294C54B26BC6
#define B1612904_E704_486C_8F8D_294C54B26BC6

// #include "lolita_lolita/lolita_core/lolita.hxx"
// #include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
// #include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"
// #include "lolita_lolita/lolita_core/lolita_core_n_2000.hxx"
// #include "lolita_lolita/lolita_core/lolita_core_n_3000.hxx"
// #include "lolita_lolita/lolita_core/lolita_core_n_4000.hxx"

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_2000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_3000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4001.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4002.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4003.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4004.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4005.hxx"

namespace lolita
{

    template<Domain t_domain>
    struct FiniteElementSet : ElementSet<FiniteElementHolder, t_domain>
    {
        
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

        template<ElementType t_ii>
        void
        setParameter(
            std::basic_string_view<Character> domain,
            std::basic_string<Character> && parameter_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto lab = std::basic_string<Character>(parameter_label);
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
                        element->template setParameter(lab, std::forward<std::function<Real(Point const &)>>(function));
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
        std::shared_ptr<Vector<Real>>
        setDegreeOfFreedom(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label
        )
        {
            auto dof = std::make_shared<Vector<Real>>();
            auto lab = std::basic_string<Character>(label);
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->template setDegreeOfFreedom<t_field, t_basis>(lab, dof);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
            dof->setZero();
            return dof;
        }

        template<ElementType t_ii, Field t_field, Basis t_basis>
        void
        addDegreeOfFreedomCoefficients(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label,
            lolita::algebra::View<Vector<Real> const> const & vector
        )
        {
            auto lab = std::basic_string<Character>(label);
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->template addDegreeOfFreedomCoefficients<t_field, t_basis>(lab, vector);
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
        setDegreeOfFreedomCoefficients(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label,
            lolita::algebra::View<Vector<Real> const> const & vector
        )
        {
            auto lab = std::basic_string<Character>(label);
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->template setDegreeOfFreedomCoefficients<t_field, t_basis>(lab, vector);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
        }

        template<ElementType t_ii, auto... t_args>
        void
        updateUnknown(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label,
            std::unique_ptr<System> const & system
        )
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
                        element->template updateUnknown<t_args...>(label, system);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
        }

        template<ElementType t_ii, auto... t_args>
        void
        updateBinding(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label,
            std::unique_ptr<System> const & system
        )
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
                        element->template updateBinding<t_args...>(label, system);
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
                        // element->setLoad(lab, load, row, col);
                        element->setLoad(lab, load);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
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
                        // element->setLoad(lab, load, row, col);
                        element->setLoad(lab, load);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
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
                        // element->setConstraint(lab, load, row, col);
                        element->setConstraint(lab, load);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
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
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        // element->setConstraint(lab, load, row, col);
                        element->setConstraint(lab, load);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
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
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->template setBehavior<t_quadrature>(behavior);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
            return behavior;
        }

        template<ElementType t_ii, Quadrature t_quadrature>
        void
        setBehavior(
            std::basic_string_view<Character> domain,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->template setBehavior<t_quadrature>(behavior);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
        }

        template<ElementType t_ii, auto t_arg, auto t_discretization>
        void
        setStrainOperators(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> finite_element_label,
            std::basic_string_view<Character> finite_element_label2
        )
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->template setStrainOperators<t_arg, t_discretization>(finite_element_label, finite_element_label2);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
        }

        template<ElementType t_ii, auto t_arg, auto t_discretization>
        void
        setStrainValues(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> strain_label
        )
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->template setStrainValues<t_arg, t_discretization>(behavior_label, strain_label);
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
        integrate(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label
        )
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->integrate(behavior_label);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
        }

        template<ElementType t_ii, auto t_arg, auto t_discretization>
        void
        setElementOperators(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> label
        )
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->template setElementOperator<t_arg, t_discretization>(label);
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
        setMaterialProperty(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->setMaterialProperty(behavior_label, material_property_label, std::forward<std::function<Real(Point const &)>>(function));
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
        setExternalVariable(
            std::basic_string_view<Character> domain,
            std::basic_string_view<Character> behavior_label,
            std::basic_string_view<Character> material_property_label,
            std::function<Real(Point const &)> && function
        )
        {
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->setExternalVariable(behavior_label, material_property_label, std::forward<std::function<Real(Point const &)>>(function));
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
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
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->template assembleUnknownBlock<t_finite_element_method, t_discretization>(behavior_label, degree_of_freedom_label, system);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
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
            auto activate_elements = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_i = t_ii.getDim();
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto const & element : this->template getElements<t_i, t_j>())
                {
                    if (element->isIn(domain))
                    {
                        element->template assembleBindingBlock<t_finite_element_method, t_discretization>(binding_label, unknown_label, constraint_label, system);
                    }
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }
            }; 
            activate_elements(activate_elements);
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
            auto print_element_outer_neighbors = [&] <Element t_element, Integer t_i = 0, Integer t_j = 0> (
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
            auto nodal_indices = std::vector<Integer>(getNumElements<0>(), 0);
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
                            nodal_indices[node->getTag()] += 1;
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
            for (auto const & node : this->template getElements<0, 0>())
            {
                if (node->degrees_of_freedom_.contains(std::string(unknown_label)))
                {
                    auto point = node->getReferenceCoordinates();
                    nodal_values[node->getTag()] += node->template getUnknownValue<t_args...>(unknown_label, point, row, col);
                    nodal_indices[node->getTag()] += 1;
                }
                nodal_values[node->getTag()] /= nodal_indices[node->getTag()];
            }
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
            auto set_quadrature_values = [&] <Integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto tag = 0;
                for (auto const & element : this->template getElements<t_coordinate, t_j>())
                {
                    if (element->degrees_of_freedom_.contains(std::string(unknown_label)))
                    {
                        if (element->quadrature_.contains(std::string(quadrature_label)))
                        {
                            for (auto const & integration_point : element->quadrature_.at(std::string(quadrature_label)).ips_)
                            {
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

#endif /* B1612904_E704_486C_8F8D_294C54B26BC6 */
