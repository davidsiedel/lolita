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

        void
        addOutput(
            std::basic_string<Character> && file_path
        )
        {
            auto outfile = std::ofstream();
            outfile.open(std::forward<std::basic_string<Character>>(file_path), std::ios_base::app);
            std::cout << "ICI !!!\n";
            // outfile << "\n Data";
        }

    // def create_output(self, res_folder_path: str):
    //     res_file_path = os.path.join(res_folder_path, "output.msh")
    //     with open(res_file_path, "w") as res_output_file:
    //         res_output_file.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n")
    //         # res_output_file.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n")
    //         nnodes = self.mesh.number_of_vertices_in_mesh + self.mesh.number_of_cell_quadrature_points_in_mesh
    //         res_output_file.write("{}\n".format(nnodes))
    //         # res_output_file.write("1 {} 1 {}\n".format(nnodes, nnodes))
    //         for v_count in range(self.mesh.number_of_vertices_in_mesh):
    //             vertex_fill = np.zeros((3,), dtype=real)
    //             vertex_fill[:len(self.mesh.vertices[:,v_count])] = self.mesh.vertices[:,v_count]
    //             res_output_file.write("{} {} {} {}\n".format(v_count + 1, vertex_fill[0], vertex_fill[1], vertex_fill[2]))
    //         q_count = self.mesh.number_of_vertices_in_mesh
    //         for element in self.elements:
    //             cell_quadrature_size = element.cell.get_quadrature_size(
    //                 # element.finite_element.construction_integration_order
    //                 element.finite_element.computation_integration_order
    //             )
    //             cell_quadrature_points = element.cell.get_quadrature_points(
    //                 # element.finite_element.construction_integration_order
    //                 element.finite_element.computation_integration_order
    //             )
    //             for qc in range(cell_quadrature_size):
    //                 x_q_c = cell_quadrature_points[:, qc]
    //                 qp_fill = np.zeros((3,), dtype=real)
    //                 qp_fill[:len(x_q_c)] = x_q_c
    //                 res_output_file.write("{} {} {} {}\n".format(q_count + 1, qp_fill[0], qp_fill[1], qp_fill[2]))
    //                 q_count += 1
    //         res_output_file.write("$EndNodes\n")
    //         res_output_file.write("$Elements\n")
    //         n_elems = nnodes + len(self.mesh.faces_vertices_connectivity) + len(self.mesh.cells_vertices_connectivity)
    //         res_output_file.write("{}\n".format(n_elems))
    //         elem_count = 1
    //         for v_count in range(self.mesh.number_of_vertices_in_mesh):
    //             res_output_file.write("{} 15 2 0 0 {}\n".format(elem_count, elem_count))
    //             elem_count += 1
    //         # q_count = self.mesh.number_of_vertices_in_mesh
    //         for element in self.elements:
    //             cell_quadrature_size = element.cell.get_quadrature_size(
    //                 # element.finite_element.construction_integration_order
    //                 element.finite_element.computation_integration_order
    //             )
    //             cell_quadrature_points = element.cell.get_quadrature_points(
    //                 # element.finite_element.construction_integration_order
    //                 element.finite_element.computation_integration_order
    //             )
    //             for qc in range(cell_quadrature_size):
    //                 x_q_c = cell_quadrature_points[:, qc]
    //                 qp_fill = np.zeros((3,), dtype=real)
    //                 qp_fill[:len(x_q_c)] = x_q_c
    //                 res_output_file.write("{} 15 2 1 1 {}\n".format(elem_count, elem_count))
    //                 elem_count += 1
    //                 # res_output_file.write("{} {} {} {}\n".format(q_count + 1, qp_fill[0], qp_fill[1], qp_fill[2]))
    //                 # q_count += 1
    //         for face_connectivity, face_shape in zip(self.mesh.faces_vertices_connectivity, self.mesh.faces_shape_types):
    //             elem_tag = get_element_tag(face_shape)
    //             res_output_file.write("{} {} 2 0 0 ".format(elem_count, elem_tag))
    //             for i_loc, coord in enumerate(face_connectivity):
    //                 if i_loc != len(face_connectivity) - 1:
    //                     res_output_file.write("{} ".format(coord + 1))
    //                 else:
    //                     res_output_file.write("{}\n".format(coord + 1))
    //             elem_count += 1
    //         for cell_connectivity, cell_shape in zip(self.mesh.cells_vertices_connectivity, self.mesh.cells_shape_types):
    //             elem_tag = get_element_tag(cell_shape)
    //             res_output_file.write("{} {} 2 0 0 ".format(elem_count, elem_tag))
    //             for i_loc, coord in enumerate(cell_connectivity):
    //                 if i_loc != len(cell_connectivity) - 1:
    //                     res_output_file.write("{} ".format(coord + 1))
    //                 else:
    //                     res_output_file.write("{}\n".format(coord + 1))
    //             elem_count += 1
    //         res_output_file.write("$EndElements\n")

    };
    
} // namespace lolita

#endif /* B1612904_E704_486C_8F8D_294C54B26BC6 */
