#ifndef CC79BDC5_49DB_4A81_8A93_18ABD6551AF1
#define CC79BDC5_49DB_4A81_8A93_18ABD6551AF1

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_00.hxx"
#include "lolita/lolita_core_n_0.hxx"
#include "lolita/lolita_core_n_001.hxx"

namespace lolita2::geometry
{

    template<template<Element, Domain, auto...> typename t_FiniteElement, Element t_element, Domain t_domain, auto... t_args>
    struct FiniteElementGeometry;

    template<Element t_element, Domain t_domain, auto t_arg>
    struct FiniteElement;
    
    template<Element t_element, Quadrature t_quadrature>
    struct ElementQuadratureRuleTraits;
    
    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.isNode() || !t_element.isNode())
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {
        
        lolita::integer static constexpr dim_ = 1;
        
        lolita::integer static constexpr size_ = 1;

        static constexpr
        lolita::integer
        size()
        {
            return 1;
        }

        std::array<std::array<lolita::real, 3>, dim_> static constexpr reference_points_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };
        
        std::array<lolita::real, dim_> static constexpr reference_weights_ = {
                +1.0000000000000000
        };

    };
        
    template<Basis t_basis>
    struct FiniteElementBasisTraits;

    template<Basis t_basis>
    requires(t_basis.isMonomial())
    struct FiniteElementBasisTraits<t_basis>
    {

        template<lolita::integer t_dim>
        static constexpr
        lolita::integer
        size()
        {
            return lolita::numerics::binomial(t_dim + t_basis.ord_, t_dim);
        }

        template<Element t_element>
        static constexpr
        lolita::integer
        size()
        {
            return lolita::numerics::binomial(t_element.dim_ + t_basis.ord_, t_element.dim_);
        }
        
        template<template<Element, Domain, auto...> typename t_FiniteElement, Element t_element, Domain t_domain, auto... t_args>
        struct Implementation : FiniteElementGeometry<t_FiniteElement, t_element, t_domain, t_args...>
        {

        private:

            lolita::integer static constexpr size_ = size<t_element>();
        
            static constexpr
            std::array<std::array<lolita::integer, 3>, size_>
            getExponents()
            {
                auto exponents = std::array<std::array<lolita::integer, 3>, size_>();
                auto row = lolita::integer(0);
                if constexpr (t_element.dim_ == 0)
                {
                    exponents[row][0] = 0;
                    exponents[row][1] = 0;
                    exponents[row][2] = 0;
                }
                else if constexpr (t_element.dim_ == 1)
                {
                    for (lolita::integer i = 0; i < t_basis.ord_ + 1; ++i)
                    {
                        exponents[row][0] = i;
                        exponents[row][1] = 0;
                        exponents[row][2] = 0;
                        row += 1;
                    }
                }
                else if constexpr (t_element.dim_ == 2)
                {
                    for (lolita::integer i = 0; i < t_basis.ord_ + 1; ++i)
                    {
                        for (lolita::integer j = 0; j < i + 1; ++j)
                        {
                            exponents[row][0] = i - j;
                            exponents[row][1] = j;
                            exponents[row][2] = 0;
                            row += 1;
                        }
                    }
                }
                else if constexpr (t_element.dim_ == 3)
                {
                    for (lolita::integer i = 0; i < t_basis.ord_ + 1; ++i)
                    {
                        for (lolita::integer j = 0; j < i + 1; ++j)
                        {
                            for (lolita::integer k = 0; k < i + 1; ++k)
                            {
                                if (j + k < i + 1)
                                {
                                    exponents[row][0] = i - (j + k);
                                    exponents[row][1] = k;
                                    exponents[row][2] = j;
                                    row += 1;
                                }
                            }
                        }
                    }
                }
                return exponents;
            }
            
            std::array<std::array<lolita::integer, 3>, size_> static constexpr exponents_ = getExponents();

        public:
        
            lolita::matrix::Vector<lolita::real, size_>
            getBasisEvaluation(
                Point const & point
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<lolita::real, size_>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (lolita::integer i = 0; i < size_; ++i)
                {
                    auto value = lolita::real(1);
                    for (lolita::integer j = 0; j < t_element.dim_; ++j)
                    {
                        auto dist = this->getRiemannianDistance(centroid, point, j);
                        value *= std::pow(2.0 * dist / diameters(j), exponents_[i][j]);
                    }
                    basis_vector_values(i) = value;
                }
                return basis_vector_values;
            }
            
            lolita::matrix::Vector<lolita::real, size_>
            getBasisDerivative(
                Point const & point,
                lolita::integer derivative_direction
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<lolita::real, size_>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (lolita::integer i = 0; i < size_; ++i)
                {
                    auto value = lolita::real(1);
                    for (lolita::integer j = 0; j < t_element.dim_; ++j)
                    {
                        if (j != derivative_direction)
                        {
                            auto dist = this->getRiemannianDistance(centroid, point, j);
                            value *= std::pow(2.0 * (dist) / diameters(j), exponents_[i][j]);
                        }
                        else
                        {
                            if (exponents_[i][j] > 0)
                            {
                                auto c = 2.0 * exponents_[i][j] / diameters(j);
                                auto dist = this->getRiemannianDistance(centroid, point, j);
                                value *= c * std::pow(2.0 * (dist) / diameters(j), exponents_[i][j] - 1);
                            }
                            else
                            {
                                value *= 0.0;
                            }
                        }
                    }
                    basis_vector_values(i) = value;
                }
                return basis_vector_values;
            }

        };

    };

    template<Element t_element, Domain t_domain>
    struct FiniteElementDegreeOfFreedom
    {

        template<Field t_field, Basis t_basis>
        static constexpr
        lolita::integer
        getSize()
        {
            auto constexpr t_field_size = FieldTraits<t_field>::template size<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template size<t_element>();
            return t_field_size * t_basis_size;
        }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, getSize<t_field, t_basis>()> const>
        getCoefficients()
        const
        {
            auto const & data = degree_of_freedom_->coefficients_.data() + index_;
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, getSize<t_field, t_basis>()> const>(data);
        }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, getSize<t_field, t_basis>()>>
        getCoefficients()
        {
            auto const & data = degree_of_freedom_->coefficients_.data() + index_;
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, getSize<t_field, t_basis>()>>(data);
        }

        lolita::integer index_;

        std::shared_ptr<DegreeOfFreedom> degree_of_freedom_;

    };

    template<Element t_element, Domain t_domain>
    struct FiniteElementUnknown
    {

        template<Field t_field, Basis t_basis>
        static constexpr
        lolita::integer
        getSize()
        {
            auto constexpr t_field_size = FieldTraits<t_field>::template size<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template size<t_element>();
            return t_field_size * t_basis_size;
        }

        template<Field t_field, Basis t_basis>
        void
        setUnknown(
            std::shared_ptr<DegreeOfFreedom> & degree_of_freedom
        )
        {
            auto constexpr t_size = FiniteElementDegreeOfFreedom<t_element, t_domain>::template getSize<t_field, t_basis>();
            unknown_ = std::make_shared<FiniteElementDegreeOfFreedom<t_element, t_domain>>(FiniteElementDegreeOfFreedom<t_element, t_domain>());
            unknown_->index_ = degree_of_freedom->coefficients_.size();
            unknown_->degree_of_freedom_ = degree_of_freedom;
            degree_of_freedom->coefficients_.resize(degree_of_freedom->coefficients_.size() + t_size);
        }

        template<Field t_field, Basis t_basis>
        void
        setBinding(
            std::shared_ptr<DegreeOfFreedom> & degree_of_freedom
        )
        {
            auto constexpr t_size = FiniteElementDegreeOfFreedom<t_element, t_domain>::template getSize<t_field, t_basis>();
            binding_ = std::make_shared<FiniteElementDegreeOfFreedom<t_element, t_domain>>(FiniteElementDegreeOfFreedom<t_element, t_domain>());
            binding_->index_ = degree_of_freedom->coefficients_.size();
            binding_->degree_of_freedom_ = degree_of_freedom;
            degree_of_freedom->coefficients_.resize(degree_of_freedom->coefficients_.size() + t_size);
        }

        std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> unknown_;

        std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> binding_;

    };
    
    template<Domain t_domain, BehaviorConcept auto t_behavior>
    struct IntegrationPoint
    {

        using t_BehaviorTraits = BehaviorTraits<t_behavior>;

        void
        initialize(
            std::shared_ptr<mgis::behaviour::Behaviour> behaviour
        )
        {
            behaviour_ = behaviour;
            material_point_ = std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behaviour_));
        }

        void
        integrate()
        {
            auto behaviour_view = mgis::behaviour::make_view(* material_point_);
            auto result = mgis::behaviour::integrate(behaviour_view, * behaviour_);
        }
        
        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, t_BehaviorTraits::template getGeneralizedStrainSize<t_domain>()> const>
        getGeneralizedStrain()
        const
        {
            auto constexpr size = t_BehaviorTraits::template getGeneralizedStrainSize<t_domain>();
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, size> const>(material_point_->s1.gradients.data());
        }
        
        Point coordinates_;
        
        lolita::real weight_;
        
        std::shared_ptr<mgis::behaviour::Behaviour> behaviour_;
        
        std::unique_ptr<mgis::behaviour::BehaviourData> material_point_;

    };
    
    template<Element t_element, Domain t_domain, auto t_finite_element_method>
    struct QuadraturePoint
    {

        // using t_FiniteELementQuadrature = ElementQuadratureRuleTraits<t_element, t_finite_element_method.getQuadrature()>;

        using t_FiniteElementMethodTraits = FiniteElementMethodTraits<t_finite_element_method>;
        
        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, t_FiniteElementMethodTraits::template getGeneralizedStrainSize<t_domain>()> const>
        getGeneralizedStrain()
        const
        {
            auto constexpr size = t_FiniteElementMethodTraits::template getGeneralizedStrainSize<t_domain>();
            auto constexpr offset = t_FiniteElementMethodTraits::template getGeneralizedStrainOffset<t_domain>();
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, size> const>(integration_point_.material_point_->s1.gradients.data() + offset);
        }
        
        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, t_FiniteElementMethodTraits::template getGeneralizedStrainSize<t_domain>()>>
        getGeneralizedStrain()
        {
            auto constexpr size = t_FiniteElementMethodTraits::template getGeneralizedStrainSize<t_domain>();
            auto constexpr offset = t_FiniteElementMethodTraits::template getGeneralizedStrainOffset<t_domain>();
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, size>>(integration_point_.material_point_->s1.gradients.data() + offset);
        }
        
        IntegrationPoint<t_domain, t_finite_element_method.getBehavior()> integration_point_;
        
        lolita::matrix::Matrix<lolita::real> generalized_strain_mapping_;

    };

    template<Element t_element, Domain t_domain, auto t_arg>
    struct FiniteElement : FiniteElementGeometry<FiniteElement, t_element, t_domain, t_arg>
    {

        struct Data
        {

            std::vector<lolita::matrix::Matrix<lolita::real>> element_operators_;

            std::vector<std::unique_ptr<FiniteElementUnknown<t_element, t_domain>>> element_unknowns_;

            std::vector<std::unique_ptr<QuadraturePoint<t_element, t_domain, t_arg>>> quadrature_points_;

        };

        void
        activate()
        {
            data_ = std::make_unique<Data>(Data());
        }

        std::unique_ptr<Data> data_;

    };

    // template<auto t_finite_element_method>
    // requires(t_finite_element_method.isHdg())
    // struct FiniteElementMethodTraits<t_finite_element_method> : FiniteElementMethodTraitsBase<t_finite_element_method>
    // {

    //     template<Element t_element, Domain t_domain>
    //     static constexpr
    //     lolita::integer
    //     getNumCellUnknowns()
    //     {
    //         auto constexpr field_size = FieldTraits<t_finite_element_method.getGeneralizedStrain().getField()>::template size<t_domain>();
    //         auto constexpr basis_size = FiniteElementBasisTraits<t_finite_element_method.discretization_.cell_basis_>::template size<t_element>();
    //         return field_size * basis_size;
    //     }

    //     template<Element t_element, Domain t_domain>
    //     static constexpr
    //     lolita::integer
    //     getNumFaceUnknowns()
    //     {
    //         auto constexpr field_size = FieldTraits<t_finite_element_method.getGeneralizedStrain().getField()>::template size<t_domain>();
    //         auto constexpr basis_size = FiniteElementBasisTraits<t_finite_element_method.discretization_.face_basis_>::template size<t_element>();
    //         return field_size * basis_size;
    //     }

    //     template<Element t_element, Domain t_domain>
    //     static constexpr
    //     lolita::integer
    //     getNumElementUnknowns()
    //     {
    //         auto num_element_unknowns = getNumCellUnknowns<t_element, t_domain>();
    //         auto set_num_faces_unknowns = [&] <lolita::integer t_i = 0> (
    //             auto & self
    //         )
    //         constexpr mutable
    //         {
    //             auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
    //             auto constexpr t_num_inner_neighbors = ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0, t_i>();
    //             num_element_unknowns += getNumFaceUnknowns<t_inner_neighbor, t_domain>() * t_num_inner_neighbors;
    //             if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
    //             {
    //                 self.template operator ()<t_i + 1>(self);
    //             }
    //         };
    //         set_num_faces_unknowns(set_num_faces_unknowns);
    //         return num_element_unknowns;
    //     }

    //     template<Element t_element, Domain t_domain>
    //     using Stabilization = lolita::matrix::Matrix<lolita::real, getNumElementUnknowns<t_element, t_domain>(), getNumElementUnknowns<t_element, t_domain>()>;

    //     // template<Element t_element, Domain t_domain>
    //     // struct Cell
    //     // {

    //     //     using t_ElementQuadratureTraits = ElementQuadratureRuleTraits<t_element, t_finite_element_method.getQuadrature()>;

    //     //     using t_Stabilization = lolita::matrix::Matrix<lolita::real, getNumElementUnknowns<t_element, t_domain>(), getNumElementUnknowns<t_element, t_domain>()>;

    //     //     lolita::matrix::Matrix<lolita::real, getNumElementUnknowns<t_element, t_domain>(), getNumElementUnknowns<t_element, t_domain>()> stabilization_;

    //     //     std::array<QuadraturePoint<t_element, t_domain, t_finite_element_method>, t_ElementQuadratureTraits::size()> quadrature_points_;
        
    //     //     template<template<Element, Domain, auto...> typename t_FiniteElement>
    //     //     struct Implementation : FiniteElementGeometry<t_FiniteElement, t_element, t_domain, t_finite_element_method>
    //     //     {

    //     //     };

    //     // };

    // };

}

#endif /* CC79BDC5_49DB_4A81_8A93_18ABD6551AF1 */
