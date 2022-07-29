#ifndef ECF0EB38_E45D_4E38_ACEB_E40188311824
#define ECF0EB38_E45D_4E38_ACEB_E40188311824

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_000.hxx"
#include "lolita/lolita_core_n_001.hxx"
#include "lolita/lolita_core_n_002.hxx"
#include "lolita/lolita_core_n_003.hxx"

namespace lolita2::geometry
{

    struct Load
    {

        auto static constexpr zero = [] (auto const &, auto const &) constexpr { return lolita::real(0); };

        Load()
        :
        loading_(std::make_shared<Loading>(zero))
        {}

        Load(
            Loading const & loading
        )
        :
        loading_(std::make_shared<Loading>(loading))
        {}

        Load(
            Loading && loading
        )
        :
        loading_(std::make_shared<Loading>(std::forward<Loading>(loading)))
        {}
        
        lolita::real
        getImposedValue(
            lolita2::Point const & point,
            lolita::real const & time
        )
        const
        {
            return loading_->operator ()(point, time);
        }

        std::shared_ptr<Loading> loading_;

    };

    template<Field t_field, Domain t_domain>
    struct FieldLoad
    {

        static constexpr
        lolita::integer
        getRows()
        {
            return FieldTraits<t_field>::template shape<t_domain>().rows();
        }

        static constexpr
        lolita::integer
        getCols()
        {
            return FieldTraits<t_field>::template shape<t_domain>().cols();
        }

        std::array<std::array<Load, getRows()>, getCols()> loads_;

    };

    template<Element t_element, Domain t_domain, auto t_arg>
    struct FiniteElement;
        
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
        
        template<Element t_element, Domain t_domain, auto t_arg>
        struct Implementation : FiniteElement<t_element, t_domain, t_arg>
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

    template<auto t_finite_element_method>
    struct DiscretizationTraits;

    template<auto t_finite_element_method>
    requires(t_finite_element_method.getDiscretization().isHybridDiscontinuousGalerkin())
    struct DiscretizationTraits<t_finite_element_method>
    {

        template<Element t_element, Domain t_domain>
        static constexpr
        lolita::integer
        getNumCellUnknowns()
        {
            auto constexpr field_size = FieldTraits<t_finite_element_method.getGeneralizedStrain().getField()>::template size<t_domain>();
            auto constexpr basis_size = FiniteElementBasisTraits<t_finite_element_method.discretization_.cell_basis_>::template size<t_element>();
            return field_size * basis_size;
        }

        template<Element t_element, Domain t_domain>
        static constexpr
        lolita::integer
        getNumFaceUnknowns()
        {
            auto constexpr field_size = FieldTraits<t_finite_element_method.getGeneralizedStrain().getField()>::template size<t_domain>();
            auto constexpr basis_size = FiniteElementBasisTraits<t_finite_element_method.discretization_.face_basis_>::template size<t_element>();
            return field_size * basis_size;
        }

        template<Element t_element, Domain t_domain>
        static constexpr
        lolita::integer
        getNumCols()
        {
            auto num_element_unknowns = getNumCellUnknowns<t_element, t_domain>();
            auto set_num_faces_unknowns = [&] <lolita::integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
                auto constexpr t_num_inner_neighbors = ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0, t_i>();
                num_element_unknowns += getNumFaceUnknowns<t_inner_neighbor, t_domain>() * t_num_inner_neighbors;
                if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            set_num_faces_unknowns(set_num_faces_unknowns);
            return num_element_unknowns;
        }

        template<Element t_element, Domain t_domain>
        struct Implementation : FiniteElement<t_element, t_domain, t_finite_element_method>
        {

        private:

            template<Basis t_gradient_basis>
            static constexpr
            lolita::integer
            getScalarBasisSize()
            {
                return FiniteElementBasisTraits<t_gradient_basis>::template size<t_element>();
            }

        public:

            template<Basis t_gradient_basis>
            lolita::matrix::Matrix<lolita::real, getScalarBasisSize<t_gradient_basis>(), getScalarBasisSize<t_gradient_basis>()>
            getLhs()
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(2 * t_gradient_basis.getOrder());
                auto lhs = lolita::matrix::Matrix<lolita::real, getScalarBasisSize<t_gradient_basis>(), getScalarBasisSize<t_gradient_basis>()>();
                lhs.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
                {
                    auto weight = this->template getReferenceQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto vec = this->template getBasisEvaluation<t_gradient_basis>(point);
                    lhs += weight * vec * vec.transpose();
                }
                std::cout << "coucou !" << std::endl;
                return lhs;
            }

        };
        
        // template<Mapping t_mapping>
        // struct MappingFactory;
        
        // template<Mapping t_mapping>
        // requires(t_mapping.isIdentity() || !t_mapping.isIdentity)
        // struct MappingFactory<t_mapping>
        // {
                
        //         template<Basis t_gradient_basis>
        //         lolita::matrix::Matrix<lolita::real, FiniteElementBasisTraits<t_gradient_basis>::template size<t_element>()>
        //         getLhs()
        //         const
        //         {
        //             auto constexpr t_quadrature = Quadrature::gauss(2 * t_gradient_basis.getOrder());
        //             auto lhs = lolita::matrix::Matrix<lolita::real, FiniteElementBasisTraits<t_gradient_basis>::template size<t_element>()>();
        //             lhs.setZero();
        //             for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
        //             {
        //                 auto weight = this->template getReferenceQuadratureWeight<t_quadrature>(i);
        //                 auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
        //                 auto vec = this->template getBasisEvaluation<t_gradient_basis>(point);
        //                 lhs += weight * vec * vec.transpose();
        //             }
        //             return lhs;
        //         }

        //     };

        // };

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
    
    template<Domain t_domain, BehaviorConcept auto t_behavior>
    struct IntegrationPoint
    {

        using t_BehaviorTraits = BehaviorTraits<t_behavior>;

        void
        setBehavior(
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

        using t_FiniteElementMethodTraits = FiniteElementMethodTraits<t_finite_element_method>;
        
        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, t_FiniteElementMethodTraits::template getGeneralizedStrainSize<t_domain>()> const>
        getGeneralizedStrain()
        const
        {
            auto constexpr size = t_FiniteElementMethodTraits::template getGeneralizedStrainSize<t_domain>();
            auto constexpr offset = t_FiniteElementMethodTraits::template getGeneralizedStrainOffset<t_domain>();
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, size> const>(integration_point_->material_point_->s1.gradients.data() + offset);
        }
        
        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, t_FiniteElementMethodTraits::template getGeneralizedStrainSize<t_domain>()>>
        getGeneralizedStrain()
        {
            auto constexpr size = t_FiniteElementMethodTraits::template getGeneralizedStrainSize<t_domain>();
            auto constexpr offset = t_FiniteElementMethodTraits::template getGeneralizedStrainOffset<t_domain>();
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, size>>(integration_point_->material_point_->s1.gradients.data() + offset);
        }
        
        std::shared_ptr<IntegrationPoint<t_domain, t_finite_element_method.getBehavior()>> integration_point_;
        
        lolita::matrix::Matrix<lolita::real> generalized_strain_mapping_;

    };

    template<Element t_element, Domain t_domain, auto t_arg>
    struct FiniteElement : FiniteElementGeometry<FiniteElement, t_element, t_domain, t_arg>
    {

        struct Data
        {

            std::vector<lolita::matrix::Matrix<lolita::real>> element_operators_;

            std::vector<std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>>> element_degrees_of_freedom_;

            std::vector<std::unique_ptr<QuadraturePoint<t_element, t_domain, t_arg>>> quadrature_points_;

        };

        template<Field t_field, Basis t_basis>
        void
        setDegreeOfFreedom(
            std::shared_ptr<DegreeOfFreedom> & degree_of_freedom
        )
        {
            auto constexpr t_size = FiniteElementDegreeOfFreedom<t_element, t_domain>::template getSize<t_field, t_basis>();
            auto element_dof = std::make_unique<FiniteElementDegreeOfFreedom<t_element, t_domain>>(FiniteElementDegreeOfFreedom<t_element, t_domain>());
            element_dof->index_ = degree_of_freedom->coefficients_.size();
            element_dof->degree_of_freedom_ = degree_of_freedom;
            degree_of_freedom->coefficients_.resize(degree_of_freedom->coefficients_.size() + t_size);
            data_->element_degrees_of_freedom_.push_back(std::move(element_dof));
        }
        
        template<Basis t_basis>
        lolita::matrix::Vector<lolita::real, FiniteElementBasisTraits<t_basis>::template size<t_element>()>
        getBasisEvaluation(
            Point const & point
        )
        const
        {
            using t_FiniteElementBasisTraits = FiniteElementBasisTraits<t_basis>;
            using t_Implementation = typename t_FiniteElementBasisTraits::template Implementation<t_element, t_domain, t_arg>;
            return static_cast<t_Implementation const *>(this)->getBasisEvaluation(point);
        }
        
        template<Basis t_basis>
        lolita::matrix::Vector<lolita::real, FiniteElementBasisTraits<t_basis>::template size<t_element>()>
        getBasisDerivative(
            Point const & point,
            lolita::integer derivative_direction
        )
        const
        {
            using t_FiniteElementBasisTraits = FiniteElementBasisTraits<t_basis>;
            using t_Implementation = typename t_FiniteElementBasisTraits::template Implementation<t_element, t_domain, t_arg>;
            return static_cast<t_Implementation const *>(this)->getBasisDerivative(point, derivative_direction);
        }

        template<Basis t_gradient_basis>
        auto
        getLhs()
        const
        {
            using t_DiscretizationTraits = DiscretizationTraits<t_arg>;
            using t_Implementation = typename t_DiscretizationTraits::template Implementation<t_element, t_domain>;
            return static_cast<t_Implementation const *>(this)->template getLhs<t_gradient_basis>();
        }

        void
        activate()
        {
            data_ = std::make_unique<Data>(Data());
            auto constexpr bas = Basis::monomial(1);
            this->template getLhs<bas>();
        }

        lolita::boolean
        isActivated()
        const
        {
            return * data_ == nullptr;
        }

        std::unique_ptr<Data> data_;

    };

    template<Element t_element, Domain t_domain, auto... t_args>
    struct FiniteElementHolder : FiniteElementGeometry<FiniteElementHolder, t_element, t_domain, t_args...>
    {

        using t_FiniteElements = std::tuple<std::shared_ptr<FiniteElement<t_element, t_domain, t_args>>...>;

        template<lolita::integer t_i>
        using t_FiniteElement = typename std::tuple_element_t<t_i, t_FiniteElements>::element_type;
          
        template<lolita::integer t_i>
        std::tuple_element_t<t_i, t_FiniteElements> const &
        getFiniteElement()
        const
        {
            return std::get<t_i>(finite_elements_);
        }
        
        template<lolita::integer t_i>
        std::tuple_element_t<t_i, t_FiniteElements> &
        getFiniteElement()
        {
            return std::get<t_i>(finite_elements_);
        }

        t_FiniteElements finite_elements_;

    };
    
}


#endif /* ECF0EB38_E45D_4E38_ACEB_E40188311824 */
