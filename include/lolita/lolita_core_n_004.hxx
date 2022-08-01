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

        Load(
            Loading const & loading,
            lolita::integer row,
            lolita::integer col
        )
        :
        row_(row),
        col_(col),
        loading_(loading)
        {}

        Load(
            Loading && loading,
            lolita::integer row,
            lolita::integer col
        )
        :
        row_(row),
        col_(col),
        loading_(std::forward<Loading>(loading))
        {}
        
        lolita::real
        getImposedValue(
            Point const & point,
            lolita::real const & time
        )
        const
        {
            return loading_(point, time);
        }

        Loading loading_;

        lolita::integer row_;

        lolita::integer col_;

    };

    template<Domain t_domain, Field t_field>
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

        void
        setLoad(
            lolita::integer row,
            lolita::integer col,
            std::shared_ptr<Loading> const & loading
        )
        {
            loads_[col][row] = loading;
        }

        void
        setLoad(
            lolita::integer row,
            lolita::integer col,
            std::shared_ptr<Loading> && loading
        )
        {
            loads_[col][row] = std::forward<std::shared_ptr<Loading>>(loading);
        }
        
        lolita::real
        getImposedValue(
            lolita::integer row,
            lolita::integer col,
            Point const & point,
            lolita::real const & time
        )
        const
        {
            if (loads_[col][row] == nullptr)
            {
                return 0.0;
            }
            else
            {
                return loads_[col][row]->operator ()(point, time);
            }
        }

    private:

        std::array<std::array<std::shared_ptr<Loading>, getRows()>, getCols()> loads_;

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

        template<Element t_element>
        static constexpr
        lolita::integer
        getSize()
        {
            return lolita::numerics::binomial(t_element.dim_ + t_basis.ord_, t_element.dim_);
        }
        
        template<Element t_element, Domain t_domain, auto t_arg>
        struct Implementation : FiniteElement<t_element, t_domain, t_arg>
        {

        private:
        
            static constexpr
            std::array<std::array<lolita::integer, 3>, getSize<t_element>()>
            getExponents()
            {
                auto exponents = std::array<std::array<lolita::integer, 3>, getSize<t_element>()>();
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
            
            std::array<std::array<lolita::integer, 3>, getSize<t_element>()> static constexpr exponents_ = getExponents();

        public:
        
            lolita::matrix::Vector<lolita::real, getSize<t_element>()>
            getBasisEvaluation(
                Point const & point
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<lolita::real, getSize<t_element>()>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (lolita::integer i = 0; i < getSize<t_element>(); ++i)
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
            
            lolita::matrix::Vector<lolita::real, getSize<t_element>()>
            getBasisDerivative(
                Point const & point,
                lolita::integer derivative_direction
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<lolita::real, getSize<t_element>()>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (lolita::integer i = 0; i < getSize<t_element>(); ++i)
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

            Field static constexpr t_field = t_finite_element_method.getField();

            template<Basis t_basis>
            static constexpr
            lolita::integer
            getBasisSize()
            {
                return FiniteElementBasisTraits<t_basis>::template size<t_element>();
            }

        public:

            template<Basis t_gradient_basis>
            lolita::matrix::Matrix<lolita::real, getBasisSize<t_gradient_basis>(), getBasisSize<t_gradient_basis>()>
            getLhs()
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(2 * t_gradient_basis.getOrder());
                auto lhs = lolita::matrix::Matrix<lolita::real, getBasisSize<t_gradient_basis>(), getBasisSize<t_gradient_basis>()>();
                lhs.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
                {
                    auto weight = this->template getReferenceQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto vec = this->template getBasisEvaluation<t_gradient_basis>(point);
                    lhs += weight * vec * vec.transpose();
                }
                // std::cout << "coucou !" << std::endl;
                return lhs;
            }

            template<Basis t_gradient_basis>
            lolita::matrix::Matrix<lolita::real, getBasisSize<t_gradient_basis>(), getBasisSize<t_gradient_basis>()>
            getGradientLhs()
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(2 * t_gradient_basis.getOrder());
                auto lhs = lolita::matrix::Matrix<lolita::real, getBasisSize<t_gradient_basis>(), getBasisSize<t_gradient_basis>()>();
                lhs.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
                {
                    auto weight = this->template getReferenceQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto vec = this->template getBasisEvaluation<t_gradient_basis>(point);
                    lhs += weight * vec * vec.transpose();
                }
                // std::cout << "coucou !" << std::endl;
                return lhs;
            }

            template<Basis t_gradient_basis>
            lolita::matrix::Matrix<lolita::real, getBasisSize<t_gradient_basis>(), getNumCols<t_element, t_domain>()>
            getGradientRhs()
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(2 * t_gradient_basis.getOrder());
                auto rhs = lolita::matrix::Matrix<lolita::real, getBasisSize<t_gradient_basis>(), getNumCols<t_element, t_domain>()>();
                rhs.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
                {
                    auto weight = this->template getReferenceQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto vec = this->template getBasisEvaluation<t_gradient_basis>(point);
                    // rhs += weight * vec * vec.transpose();
                }
                // std::cout << "coucou !" << std::endl;
                return rhs;
            }

            // template<Mapping t_mapping, Basis t_gradient_basis>
            // lolita::matrix::Matrix<lolita::real, getMappingSize<t_mapping>(), getNumCols<t_element, t_domain>()>
            // getGradient(
            //     lolita::integer i
            // )
            // const
            // requires(t_mapping.isGradient())
            // {
            //     auto rhs = lolita::matrix::Matrix<lolita::real, getMappingSize<t_mapping>(), getNumCols<t_element, t_domain>()>();
            //     rhs.setZero();
            //     return rhs;
            // }

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

        // template<Basis t_basis>
        // lolita::matrix::Span<lolita::matrix::Vector<lolita::real, getSize<t_field, t_basis>()> const>
        // getCoefficients()
        // const
        // {
        //     auto const & data = degree_of_freedom_->coefficients_.data() + index_;
        //     return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, getSize<t_field, t_basis>()> const>(data);
        // }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>
        getCoefficients(
            lolita::integer row,
            lolita::integer col
        )
        {
            auto constexpr t_field_shape = FieldTraits<t_field>::template shape<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
            auto const & data = degree_of_freedom_->coefficients_.data() + index_ + (t_field_shape.cols() * row  + col) * t_basis_size;
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>(data);
        }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()> const>
        getCoefficients(
            lolita::integer row,
            lolita::integer col
        )
        const
        {
            auto constexpr t_field_shape = FieldTraits<t_field>::template shape<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
            auto const & data = degree_of_freedom_->coefficients_.data() + index_ + (t_field_shape.cols() * row  + col) * t_basis_size;
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>(data);
        }
        
        std::shared_ptr<DegreeOfFreedom> const &
        getDegreeOfFreedom()
        const
        {
            return degree_of_freedom_;
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

            std::unique_ptr<FieldLoad<t_domain, t_arg.getField()>> loads_;

        };

        void
        setLoad(
            lolita::integer row,
            lolita::integer col,
            std::shared_ptr<Loading> const & load
        )
        {
            if (data_->loads_ == nullptr)
            {
                data_->loads_ = std::make_unique<FieldLoad<t_domain, t_arg.getField()>>(FieldLoad<t_domain, t_arg.getField()>());
            }
            data_->loads_->setLoad(row, col, load);
        }

        lolita::boolean
        hasDegreeOfFreedom(
            std::basic_string_view<lolita::character> label
        )
        const
        {
            auto has_dof = [&] (
                std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> const & dof
            )
            {
                return dof->getDegreeOfFreedom()->getLabel() == label;
            };
            auto iter = std::find_if(data_->element_degrees_of_freedom_.begin(), data_->element_degrees_of_freedom_.end(), has_dof);
            return iter != data_->element_degrees_of_freedom_.end();
        }

        template<Field t_field, Basis t_basis>
        void
        setDegreeOfFreedom(
            std::shared_ptr<DegreeOfFreedom> & degree_of_freedom
        )
        {
            auto constexpr t_size = FiniteElementDegreeOfFreedom<t_element, t_domain>::template getSize<t_field, t_basis>();
            auto has_dof = [&] (
                std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> const & dof
            )
            {
                return dof->getDegreeOfFreedom() == degree_of_freedom;
            };
            auto iter = std::find_if(data_->element_degrees_of_freedom_.begin(), data_->element_degrees_of_freedom_.end(), has_dof);
            if (iter == data_->element_degrees_of_freedom_.end())
            {
                auto element_dof = std::make_unique<FiniteElementDegreeOfFreedom<t_element, t_domain>>(FiniteElementDegreeOfFreedom<t_element, t_domain>());
                element_dof->index_ = degree_of_freedom->coefficients_.size();
                element_dof->degree_of_freedom_ = degree_of_freedom;
                degree_of_freedom->coefficients_.resize(degree_of_freedom->coefficients_.size() + t_size);
                data_->element_degrees_of_freedom_.push_back(std::move(element_dof));
            }
        }

        std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> const &
        getDegreeOfFreedom(
            std::basic_string_view<lolita::character> label
        )
        const
        {
            auto has_dof = [&] (
                std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> const & ld
            )
            {
                return ld->getDegreeOfFreedom()->getLabel() == label;
            };
            auto iter = std::find_if(data_->element_degrees_of_freedom_.begin(), data_->element_degrees_of_freedom_.end(), has_dof);
            if (iter != data_->element_degrees_of_freedom_.end())
            {
                return * iter;
            }
            else
            {
                throw std::runtime_error("NO");
            }
        }

        std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> &
        getDegreeOfFreedom(
            std::basic_string_view<lolita::character> label
        )
        {
            auto has_dof = [&] (
                std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> const & ld
            )
            {
                return ld->getDegreeOfFreedom()->getLabel() == label;
            };
            auto iter = std::find_if(data_->element_degrees_of_freedom_.begin(), data_->element_degrees_of_freedom_.end(), has_dof);
            if (iter != data_->element_degrees_of_freedom_.end())
            {
                return * iter;
            }
            else
            {
                throw std::runtime_error("NO");
            }
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

        template<Basis t_gradient_basis>
        auto
        getGradientRhs()
        const
        {
            using t_DiscretizationTraits = DiscretizationTraits<t_arg>;
            using t_Implementation = typename t_DiscretizationTraits::template Implementation<t_element, t_domain>;
            return static_cast<t_Implementation const *>(this)->template getGradientRhs<t_gradient_basis>();
        }

        // template<Mapping t_mapping, Basis t_gradient_basis>
        // lolita::matrix::Matrix<lolita::real, DiscretizationTraits<t_finite_element_method>::template getMappingSize<t_mapping>(), getNumCols<t_element, t_domain>()>
        // getGradient(
        //     lolita::integer i
        // )
        // const
        // {
        //     using t_DiscretizationTraits = DiscretizationTraits<t_arg>;
        //     using t_Implementation = typename t_DiscretizationTraits::template Implementation<t_element, t_domain>;
        //     return static_cast<t_Implementation const *>(this)->template getGradientRhs<t_mapping, t_gradient_basis>();
        // }

        void
        activate()
        {
            data_ = std::make_unique<Data>(Data());
            auto constexpr bas = Basis::monomial(1);
            this->template getLhs<bas>();
            this->template getGradientRhs<bas>();
        }

        lolita::boolean
        isActivated()
        const
        {
            return data_ == nullptr;
        }

        std::unique_ptr<Data> data_;

    };

    template<Element t_element, Domain t_domain, auto... t_args>
    struct FiniteElementHolder : FiniteElementGeometry<FiniteElementHolder, t_element, t_domain, t_args...>
    {

        using t_FiniteElements = std::tuple<std::shared_ptr<FiniteElement<t_element, t_domain, t_args>>...>;

        template<lolita::integer t_i>
        using t_FiniteElement = typename std::tuple_element_t<t_i, t_FiniteElements>::element_type;

        template<FiniteElementMethodConcept auto t_arg>
        static constexpr
        lolita::integer
        getArgIndex()
        {
            auto index = lolita::integer(0);
            auto found = lolita::boolean(false);
            auto set_index = [&] (
                FiniteElementMethodConcept auto const & arg
            )
            constexpr
            {
                if constexpr (std::is_same_v<std::decay_t<decltype(t_arg)>, std::decay_t<decltype(arg)>>)
                {
                    if (t_arg == arg)
                    {
                        found = true;
                    }
                }
                if (!found)
                {
                    index += 1;
                }
            };
            (set_index(t_args), ...);
            return index;
        }
          
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
        
        template<FiniteElementMethodConcept auto t_arg>
        std::tuple_element_t<getArgIndex<t_arg>(), t_FiniteElements> const &
        getFiniteElement()
        const
        {
            return std::get<getArgIndex<t_arg>()>(finite_elements_);
        }
        
        template<FiniteElementMethodConcept auto t_arg>
        std::tuple_element_t<getArgIndex<t_arg>(), t_FiniteElements> &
        getFiniteElement()
        {
            return std::get<getArgIndex<t_arg>()>(finite_elements_);
        }

        t_FiniteElements finite_elements_;

    };
    
}


#endif /* ECF0EB38_E45D_4E38_ACEB_E40188311824 */
