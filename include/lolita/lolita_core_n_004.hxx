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

    template<auto... t_args>
    using RealMatrix = lolita::matrix::Matrix<lolita::real, t_args...>;

    template<auto... t_args>
    using RealVector = lolita::matrix::Vector<lolita::real, t_args...>;

    template<Element t_element, Domain t_domain, auto t_finite_element_method>
    struct FiniteElementTraits;

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

    template<Element t_element, Domain t_domain, auto t_finite_element_method>
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

        template<lolita::integer t_dim>
        static constexpr
        lolita::integer
        getSize()
        {
            return lolita::numerics::binomial(t_dim + t_basis.ord_, t_dim);
        }

        template<Element t_element>
        static constexpr
        lolita::integer
        getSize()
        {
            return lolita::numerics::binomial(t_element.dim_ + t_basis.ord_, t_element.dim_);
        }
        
        template<Element t_element, Domain t_domain, auto t_finite_element_method>
        struct Implementation : FiniteElement<t_element, t_domain, t_finite_element_method>
        {

            static constexpr
            lolita::integer
            getSize()
            {
                return FiniteElementBasisTraits::template getSize<t_element>();
            }

        private:
        
            static constexpr
            std::array<std::array<lolita::integer, 3>, getSize()>
            getExponents()
            {
                auto exponents = std::array<std::array<lolita::integer, 3>, getSize()>();
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
            
            std::array<std::array<lolita::integer, 3>, getSize()> static constexpr exponents_ = getExponents();

        public:
        
            lolita::matrix::Vector<lolita::real, getSize()>
            getBasisEvaluation(
                Point const & point
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<lolita::real, getSize()>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (lolita::integer i = 0; i < getSize(); ++i)
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
            
            lolita::matrix::Vector<lolita::real, getSize()>
            getBasisDerivative(
                Point const & point,
                lolita::integer derivative_direction
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<lolita::real, getSize()>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (lolita::integer i = 0; i < getSize(); ++i)
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

    template<auto t_discretization>
    struct HybridDiscontinuousGalerkinTraits
    {

        static constexpr
        Basis
        getCellBasis()
        {
            return t_discretization.cell_basis_;
        }

        static constexpr
        Basis
        getFaceBasis()
        {
            return t_discretization.face_basis_;
        }

        static constexpr
        Basis
        getGradBasis()
        {
            return t_discretization.grad_basis_;
        }

        template<Element t_element>
        static constexpr
        lolita::integer
        getCellBasisSize()
        {
            return FiniteElementBasisTraits<t_discretization.cell_basis_>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        lolita::integer
        getFaceBasisSize()
        {
            return FiniteElementBasisTraits<t_discretization.face_basis_>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        lolita::integer
        getGradBasisSize()
        {
            return FiniteElementBasisTraits<t_discretization.grad_basis_>::template getSize<t_element>();
        }

        template<Element t_element, Domain t_domain, Field t_field>
        static constexpr
        lolita::integer
        getNumCellUnknowns()
        {
            auto constexpr field_size = FieldTraits<t_field>::template size<t_domain>();
            auto constexpr basis_size = FiniteElementBasisTraits<t_discretization.cell_basis_>::template getSize<t_element>();
            return field_size * basis_size;
        }

        template<Element t_element, Domain t_domain, Field t_field>
        static constexpr
        lolita::integer
        getNumFaceUnknowns()
        {
            auto constexpr field_size = FieldTraits<t_field>::template size<t_domain>();
            auto constexpr basis_size = FiniteElementBasisTraits<t_discretization.face_basis_>::template getSize<t_element>();
            return field_size * basis_size;
        }

        template<Element t_element, Domain t_domain, Field t_field>
        static constexpr
        lolita::integer
        getNumElementUnknowns()
        {
            auto num_element_unknowns = getNumCellUnknowns<t_element, t_domain, t_field>();
            auto set_num_faces_unknowns = [&] <lolita::integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
                auto constexpr t_num_inner_neighbors = ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0, t_i>();
                num_element_unknowns += getNumFaceUnknowns<t_inner_neighbor, t_domain, t_field>() * t_num_inner_neighbors;
                if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            set_num_faces_unknowns(set_num_faces_unknowns);
            return num_element_unknowns;
        }

        template<Element t_element, Domain t_domain, auto t_finite_element_method>
        struct Implementation : FiniteElement<t_element, t_domain, t_finite_element_method>
        {

        private:

            using t_FiniteElementTraits = FiniteElementTraits<t_element, t_domain, t_finite_element_method>;

            template<Mapping t_mapping>
            static constexpr
            lolita::integer
            getMappingSize()
            {
                return MappingTraits<t_mapping>::template size<t_domain, t_finite_element_method.getField()>();
            }

        public:

            static constexpr
            lolita::integer
            getNumElementUnknowns()
            {
                return HybridDiscontinuousGalerkinTraits::template getNumElementUnknowns<t_element, t_domain, t_finite_element_method.getField()>();
            }
            
            RealMatrix<getGradBasisSize<t_element>(), getGradBasisSize<t_element>()>
            getGradientLhs()
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(2 * getGradBasis().getOrder());
                auto lhs = RealMatrix<getGradBasisSize<t_element>(), getGradBasisSize<t_element>()>();
                lhs.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
                {
                    auto weight = this->template getReferenceQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto test_vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    lhs += weight * test_vector * test_vector.transpose();
                }
                return lhs.llt().solve(decltype(lhs)::Identity());
            }

            RealMatrix<getGradBasisSize<t_element>(), getNumElementUnknowns()>
            getGradientRhs(
                lolita::integer row,
                lolita::integer col
            )
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(2 * getGradBasis().getOrder());
                auto constexpr t_field_size = FieldTraits<t_finite_element_method.getField()>::template size<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = RealMatrix<getGradBasisSize<t_element>(), getNumElementUnknowns()>();
                rhs.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
                {
                    auto weight = this->template getReferenceQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto test_vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    auto left_vector = this->template getBasisDerivative<getCellBasis()>(point, col);
                    auto col_offset = row * getCellBasisSize<t_element>();
                    auto block = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset);
                    block += weight * test_vector * left_vector.transpose();
                }
                auto set_faces_blocks = [&] <lolita::integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0, t_i>();
                    auto face_index = lolita::integer(0);
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        for (auto i = 0; i < ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::size(); i++)
                        {
                            auto weight = face->template getReferenceQuadratureWeight<t_quadrature>(i);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(i);
                            auto test_vector = this->template getBasisEvaluation<getGradBasis()>(point);
                            auto left_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto left_cell_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                            auto normal_vector = face->getNormalVector(point)(col) * this->template getInnerNeighborOrientation<0, t_i>(face);
                            auto face_col_offset = face_offset + (t_field_size * face_index + row) * getFaceBasisSize<t_inner_neighbor>();
                            auto face_block = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset);
                            face_block += weight * test_vector * left_face_vector.transpose();
                            auto cell_col_offset = row * getCellBasisSize<t_element>();
                            auto cell_block = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, cell_col_offset);
                            cell_block -= weight * test_vector * left_cell_vector.transpose();
                            face_index ++;
                        }
                    }
                    face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>() * t_num_inner_neighbors;
                    if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator()<t_i + 1>(self);
                    }
                };
                set_faces_blocks(set_faces_blocks);
                return rhs;
            }

            template<Mapping t_mapping>
            RealMatrix<getMappingSize<t_mapping>(), getNumElementUnknowns()>
            getMapping(
                Point const & point
            )
            const
            {
                auto grd = RealMatrix<getMappingSize<t_mapping>(), getNumElementUnknowns()>();
                grd.setZero();
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain, t_finite_element_method.getField()>())
                {
                    if (mapping_value.value() > 1.e-12)
                    {
                        auto rhs = getGradientRhs(mapping_value.row(), mapping_value.col());
                        auto line = grd.template block<1, getNumElementUnknowns()>(mapping_value.rank(), 0);
                        line = mapping_value.value() * this->template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
                    }
                    else
                    {
                        auto line = grd.template block<1, getNumElementUnknowns()>(mapping_value.rank(), 0);
                        line = RealMatrix<1, getNumElementUnknowns()>::Zero();
                    }
                }
                return grd;
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
    
    template<auto t_finite_element_method>
    struct DiscretizationTraits;

    template<Element t_element, Domain t_domain, auto t_finite_element_method>
    struct FiniteElementTraits
    {
        
        static constexpr
        Field const &
        getField()
        {
            return t_finite_element_method.getField();
        }

        template<Basis t_basis>
        static constexpr
        lolita::integer
        getBasisSize()
        {
            return FiniteElementBasisTraits<t_basis>::template size<t_element>();
        }
        
        static constexpr
        lolita::integer
        getGeneralizedStrainSize()
        {
            return GeneralizedStrainTraits<t_finite_element_method.getGeneralizedStrain()>::template size<t_domain>();
        }
        
        static constexpr
        lolita::integer
        getGeneralizedStrainOffset()
        {
            auto constexpr t_finite_element_generalized_strain = t_finite_element_method.getGeneralizedStrain();
            auto offset = lolita::integer(0);
            auto is_set = false;
            auto set_offset = [&] <lolita::integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr t_generalized_strain = t_finite_element_method.getBehavior().template getGeneralizedStrain<t_i>();
                if constexpr (std::is_same_v<std::decay_t<decltype(t_generalized_strain)>, std::decay_t<decltype(t_finite_element_generalized_strain)>>)
                {
                    if constexpr (t_generalized_strain == t_finite_element_generalized_strain)
                    {
                        is_set = true;
                    }
                }
                if (!is_set)
                {
                    offset += GeneralizedStrainTraits<t_generalized_strain>::template size<t_domain>();
                }
                if constexpr (t_i < t_finite_element_method.getBehavior().getNumGeneralizedStrains() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            set_offset(set_offset);
            return offset;
        }

        template<Mapping t_mapping>
        static constexpr
        lolita::integer
        getMappingSize()
        {
            return MappingTraits<t_mapping>::template size<t_domain, t_finite_element_method.getField()>();
        }

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

        using t_FiniteElementTraits = FiniteElementTraits<t_element, t_domain, t_finite_element_method>;
        
        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::getGeneralizedStrainSize()> const>
        getGeneralizedStrain()
        const
        {
            auto constexpr size = t_FiniteElementTraits::getGeneralizedStrainSize();
            auto constexpr offset = t_FiniteElementTraits::getGeneralizedStrainOffset();
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, size> const>(integration_point_->material_point_->s1.gradients.data() + offset);
        }
        
        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::getGeneralizedStrainSize()>>
        getGeneralizedStrain()
        {
            auto constexpr size = t_FiniteElementTraits::getGeneralizedStrainSize();
            auto constexpr offset = t_FiniteElementTraits::getGeneralizedStrainOffset();
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, size>>(integration_point_->material_point_->s1.gradients.data() + offset);
        }
        
        std::shared_ptr<IntegrationPoint<t_domain, t_finite_element_method.getBehavior()>> integration_point_;
        
        std::vector<lolita::matrix::Matrix<lolita::real>> quadrature_operators_;
        
        lolita::matrix::Matrix<lolita::real> generalized_strain_mapping_;

    };

    template<Element t_element, Domain t_domain, auto t_finite_element_method>
    struct FiniteElement : FiniteElementGeometry<FiniteElement, t_element, t_domain, t_finite_element_method>
    {

        using t_FiniteElementTraits = FiniteElementTraits<t_element, t_domain, t_finite_element_method>;

        template<Basis t_basis>
        using t_Basis = typename FiniteElementBasisTraits<t_basis>::template Implementation<t_element, t_domain, t_finite_element_method>;

        template<auto t_discretization>
        using t_Disc = typename HybridDiscontinuousGalerkinTraits<t_discretization>::template Implementation<t_element, t_domain, t_finite_element_method>;

        struct Data
        {

            std::vector<lolita::matrix::Matrix<lolita::real>> element_operators_;

            std::vector<std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>>> element_degrees_of_freedom_;

            std::vector<std::unique_ptr<QuadraturePoint<t_element, t_domain, t_finite_element_method>>> quadrature_points_;

            std::unique_ptr<FieldLoad<t_domain, t_finite_element_method.getField()>> loads_;

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
                data_->loads_ = std::make_unique<FieldLoad<t_domain, t_finite_element_method.getField()>>(FieldLoad<t_domain, t_finite_element_method.getField()>());
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
        RealVector<t_Basis<t_basis>::getSize()>
        getBasisEvaluation(
            Point const & point
        )
        const
        {
            return static_cast<t_Basis<t_basis> const *>(this)->getBasisEvaluation(point);
        }
        
        template<Basis t_basis>
        RealVector<t_Basis<t_basis>::getSize()>
        getBasisDerivative(
            Point const & point,
            lolita::integer derivative_direction
        )
        const
        {
            return static_cast<t_Basis<t_basis> const *>(this)->getBasisDerivative(point, derivative_direction);
        }

        template<Mapping t_mapping, auto t_discretization>
        RealMatrix<t_FiniteElementTraits::template getMappingSize<t_mapping>(), t_Disc<t_discretization>::getNumElementUnknowns()>
        getMapping(
            Point const & point
        )
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getMapping<t_mapping>(point);
        }

        void
        activate()
        {
            data_ = std::make_unique<Data>(Data());
            auto constexpr cell_basis = lolita2::Basis::monomial(1);
            auto constexpr face_basis = lolita2::Basis::monomial(1);
            auto point = Point();
            point.setZero();
            auto constexpr hdg = lolita2::HybridDiscontinuousGalerkin(cell_basis, face_basis, lolita2::HybridDiscontinuousGalerkin::Stabilization::Hdg);
            // auto constexpr bas = Basis::monomial(1);
            this->template getBasisEvaluation<cell_basis>(point);
            if constexpr (t_element.isSub(t_domain, 0))
            {
                auto mapp = this->template getMapping<Mapping::gradient(), hdg>(point);
                mapp.setZero();
                std::cout << "mapp : " << std::endl;
                std::cout << mapp << std::endl;
            }
            // this->template getGradientRhs<bas>();
        }

        template<auto t_discretization, Quadrature t_quadrature>
        void
        make()
        {
            using MappingOperator = RealMatrix<t_FiniteElementTraits::getGeneralizedStrainSize(), t_Disc<t_discretization>::getNumElementUnknowns()>;
            data_ = std::make_unique<Data>();
            for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
            {
                auto reference_weight = this->template getReferenceQuadratureWeight<t_quadrature>(i);
                auto reference_point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                auto current_weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                auto current_point = this->template getCurrentQuadraturePoint<t_quadrature>(i);
                auto qp = std::make_unique<QuadraturePoint<t_element, t_domain, t_finite_element_method>>();
                qp->generalized_strain_mapping_ = MappingOperator::Zero();
                auto set_quadrature_mapping = [&] <lolita::integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_mapping = t_finite_element_method.template getMapping<t_i>();
                    auto constexpr t_block_cols = t_Disc<t_discretization>::getNumElementUnknowns();
                    auto constexpr t_block_rows = FiniteElementMethodTraits<t_finite_element_method>::template getMappingSize<t_domain, t_mapping>();
                    auto constexpr t_block_offset = FiniteElementMethodTraits<t_finite_element_method>::template getMappingOffset<t_domain, t_mapping>();
                    std::cout << "*** t_i : " << t_i << std::endl;
                    std::cout << "mapping offset : " << t_block_offset << std::endl;
                    std::cout << "t_block_rows : " << t_block_rows << std::endl;
                    std::cout << "t_block_cols : " << t_block_cols << std::endl;
                    auto mapping_block = qp->generalized_strain_mapping_.template block<t_block_rows, t_block_cols>(t_block_offset, 0);
                    mapping_block = getMapping<t_mapping, t_discretization>(reference_point);
                    if constexpr (t_i < t_finite_element_method.getGeneralizedStrain().getNumMappings() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_quadrature_mapping(set_quadrature_mapping);
                data_->quadrature_points_.push_back(std::move(qp));
            }
            
        }

        lolita::boolean
        isActivated()
        const
        {
            return data_ == nullptr;
        }

        std::unique_ptr<Data> data_;

    };

    
    // template<Domain t_domain, BehaviorConcept auto t_behavior>
    // struct IntegrationPoints
    // {

    //     template<Quadrature t_quadrature>
    //     void
    //     makeBehavior()
    //     {
    //         for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
    //         {
    //             auto ip = std::make_shared<t_IntegrationPoint<getArgIndex<t_behavior>()>>;
    //             getIntegrationPoints<t_behavior>().push_back(ip);
    //         }
            
    //     }

    //     std::vector<std::shared_ptr<IntegrationPoint<t_domain, t_behavior>>>;

    // };

    template<Element t_element, Domain t_domain, auto... t_finite_element_methods>
    struct FiniteElementHolder : FiniteElementGeometry<FiniteElementHolder, t_element, t_domain, t_finite_element_methods...>
    {

        using t_IntegrationPoints = lolita::utility::unique_tuple_t<std::tuple<std::vector<std::shared_ptr<IntegrationPoint<t_domain, t_finite_element_methods.getBehavior()>>>...>>;

        // template<typename T, typename... U>
        // static constexpr
        // auto
        // data(
        //     T const & arg
        //     std::tuple<U...> const & tpl
        // )
        // {
        //     if constexpr ((std::is_same_v<T, U> || ...))
        //     {
        //         return tpl;
        //     }
        //     else
        //     {
        //         return std::make_tuple(std::get<U>(tpl)..., arg);
        //     }
        // }

        // data()

        using t_FiniteElements = std::tuple<std::shared_ptr<FiniteElement<t_element, t_domain, t_finite_element_methods>>...>;

        // using t_IntegrationPoints = std::tuple<std::vector<std::shared_ptr<IntegrationPoint<t_domain, t_finite_element_methods.getBehavior()>>>...>;

        template<lolita::integer t_i>
        using t_FiniteElement = typename std::tuple_element_t<t_i, std::tuple<FiniteElement<t_element, t_domain, t_finite_element_methods>...>>;

        template<lolita::integer t_i>
        using t_IntegrationPoint = typename std::tuple_element_t<t_i, lolita::utility::unique_tuple_t<std::tuple<IntegrationPoint<t_domain, t_finite_element_methods.getBehavior()>...>>>;

        // template<lolita::integer t_i>
        // using t_FiniteElement = typename std::tuple_element_t<t_i, std::tuple<FiniteElement<t_element, t_domain, t_finite_element_methods>...>>;

        // template<BehaviorConcept auto t_finite_element_method>
        // using t_IntegrationPoint = typename std::tuple_element_t<t_i, std::tuple<IntegrationPoint<t_domain, t_finite_element_methods.getBehavior()>...>>;

        template<BehaviorConcept auto t_finite_element_method>
        static constexpr
        lolita::integer
        getArgIndex()
        {
            auto index = lolita::integer(0);
            auto found = lolita::boolean(false);
            auto set_index = [&] (
                BehaviorConcept auto const & arg
            )
            constexpr
            {
                if constexpr (std::is_same_v<std::decay_t<decltype(t_finite_element_method)>, std::decay_t<decltype(arg)>>)
                {
                    if (t_finite_element_method == arg)
                    {
                        found = true;
                    }
                }
                if (!found)
                {
                    index += 1;
                }
            };
            (set_index(t_finite_element_methods.getBehavior()), ...);
            return index;
        }

        template<FiniteElementMethodConcept auto t_finite_element_method>
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
                if constexpr (std::is_same_v<std::decay_t<decltype(t_finite_element_method)>, std::decay_t<decltype(arg)>>)
                {
                    if (t_finite_element_method == arg)
                    {
                        found = true;
                    }
                }
                if (!found)
                {
                    index += 1;
                }
            };
            (set_index(t_finite_element_methods), ...);
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
        
        template<FiniteElementMethodConcept auto t_finite_element_method>
        std::tuple_element_t<getArgIndex<t_finite_element_method>(), t_FiniteElements> const &
        getFiniteElement()
        const
        {
            return std::get<getArgIndex<t_finite_element_method>()>(finite_elements_);
        }
        
        template<FiniteElementMethodConcept auto t_finite_element_method>
        std::tuple_element_t<getArgIndex<t_finite_element_method>(), t_FiniteElements> &
        getFiniteElement()
        {
            return std::get<getArgIndex<t_finite_element_method>()>(finite_elements_);
        }
        
        template<lolita::integer t_i>
        std::tuple_element_t<t_i, t_IntegrationPoints> const &
        getIntegrationPoints()
        const
        {
            return std::get<t_i>(integrations_points_);
        }
        
        template<lolita::integer t_i>
        std::tuple_element_t<t_i, t_IntegrationPoints> &
        getIntegrationPoints()
        {
            return std::get<t_i>(integrations_points_);
        }
        
        template<BehaviorConcept auto t_finite_element_method>
        std::tuple_element_t<getArgIndex<t_finite_element_method>(), t_IntegrationPoints> const &
        getIntegrationPoints()
        const
        {
            return std::get<getArgIndex<t_finite_element_method>()>(integrations_points_);
        }
        
        template<BehaviorConcept auto t_finite_element_method>
        std::tuple_element_t<getArgIndex<t_finite_element_method>(), t_IntegrationPoints> &
        getIntegrationPoints()
        {
            return std::get<getArgIndex<t_finite_element_method>()>(integrations_points_);
        }

        // template<BehaviorConcept auto t_behavior, Quadrature t_quadrature>
        // void
        // makeBehavior()
        // {
        //     for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
        //     {
        //         getIntegrationPoints<t_behavior>().push_back(std::make_shared<t_IntegrationPoint<getArgIndex<t_behavior>()>>());
        //     }
        //     auto jklm = [&] <lolita::integer t_i = 0> (
        //         auto & self
        //     )
        //     constexpr mutable
        //     {
        //         auto constexpr t_finite_element_method = std::get<t_i>(std::make_tuple(t_finite_element_methods...));
        //         // using t_Behavior = std::decay_t<decltype(t_behavior)>;
        //         // using t_FiniteElementBehavior = std::tuple_element_t<t_i, std::tuple<std::decay_t<decltype(t_finite_element_method)>...>>;
        //         if constexpr (std::is_same_v<std::decay_t<decltype(t_behavior)>, std::decay_t<decltype(t_finite_element_method.getBehavior())>>)
        //         {
        //             getElement<t_i>()->data_->quadrature_points_ = std::vector<std::unique_ptr<QuadraturePoint<t_element, t_domain, t_finite_element_method>>>();
        //             for (size_t i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
        //             {
        //                 getElement<t_i>()->data_->quadrature_points_[i] = std::make_unique<QuadraturePoint<t_element, t_domain, t_finite_element_method>>();
        //                 getElement<t_i>()->data_->quadrature_points_[i] = std::make_unique<QuadraturePoint<t_element, t_domain, t_finite_element_method>>();
        //             }
                    
        //             getElement<t_i>()->data_->quadrature_points_ = std::make_shared<QuadraturePoint<t_element, t_domain, ttt>>();
        //             getElement<t_i>()->data_->quadrature_points_ = getIntegrationPoints<t_behavior>();
        //         }
        //     };
        // }

        t_FiniteElements finite_elements_;

        t_IntegrationPoints integrations_points_;

    };
    
}


#endif /* ECF0EB38_E45D_4E38_ACEB_E40188311824 */
