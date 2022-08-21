#ifndef AB07387F_F0B7_4F5A_875F_1CC5F0206148
#define AB07387F_F0B7_4F5A_875F_1CC5F0206148

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_2000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_3000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4001.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4002.hxx"

namespace lolita
{

    template<Basis t_basis>
    requires(t_basis.isMonomial())
    struct FiniteElementBasisTraits<t_basis>
    {

        template<Integer t_dim>
        static constexpr
        Integer
        getSize()
        {
            return lolita::numerics::binomial(t_dim + t_basis.ord_, t_dim);
        }

        template<Element t_element>
        static constexpr
        Integer
        getSize()
        {
            return lolita::numerics::binomial(t_element.dim_ + t_basis.ord_, t_element.dim_);
        }

        template<Element t_element, Field t_field>
        static constexpr
        Integer
        getSize()
        {
            return lolita::numerics::binomial(t_element.dim_ + t_basis.ord_, t_element.dim_) * FieldTraits<t_field>::getSize();
        }
        
        template<Element t_element, Domain t_domain>
        struct Implementation : FiniteElementHolder<t_element, t_domain>
        {

            static constexpr
            Integer
            getSize()
            {
                return FiniteElementBasisTraits::template getSize<t_element>();
            }
            
            template<Field t_field>
            static constexpr
            Integer
            getSize()
            {
                return FiniteElementBasisTraits::template getSize<t_element, t_field>();
            }

        private:
        
            static constexpr
            std::array<std::array<Integer, 3>, getSize()>
            getExponents()
            {
                auto exponents = std::array<std::array<Integer, 3>, getSize()>();
                auto row = Integer(0);
                if constexpr (t_element.dim_ == 0)
                {
                    exponents[row][0] = 0;
                    exponents[row][1] = 0;
                    exponents[row][2] = 0;
                }
                else if constexpr (t_element.dim_ == 1)
                {
                    for (auto i = 0; i < t_basis.ord_ + 1; ++i)
                    {
                        exponents[row][0] = i;
                        exponents[row][1] = 0;
                        exponents[row][2] = 0;
                        row += 1;
                    }
                }
                else if constexpr (t_element.dim_ == 2)
                {
                    for (auto i = 0; i < t_basis.ord_ + 1; ++i)
                    {
                        for (auto j = 0; j < i + 1; ++j)
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
                    for (auto i = 0; i < t_basis.ord_ + 1; ++i)
                    {
                        for (auto j = 0; j < i + 1; ++j)
                        {
                            for (auto k = 0; k < i + 1; ++k)
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
            
            std::array<std::array<Integer, 3>, getSize()> static constexpr exponents_ = getExponents();

        public:
        
            lolita::algebra::Vector<Real, getSize()>
            getBasisEvaluation(
                Point const & point
            )
            const
            {
                auto basis_vector_values = lolita::algebra::Vector<Real, getSize()>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (auto i = 0; i < getSize(); ++i)
                {
                    auto value = Real(1);
                    for (auto j = 0; j < t_element.dim_; ++j)
                    {
                        auto dist = this->getRiemannianDistance(centroid, point, j);
                        value *= std::pow(2.0 * dist / diameters(j), exponents_[i][j]);
                    }
                    basis_vector_values(i) = value;
                }
                return basis_vector_values;
            }
            
            lolita::algebra::Vector<Real, getSize()>
            getBasisDerivative(
                Point const & point,
                Integer derivative_direction
            )
            const
            {
                auto basis_vector_values = lolita::algebra::Vector<Real, getSize()>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (auto i = 0; i < getSize(); ++i)
                {
                    auto value = Real(1);
                    for (auto j = 0; j < t_element.dim_; ++j)
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
        Integer
        getCellBasisSize()
        {
            return FiniteElementBasisTraits<t_discretization.cell_basis_>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        Integer
        getFaceBasisSize()
        {
            return FiniteElementBasisTraits<t_discretization.face_basis_>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        Integer
        getGradBasisSize()
        {
            return FiniteElementBasisTraits<t_discretization.grad_basis_>::template getSize<t_element>();
        }

        template<Element t_element, Domain t_domain, Field t_field>
        static constexpr
        Integer
        getNumCellUnknowns()
        {
            auto constexpr field_size = FieldTraits<t_field>::template getSize<t_domain>();
            auto constexpr basis_size = FiniteElementBasisTraits<t_discretization.cell_basis_>::template getSize<t_element>();
            return field_size * basis_size;
        }

        template<Element t_element, Domain t_domain, Field t_field>
        static constexpr
        Integer
        getNumFaceUnknowns()
        {
            auto constexpr field_size = FieldTraits<t_field>::template getSize<t_domain>();
            auto constexpr basis_size = FiniteElementBasisTraits<t_discretization.face_basis_>::template getSize<t_element>();
            return field_size * basis_size;
        }

        template<Element t_element, Domain t_domain, Field t_field>
        static constexpr
        Integer
        getNumElementUnknowns()
        {
            auto num_element_unknowns = getNumCellUnknowns<t_element, t_domain, t_field>();
            auto set_num_faces_unknowns = [&] <Integer t_i = 0> (
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

        template<Element t_element, Domain t_domain>
        struct Implementation : FiniteElementHolder<t_element, t_domain>
        {

        private:

            // template<Field t_field>
            // static constexpr
            // Integer
            // getNumCellUnknowns()
            // {
            //     return HybridDiscontinuousGalerkinTraits::template getNumCellUnknowns<t_element, t_domain, t_field>();
            // }

            // template<Field t_field>
            // static constexpr
            // Integer
            // getNumFaceUnknowns()
            // {
            //     return HybridDiscontinuousGalerkinTraits::template getNumFaceUnknowns<t_element, t_domain, t_field>();
            // }

            template<Mapping t_mapping, Field t_field>
            static constexpr
            Integer
            getMappingSize()
            {
                return MappingTraits<t_mapping>::template getSize<t_domain, t_field>();
            }

        public:

            template<Field t_field>
            static constexpr
            Integer
            getNumElementUnknowns()
            {
                return HybridDiscontinuousGalerkinTraits::template getNumElementUnknowns<t_element, t_domain, t_field>();
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

            template<Field t_field>
            RealMatrix<getGradBasisSize<t_element>(), getNumElementUnknowns<t_field>()>
            getGradientRhs(
                Integer row,
                Integer col
            )
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(2 * getGradBasis().getOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = RealMatrix<getGradBasisSize<t_element>(), getNumElementUnknowns<t_field>()>();
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
                auto set_faces_blocks = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0, t_i>();
                    auto face_index = Integer(0);
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

            template<Field t_field>
            RealMatrix<getGradBasisSize<t_element>(), getNumElementUnknowns<t_field>()>
            getSymmetricGradientRhs(
                Integer row,
                Integer col
            )
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(2 * getGradBasis().getOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = RealMatrix<getGradBasisSize<t_element>(), getNumElementUnknowns<t_field>()>();
                rhs.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::size(); i++)
                {
                    auto weight = this->template getReferenceQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto test_vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    auto left_vector = this->template getBasisDerivative<getCellBasis()>(point, col);
                    auto left_vector2 = this->template getBasisDerivative<getCellBasis()>(point, row);
                    auto col_offset = row * getCellBasisSize<t_element>();
                    auto col_offset2 = col * getCellBasisSize<t_element>();
                    auto block = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset);
                    auto block2 = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset2);
                    block += (1./2.) * weight * test_vector * left_vector.transpose();
                    block2 += (1./2.) * weight * test_vector * left_vector2.transpose();
                }
                auto set_faces_blocks = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0, t_i>();
                    auto face_index = Integer(0);
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

            template<Mapping t_mapping, Field t_field>
            Matrix<Real, getMappingSize<t_mapping, t_field>(), getNumElementUnknowns<t_field>()>
            getMapping(
                Point const & point
            )
            const
            requires(t_mapping.isGradient() || t_mapping.isLargeStrain())
            {
                auto mapping = Matrix<Real, getMappingSize<t_mapping, t_field>(), getNumElementUnknowns<t_field>()>();
                mapping.setZero();
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain, t_field>())
                {
                    auto rhs = getGradientRhs<t_field>(mapping_value.row(), mapping_value.col());
                    auto line = mapping.template block<1, getNumElementUnknowns<t_field>()>(mapping_value.rank(), 0);
                    line = mapping_value.value() * this->template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
                }
                return mapping;
            }

            template<Mapping t_mapping, Field t_field>
            Matrix<Real, getMappingSize<t_mapping, t_field>(), getNumElementUnknowns<t_field>()>
            getMapping(
                Point const & point
            )
            const
            requires(t_mapping.isSmallStrain())
            {
                auto mapping = Matrix<Real, getMappingSize<t_mapping, t_field>(), getNumElementUnknowns<t_field>()>();
                mapping.setZero();
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain, t_field>())
                {
                    auto rhs = getSymmetricGradientRhs<t_field>(mapping_value.row(), mapping_value.col());
                    auto line = mapping.template block<1, getNumElementUnknowns<t_field>()>(mapping_value.rank(), 0);
                    line = mapping_value.value() * this->template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
                }
                return mapping;
            }

            template<Mapping t_mapping, Field t_field>
            Matrix<Real, getMappingSize<t_mapping, t_field>(), getNumElementUnknowns<t_field>()>
            getMapping(
                Point const & point
            )
            const
            requires(t_mapping.isIdentity())
            {
                auto mapping = Matrix<Real, getMappingSize<t_mapping, t_field>(), getNumElementUnknowns<t_field>()>();
                mapping.setZero();
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain, t_field>())
                {
                    auto left_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                    auto col_offset = mapping_value.row() * getCellBasisSize<t_element>();
                    auto line = mapping.template block<1, getCellBasisSize<t_element>()>(mapping_value.rank(), col_offset);
                    line = mapping_value.value() * this->template getBasisEvaluation<getCellBasis()>(point);
                }
                return mapping;
            }

            template<Field t_field>
            Matrix<Real, getNumElementUnknowns<t_field>(), getNumElementUnknowns<t_field>()>
            getStabilization()
            const
            {
                auto stabilization = Matrix<Real, getNumElementUnknowns<t_field>(), getNumElementUnknowns<t_field>()>();
                stabilization.setIdentity();
                return stabilization;
            }

            template<Field t_field>
            Vector<Real, getNumElementUnknowns<t_field>()>
            getUnknowns(
                std::basic_string_view<Character> label
            )
            const
            {
                auto offset = 0;
                auto unknown = Vector<Real, getNumElementUnknowns<t_field>()>();
                auto const & cell_dof = this->degrees_of_freedom_.at(std::string(label));
                auto cell_block = unknown.template segment<cell_dof.template getSize<t_field, getCellBasis()>()>(offset);
                cell_block = cell_dof.template getCoefficients<t_field, getCellBasis()>();
                offset += cell_dof.template getSize<t_field, getCellBasis()>();
                auto set_faces_unknowns = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto const & face_dof = face->degrees_of_freedom_.at(std::string(label));
                        auto face_block = unknown.template segment<face_dof.template getSize<t_field, getFaceBasis()>()>(offset);
                        face_block = face_dof.template getCoefficients<t_field, getFaceBasis()>();
                        offset += face_dof.template getSize<t_field, getFaceBasis()>();
                    }
                    if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_faces_unknowns(set_faces_unknowns);
                return unknown;
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            Vector<Real, getNumElementUnknowns<t_finite_element_method.getField()>()>
            getInternalForces(
                std::basic_string_view<Character> behavior_label,
                std::basic_string_view<Character> degree_of_freedom_label
            )
            const
            {
                auto constexpr strain_operator_num_rows = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
                auto constexpr strain_operator_num_cols = getNumElementUnknowns<t_finite_element_method.getField()>();
                using StrainOperatorView = lolita::algebra::View<Matrix<Real, strain_operator_num_rows, strain_operator_num_cols> const>;
                using StrainView = lolita::algebra::View<Vector<Real, strain_operator_num_rows>>;
                auto unknown = this->template getUnknowns<t_finite_element_method.getField()>(degree_of_freedom_label);
                auto internal_forces = Vector<Real, strain_operator_num_cols>();
                internal_forces.setZero();
                for (auto & ip : this->quadrature_.at(std::string(behavior_label)).ips_)
                {
                    auto strain_operator_view = StrainOperatorView(ip.ops_.at(std::string(degree_of_freedom_label)).data());
                    auto stress_view = StrainView(ip.behavior_data_->s1.thermodynamic_forces.data());
                    internal_forces += ip.weight_ * strain_operator_view.transpose() * stress_view;
                }
                internal_forces += this->operators_.at("Stabilization") * unknown;
                return internal_forces;
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            Vector<Real, getNumElementUnknowns<t_finite_element_method.getField()>()>
            getExternalForces(
                std::basic_string_view<Character> behavior_label,
                std::basic_string_view<Character> degree_of_freedom_label
            )
            const
            {
                auto constexpr strain_operator_num_rows = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
                auto constexpr strain_operator_num_cols = getNumElementUnknowns<t_finite_element_method.getField()>();
                auto external_forces = Vector<Real, strain_operator_num_cols>();
                external_forces.setZero();
                return external_forces;
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            Matrix<Real, getNumElementUnknowns<t_finite_element_method.getField()>(), getNumElementUnknowns<t_finite_element_method.getField()>()>
            getJacobianMatrix(
                std::basic_string_view<Character> behavior_label,
                std::basic_string_view<Character> degree_of_freedom_label
            )
            const
            {
                auto constexpr strain_operator_num_rows = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
                auto constexpr strain_operator_num_cols = getNumElementUnknowns<t_finite_element_method.getField()>();
                using JacobianOperator = Matrix<Real, strain_operator_num_cols, strain_operator_num_cols>;
                using StrainOperatorView = lolita::algebra::View<Matrix<Real, strain_operator_num_rows, strain_operator_num_cols> const>;
                auto jacobian = Matrix<Real, strain_operator_num_cols, strain_operator_num_cols>();
                jacobian.setZero();
                for (auto & ip : this->quadrature_.at(std::string(behavior_label)).ips_)
                {
                    auto strain_operator_view = StrainOperatorView(ip.ops_.at(std::string(degree_of_freedom_label)).data());
                    auto jacob = ip.template getJacobian<t_finite_element_method>();
                    // auto jacob = Matrix<Real, strain_operator_num_rows, strain_operator_num_rows>();
                    jacobian += ip.weight_ * strain_operator_view.transpose() * jacob * strain_operator_view;
                }
                jacobian += this->operators_.at("Stabilization");
                return jacobian;
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            void
            assemble(
                std::basic_string_view<Character> behavior_label,
                std::basic_string_view<Character> degree_of_freedom_label,
                std::unique_ptr<System> const & system
            )
            const
            {
                auto constexpr num_cell_unknowns = getNumCellUnknowns<t_element, t_domain, t_finite_element_method.getField()>();
                auto constexpr strain_operator_num_cols = getNumElementUnknowns<t_finite_element_method.getField()>();
                auto constexpr num_face_unknowns = getNumElementUnknowns<t_finite_element_method.getField()>() - num_cell_unknowns;
                auto internal_forces = getInternalForces<t_finite_element_method>(behavior_label, degree_of_freedom_label);
                auto external_forces = getExternalForces<t_finite_element_method>(behavior_label, degree_of_freedom_label);
                auto jac = getJacobianMatrix<t_finite_element_method>(behavior_label, degree_of_freedom_label);
                auto residual = internal_forces - external_forces;
                auto k_tt = jac.template block<num_cell_unknowns, num_cell_unknowns>(0, 0);
                auto k_tf = jac.template block<num_cell_unknowns, num_face_unknowns>(0, num_cell_unknowns);
                auto k_ft = jac.template block<num_face_unknowns, num_cell_unknowns>(num_cell_unknowns, 0);
                auto k_ff = jac.template block<num_face_unknowns, num_face_unknowns>(num_cell_unknowns, num_cell_unknowns);
                auto r_t = residual.template segment<num_cell_unknowns>(0);
                auto r_f = residual.template segment<num_face_unknowns>(num_cell_unknowns);
                // auto k_tt_inv = k_tt.llt().solve(decltype(k_tt)::Identity());
                //
                //
                auto constexpr t_field = t_finite_element_method.getField();
                auto offset_i = 0;
                auto set_faces_unknowns = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    for (auto const & face_i : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto const & face_i_dof = face_i->degrees_of_freedom_.at(std::string(degree_of_freedom_label));
                        auto size_i = face_i_dof.template getSize<t_field, getFaceBasis()>();
                        auto index_i = face_i_dof.getTag();
                        for (auto i = index_i; i < index_i + size_i; i++)
                        {
                            system->addRhsValue(i, r_f(offset_i));
                            auto offset_j = 0;
                            auto set_faces_unknowns2 = [&] <Integer t_j = 0> (
                                auto & self2
                            )
                            constexpr mutable
                            {
                                for (auto const & face_j : this->template getInnerNeighbors<0, t_j>())
                                {
                                    auto const & face_j_dof = face_j->degrees_of_freedom_.at(std::string(degree_of_freedom_label));
                                    auto size_j = face_j_dof.template getSize<t_field, getFaceBasis()>();
                                    auto index_j = face_j_dof.getTag();
                                    for (auto j = index_j; j < index_j + size_j; j++)
                                    {
                                        system->addLhsValue(i, j, k_ff(offset_i, offset_j));
                                        offset_j ++;
                                    }
                                }
                                if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                                {
                                    self2.template operator ()<t_j + 1>(self2);
                                }
                            };
                            set_faces_unknowns2(set_faces_unknowns2);
                            offset_i ++;
                        }
                    }
                    if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_faces_unknowns(set_faces_unknowns);
            }

            // template<FiniteElementMethodConcept auto t_finite_element_method>
            // void
            // makeLoad(
            //     std::basic_string_view<Character> behavior_label,
            //     std::basic_string_view<Character> degree_of_freedom_label,
            //     std::function<Real(Point const &)> && function,
            //     std::unique_ptr<System> const & system
            // )
            // const
            // {
            //     auto tangent_tt = this->operators_.at("KTT");
            //     auto tangent_tt = this->operators_.at("KFT");
            //     auto tangent_tt = this->operators_.at("KTF");
            //     auto constexpr num_cell_unknowns = getNumCellUnknowns<t_element, t_domain, t_finite_element_method.getField()>();
            //     auto constexpr strain_operator_num_cols = getNumElementUnknowns<t_finite_element_method.getField()>();
            //     auto constexpr num_face_unknowns = getNumElementUnknowns<t_finite_element_method.getField()>() - num_cell_unknowns;
            //     auto internal_forces = getInternalForces<t_finite_element_method>(behavior_label, degree_of_freedom_label);
            //     auto external_forces = getExternalForces<t_finite_element_method>(behavior_label, degree_of_freedom_label);
            //     auto jac = getJacobianMatrix<t_finite_element_method>(behavior_label, degree_of_freedom_label);
            //     auto residual = internal_forces - external_forces;
            //     auto k_tt = jac.template block<num_cell_unknowns, num_cell_unknowns>(0, 0);
            //     auto k_tf = jac.template block<num_cell_unknowns, num_face_unknowns>(0, num_cell_unknowns);
            //     auto k_ft = jac.template block<num_face_unknowns, num_cell_unknowns>(num_cell_unknowns, 0);
            //     auto k_ff = jac.template block<num_face_unknowns, num_face_unknowns>(num_cell_unknowns, num_cell_unknowns);
            //     auto r_t = residual.template segment<num_cell_unknowns>(0);
            //     auto r_f = residual.template segment<num_face_unknowns>(num_cell_unknowns);
            //     // auto k_tt_inv = k_tt.llt().solve(decltype(k_tt)::Identity());
            //     //
            //     //
            //     auto constexpr t_field = t_finite_element_method.getField();
            //     auto offset_i = 0;
            //     auto set_faces_unknowns = [&] <Integer t_i = 0> (
            //         auto & self
            //     )
            //     constexpr mutable
            //     {
            //         for (auto const & face_i : this->template getInnerNeighbors<0, t_i>())
            //         {
            //             auto const & face_i_dof = face_i->degrees_of_freedom_.at(std::string(degree_of_freedom_label));
            //             auto size_i = face_i_dof.template getSize<t_field, getFaceBasis()>();
            //             auto index_i = face_i_dof.getTag();
            //             for (auto i = index_i; i < index_i + size_i; i++)
            //             {
            //                 system->addRhsValue(i, r_f(offset_i));
            //                 auto offset_j = 0;
            //                 auto set_faces_unknowns2 = [&] <Integer t_j = 0> (
            //                     auto & self2
            //                 )
            //                 constexpr mutable
            //                 {
            //                     for (auto const & face_j : this->template getInnerNeighbors<0, t_j>())
            //                     {
            //                         auto const & face_j_dof = face_j->degrees_of_freedom_.at(std::string(degree_of_freedom_label));
            //                         auto size_j = face_j_dof.template getSize<t_field, getFaceBasis()>();
            //                         auto index_j = face_j_dof.getTag();
            //                         for (auto j = index_j; j < index_j + size_j; j++)
            //                         {
            //                             system->addLhsValue(i, j, k_ff(offset_i, offset_j));
            //                             offset_j ++;
            //                         }
            //                     }
            //                     if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
            //                     {
            //                         self2.template operator ()<t_j + 1>(self2);
            //                     }
            //                 };
            //                 set_faces_unknowns2(set_faces_unknowns2);
            //                 offset_i ++;
            //             }
            //         }
            //         if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
            //         {
            //             self.template operator ()<t_i + 1>(self);
            //         }
            //     };
            //     set_faces_unknowns(set_faces_unknowns);
            // }

        };

    };

    // struct QuadraturePoint
    // {

    //     std::shared_ptr<mgis::behaviour::BehaviourData> const &
    //     getBehaviourData()
    //     const
    //     {
    //         return behavior_data_;
    //     }

    //     std::shared_ptr<mgis::behaviour::BehaviourData> &
    //     getBehaviourData()
    //     {
    //         return behavior_data_;
    //     }

    //     void
    //     integrate(
    //         Integer & res
    //     )
    //     {
    //         res = mgis::behaviour::integrate(* std::make_unique<mgis::behaviour::BehaviourDataView>(mgis::behaviour::make_view(* behavior_data_)), * behavior_);
    //     }
        
    //     std::shared_ptr<mgis::behaviour::Behaviour> behavior_;
        
    //     std::shared_ptr<mgis::behaviour::BehaviourData> behavior_data_;

    // };

    // struct QuadraturePoint2 : QuadraturePoint
    // {

    //     std::vector<lolita::algebra::Matrix<Real>> element_operators_;

    // };

    // struct QuadratureOperator
    // {

    //     inline
    //     std::basic_string_view<Character>
    //     getLabel()
    //     const
    //     {
    //         return label_;
    //     }

    //     inline
    //     lolita::algebra::Matrix<Real>
    //     getOperator(
    //         Integer i
    //     )
    //     const
    //     {
    //         return operators_[i];
    //     }

    //     std::basic_string_view<Character> label_;

    //     std::vector<lolita::algebra::Matrix<Real>> operators_;

    // };
    
    // template<Domain t_domain>
    // struct IntegrationPoint
    // {

    //     IntegrationPoint(
    //         Point const & coordinates,
    //         Real const & weight,
    //         std::shared_ptr<mgis::behaviour::Behaviour> const & behaviour
    //     )
    //     :
    //     coordinates_(coordinates),
    //     weight_(weight),
    //     behavior_data_(std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behaviour)))
    //     {}

    //     std::unique_ptr<mgis::behaviour::BehaviourData> const &
    //     getMaterialPoint()
    //     const
    //     {
    //         return behavior_data_;
    //     }

    //     std::unique_ptr<mgis::behaviour::BehaviourData> &
    //     getMaterialPoint()
    //     {
    //         return behavior_data_;
    //     }
        
    //     template<BehaviorConcept auto t_behavior>
    //     lolita::algebra::Span<lolita::algebra::Vector<Real, BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>()> const>
    //     getGeneralizedStrain()
    //     const
    //     {
    //         auto constexpr size = BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>();
    //         return lolita::algebra::Span<lolita::algebra::Vector<Real, size> const>(behavior_data_->s1.gradients.data());
    //     }
        
    //     template<auto t_finite_element_method>
    //     lolita::algebra::Span<RealVector<FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()> const>
    //     getGeneralizedStrain()
    //     const
    //     {
    //         auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
    //         auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
    //         return lolita::algebra::Span<lolita::algebra::Vector<Real, size> const>(behavior_data_->s1.gradients.data() + offset);
    //     }
        
    //     template<auto t_finite_element_method>
    //     lolita::algebra::Span<RealVector<FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()>>
    //     getGeneralizedStrain()
    //     {
    //         auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
    //         auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
    //         return lolita::algebra::Span<lolita::algebra::Vector<Real, size>>(behavior_data_->s1.gradients.data() + offset);
    //     }
        
    //     Point coordinates_;
        
    //     Real weight_;
        
    //     std::unique_ptr<mgis::behaviour::BehaviourData> behavior_data_;

    // };

    // template<Domain t_domain>
    // struct ElementIntegrationPoints
    // {

    //     ElementIntegrationPoints(
    //         std::shared_ptr<mgis::behaviour::Behaviour> const & behaviour
    //     )
    //     :
    //     behavior_(behaviour)
    //     {}

    //     Integer
    //     getSize()
    //     const
    //     {
    //         return integration_points_.size();
    //     }

    //     std::basic_string_view<Character>
    //     getLabel()
    //     const
    //     {
    //         return behavior_->behaviour;
    //     }

    //     void
    //     integrate()
    //     {
    //         for (auto const & integration_point : integration_points_)
    //         {
    //             auto behaviour_view = mgis::behaviour::make_view(integration_point->behavior_data_);
    //             auto result = mgis::behaviour::integrate(behaviour_view, * behavior_);
    //         }
    //     }
        
    //     std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

    //     std::vector<IntegrationPoint<t_domain>> integration_points_;

    // };
    
} // namespace lolita

#endif /* AB07387F_F0B7_4F5A_875F_1CC5F0206148 */
