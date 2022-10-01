#ifndef B536F722_1A61_4366_8C7E_4B005041A1BF
#define B536F722_1A61_4366_8C7E_4B005041A1BF

#include "lolita.hxx"
#include "core/000_physics_traits.hxx"
#include "core/100_geometry.hxx"
#include "core/linear_system.hxx"
#include "core/200_quadrature.hxx"
#include "core/201_finite_element_dof.hxx"
#include "core/202_finite_element_frm.hxx"
#include "core/300_finite_element.hxx"
#include "core/400_finite_element_basis.hxx"

namespace lolita
{

    template<auto t_discretization>
    struct DiscretizationTraits
    {

        static constexpr
        Basis
        getCellBasis()
        {
            return t_discretization.getCellBasis();
        }

        static constexpr
        Basis
        getFaceBasis()
        {
            return t_discretization.getFaceBasis();
        }

        static constexpr
        Basis
        getGradBasis()
        {
            return t_discretization.getGradBasis();
        }

        static constexpr
        Basis
        getPotentialBasis()
        {
            return t_discretization.getFaceBasis().toOrder(t_discretization.getFaceBasis().getOrd() + 1);
        }

        template<Element t_element>
        static constexpr
        Integer
        getCellBasisSize()
        {
            return BasisTraits<getCellBasis()>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        Integer
        getFaceBasisSize()
        {
            return BasisTraits<getFaceBasis()>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        Integer
        getGradBasisSize()
        {
            return BasisTraits<getGradBasis()>::template getSize<t_element>();
        }

        template<Element t_element>
        static constexpr
        Integer
        getPotentialBasisSize()
        {
            return BasisTraits<getPotentialBasis()>::template getSize<t_element>();
        }

        template<Element t_element, Domain t_domain, FieldConcept auto t_field>
        static constexpr
        Integer
        getNumCellUnknowns()
        requires(t_element.isSub(t_domain, 0))
        {
            return FieldTraits<t_field>::template getSize<t_domain>() * BasisTraits<t_discretization.getCellBasis()>::template getSize<t_element>();
        }

        template<Element t_element, Domain t_domain, FieldConcept auto t_field>
        static constexpr
        Integer
        getNumFaceUnknowns()
        requires(t_element.isSub(t_domain, 1))
        {
            return FieldTraits<t_field>::template getSize<t_domain>() * BasisTraits<t_discretization.getFaceBasis()>::template getSize<t_element>();
        }

        template<Element t_element, Domain t_domain, FieldConcept auto t_field>
        static constexpr
        Integer
        getNumElementUnknowns()
        requires(t_element.isSub(t_domain, 0))
        {
            auto num_element_unknowns = getNumCellUnknowns<t_element, t_domain, t_field>();
            auto set_num_faces_unknowns = [&] <Integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                auto constexpr num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                num_element_unknowns += getNumFaceUnknowns<inner_neighbor, t_domain, t_field>() * num_inner_neighbors;
                if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            set_num_faces_unknowns(set_num_faces_unknowns);
            return num_element_unknowns;
        }

        /**
         * @brief 
         * 
         * @tparam t_element 
         * @tparam t_domain 
         * @tparam t_field 
         */
        template<Element t_element, Domain t_domain, FieldConcept auto t_field>
        static constexpr
        Integer
        getDiscreteFieldDegreeOfFreedomSize()
        requires(t_element.isSub(t_domain, 0))
        {
            auto num_element_unknowns = getNumCellUnknowns<t_element, t_domain, t_field>();
            auto set_num_faces_unknowns = [&] <Integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                auto constexpr num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                num_element_unknowns += getNumFaceUnknowns<inner_neighbor, t_domain, t_field>() * num_inner_neighbors;
                if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            set_num_faces_unknowns(set_num_faces_unknowns);
            return num_element_unknowns;
        }

        /**
         * @brief 
         * 
         * @tparam t_element 
         * @tparam t_domain 
         * @tparam t_field 
         */
        template<Element t_element, Domain t_domain, FieldConcept auto t_field>
        static constexpr
        Integer
        getDiscreteFieldDegreeOfFreedomSize()
        requires(t_element.isSub(t_domain, 1))
        {
            return getNumFaceUnknowns<t_element, t_domain, t_field>();
        }

        template<Element t_element, Domain t_domain>
        struct Implementation : FiniteElement<t_element, t_domain>
        {

        private:

            using Base = FiniteElement<t_element, t_domain>;

            template<MappingConcept auto t_mapping, FieldConcept auto t_field>
            static constexpr
            Integer
            getMappingSize()
            {
                return MappingTraits<t_mapping>::template getSize<t_domain, t_field>();
            }

            template<MappingConcept auto t_mapping>
            static constexpr
            Integer
            getMappingSize()
            {
                return MappingTraits<t_mapping>::template getSize<t_domain>();
            }

            template<FieldConcept auto t_field>
            static constexpr
            Integer
            getNumCellUnknowns()
            requires(t_element.isSub(t_domain, 0))
            {
                return DiscretizationTraits::template getNumCellUnknowns<t_element, t_domain, t_field>();
            }

            template<FieldConcept auto t_field>
            static constexpr
            Integer
            getNumFaceUnknowns()
            requires(t_element.isSub(t_domain, 1))
            {
                return DiscretizationTraits::template getNumFaceUnknowns<t_element, t_domain, t_field>();
            }

            template<FieldConcept auto t_field>
            static constexpr
            Integer
            getNumElementUnknowns()
            requires(t_element.isSub(t_domain, 0))
            {
                return DiscretizationTraits::template getNumElementUnknowns<t_element, t_domain, t_field>();
            }

        public:

            template<FieldConcept auto t_field>
            static constexpr
            Integer
            getDiscreteFieldDegreeOfFreedomSize()
            {
                return DiscretizationTraits::template getDiscreteFieldDegreeOfFreedomSize<t_element, t_domain, t_field>();
            }

            template<FieldConcept auto t_field>
            constexpr
            Integer
            getDiscreteFieldDegreeOfFreedomNumCoefficients()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                return DiscretizationTraits::template getNumElementUnknowns<t_element, t_domain, t_field>();
            }

            template<FieldConcept auto t_field>
            constexpr
            Integer
            getDiscreteFieldDegreeOfFreedomNumCoefficients()
            const
            requires(t_element.isSub(t_domain, 1))
            {
                return DiscretizationTraits::template getNumFaceUnknowns<t_element, t_domain, t_field>();
            }

            static constexpr
            Integer
            getGradientConstructionQuadratureOrder()
            {
                if constexpr (t_domain.isAxiSymmetric())
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1, getGradBasis().getOrd()) + 1;
                }
                else
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1, getGradBasis().getOrd());
                }
            }

            static constexpr
            Integer
            getPotentialConstructionQuadratureOrder()
            {
                if constexpr (t_domain.isAxiSymmetric())
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1) + 1;
                }
                else
                {
                    return 2 * numerics::max(getCellBasis().getOrd(), getFaceBasis().getOrd() + 1);
                }
            }
            
            DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getPotentialBasisSize<t_element>() - 1>
            getPotentialLhs()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getPotentialConstructionQuadratureOrder());
                auto lhs = DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getPotentialBasisSize<t_element>() - 1>();
                lhs.setZero();
                for (auto i_component = 0; i_component < t_domain.getDim(); i_component++)
                {
                    for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                        auto vector = this->template getBasisDerivative<getPotentialBasis()>(point, i_component);
                        auto vector_j = vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                        lhs += weight * vector_j * vector_j.transpose();
                    }
                }
                return lhs.llt().solve(decltype(lhs)::Identity());
            }

            template<FieldConcept auto t_field>
            DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getDiscreteFieldDegreeOfFreedomSize<t_field>()>
            getPotentialRhs(
                Integer row
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getPotentialConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto rhs = DenseMatrix<Real, getPotentialBasisSize<t_element>() - 1, getDiscreteFieldDegreeOfFreedomSize<t_field>()>();
                rhs.setZero();
                for (auto i_component = 0; i_component < t_domain.getDim(); i_component++)
                {
                    auto face_offset = t_field_size * getCellBasisSize<t_element>();
                    for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                        auto row_cell_vector = this->template getBasisDerivative<getPotentialBasis()>(point, i_component);
                        auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                        auto col_cell_vector = this->template getBasisDerivative<getCellBasis()>(point, i_component);
                        auto col_offset = row * getCellBasisSize<t_element>();
                        auto block = rhs.template block<getPotentialBasisSize<t_element>() - 1, getCellBasisSize<t_element>()>(0, col_offset);
                        block += weight * row_cell_vector_j * col_cell_vector.transpose();
                    }
                    auto set_faces_blocks = [&] <Integer t_i = 0> (
                        auto & self
                    )
                    constexpr mutable
                    {
                        auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                        auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                        auto i_face = 0;
                        for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                        {
                            auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                            auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                            for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                            {
                                auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                                auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                                auto argf = i;
                                //
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                                auto normal_vector = face->getNormalVector(point);
                                //
                                auto row_cell_vector = this->template getBasisDerivative<getPotentialBasis()>(i_p, i_component);
                                auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                                auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                                auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                                auto normal_vector_component = face_orientation * normal_vector(i_component);
                                auto face_col_offset = face_offset + row * getFaceBasisSize<t_inner_neighbor>();
                                auto face_block = rhs.template block<getPotentialBasisSize<t_element>() - 1, getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset);
                                face_block += normal_vector_component * weight * row_cell_vector_j * col_face_vector.transpose();
                                //
                                auto cell_col_offset = row * getCellBasisSize<t_element>();
                                auto cell_block = rhs.template block<getPotentialBasisSize<t_element>() - 1, getCellBasisSize<t_element>()>(0, cell_col_offset);
                                cell_block -= normal_vector_component * weight * row_cell_vector_j * col_cell_vector.transpose();
                            }
                            face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                            i_face ++;
                        }
                        if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                        {
                            self.template operator()<t_i + 1>(self);
                        }
                    };
                    set_faces_blocks(set_faces_blocks);
                }
                return rhs;
            }

            template<FieldConcept auto t_field>
            DenseMatrix<Real, getDiscreteFieldDegreeOfFreedomSize<t_field>(), getDiscreteFieldDegreeOfFreedomSize<t_field>()>
            getStabilization()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getPotentialConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                //
                auto S_T_A_B_I_L_I_Z_A_T_I_O_N = DenseMatrix<Real, getDiscreteFieldDegreeOfFreedomSize<t_field>(), getDiscreteFieldDegreeOfFreedomSize<t_field>()>();
                auto potential_operator = DenseMatrix<Real, getPotentialBasisSize<t_element>(), getDiscreteFieldDegreeOfFreedomSize<t_field>()>();
                auto potential_lhs = getPotentialLhs();
                auto denom = 1.0 / QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize();
                S_T_A_B_I_L_I_Z_A_T_I_O_N.setZero();
                for (auto i_component = 0; i_component < FieldTraits<t_field>::template getSize<t_domain>(); i_component++)
                {
                    auto face_offset = t_field_size * getCellBasisSize<t_element>();
                    potential_operator.setZero();
                    auto oppp = potential_operator.template block<getPotentialBasisSize<t_element>() - 1, getDiscreteFieldDegreeOfFreedomSize<t_field>()>(1, 0);
                    oppp = potential_lhs * getPotentialRhs<t_field>(i_component);
                    for (auto i_quadrature = 0; i_quadrature < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i_quadrature++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                        auto row_cell_vector = this->template getBasisEvaluation<getPotentialBasis()>(point);
                        auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                        potential_operator.template block<1, getDiscreteFieldDegreeOfFreedomSize<t_field>()>(0, 0) -= denom * weight * row_cell_vector_j.transpose() * oppp;
                        potential_operator.template block<1, getCellBasisSize<t_element>()>(0, i_component * getCellBasisSize<t_element>()) +=
                        denom * weight * this->template getBasisEvaluation<getCellBasis()>(point).transpose();
                    }
                    //
                    auto cell_projection_lhs = DenseMatrix<Real, getCellBasisSize<t_element>(), getCellBasisSize<t_element>()>();
                    auto cell_projection_rhs = DenseMatrix<Real, getCellBasisSize<t_element>(), getDiscreteFieldDegreeOfFreedomSize<t_field>()>();
                    cell_projection_lhs.setZero();
                    cell_projection_rhs.setZero();
                    for (auto i_quadrature = 0; i_quadrature < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i_quadrature++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                        auto row_cell_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                        cell_projection_lhs += weight * row_cell_vector * row_cell_vector.transpose();
                        cell_projection_rhs += weight * row_cell_vector * this->template getBasisEvaluation<getPotentialBasis()>(point).transpose() * potential_operator;
                    }
                    auto proj = cell_projection_lhs.llt().solve(decltype(cell_projection_lhs)::Identity()) * cell_projection_rhs;
                    //
                    auto set_faces_blocks = [&] <Integer t_i = 0> (
                        auto & self
                    )
                    constexpr mutable
                    {
                        auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                        auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                        auto i_face = 0;
                        for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                        {
                            auto constexpr quadrature_size = QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize();
                            auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                            auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                            auto face_projection_lhs = DenseMatrix<Real, getFaceBasisSize<t_inner_neighbor>(), getFaceBasisSize<t_inner_neighbor>()>();
                            auto face_projection_rhs = DenseMatrix<Real, getFaceBasisSize<t_inner_neighbor>(), getDiscreteFieldDegreeOfFreedomSize<t_field>()>();
                            face_projection_lhs.setZero();
                            face_projection_rhs.setZero();
                            for (auto i_quadrature = 0; i_quadrature < quadrature_size; i_quadrature++)
                            {
                                auto argh = face_sign == 1 ? i_quadrature : quadrature_size - (i_quadrature + 1);
                                auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                                auto argf = i_quadrature;
                                //
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                                auto normal_vector = face->getNormalVector(point);
                                auto row_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                                auto col_potential_vector = this->template getBasisEvaluation<getPotentialBasis()>(i_p);
                                auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                                //
                                face_projection_lhs += weight * row_face_vector * row_face_vector.transpose();
                                face_projection_rhs -= weight * row_face_vector * col_potential_vector.transpose() * potential_operator;
                                face_projection_rhs += weight * row_face_vector * col_cell_vector.transpose() * proj;
                                //
                                auto cell_col_offset = i_component * getCellBasisSize<t_element>();
                                auto cell_block = face_projection_rhs.template block<getFaceBasisSize<t_inner_neighbor>(), getCellBasisSize<t_element>()>(0, cell_col_offset);
                                cell_block -= weight * row_face_vector * col_cell_vector.transpose();
                                //
                                auto face_col_offset = face_offset + i_component * getFaceBasisSize<t_inner_neighbor>();
                                auto face_block = face_projection_rhs.template block<getFaceBasisSize<t_inner_neighbor>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset);
                                face_block += weight * row_face_vector * row_face_vector.transpose();
                            }
                            auto fa_ce_proj = face_projection_lhs.llt().solve(decltype(face_projection_lhs)::Identity()) * face_projection_rhs;
                            for (auto i_quadrature = 0; i_quadrature < quadrature_size; i_quadrature++)
                            {
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                                auto row_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                                auto V_E_C_T_O_R = row_face_vector.transpose() * fa_ce_proj;
                                auto factor = (1.0 / face->getLocalFrameDiameters().norm());
                                S_T_A_B_I_L_I_Z_A_T_I_O_N += factor * weight * V_E_C_T_O_R.transpose() * V_E_C_T_O_R;
                            }
                            face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                            i_face ++;
                        }
                        if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                        {
                            self.template operator()<t_i + 1>(self);
                        }
                    };
                    set_faces_blocks(set_faces_blocks);
                }
                return S_T_A_B_I_L_I_Z_A_T_I_O_N;
            }
            
            DenseMatrix<Real, getGradBasisSize<t_element>(), getGradBasisSize<t_element>()>
            getGradientLhs()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getGradientConstructionQuadratureOrder());
                auto lhs = DenseMatrix<Real, getGradBasisSize<t_element>(), getGradBasisSize<t_element>()>();
                lhs.setZero();
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    lhs += weight * vector * vector.transpose();
                }
                return lhs.llt().solve(decltype(lhs)::Identity());
            }

            template<FieldConcept auto t_field>
            DenseMatrix<Real, getGradBasisSize<t_element>(), getDiscreteFieldDegreeOfFreedomSize<t_field>()>
            getGradientRhs(
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getGradientConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = DenseMatrix<Real, getGradBasisSize<t_element>(), getDiscreteFieldDegreeOfFreedomSize<t_field>()>();
                rhs.setZero();
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    auto col_cell_vector = this->template getBasisDerivative<getCellBasis()>(point, col);
                    auto col_offset = row * getCellBasisSize<t_element>();
                    auto block = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset);
                    block += weight * row_cell_vector * col_cell_vector.transpose();
                }
                auto set_faces_blocks = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                    auto i_face = 0;
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                        auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                        for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                        {
                            auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                            auto argf = i;
                            auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            //
                            auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                            auto normal_vector = face->getNormalVector(point);
                            auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(i_p);
                            auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                            auto normal_vector_component = face_orientation * normal_vector(col);
                            auto face_col_offset = face_offset + row * getFaceBasisSize<t_inner_neighbor>();
                            auto face_block = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset);
                            face_block += normal_vector_component * weight * row_cell_vector * col_face_vector.transpose();
                            auto cell_col_offset = row * getCellBasisSize<t_element>();
                            auto cell_block = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, cell_col_offset);
                            cell_block -= normal_vector_component * weight * row_cell_vector * col_cell_vector.transpose();
                        }
                        face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                        i_face ++;
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator()<t_i + 1>(self);
                    }
                };
                set_faces_blocks(set_faces_blocks);
                return rhs;
            }

            template<FieldConcept auto t_field>
            DenseMatrix<Real, getGradBasisSize<t_element>(), getDiscreteFieldDegreeOfFreedomSize<t_field>()>
            getSymmetricGradientRhs(
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_quadrature = Quadrature("Gauss", getGradientConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = DenseMatrix<Real, getGradBasisSize<t_element>(), getDiscreteFieldDegreeOfFreedomSize<t_field>()>();
                rhs.setZero();
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    auto col_cell_vector_i = this->template getBasisDerivative<getCellBasis()>(point, col);
                    auto col_cell_vector_j = this->template getBasisDerivative<getCellBasis()>(point, row);
                    auto col_offset_i = row * getCellBasisSize<t_element>();
                    auto block_i = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset_i);
                    block_i += (1.0 / 2.0) * weight * row_cell_vector * col_cell_vector_i.transpose();
                    auto col_offset_j = col * getCellBasisSize<t_element>();
                    auto block_j = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, col_offset_j);
                    block_j += (1.0 / 2.0) * weight * row_cell_vector * col_cell_vector_j.transpose();
                }
                auto set_faces_blocks = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<0, t_i>();
                    auto i_face = 0;
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                        auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                        for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize(); i++)
                        {
                            auto argh = face_sign == 1 ? i : QuadratureTraits<t_quadrature>::template Rule<t_inner_neighbor>::getSize() - (i + 1);
                            auto argf = i;
                            auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                            auto normal_vector = face->getNormalVector(point);
                            auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(i_p);
                            auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                            auto normal_vector_component_i = face_orientation * normal_vector(col);
                            auto normal_vector_component_j = face_orientation * normal_vector(row);
                            auto face_col_offset_i = face_offset + row * getFaceBasisSize<t_inner_neighbor>();
                            auto face_block_i = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset_i);
                            face_block_i += (1.0 / 2.0) * normal_vector_component_i * weight * row_cell_vector * col_face_vector.transpose();
                            //
                            auto face_col_offset_j = face_offset + col * getFaceBasisSize<t_inner_neighbor>();
                            auto face_block_j = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset_j);
                            face_block_j += (1.0 / 2.0) * normal_vector_component_j * weight * row_cell_vector * col_face_vector.transpose();
                            //
                            auto cell_col_offset_i = row * getCellBasisSize<t_element>();
                            auto cell_block_i = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, cell_col_offset_i);
                            cell_block_i -= (1.0 / 2.0) * normal_vector_component_i * weight * row_cell_vector * col_cell_vector.transpose();
                            //
                            auto cell_col_offset_j = col * getCellBasisSize<t_element>();
                            auto cell_block_j = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, cell_col_offset_j);
                            cell_block_j -= (1.0 / 2.0) * normal_vector_component_j * weight * row_cell_vector * col_cell_vector.transpose();
                        }
                        face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                        i_face ++;
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator()<t_i + 1>(self);
                    }
                };
                set_faces_blocks(set_faces_blocks);
                return rhs;
            }
            
            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getDiscreteFieldDegreeOfFreedomSize<t_mapping.getField()>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.isGradient() || t_mapping.isLargeStrain())
            {
                auto mapping = DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getDiscreteFieldDegreeOfFreedomSize<t_mapping.getField()>()>();
                mapping.setZero();
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain>())
                {
                    auto rhs = getGradientRhs<t_mapping.getField()>(mapping_value.row(), mapping_value.col());
                    auto line = mapping.template block<1, getDiscreteFieldDegreeOfFreedomSize<t_mapping.getField()>()>(mapping_value.rank(), 0);
                    line = mapping_value.value() * this->template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
                }
                return mapping;
            }

            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getDiscreteFieldDegreeOfFreedomSize<t_mapping.getField()>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.isSmallStrain())
            {
                auto mapping = DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getDiscreteFieldDegreeOfFreedomSize<t_mapping.getField()>()>();
                mapping.setZero();
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain>())
                {
                    auto rhs = getSymmetricGradientRhs<t_mapping.getField()>(mapping_value.row(), mapping_value.col());
                    auto line = mapping.template block<1, getDiscreteFieldDegreeOfFreedomSize<t_mapping.getField()>()>(mapping_value.rank(), 0);
                    line = mapping_value.value() * this->template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
                }
                return mapping;
            }

            template<MappingConcept auto t_mapping>
            DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getDiscreteFieldDegreeOfFreedomSize<t_mapping.getField()>()>
            getMapping(
                PointConcept auto const & point
            )
            const
            requires(t_mapping.isIdentity())
            {
                auto mapping = DenseMatrix<Real, MappingTraits<t_mapping>::template getSize<t_domain>(), getDiscreteFieldDegreeOfFreedomSize<t_mapping.getField()>()>();
                mapping.setZero();
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain>())
                {
                    auto left_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                    auto col_offset = mapping_value.row() * getCellBasisSize<t_element>();
                    auto line = mapping.template block<1, getCellBasisSize<t_element>()>(mapping_value.rank(), col_offset);
                    line = mapping_value.value() * this->template getBasisEvaluation<getCellBasis()>(point);
                }
                return mapping;
            }

            // template<FieldConcept auto t_field>
            // DenseVector<Real, getDiscreteFieldDegreeOfFreedomSize<t_field>()>
            // getUnknowns()
            // const
            // {
            //     auto offset = 0;
            //     auto unknown = DenseVector<Real, getDiscreteFieldDegreeOfFreedomSize<t_field>()>();
            //     auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
            //     auto cell_block = unknown.template segment<cell_dof.template getSize<t_element, t_field, getCellBasis()>()>(offset);
            //     cell_block = cell_dof.template getCoefficients<t_element, t_field, getCellBasis()>();
            //     offset += cell_dof.template getSize<t_element, t_field, getCellBasis()>();
            //     auto set_faces_unknowns = [&] <Integer t_i = 0> (
            //         auto & self
            //     )
            //     constexpr mutable
            //     {
            //         auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
            //         for (auto const & face : this->template getInnerNeighbors<0, t_i>())
            //         {
            //             auto const & face_dof = face->template getDiscreteField<t_field>().getDegreeOfFreedom();
            //             auto face_block = unknown.template segment<face_dof.template getSize<t_inner_neighbor, t_field, getFaceBasis()>()>(offset);
            //             face_block = face_dof.template getCoefficients<t_inner_neighbor, t_field, getFaceBasis()>();
            //             offset += face_dof.template getSize<t_inner_neighbor, t_field, getFaceBasis()>();
            //         }
            //         if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
            //         {
            //             self.template operator ()<t_i + 1>(self);
            //         }
            //     };
            //     set_faces_unknowns(set_faces_unknowns);
            //     return unknown;
            // }

            template<FieldConcept auto t_field>
            DenseVector<Real, Base::template getDiscreteFieldDegreeOfFreedomSize<t_field, t_discretization>()>
            getDiscreteFieldDegreeOfFreedomCoefficients()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto offset = 0;
                auto unknown = DenseVector<Real, Base::template getDiscreteFieldDegreeOfFreedomSize<t_field, t_discretization>()>();
                auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
                auto cell_block = unknown.template segment<cell_dof.template getSize<t_element, t_field, getCellBasis()>()>(offset);
                cell_block = cell_dof.template getCoefficients<t_element, t_field, getCellBasis()>();
                offset += cell_dof.template getSize<t_element, t_field, getCellBasis()>();
                auto set_faces_unknowns = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<0, t_i>();
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto const & face_dof = face->template getDiscreteField<t_field>().getDegreeOfFreedom();
                        auto face_block = unknown.template segment<face_dof.template getSize<t_inner_neighbor, t_field, getFaceBasis()>()>(offset);
                        face_block = face_dof.template getCoefficients<t_inner_neighbor, t_field, getFaceBasis()>();
                        offset += face_dof.template getSize<t_inner_neighbor, t_field, getFaceBasis()>();
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_faces_unknowns(set_faces_unknowns);
                return unknown;
            }

            template<FieldConcept auto t_field>
            void
            addDiscreteFieldDegreeOfFreedom()
            requires(t_element.isSub(t_domain, 0))
            {
                static_cast<Base *>(this)->template addDiscreteFieldDegreeOfFreedom<t_field, getCellBasis()>();
            }

            template<FieldConcept auto t_field, Strategy t_s>
            void
            addDiscreteFieldDegreeOfFreedom(
                std::unique_ptr<LinearSystem<t_s>> const & linear_system
            )
            requires(t_element.isSub(t_domain, 1))
            {
                static_cast<Base *>(this)->template addDiscreteFieldDegreeOfFreedom<t_field, getFaceBasis()>(linear_system);
            }

            template<FieldConcept auto t_field, Label t_label>
            void
            addDiscreteFieldOperator()
            requires(t_label == "Stabilization")
            {
                this->template getDiscreteField<t_field>().addMatrix(t_label, getStabilization<t_field>());
            }

            template<FieldConcept auto t_field>
            void
            upgradeDiscreteFieldDegreeOfFreedom(
                DenseVectorConcept<Real> auto const & increment
            )
            {}

            template<FieldConcept auto t_field>
            void
            upgradeDiscreteFieldDegreeOfFreedom(
                DenseVectorConcept<Real> auto const & increment
            )
            requires(t_element.isSub(t_domain, 0))
            {
                auto constexpr num_cell_unknowns = getNumCellUnknowns<t_field>();
                auto constexpr strain_operator_num_cols = getDiscreteFieldDegreeOfFreedomSize<t_field>();
                auto constexpr num_face_unknowns = getDiscreteFieldDegreeOfFreedomSize<t_field>() - num_cell_unknowns;
                auto faces_correction = DenseVector<Real, num_face_unknowns>();
                faces_correction.setZero();
                auto faces_correction_offset = 0;
                auto set_faces_increment = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto const & face_dof = face->template getDiscreteField<t_field>().getDegreeOfFreedom();
                        auto lhs = faces_correction.template segment<face_dof.template getSize<t_field, getFaceBasis()>()>(faces_correction_offset);
                        auto rhs = increment.template segment<face_dof.template getSize<t_field, getFaceBasis()>()>(face_dof.getOffset());
                        lhs = rhs;
                        faces_correction_offset += face_dof.template getSize<t_field, getFaceBasis()>();
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_faces_increment(set_faces_increment);
                auto k_tt_inv = this->template getDiscreteField<t_field>().getMatrix("KTT");
                auto k_tf = this->template getDiscreteField<t_field>().getMatrix("KTF");
                auto r_t = this->template getDiscreteField<t_field>().getVector("RT");
                auto cell_corr = k_tt_inv * (r_t - k_tf * faces_correction);
                this->template getDiscreteField<t_field>().getDegreeOfFreedom().addCoefficients(cell_corr);
            }

            template<FieldConcept auto t_field>
            void
            upgradeDiscreteFieldDegreeOfFreedom(
                DenseVectorConcept<Real> auto const & increment
            )
            requires(t_element.isSub(t_domain, 1))
            {
                static_cast<Base *>(this)->template upgradeDiscreteFieldDegreeOfFreedom<t_field, getFaceBasis()>(increment);
            }

            template<FieldConcept auto t_field>
            Real
            getDiscreteFieldDegreeOfFreedomValue(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            {
                return 0.0;
            }

            template<FieldConcept auto t_field>
            Real
            getDiscreteFieldDegreeOfFreedomValue(
                PointConcept auto const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                return static_cast<Base *>(this)->template getDiscreteFieldDegreeOfFreedomValue<t_field, getCellBasis()>(point, row, col);
            }

            template<FieldConcept auto t_field>
            Real
            getDiscreteFieldDegreeOfFreedomValue(
                Point const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                return static_cast<Base *>(this)->template getDiscreteFieldDegreeOfFreedomValue<t_field, getFaceBasis()>(point, row, col);
            }

            template<FieldConcept auto t_field, Quadrature t_quadrature>
            Real
            getDiscreteFieldDegreeOfFreedomIntegrationValue(
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                this->template getDiscreteFieldDegreeOfFreedomValue<t_field>(row, col);
                auto value = Real(0);
                for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
                {
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    value += weight * this->template getDiscreteFieldDegreeOfFreedomValue<t_field>(point, row, col);
                }
                return value;
            }

            template<PotentialConcept auto t_behavior>
            void
            hhh()
            const
            {

            }

            template<PotentialConcept auto t_behavior, MappingConcept auto t_strain>
            DenseVector<Real, getDiscreteFieldDegreeOfFreedomSize<t_strain.getField()>()>
            getInternalForces()
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto constexpr t_field = t_strain.getField();
                auto constexpr strain_operator_num_rows = MappingTraits<t_strain>::template getSize<t_domain>();
                auto constexpr strain_operator_num_cols = getDiscreteFieldDegreeOfFreedomSize<t_field>();
                // using StrainOperatorView = algebra::View<DenseMatrix<Real, strain_operator_num_rows, strain_operator_num_cols> const>;
                using StrainView = algebra::View<DenseVector<Real, strain_operator_num_rows>>;
                auto const unknown = this->template getDiscreteFieldDegreeOfFreedomCoefficients<t_field>();
                auto internal_forces = DenseVector<Real, strain_operator_num_cols>();
                internal_forces.setZero();
                for (auto const & ip : this->template getFormulation<t_behavior>().getIntegrationPoints())
                {
                    // auto strain_operator_view = StrainOperatorView(ip.template getStrainOperator<t_strain>().data());
                    auto stress_view = StrainView(ip.behavior_data_->s1.thermodynamic_forces.data());
                    internal_forces += ip.getCurrentWeight() * ip.template getStrainOperator<t_strain>().transpose() * stress_view;
                }
                auto const & s_o = this->template getDiscreteField<t_field>().getMatrix("Stabilization");
                auto const & s_p = this->template getDiscreteField<t_field>().getScalar("Stabilization");
                internal_forces += s_p * s_o * unknown;
                return internal_forces;
            }

            template<PotentialConcept auto t_behavior, MappingConcept auto t_strain>
            void
            setInternalForces(
                DenseVector<Real> & internal_forces,
                Integer & offset
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto constexpr t_field = t_strain.getField();
                auto constexpr strain_operator_num_rows = MappingTraits<t_strain>::template getSize<t_domain>();
                auto constexpr strain_operator_num_cols = getDiscreteFieldDegreeOfFreedomSize<t_field>();
                using StrainView = algebra::View<DenseVector<Real, strain_operator_num_rows>>;
                auto internal_forces_block = internal_forces.template segment<getDiscreteFieldDegreeOfFreedomSize<t_field>()>(offset);
                auto const unknown = this->template getDiscreteFieldDegreeOfFreedomCoefficients<t_field>();
                for (auto const & ip : this->template getFormulation<t_behavior>().getIntegrationPoints())
                {
                    auto stress_view = StrainView(ip.behavior_data_->s1.thermodynamic_forces.data());
                    internal_forces_block += ip.getCurrentWeight() * ip.template getStrainOperator<t_strain>().transpose() * stress_view;
                }
                offset += getDiscreteFieldDegreeOfFreedomSize<t_field>();
            }

            template<PotentialConcept auto t_behavior, MappingConcept auto t_strain>
            void
            setStabilizationForces(
                DenseVector<Real> & internal_forces,
                Integer & offset
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto constexpr t_field = t_strain.getField();
                auto constexpr strain_operator_num_cols = getDiscreteFieldDegreeOfFreedomSize<t_field>();
                auto internal_forces_block = internal_forces.template segment<strain_operator_num_cols>(offset);
                auto const unknown = this->template getDiscreteFieldDegreeOfFreedomCoefficients<t_field>();
                auto const & s_o = this->template getDiscreteField<t_field>().getMatrix("Stabilization");
                auto const & s_p = this->template getDiscreteField<t_field>().getScalar("Stabilization");
                internal_forces_block += s_p * s_o * unknown;
            }

        };

    };
    
} // namespace lolita


#endif /* B536F722_1A61_4366_8C7E_4B005041A1BF */
