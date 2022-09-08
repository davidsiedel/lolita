#ifndef C17BAFB5_3EEB_4E27_8A93_F8EEC5AAF622
#define C17BAFB5_3EEB_4E27_8A93_F8EEC5AAF622

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_2000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_3000.hxx"

namespace lolita
{

    template<Element t_element, Domain t_domain>
    struct FiniteElementHolder;
        
    template<Basis t_basis>
    struct FiniteElementBasisTraits;

    template<auto t_discretization>
    struct HybridDiscontinuousGalerkinTraits;

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
        
        template<Element t_element, Domain t_domain>
        struct Implementation : FiniteElementHolder<t_element, t_domain>
        {

            static constexpr
            Integer
            getSize()
            {
                return FiniteElementBasisTraits::template getSize<t_element>();
            }
            
            // template<Field t_field>
            // static constexpr
            // Integer
            // getSize()
            // {
            //     return FiniteElementBasisTraits::template getSize<t_element, t_field>();
            // }

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
        
            // Real
            // getLocalFrameDistance1(
            //     Point const & first_point,
            //     Point const & second_point,
            //     Integer kkk
            // )
            // const
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     // auto first_point_mapping = Point();
            //     // auto second_point_mapping = Point();
            //     // auto const & current_coordinates = this->getCurrentCoordinates();
            //     // first_point_mapping.setZero();
            //     // second_point_mapping.setZero();
            //     // for (auto i = 0; i < t_domain.getDim(); ++i)
            //     // {
            //     //     first_point_mapping(i) = FiniteElementHolder::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
            //     //     second_point_mapping(i) = FiniteElementHolder::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
            //     // }
            //     // return (second_point_mapping - first_point_mapping)(kkk);
            //     return (second_point - first_point)(kkk);
            // }
            
            // Real
            // getLocalFrameDistance1(
            //     Point const & first_point,
            //     Point const & second_point,
            //     Integer kkk
            // )
            // const
            // requires(t_element.isSub(t_domain, 1))
            // {
            //     auto rotation_matrix = this->getRotationMatrix(this->getReferenceCentroid());
            //     // auto first_point_mapping = Point();
            //     // auto second_point_mapping = Point();
            //     // auto const & current_coordinates = this->getCurrentCoordinates();
            //     // first_point_mapping.setZero();
            //     // second_point_mapping.setZero();
            //     // for (auto i = 0; i < t_domain.getDim(); ++i)
            //     // {
            //     //     first_point_mapping(i) = FiniteElementHolder::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
            //     //     second_point_mapping(i) = FiniteElementHolder::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
            //     // }
            //     // // return (rotation_matrix * (second_point_mapping - first_point_mapping))(kkk);
            //     return (rotation_matrix * (second_point - first_point))(kkk);
            // }

        public:
        
            lolita::algebra::Vector<Real, getSize()>
            getBasisEvaluation(
                Point const & point
            )
            const
            {
                auto basis_vector_values = lolita::algebra::Vector<Real, getSize()>();
                auto const centroid = this->getReferenceCentroid();
                // auto const centroid = this->getCurrentCentroid();
                // auto const diameters = this->getCurrentDiameters();
                auto const diameters = this->getLocalFrameDiameters();
                for (auto i = 0; i < getSize(); ++i)
                {
                    auto value = Real(1);
                    for (auto j = 0; j < t_element.dim_; ++j)
                    {
                        // auto dist = this->getRiemannianDistance(centroid, point, j);
                        auto dist = this->getLocalFrameDistance(centroid, point, j);
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
                // auto const centroid = this->getCurrentCentroid();
                // auto const diameters = this->getCurrentDiameters();
                auto const diameters = this->getLocalFrameDiameters();
                for (auto i = 0; i < getSize(); ++i)
                {
                    auto value = Real(1);
                    for (auto j = 0; j < t_element.dim_; ++j)
                    {
                        if (j != derivative_direction)
                        {
                            // auto dist = this->getRiemannianDistance(centroid, point, j);
                            auto dist = this->getLocalFrameDistance(centroid, point, j);
                            value *= std::pow(2.0 * (dist) / diameters(j), exponents_[i][j]);
                        }
                        else
                        {
                            if (exponents_[i][j] > 0)
                            {
                                auto c = 2.0 * exponents_[i][j] / diameters(j);
                                // auto dist = this->getRiemannianDistance(centroid, point, j);
                                auto dist = this->getLocalFrameDistance(centroid, point, j);
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

        private:
            
            std::array<std::array<Integer, 3>, getSize()> static constexpr exponents_ = getExponents();

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

        static constexpr
        Basis
        getPotentialBasis()
        {
            return Basis::monomial(t_discretization.face_basis_.getOrd() + 1);
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

        template<Element t_element>
        static constexpr
        Integer
        getPotentialBasisSize()
        {
            return FiniteElementBasisTraits<getPotentialBasis()>::template getSize<t_element>();
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
            getNumCellUnknowns()
            requires(t_element.isSub(t_domain, 0))
            {
                // auto constexpr field_size = FieldTraits<t_field>::template getSize<t_domain>();
                // auto constexpr basis_size = FiniteElementBasisTraits<t_discretization.cell_basis_>::template getSize<t_element>();
                // return field_size * basis_size;
                return HybridDiscontinuousGalerkinTraits::template getNumCellUnknowns<t_element, t_domain, t_field>();
            }

            template<Field t_field>
            static constexpr
            Integer
            getNumFaceUnknowns()
            requires(t_element.isSub(t_domain, 1))
            {
                // auto constexpr field_size = FieldTraits<t_field>::template getSize<t_domain>();
                // auto constexpr basis_size = FiniteElementBasisTraits<t_discretization.face_basis_>::template getSize<t_element>();
                // return field_size * basis_size;
                return HybridDiscontinuousGalerkinTraits::template getNumFaceUnknowns<t_element, t_domain, t_field>();
            }

            template<Field t_field>
            static constexpr
            Integer
            getNumElementUnknowns()
            requires(t_element.isSub(t_domain, 0))
            {
                return HybridDiscontinuousGalerkinTraits::template getNumElementUnknowns<t_element, t_domain, t_field>();
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
            
            RealMatrix<getGradBasisSize<t_element>(), getGradBasisSize<t_element>()>
            getGradientLhs()
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto lhs = RealMatrix<getGradBasisSize<t_element>(), getGradBasisSize<t_element>()>();
                lhs.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    // auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i);
                    auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                    auto vector = this->template getBasisEvaluation<getGradBasis()>(point);
                    lhs += weight * vector * vector.transpose();
                }
                return lhs.llt().solve(decltype(lhs)::Identity());
            }
            
            RealMatrix<getPotentialBasisSize<t_element>() - 1, getPotentialBasisSize<t_element>() - 1>
            getPotentialLhs()
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto lhs = RealMatrix<getPotentialBasisSize<t_element>() - 1, getPotentialBasisSize<t_element>() - 1>();
                lhs.setZero();
                for (auto i_component = 0; i_component < t_domain.getDim(); i_component++)
                {
                    for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                        // auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i);
                        auto vector = this->template getBasisDerivative<getPotentialBasis()>(point, i_component);
                        auto vector_j = vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                        lhs += weight * vector_j * vector_j.transpose();
                    }
                }
                return lhs.llt().solve(decltype(lhs)::Identity());
            }

            template<Field t_field>
            RealMatrix<getPotentialBasisSize<t_element>() - 1, getNumElementUnknowns<t_field>()>
            getPotentialRhs(
                Integer row
            )
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto rhs = RealMatrix<getPotentialBasisSize<t_element>() - 1, getNumElementUnknowns<t_field>()>();
                rhs.setZero();
                for (auto i_component = 0; i_component < t_domain.getDim(); i_component++)
                {
                    auto face_offset = t_field_size * getCellBasisSize<t_element>();
                    for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                        // auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i);
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
                        auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
                        auto constexpr t_num_inner_neighbors = ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0, t_i>();
                        auto i_face = 0;
                        for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                        {
                            auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                            auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                            for (auto i = 0; i < ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize(); i++)
                            {
                                auto argh = face_sign == 1 ? i : ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize() - (i + 1);
                                auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                                auto argf = i;
                                //
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                                // auto point = face->template getCurrentQuadraturePoint<t_quadrature>(i);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                                // auto normal_vector = face->getNormalVector(face->template getReferenceQuadraturePoint<t_quadrature>(i));
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
                        if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                        {
                            self.template operator()<t_i + 1>(self);
                        }
                    };
                    set_faces_blocks(set_faces_blocks);
                }
                return rhs;
            }

            template<Field t_field>
            Matrix<Real, getNumElementUnknowns<t_field>(), getNumElementUnknowns<t_field>()>
            getStabilization()
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(getPotentialConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                //
                auto S_T_A_B_I_L_I_Z_A_T_I_O_N = Matrix<Real, getNumElementUnknowns<t_field>(), getNumElementUnknowns<t_field>()>();
                auto potential_operator = Matrix<Real, getPotentialBasisSize<t_element>(), getNumElementUnknowns<t_field>()>();
                auto potential_lhs = getPotentialLhs();
                auto denom = 1.0 / ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize();
                S_T_A_B_I_L_I_Z_A_T_I_O_N.setZero();
                for (auto i_component = 0; i_component < FieldTraits<t_field>::template getSize<t_domain>(); i_component++)
                {
                    auto face_offset = t_field_size * getCellBasisSize<t_element>();
                    potential_operator.setZero();
                    auto oppp = potential_operator.template block<getPotentialBasisSize<t_element>() - 1, getNumElementUnknowns<t_field>()>(1, 0);
                    oppp = potential_lhs * getPotentialRhs<t_field>(i_component);
                    for (auto i_quadrature = 0; i_quadrature < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i_quadrature++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                        // auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i_quadrature);
                        auto point = this->template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                        auto row_cell_vector = this->template getBasisEvaluation<getPotentialBasis()>(point);
                        auto row_cell_vector_j = row_cell_vector.template segment<getPotentialBasisSize<t_element>() - 1>(1);
                        potential_operator.template block<1, getNumElementUnknowns<t_field>()>(0, 0) -= denom * weight * row_cell_vector_j.transpose() * oppp;
                        potential_operator.template block<1, getCellBasisSize<t_element>()>(0, i_component * getCellBasisSize<t_element>()) +=
                        denom * weight * this->template getBasisEvaluation<getCellBasis()>(point).transpose();
                    }
                    //
                    auto cell_projection_lhs = Matrix<Real, getCellBasisSize<t_element>(), getCellBasisSize<t_element>()>();
                    auto cell_projection_rhs = Matrix<Real, getCellBasisSize<t_element>(), getNumElementUnknowns<t_field>()>();
                    cell_projection_lhs.setZero();
                    cell_projection_rhs.setZero();
                    for (auto i_quadrature = 0; i_quadrature < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i_quadrature++)
                    {
                        auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                        // auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i_quadrature);
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
                        auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
                        auto constexpr t_num_inner_neighbors = ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0, t_i>();
                        auto i_face = 0;
                        for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                        {
                            auto constexpr quadrature_size = ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize();
                            auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                            auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                            auto face_projection_lhs = Matrix<Real, getFaceBasisSize<t_inner_neighbor>(), getFaceBasisSize<t_inner_neighbor>()>();
                            auto face_projection_rhs = Matrix<Real, getFaceBasisSize<t_inner_neighbor>(), getNumElementUnknowns<t_field>()>();
                            face_projection_lhs.setZero();
                            face_projection_rhs.setZero();
                            for (auto i_quadrature = 0; i_quadrature < quadrature_size; i_quadrature++)
                            {
                                auto argh = face_sign == 1 ? i_quadrature : quadrature_size - (i_quadrature + 1);
                                auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                                auto argf = i_quadrature;
                                //
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                                // auto point = face->template getCurrentQuadraturePoint<t_quadrature>(i);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                                // auto normal_vector = face->getNormalVector(face->template getReferenceQuadraturePoint<t_quadrature>(i));
                                auto normal_vector = face->getNormalVector(point);
                                //
                                // //
                                // auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                                // //
                                // // auto face_rotation_matrix = face->getRotationMatrix(face->getReferenceCentroid());
                                // // auto hhh = face_rotation_matrix * (face->getCurrentCentroid() - this->getCurrentCentroid());
                                // // auto face_orientation = hhh(t_domain.getDim() - 1) > 0 ? 1 : -1;
                                // //
                                // auto argh = face_orientation == 1 ? i_quadrature : ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize() - (i_quadrature + 1);
                                // auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                                // //
                                // auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                                // auto point = face->template getCurrentQuadraturePoint<t_quadrature>(i_quadrature);
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
                            for (auto i_quadrature = 0; i_quadrature < ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize(); i_quadrature++)
                            {
                                auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(i_quadrature);
                                // auto point = face->template getCurrentQuadraturePoint<t_quadrature>(i_quadrature);
                                auto point = face->template getReferenceQuadraturePoint<t_quadrature>(i_quadrature);
                                auto row_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                                auto V_E_C_T_O_R = row_face_vector.transpose() * fa_ce_proj;
                                auto factor = (1.0 / face->getLocalFrameDiameters().norm());
                                S_T_A_B_I_L_I_Z_A_T_I_O_N += factor * weight * V_E_C_T_O_R.transpose() * V_E_C_T_O_R;
                            }
                            face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                            i_face ++;
                        }
                        if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                        {
                            self.template operator()<t_i + 1>(self);
                        }
                    };
                    set_faces_blocks(set_faces_blocks);
                }
                return S_T_A_B_I_L_I_Z_A_T_I_O_N;
            }

            template<Field t_field>
            RealMatrix<getGradBasisSize<t_element>(), getNumElementUnknowns<t_field>()>
            getGradientRhs(
                Integer row,
                Integer col
            )
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = RealMatrix<getGradBasisSize<t_element>(), getNumElementUnknowns<t_field>()>();
                rhs.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    // auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i);
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
                    auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0, t_i>();
                    auto i_face = 0;
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                        auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                        for (auto i = 0; i < ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize(); i++)
                        {
                            auto argh = face_sign == 1 ? i : ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize() - (i + 1);
                            auto argf = i;
                            auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            //
                            auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                            // auto point = face->template getCurrentQuadraturePoint<t_quadrature>(i);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                            // auto normal_vector = face->getNormalVector(face->template getReferenceQuadraturePoint<t_quadrature>(i));
                            auto normal_vector = face->getNormalVector(point);
                            //
                            //
                            //
                            // auto face_rotation_matrix = face->getRotationMatrix(face->getReferenceCentroid());
                            // auto hhh = face_rotation_matrix * (face->getCurrentCentroid() - this->getCurrentCentroid());
                            // auto face_orientation = hhh(t_domain.getDim() - 1) > 0 ? 1 : -1;
                            //
                            // auto argh = face_orientation == 1 ? i : ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize() - (i + 1);
                            // auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            //
                            // auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(i);
                            // auto point = face->template getCurrentQuadraturePoint<t_quadrature>(i);
                            //
                            //
                            auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(i_p);
                            auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                            // auto normal_vector = face->getNormalVector(face->template getReferenceQuadraturePoint<t_quadrature>(i));
                            auto normal_vector_component = face_orientation * normal_vector(col);
                            auto face_col_offset = face_offset + row * getFaceBasisSize<t_inner_neighbor>();
                            auto face_block = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset);
                            face_block += normal_vector_component * weight * row_cell_vector * col_face_vector.transpose();
                            //
                            auto cell_col_offset = row * getCellBasisSize<t_element>();
                            auto cell_block = rhs.template block<getGradBasisSize<t_element>(), getCellBasisSize<t_element>()>(0, cell_col_offset);
                            cell_block -= normal_vector_component * weight * row_cell_vector * col_cell_vector.transpose();
                            // PRINTS
                            // std::cout << "i_face ------------------ " << i << " \n";
                            // std::cout << i_face << "\n";
                            // std::cout << "------------------ " << "\n";
                            // std::cout << "face_orientation ------------------ " << i << " \n";
                            // std::cout << face_orientation << "\n";
                            // std::cout << "------------------ " << "\n";
                            // std::cout << "normal_vector ------------------ " << i << " \n";
                            // std::cout << normal_vector << "\n";
                            // std::cout << "------------------ " << "\n";
                            // std::cout << "col_face_vector ------------------ " << i << " \n";
                            // std::cout << col_face_vector << "\n";
                            // std::cout << "------------------ " << "\n";
                        }
                        face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                        i_face ++;
                    }
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
                auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = RealMatrix<getGradBasisSize<t_element>(), getNumElementUnknowns<t_field>()>();
                rhs.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    // auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i);
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
                    auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0, t_i>();
                    auto i_face = 0;
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                        auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                        for (auto i = 0; i < ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize(); i++)
                        {
                            auto argh = face_sign == 1 ? i : ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize() - (i + 1);
                            auto argf = i;
                            auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            //
                            auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                            // auto point = face->template getCurrentQuadraturePoint<t_quadrature>(i);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                            // auto normal_vector = face->getNormalVector(face->template getReferenceQuadraturePoint<t_quadrature>(i));
                            auto normal_vector = face->getNormalVector(point);
                            //
                            // auto argh = face_orientation == 1 ? i : ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize() - (i + 1);
                            // auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            //
                            // auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(i);
                            // auto point = face->template getCurrentQuadraturePoint<t_quadrature>(i);
                            //
                            //
                            auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(i_p);
                            auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                            // auto normal_vector = face->getNormalVector(face->template getReferenceQuadraturePoint<t_quadrature>(i));
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
            getSymmetricGradientRhsDEBUG(
                Integer row,
                Integer col
            )
            const
            {
                auto constexpr t_quadrature = Quadrature::gauss(getGradientConstructionQuadratureOrder());
                auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
                auto face_offset = t_field_size * getCellBasisSize<t_element>();
                auto rhs = RealMatrix<getGradBasisSize<t_element>(), getNumElementUnknowns<t_field>()>();
                rhs.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i++)
                {
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    // auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i);
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
                std::cout << "----- rhs" << row << " " << col << std::endl;
                std::cout << mat2str(rhs) << std::endl;
                std::cout << "---------------------------------------------------------------> element " << std::endl;
                std::cout << mat2str(this->getCurrentCoordinates()) << std::endl;
                std::cout << "------> centroid " << std::endl;
                std::cout << mat2str(this->getCurrentCentroid()) << std::endl;
                auto set_faces_blocks = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<0, t_i>();
                    auto constexpr t_num_inner_neighbors = ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0, t_i>();
                    auto i_face = 0;
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto face_orientation = this->template getInnerNeighborOrientation<0, t_i>(i_face);
                        auto face_sign = this->template getInnerNeighborSign<0, t_i>(i_face);
                        std::cout << "----------------------> face " << i_face << " : " << std::endl;
                        std::cout << mat2str(face->getCurrentCoordinates()) << std::endl;
                        std::cout << "------> face_centroid " << std::endl;
                        std::cout << mat2str(face->getCurrentCentroid()) << std::endl;
                        std::cout << "------> face_orientation " << std::endl;
                        std::cout << face_orientation << std::endl;
                        std::cout << "------> face_rotation " << std::endl;
                        std::cout << mat2str(face->getRotationMatrix(face->getReferenceCentroid())) << std::endl;
                        for (auto i = 0; i < ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize(); i++)
                        {
                            auto argh = face_sign == 1 ? i : ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize() - (i + 1);
                            auto argf = i;
                            // argh = i;
                            auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            //
                            auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(argf);
                            // auto point = face->template getCurrentQuadraturePoint<t_quadrature>(i);
                            auto point = face->template getReferenceQuadraturePoint<t_quadrature>(argf);
                            // auto normal_vector = face->getNormalVector(face->template getReferenceQuadraturePoint<t_quadrature>(i));
                            auto normal_vector = face->getNormalVector(point);
                            //
                            //
                            // std::cout << "------> face_orientation : " << face_orientation << std::endl;
                            //
                            // auto face_rotation_matrix = face->getRotationMatrix(face->getReferenceCentroid());
                            // auto hhh = face_rotation_matrix * (face->getCurrentCentroid() - this->getCurrentCentroid());
                            // auto face_orientation = hhh(t_domain.getDim() - 1) > 0 ? 1 : -1;
                            //
                            // auto argh = face_orientation == 1 ? i : ElementQuadratureRuleTraits<t_inner_neighbor, t_quadrature>::getSize() - (i + 1);
                            // // auto i_p = this->template getInnerNeighborReferenceQuadraturePoint<t_quadrature, 0, t_i>(i_face, argh);
                            // auto i_p2 = face->template getCurrentQuadraturePoint<t_quadrature>(argh);
                            // //
                            // auto weight = face->template getCurrentQuadratureWeight<t_quadrature>(i);
                            // auto point = face->template getCurrentQuadraturePoint<t_quadrature>(i);
                            std::cout << "------> weight " << i << " : " << std::endl;
                            std::cout << weight << std::endl;
                            std::cout << "------> point " << i << " : " << std::endl;
                            std::cout << mat2str(face->template getCurrentQuadraturePoint<t_quadrature>(i)) << std::endl;
                            //
                            //
                            auto row_cell_vector = this->template getBasisEvaluation<getGradBasis()>(i_p);
                            auto col_face_vector = face->template getBasisEvaluation<getFaceBasis()>(point);
                            auto col_cell_vector = this->template getBasisEvaluation<getCellBasis()>(i_p);
                            // auto normal_vector = face->getNormalVector(face->template getReferenceQuadraturePoint<t_quadrature>(i));
                            std::cout << "------> row_cell_vector " << i << " : " << std::endl;
                            std::cout << mat2str(row_cell_vector) << std::endl;
                            std::cout << "------> col_face_vector " << i << " : " << std::endl;
                            std::cout << mat2str(col_face_vector) << std::endl;
                            std::cout << "------> col_cell_vector " << i << " : " << std::endl;
                            std::cout << mat2str(col_cell_vector) << std::endl;
                            std::cout << "------> normal_vector " << i << " : " << std::endl;
                            std::cout << mat2str(normal_vector) << std::endl;
                            auto normal_vector_component_i = face_orientation * normal_vector(col);
                            auto normal_vector_component_j = face_orientation * normal_vector(row);
                            //
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
                        auto face_col_offset_i = face_offset + row * getFaceBasisSize<t_inner_neighbor>();
                        auto face_col_offset_j = face_offset + col * getFaceBasisSize<t_inner_neighbor>();
                        auto face_block_i = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset_i);
                        auto face_block_j = rhs.template block<getGradBasisSize<t_element>(), getFaceBasisSize<t_inner_neighbor>()>(0, face_col_offset_j);
                        std::cout << "------> face block " << row << " : " << std::endl;
                        std::cout << mat2str(face_block_i) << std::endl;
                        std::cout << "------> face block " << col << " : " << std::endl;
                        std::cout << mat2str(face_block_j) << std::endl;
                        face_offset += t_field_size * getFaceBasisSize<t_inner_neighbor>();
                        i_face ++;
                    }
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
            Vector<Real, getNumCellUnknowns<t_field>()>
            getCellUnknowns(
                std::basic_string_view<Character> label
            )
            const
            {
                auto unknown = Vector<Real, getNumCellUnknowns<t_field>()>();
                auto const & cell_dof = this->degrees_of_freedom_.at(std::string(label));
                auto cell_block = unknown.template segment<cell_dof.template getSize<t_field, getCellBasis()>()>(0);
                cell_block = cell_dof.template getCoefficients<t_field, getCellBasis()>();
                return unknown;
            }

            template<Field t_field>
            Vector<Real, getNumFaceUnknowns<t_field>()>
            getFaceUnknowns(
                std::basic_string_view<Character> label
            )
            const
            {
                auto unknown = Vector<Real, getNumFaceUnknowns<t_field>()>();
                auto const & face_dof = this->degrees_of_freedom_.at(std::string(label));
                auto face_block = unknown.template segment<face_dof.template getSize<t_field, getFaceBasis()>()>(0);
                face_block = face_dof.template getCoefficients<t_field, getFaceBasis()>();
                return unknown;
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
                for (auto & integration_point : this->quadrature_.at(std::string(behavior_label)).ips_)
                {
                    auto strain_operator_view = StrainOperatorView(integration_point.ops_.at(std::string(degree_of_freedom_label)).data());
                    auto stress_view = StrainView(integration_point.behavior_data_->s1.thermodynamic_forces.data());
                    // std::cout << "stress_view : " << Matrix<Real, 1, -1>(stress_view) << std::endl;
                    internal_forces += integration_point.weight_ * strain_operator_view.transpose() * stress_view;
                }
                auto stabilization_label = std::string(degree_of_freedom_label) + "Stabilization";
                internal_forces += this->parameters_.at(stabilization_label) * this->operators_.at(stabilization_label) * unknown;
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
                for (auto const & load : this->loads_)
                {
                    auto constexpr sizee = Implementation::template getBasisSize<getCellBasis()>();
                    auto offset = load.second.getFunction().getCol() * sizee;
                    auto cell_block = external_forces.template segment<sizee>(offset);
                    auto i_quadrature = 0;
                    for (auto & ip : this->quadrature_.at(std::string(behavior_label)).ips_)
                    {
                        auto const & current_point = ip.getCurrentCoordinates();
                        auto const & point = ip.getReferenceCoordinates();
                        auto vector = load.second.getFunction().getImposedValue(current_point, 0.0) * this->template getBasisEvaluation<getCellBasis()>(point);
                        cell_block += vector;
                        i_quadrature ++;
                    }
                }                
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
                auto ipcount = 0;
                for (auto & integration_point : this->quadrature_.at(std::string(behavior_label)).ips_)
                {
                    auto strain_operator_view = StrainOperatorView(integration_point.ops_.at(std::string(degree_of_freedom_label)).data());
                    auto jacob = integration_point.template getJacobian<t_finite_element_method>();
                    // std::cout << "jacob :\n" << jacob << std::endl;
                    // std::cout << "strain_operator_view :\n" << strain_operator_view << std::endl;
                    jacobian += integration_point.weight_ * strain_operator_view.transpose() * jacob * strain_operator_view;
                }
                auto stabilization_label = std::string(degree_of_freedom_label) + "Stabilization";
                jacobian += this->parameters_.at(stabilization_label) * this->operators_.at(stabilization_label);
                return jacobian;
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            static constexpr
            Integer
            getFieldSize()
            {
                return FieldTraits<t_finite_element_method.getField()>::template getSize<t_domain>();
            }

            // -------------------------------------------------------------------------------------------------------------------------------------------------
            // -------------------------------------------------------------------------------------------------------------------------------------------------
            // -------------------------------------------------------------------------------------------------------------------------------------------------
            
            // Integer
            // getElementarySystemSize(
            //     std::basic_string_view<Character> behavior_label
            // )
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     auto const & formulation = this->quadrature_.at(std::string(behavior_label));
            //     return formulation.residual_vector_.size();
            // }

            // template<FiniteElementMethodConcept auto t_finite_element_method>
            // void
            // setElementarySystem(
            //     std::basic_string_view<Character> behavior_label,
            //     std::basic_string_view<Character> degree_of_freedom_label,
            //     std::unique_ptr<System> const & system
            // )
            // requires(t_element.isSub(t_domain, 0))
            // {
            //     auto constexpr num_cell_unknowns = getNumCellUnknowns<t_finite_element_method.getField()>();
            //     auto constexpr num_face_unknowns = getNumElementUnknowns<t_finite_element_method.getField()>() - num_cell_unknowns;
            //     auto & formulation = this->quadrature_.at(std::string(behavior_label));
            //     auto internal_forces = getInternalForces<t_finite_element_method>(behavior_label, degree_of_freedom_label);
            //     auto external_forces = getExternalForces<t_finite_element_method>(behavior_label, degree_of_freedom_label);
            //     auto residual = external_forces - internal_forces;
            //     auto jacobian = getJacobianMatrix<t_finite_element_method>(behavior_label, degree_of_freedom_label);
            //     auto k_tt = jacobian.template block<num_cell_unknowns, num_cell_unknowns>(0, 0);
            //     auto k_tf = jacobian.template block<num_cell_unknowns, num_face_unknowns>(0, num_cell_unknowns);
            //     auto k_ft = jacobian.template block<num_face_unknowns, num_cell_unknowns>(num_cell_unknowns, 0);
            //     auto k_ff = jacobian.template block<num_face_unknowns, num_face_unknowns>(num_cell_unknowns, num_cell_unknowns);
            //     auto r_t = residual.template segment<num_cell_unknowns>(0);
            //     auto r_f = residual.template segment<num_face_unknowns>(num_cell_unknowns);
            //     auto k_tt_inv = k_tt.llt().solve(decltype(k_tt)::Identity());
            //     auto k_c = k_ff - k_ft * k_tt_inv * k_tf;
            //     auto r_c = r_f - k_ft * k_tt_inv * r_t;
            //     formulation.jacobian_matrix_ = k_c;
            //     formulation.residual_vector_ = r_c;
            //     formulation.auxiliary_matrix_terms_["KTT"] = k_tt_inv;
            //     formulation.auxiliary_matrix_terms_["KTF"] = k_tf;
            //     formulation.auxiliary_matrix_terms_["KFT"] = k_ft;
            //     formulation.auxiliary_matrix_terms_["KFF"] = k_ff;
            //     formulation.auxiliary_vector_terms_["RC"] = r_c;
            //     system->setNormalization(external_forces.cwiseAbs().maxCoeff());
            // }

            // template<FiniteElementMethodConcept auto t_finite_element_method>
            // void
            // assembleUnknownBlock(
            //     std::basic_string_view<Character> behavior_label,
            //     std::basic_string_view<Character> degree_of_freedom_label,
            //     std::vector<System::MatrixEntry> & lhs_values,
            //     std::unique_ptr<System> const & system
            // )
            // {
            //     auto constexpr t_field = t_finite_element_method.getField();
            //     auto const & formulation = this->quadrature_.at(std::string(behavior_label));
            //     auto global_offset = system->getUnknownOffset(degree_of_freedom_label);
            //     auto offset_i = 0;
            //     auto set_faces_unknowns = [&] <Integer t_i = 0> (
            //         auto & self
            //     )
            //     constexpr mutable
            //     {
            //         for (auto const & face_i : this->template getInnerNeighbors<0, t_i>())
            //         {
            //             auto const & face_i_dof = face_i->degrees_of_freedom_.at(std::string(degree_of_freedom_label));
            //             auto constexpr size_iii = face_i_dof.template getSize<t_field, getFaceBasis()>();
            //             auto iii = global_offset + face_i_dof.getTag();
            //             auto offset_j = 0;
            //             auto set_faces_unknowns2 = [&] <Integer t_j = 0> (
            //                 auto & self2
            //             )
            //             constexpr mutable
            //             {
            //                 for (auto const & face_j : this->template getInnerNeighbors<0, t_j>())
            //                 {
            //                     auto const & face_j_dof = face_j->degrees_of_freedom_.at(std::string(degree_of_freedom_label));
            //                     auto constexpr size_jjj = face_j_dof.template getSize<t_field, getFaceBasis()>();
            //                     auto jjj = global_offset + face_j_dof.getTag();
            //                     system->addLhsValues(iii, jjj, formulation.jacobian_matrix_.template block<size_iii, size_jjj>(offset_i, offset_j));
            //                     offset_j += size_jjj;
            //                 }
            //                 if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
            //                 {
            //                     self2.template operator ()<t_j + 1>(self2);
            //                 }
            //             };
            //             set_faces_unknowns2(set_faces_unknowns2);
            //             offset_i += size_iii;
            //         }
            //         if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
            //         {
            //             self.template operator ()<t_i + 1>(self);
            //         }
            //     };
            //     set_faces_unknowns(set_faces_unknowns);
            //     // <- TEST
            // }

            // -------------------------------------------------------------------------------------------------------------------------------------------------
            // -------------------------------------------------------------------------------------------------------------------------------------------------
            // -------------------------------------------------------------------------------------------------------------------------------------------------

            template<FiniteElementMethodConcept auto t_finite_element_method>
            void
            condensateUnknown(
                std::basic_string_view<Character> behavior_label,
                std::basic_string_view<Character> degree_of_freedom_label
            )
            {
                auto constexpr num_cell_unknowns = getNumCellUnknowns<t_finite_element_method.getField()>();
                auto constexpr num_face_unknowns = getNumElementUnknowns<t_finite_element_method.getField()>() - num_cell_unknowns;
                auto jac = getJacobianMatrix<t_finite_element_method>(behavior_label, degree_of_freedom_label);
                auto k_tt = jac.template block<num_cell_unknowns, num_cell_unknowns>(0, 0);
                auto k_tf = jac.template block<num_cell_unknowns, num_face_unknowns>(0, num_cell_unknowns);
                auto k_ft = jac.template block<num_face_unknowns, num_cell_unknowns>(num_cell_unknowns, 0);
                auto k_ff = jac.template block<num_face_unknowns, num_face_unknowns>(num_cell_unknowns, num_cell_unknowns);
                auto k_tt_inv = k_tt.llt().solve(decltype(k_tt)::Identity());
                this->operators_["KTT"] = k_tt_inv;
                this->operators_["KTF"] = k_tf;
                this->operators_["KFT"] = k_ft;
                this->operators_["KFF"] = k_ff;
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            void
            assembleUnknownVector(
                std::basic_string_view<Character> behavior_label,
                std::basic_string_view<Character> degree_of_freedom_label,
                std::unique_ptr<System> const & system
            )
            {
                auto constexpr num_cell_unknowns = getNumCellUnknowns<t_finite_element_method.getField()>();
                auto constexpr strain_operator_num_cols = getNumElementUnknowns<t_finite_element_method.getField()>();
                auto constexpr num_face_unknowns = getNumElementUnknowns<t_finite_element_method.getField()>() - num_cell_unknowns;
                this->template condensateUnknown<t_finite_element_method>(behavior_label, degree_of_freedom_label);
                auto internal_forces = getInternalForces<t_finite_element_method>(behavior_label, degree_of_freedom_label);
                auto external_forces = getExternalForces<t_finite_element_method>(behavior_label, degree_of_freedom_label);
                system->setNormalization(external_forces.cwiseAbs().maxCoeff());
                // auto residual = internal_forces - external_forces;
                auto residual = external_forces - internal_forces;
                // std::cout << "element_residual : " << Matrix<Real, 1, -1>(residual) << std::endl;
                // -> TEST
                // auto jac = getJacobianMatrix<t_finite_element_method>(behavior_label, degree_of_freedom_label);
                // auto k_tt = jac.template block<num_cell_unknowns, num_cell_unknowns>(0, 0);
                // auto k_tf = jac.template block<num_cell_unknowns, num_face_unknowns>(0, num_cell_unknowns);
                // auto k_ft = jac.template block<num_face_unknowns, num_cell_unknowns>(num_cell_unknowns, 0);
                // auto k_ff = jac.template block<num_face_unknowns, num_face_unknowns>(num_cell_unknowns, num_cell_unknowns);
                // auto k_tt_inv = k_tt.llt().solve(decltype(k_tt)::Identity());
                // <- TEST
                auto r_t = residual.template segment<num_cell_unknowns>(0);
                auto r_f = residual.template segment<num_face_unknowns>(num_cell_unknowns);
                // -> TEST
                // this->operators_["KTT"] = k_tt_inv;
                // this->operators_["KTF"] = k_tf;
                // <- TEST
                this->operators_["RT"] = r_t;
                // auto k_c = k_ff - k_ft * k_tt_inv * k_tf;
                auto const & k_ft = this->operators_.at("KFT");
                auto const & k_tt_inv = this->operators_.at("KTT");
                auto r_c = r_f - k_ft * k_tt_inv * r_t;
                //
                //
                // -> DEBUG
                // std::cout << "********* cell : " << lolita::mat2str(this->getCurrentCentroid()) << std::endl;
                // std::cout << "element_residual : " << lolita::mat2str(residual) << std::endl;
                // <- DEBUG
                auto constexpr t_field = t_finite_element_method.getField();
                auto offset_i = 0;
                auto offset_count = 0;
                auto global_offset = system->getUnknownOffset(degree_of_freedom_label);
                auto set_faces_unknowns = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    for (auto const & face_i : this->template getInnerNeighbors<0, t_i>())
                    {
                        // -> DEBUG
                        // std::cout << "*** face : " << lolita::mat2str(face_i->getCurrentCentroid()) << std::endl;
                        // auto const unknown = face_i->degrees_of_freedom_.at("Displacement").template getCoefficients<lolita::Field::vector(), getFaceBasis()>();
                        // std::cout << "unknown : " << lolita::mat2str(unknown) << std::endl;
                        // if (face_i->isIn("TOP"))
                        // {
                        //     std::cout << "*** TOP" << std::endl;
                        //     auto const force = face_i->degrees_of_freedom_.at("TopForce").template getCoefficients<lolita::Field::scalar(), getFaceBasis()>();
                        //     std::cout << "force : " << lolita::mat2str(force) << std::endl;
                        // }
                        // if (face_i->isIn("BOTTOM"))
                        // {
                        //     std::cout << "*** BOTTOM" << std::endl;
                        //     auto const force = face_i->degrees_of_freedom_.at("BottomForce").template getCoefficients<lolita::Field::scalar(), getFaceBasis()>();
                        //     std::cout << "force : "  << lolita::mat2str(force) << std::endl;
                        // }
                        // if (face_i->isIn("LEFT"))
                        // {
                        //     std::cout << "*** LEFT" << std::endl;
                        //     auto const force = face_i->degrees_of_freedom_.at("LeftForce").template getCoefficients<lolita::Field::scalar(), getFaceBasis()>();
                        //     std::cout << "force : "  << lolita::mat2str(force) << std::endl;
                        // }
                        // <- DEBUG
                        auto const & face_i_dof = face_i->degrees_of_freedom_.at(std::string(degree_of_freedom_label));
                        auto constexpr size_iii = face_i_dof.template getSize<t_field, getFaceBasis()>();
                        auto iii = global_offset + face_i_dof.getTag();
                        system->addRhsValues(iii, r_c.template segment<size_iii>(offset_i));
                        offset_i += size_iii;
                        // -> TEST
                        // for (auto i = 0; i < face_i_dof.template getSize<t_field, getFaceBasis()>(); i++)
                        // {
                        //     system->addRhsValue(i + face_i_dof.getTag() + global_offset, r_c(offset_i));
                        //     offset_i ++;
                        // }
                        // <- TEST
                        // face_count ++;
                    }
                    if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_faces_unknowns(set_faces_unknowns);
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            void
            assembleUnknownBlock(
                std::basic_string_view<Character> behavior_label,
                std::basic_string_view<Character> degree_of_freedom_label,
                std::unique_ptr<System> const & system
            )
            {
                auto constexpr num_cell_unknowns = getNumCellUnknowns<t_finite_element_method.getField()>();
                auto constexpr strain_operator_num_cols = getNumElementUnknowns<t_finite_element_method.getField()>();
                auto constexpr num_face_unknowns = getNumElementUnknowns<t_finite_element_method.getField()>() - num_cell_unknowns;
                // auto internal_forces = getInternalForces<t_finite_element_method>(behavior_label, degree_of_freedom_label);
                // auto external_forces = getExternalForces<t_finite_element_method>(behavior_label, degree_of_freedom_label);
                // system->setNormalization(external_forces.cwiseAbs().maxCoeff());
                // auto residual = internal_forces - external_forces;
                // -> TEST
                // auto jac = getJacobianMatrix<t_finite_element_method>(behavior_label, degree_of_freedom_label);
                // auto k_tt = jac.template block<num_cell_unknowns, num_cell_unknowns>(0, 0);
                // auto k_tf = jac.template block<num_cell_unknowns, num_face_unknowns>(0, num_cell_unknowns);
                // auto k_ft = jac.template block<num_face_unknowns, num_cell_unknowns>(num_cell_unknowns, 0);
                // auto k_ff = jac.template block<num_face_unknowns, num_face_unknowns>(num_cell_unknowns, num_cell_unknowns);
                // auto k_tt_inv = k_tt.llt().solve(decltype(k_tt)::Identity());
                // <- TEST
                auto const & k_tt_inv = this->operators_.at("KTT");
                auto const & k_tf = this->operators_.at("KTF");
                auto const & k_ft = this->operators_.at("KFT");
                auto const & k_ff = this->operators_.at("KFF");
                // auto r_t = residual.template segment<num_cell_unknowns>(0);
                // auto r_f = residual.template segment<num_face_unknowns>(num_cell_unknowns);
                // this->operators_["KTT"] = k_tt_inv;
                // this->operators_["KTF"] = k_tf;
                // this->operators_["RT"] = r_t;
                auto k_c = k_ff - k_ft * k_tt_inv * k_tf;
                // -> DEBUG
                // std::cout << "********* cell : " << lolita::mat2str(this->getCurrentCentroid()) << std::endl;
                // std::cout << lolita::mat2str(k_c) << std::endl;
                // <- DEBUG
                // auto r_c = r_f - k_ft * k_tt_inv * r_t;
                //
                //
                // auto constexpr t_field = t_finite_element_method.getField();
                // auto offset_i = 0;
                // auto global_offset = system->getUnknownOffset(degree_of_freedom_label);
                // auto set_faces_unknowns = [&] <Integer t_i = 0> (
                //     auto & self
                // )
                // constexpr mutable
                // {
                //     for (auto const & face_i : this->template getInnerNeighbors<0, t_i>())
                //     {
                //         auto const & face_i_dof = face_i->degrees_of_freedom_.at(std::string(degree_of_freedom_label));
                //         for (auto i = 0; i < face_i_dof.template getSize<t_field, getFaceBasis()>(); i++)
                //         {
                //             // system->addRhsValue(i + face_i_dof.getTag() + global_offset, - r_c(offset_i));
                //             auto offset_j = 0;
                //             auto set_faces_unknowns2 = [&] <Integer t_j = 0> (
                //                 auto & self2
                //             )
                //             constexpr mutable
                //             {
                //                 for (auto const & face_j : this->template getInnerNeighbors<0, t_j>())
                //                 {
                //                     auto const & face_j_dof = face_j->degrees_of_freedom_.at(std::string(degree_of_freedom_label));
                //                     for (auto j = 0; j < face_j_dof.template getSize<t_field, getFaceBasis()>(); j++)
                //                     {
                //                         system->addLhsValue(i + face_i_dof.getTag() + global_offset, j + face_j_dof.getTag() + global_offset, k_c(offset_i, offset_j));
                //                         offset_j ++;
                //                     }
                //                 }
                //                 if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                //                 {
                //                     self2.template operator ()<t_j + 1>(self2);
                //                 }
                //             };
                //             set_faces_unknowns2(set_faces_unknowns2);
                //             offset_i ++;
                //         }
                //     }
                //     if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                //     {
                //         self.template operator ()<t_i + 1>(self);
                //     }
                // };
                // set_faces_unknowns(set_faces_unknowns);
                // -> TEST
                auto constexpr t_field = t_finite_element_method.getField();
                auto global_offset = system->getUnknownOffset(degree_of_freedom_label);
                auto offset_i = 0;
                auto set_faces_unknowns = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    for (auto const & face_i : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto const & face_i_dof = face_i->degrees_of_freedom_.at(std::string(degree_of_freedom_label));
                        auto constexpr size_iii = face_i_dof.template getSize<t_field, getFaceBasis()>();
                        auto iii = global_offset + face_i_dof.getTag();
                        auto offset_j = 0;
                        auto set_faces_unknowns2 = [&] <Integer t_j = 0> (
                            auto & self2
                        )
                        constexpr mutable
                        {
                            for (auto const & face_j : this->template getInnerNeighbors<0, t_j>())
                            {
                                auto const & face_j_dof = face_j->degrees_of_freedom_.at(std::string(degree_of_freedom_label));
                                auto constexpr size_jjj = face_j_dof.template getSize<t_field, getFaceBasis()>();
                                auto jjj = global_offset + face_j_dof.getTag();
                                system->addLhsValues(iii, jjj, k_c.template block<size_iii, size_jjj>(offset_i, offset_j));
                                offset_j += size_jjj;
                            }
                            if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                            {
                                self2.template operator ()<t_j + 1>(self2);
                            }
                        };
                        set_faces_unknowns2(set_faces_unknowns2);
                        offset_i += size_iii;
                    }
                    if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_faces_unknowns(set_faces_unknowns);
                // <- TEST
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            void
            assembleBindingVector(
                std::basic_string_view<Character> binding_label,
                std::basic_string_view<Character> unknown_label,
                std::basic_string_view<Character> constraint_label,
                std::unique_ptr<System> const & system,
                Real const & time
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto constexpr quadrature = Quadrature::gauss(2 * getFaceBasis().getOrd());
                auto constexpr field = t_finite_element_method.getField();
                auto const & face_unknown = this->degrees_of_freedom_.at(std::string(unknown_label));
                auto const & face_binding = this->degrees_of_freedom_.at(std::string(binding_label));
                auto const & constraint = this->constraints_.at(std::string(constraint_label)).getFunction();
                auto lag_label = std::string(binding_label) + "Lagrange";
                auto const & lagrange_parameter = this->parameters_.at(lag_label);
                // -> DEBUG
                // auto matrix = Matrix<Real, getFaceBasisSize<t_element>(), getFaceBasisSize<t_element>()>();
                // <- DEBUG
                auto binding_external_forces_vector = Vector<Real, getFaceBasisSize<t_element>()>();
                auto binding_internal_forces_vector = Vector<Real, getFaceBasisSize<t_element>()>();
                auto unknown_internal_forces_vector = Vector<Real, getFaceBasisSize<t_element>()>();
                auto const unknown_vector = face_unknown.template getCoefficients<field, getFaceBasis()>(constraint.getRow(), constraint.getCol());
                auto const binding_vector = face_binding.template getCoefficients<Field::scalar(), getFaceBasis()>(0, 0);
                // -> DEBUG
                // matrix.setZero();
                // <- DEBUG
                binding_external_forces_vector.setZero();
                binding_internal_forces_vector.setZero();
                unknown_internal_forces_vector.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, quadrature>::getSize(); i++)
                {
                    auto point = this->template getCurrentQuadraturePoint<quadrature>(i);
                    // auto point = this->template getReferenceQuadraturePoint<quadrature>(i);
                    auto point_ref = this->template getReferenceQuadraturePoint<quadrature>(i);
                    auto weight = this->template getCurrentQuadratureWeight<quadrature>(i);
                    auto basis_vector = this->template getBasisEvaluation<getFaceBasis()>(point_ref);
                    // auto diam = 1.0 / this->getCurrentDiameters().norm();
                    // auto diam = 1.0;
                    // binding_external_forces_vector += weight * constraint.getImposedValue(point, time) * basis_vector;
                    binding_external_forces_vector += weight * lagrange_parameter * constraint.getImposedValue(point, time) * basis_vector;
                    // binding_internal_forces_vector += weight * lagrange_parameter * unknown_vector.dot(basis_vector) * basis_vector;
                    // binding_internal_forces_vector += diam * weight * lagrange_parameter * unknown_vector;
                    binding_internal_forces_vector += weight * lagrange_parameter * basis_vector.transpose() * unknown_vector * basis_vector;
                    // binding_internal_forces_vector += weight * unknown_vector;
                    // unknown_internal_forces_vector += diam * weight * lagrange_parameter * binding_vector;
                    unknown_internal_forces_vector += weight * lagrange_parameter * basis_vector.transpose() * binding_vector * basis_vector;
                    // -> DEBUG
                    // matrix += weight * basis_vector * basis_vector.transpose();
                    // <- DEBUG
                }
                // -> DEBUG
                // std::cout << binding_label << " " << this->getTag() << std::endl;
                // std::cout << "--  unknown :" << std::endl;
                // std::cout << unknown_vector.transpose() << std::endl;
                // std::cout << "--  binding with coef " << lagrange_parameter << " :" << std::endl;
                // std::cout << binding_vector.transpose() << std::endl;
                // std::cout << "-- residual :" << std::endl;
                // std::cout << mat2str(binding_external_forces_vector - binding_internal_forces_vector) << std::endl;
                // std::cout << "-- internal force :" << std::endl;
                // std::cout << unknown_internal_forces_vector.transpose() << std::endl;
                // unknown_internal_forces_vector = binding_vector;
                // auto M_A_T_R_I_X = matrix.llt().solve(decltype(matrix)::Identity());
                // auto V_E_C_T_O_R_0 = M_A_T_R_I_X * binding_external_forces_vector;
                // auto face_displacement_difference = unknown_vector - V_E_C_T_O_R_0;
                // auto V_E_C_T_O_R = lagrange_parameter * (M_A_T_R_I_X * binding_external_forces_vector - unknown_vector);
                // auto V_E_C_T_O_R = lagrange_parameter * (M_A_T_R_I_X * (binding_external_forces_vector - binding_internal_forces_vector));
                // auto lagrange_external_forces = lagrange_parameter * V_E_C_T_O_R_0;
                // <- DEBUG
                // system->setNormalization(lagrange_external_forces.cwiseAbs().maxCoeff());
                system->setNormalization(binding_external_forces_vector.cwiseAbs().maxCoeff());
                auto const binding_offset = system->getBindingOffset(binding_label) + face_binding.getTag();
                auto const unknown_offset = system->getUnknownOffset(unknown_label) + face_unknown.getTag() + getFaceBasisSize<t_element>() * constraint.getRow();
                // -> TEST
                // for (auto iii = 0; iii < getFaceBasisSize<t_element>(); iii++)
                // {
                //     // system->addRhsValue(iii + binding_offset, - lagrange_parameter * face_displacement_difference(iii));
                //     // system->addRhsValue(iii + unknown_offset, - lagrange_parameter * binding_vector(iii));
                //     system->addRhsValue(iii + binding_offset, (binding_external_forces_vector - binding_internal_forces_vector)(iii));
                //     system->addRhsValue(iii + unknown_offset, - unknown_internal_forces_vector(iii));
                //     // for (auto jjj = 0; jjj < getFaceBasisSize<t_element>(); jjj++)
                //     // {
                //     //     system->addLhsValue(iii + binding_offset, jjj + unknown_offset, matrix(iii, jjj));
                //     //     system->addLhsValue(jjj + unknown_offset, iii + binding_offset, matrix(iii, jjj));
                //     // }
                // }
                // <- TEST
                system->addRhsValues(binding_offset, binding_external_forces_vector - binding_internal_forces_vector);
                system->addRhsValues(unknown_offset, - unknown_internal_forces_vector);
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            void
            assembleBindingBlock(
                std::basic_string_view<Character> binding_label,
                std::basic_string_view<Character> unknown_label,
                std::basic_string_view<Character> constraint_label,
                std::unique_ptr<System> const & system
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto constexpr quadrature = Quadrature::gauss(2 * getFaceBasis().getOrd());
                auto constexpr field = t_finite_element_method.getField();
                auto const & face_unknown = this->degrees_of_freedom_.at(std::string(unknown_label));
                auto const & face_binding = this->degrees_of_freedom_.at(std::string(binding_label));
                auto const & constraint = this->constraints_.at(std::string(constraint_label)).getFunction();
                auto matrix = Matrix<Real, getFaceBasisSize<t_element>(), getFaceBasisSize<t_element>()>();
                // auto binding_external_forces_vector = Vector<Real, getFaceBasisSize<t_element>()>();
                // auto binding_internal_forces_vector = Vector<Real, getFaceBasisSize<t_element>()>();
                // auto unknown_internal_forces_vector = Vector<Real, getFaceBasisSize<t_element>()>();
                // auto unknown_vector = face_unknown.template getCoefficients<field, getFaceBasis()>(constraint.getRow(), constraint.getCol());
                // auto binding_vector = face_binding.template getCoefficients<Field::scalar(), getFaceBasis()>(0, 0);
                auto lag_label = std::string(binding_label) + "Lagrange";
                auto const & lagrange_parameter = this->parameters_.at(lag_label);
                matrix.setZero();
                // binding_external_forces_vector.setZero();
                // binding_internal_forces_vector.setZero();
                // unknown_internal_forces_vector.setZero();
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, quadrature>::getSize(); i++)
                {
                    auto point = this->template getCurrentQuadraturePoint<quadrature>(i);
                    // auto point = this->template getReferenceQuadraturePoint<quadrature>(i);
                    auto point_ref = this->template getReferenceQuadraturePoint<quadrature>(i);
                    auto weight = this->template getCurrentQuadratureWeight<quadrature>(i);
                    auto basis_vector = this->template getBasisEvaluation<getFaceBasis()>(point_ref);
                    // binding_external_forces_vector += weight * lagrange_parameter * constraint.getImposedValue(point, 0.0) * basis_vector;
                    // binding_internal_forces_vector += weight * lagrange_parameter * unknown_vector.dot(basis_vector) * basis_vector;
                    // unknown_internal_forces_vector += weight * lagrange_parameter * binding_vector;
                    matrix += weight * lagrange_parameter * basis_vector * basis_vector.transpose();
                }
                // -> DEBUG
                // auto M_A_T_R_I_X = matrix.llt().solve(decltype(matrix)::Identity());
                // auto M_A_T_R_I_X = decltype(matrix)();
                // M_A_T_R_I_X.setIdentity();
                // M_A_T_R_I_X *= lagrange_parameter;
                // <- DEBUG
                // auto M_A_T_R_I_X = matrix.llt().solve(decltype(matrix)::Identity());
                // auto V_E_C_T_O_R = M_A_T_R_I_X * (binding_internal_forces_vector - binding_external_forces_vector);
                // auto binding_residual_vector = binding_internal_forces_vector - binding_external_forces_vector;
                // system->setNormalization(binding_external_forces_vector.cwiseAbs().maxCoeff());
                auto const binding_offset = system->getBindingOffset(binding_label) + face_binding.getTag();
                auto const unknown_offset = system->getUnknownOffset(unknown_label) + face_unknown.getTag() + getFaceBasisSize<t_element>() * constraint.getRow();
                // -> TEST
                // for (auto iii = 0; iii < getFaceBasisSize<t_element>(); iii++)
                // {
                //     // system->addRhsValue(iii + binding_offset, - binding_residual_vector(iii));
                //     // system->addRhsValue(iii + unknown_offset, - unknown_internal_forces_vector(iii));
                //     for (auto jjj = 0; jjj < getFaceBasisSize<t_element>(); jjj++)
                //     {
                //         system->addLhsValue(iii + binding_offset, jjj + unknown_offset, matrix(iii, jjj));
                //         system->addLhsValue(jjj + unknown_offset, iii + binding_offset, matrix(iii, jjj));
                //     }
                // }
                // <- TEST
                system->addLhsValues(binding_offset, unknown_offset, matrix);
                system->addLhsValues(unknown_offset, binding_offset, matrix);
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            void
            updateUnknown(
                std::basic_string_view<Character> unknown_label,
                std::unique_ptr<System> const & system
            )
            requires(t_element.isSub(t_domain, 0))
            {
                auto constexpr field = t_finite_element_method.getField();
                auto constexpr num_cell_unknowns = getNumCellUnknowns<t_finite_element_method.getField()>();
                auto constexpr strain_operator_num_cols = getNumElementUnknowns<t_finite_element_method.getField()>();
                auto constexpr num_face_unknowns = getNumElementUnknowns<t_finite_element_method.getField()>() - num_cell_unknowns;
                //
                auto faces_correction = Vector<Real, num_face_unknowns>();
                faces_correction.setZero();
                auto faces_correction_offset = 0;
                auto global_offset = system->getUnknownOffset(unknown_label);
                auto set_faces_unknowns = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    for (auto const & face : this->template getInnerNeighbors<0, t_i>())
                    {
                        auto const & face_unknown = face->degrees_of_freedom_.at(std::string(unknown_label));
                        auto constexpr face_unknown_size = face_unknown.template getSize<field, getFaceBasis()>();
                        auto face_unknown_offset = face_unknown.getTag() + global_offset;
                        auto face_correction = system->getUnknownCorrection(unknown_label).template segment<face_unknown_size>(face_unknown_offset);
                        faces_correction.template segment<face_unknown_size>(faces_correction_offset) = face_correction;
                        faces_correction_offset += face_unknown_size;
                    }
                    if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumInnerNeighbors<0>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_faces_unknowns(set_faces_unknowns);
                // std::cout << "faces_correction : " << "\n";
                // std::cout << faces_correction << "\n";
                // auto k_tt_inv = algebra::View<Matrix<Real, num_cell_unknowns, num_cell_unknowns> const>(this->operators_["KTT"].data());
                // auto k_tf = algebra::View<Matrix<Real, num_cell_unknowns, num_face_unknowns> const>(this->operators_["KTF"].data());
                // auto r_t = algebra::View<Vector<Real, num_cell_unknowns> const>(this->operators_["RT"].data());
                auto k_tt_inv = this->operators_.at("KTT").template block<num_cell_unknowns, num_cell_unknowns>(0, 0);
                auto k_tf = this->operators_.at("KTF").template block<num_cell_unknowns, num_face_unknowns>(0, 0);
                auto r_t = this->operators_.at("RT").template block<num_cell_unknowns, 1>(0, 0);
                // auto cell_corr = - (k_tt_inv * r_t + k_tt_inv * k_tf * faces_correction);
                // auto cell_corr = k_tt_inv * (r_t - k_tf * faces_correction);
                auto cell_corr = k_tt_inv * (r_t - k_tf * faces_correction);
                this->degrees_of_freedom_.at(std::string(unknown_label)).template getCoefficients<field, getCellBasis()>() += cell_corr;
                //
                // auto & face_unknown = this->degrees_of_freedom_.at(std::string(unknown_label));
                // auto unknown_vector = face_unknown.template getCoefficientsCopy<field, getFaceBasis()>();
                // face_unknown.template addCoefficients<field, getFaceBasis()>(system->getUnknownCorrection(unknown_label));
                // auto unknown_increment = face_unknown.template getCoefficients<field, getFaceBasis()>() - unknown_vector;
                // auto k_tf = this->operators_["KTF"];
                // auto r_t = this->operators_["RT"];
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            void
            updateUnknown(
                std::basic_string_view<Character> unknown_label,
                std::unique_ptr<System> const & system
            )
            requires(t_element.isSub(t_domain, 1))
            {
                auto constexpr field = t_finite_element_method.getField();
                auto & face_unknown = this->degrees_of_freedom_.at(std::string(unknown_label));
                face_unknown.template addCoefficients<field, getFaceBasis()>(system->getUnknownCorrection(unknown_label));
                // auto constexpr face_unknown_size = face_unknown.template getSize<field, getFaceBasis()>();
                // auto face_unknown_offset = face_unknown.getTag() + system->getUnknownOffset(unknown_label);
                // auto face_correction = system->getUnknownCorrection(unknown_label).template segment<face_unknown_size>(face_unknown_offset);
                // faces_correction.template segment<face_unknown_size>(faces_correction_offset) = face_correction;
                // faces_correction_offset += face_unknown_size;
                // auto constexpr field = t_finite_element_method.getField();
                // auto constexpr num_cell_unknowns = getNumCellUnknowns<t_finite_element_method.getField()>();
                // auto constexpr strain_operator_num_cols = getNumElementUnknowns<t_finite_element_method.getField()>();
                // auto constexpr num_face_unknowns = getNumElementUnknowns<t_finite_element_method.getField()>() - num_cell_unknowns;
                // auto & face_unknown = this->degrees_of_freedom_.at(std::string(unknown_label));
                // auto unknown_vector = face_unknown.template getCoefficientsCopy<field, getFaceBasis()>();
                // face_unknown.template addCoefficients<field, getFaceBasis()>(system->getUnknownCorrection(unknown_label));
                // auto unknown_increment = face_unknown.template getCoefficients<field, getFaceBasis()>() - unknown_vector;
                // auto k_tt_inv = algebra::View<Matrix<Real, num_cell_unknowns, num_cell_unknowns> const>(this->operators_["KTT"].data());
                // // auto k_tf = this->operators_["KTF"];
                // // auto r_t = this->operators_["RT"];
                // auto k_tf = algebra::View<Matrix<Real, num_cell_unknowns, num_face_unknowns> const>(this->operators_["KTF"].data());
                // auto r_t = algebra::View<Vector<Real, num_cell_unknowns> const>(this->operators_["RT"].data());
                // auto cell_corr = - (k_tt_inv * r_t + k_tt_inv * k_tf * unknown_increment);
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            Real
            getUnknownValue(
                std::basic_string_view<Character> unknown_label,
                Point const & point,
                Integer row,
                Integer col
            )
            const
            {
                return 0.0;
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            Real
            getUnknownValue(
                std::basic_string_view<Character> unknown_label,
                Point const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto constexpr field = t_finite_element_method.getField();
                auto coefficients = this->degrees_of_freedom_.at(std::string(unknown_label)).template getCoefficients<field, getCellBasis()>(row, col);
                auto basis_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                return coefficients.dot(basis_vector);
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            Real
            getUnknownValue(
                std::basic_string_view<Character> unknown_label,
                Point const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto constexpr field = t_finite_element_method.getField();
                auto coefficients = this->degrees_of_freedom_.at(std::string(unknown_label)).template getCoefficients<field, getFaceBasis()>(row, col);
                auto basis_vector = this->template getBasisEvaluation<getFaceBasis()>(point);
                return coefficients.dot(basis_vector);
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            Real
            getBindingValue(
                std::basic_string_view<Character> binding_label,
                Point const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 0))
            {
                auto constexpr field = t_finite_element_method.getField();
                auto coefficients = this->degrees_of_freedom_.at(std::string(binding_label)).template getCoefficients<field, getCellBasis()>(row, col);
                auto basis_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                return coefficients.dot(basis_vector);
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            Real
            getBindingValue(
                std::basic_string_view<Character> binding_label,
                Point const & point,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto constexpr field = t_finite_element_method.getField();
                auto coefficients = this->degrees_of_freedom_.at(std::string(binding_label)).template getCoefficients<field, getFaceBasis()>(row, col);
                auto basis_vector = this->template getBasisEvaluation<getFaceBasis()>(point);
                return coefficients.dot(basis_vector);
            }

            template<FiniteElementMethodConcept auto t_finite_element_method>
            Real
            getBindingIntegral(
                std::basic_string_view<Character> binding_label,
                Integer row,
                Integer col
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto constexpr quadrature = Quadrature::gauss(2 * getFaceBasis().getOrd());
                auto constexpr field = Field::scalar();
                auto coefficients = this->degrees_of_freedom_.at(std::string(binding_label)).template getCoefficients<field, getFaceBasis()>(row, col);
                auto value = Real(0);
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, quadrature>::getSize(); i++)
                {
                    auto point_ref = this->template getReferenceQuadraturePoint<quadrature>(i);
                    auto weight = this->template getCurrentQuadratureWeight<quadrature>(i);
                    auto basis_vector = this->template getBasisEvaluation<getFaceBasis()>(point_ref);
                    value += weight * coefficients.dot(basis_vector);
                }
                auto lag_label = std::string(binding_label) + "Lagrange";
                auto const & lagrange_parameter = this->parameters_.at(lag_label);
                // auto const & dofff = this->degrees_of_freedom_.at(std::string(binding_label));
                // auto constexpr ones_size = dofff.template getSize<field, getFaceBasis()>();
                // return coefficients.dot(Vector<Real, ones_size>::Ones()) * lagrange_parameter;
                return lagrange_parameter * value;
            }
        
            template<FiniteElementMethodConcept auto t_finite_element_method>
            Integer
            getBandWidth(
                std::basic_string_view<Character> unknown_label
            )
            const
            requires(t_element.isSub(t_domain, 1))
            {
                auto constexpr field = t_finite_element_method.getField();
                auto band_width = 0;
                auto set_faces_unknowns = [&] <Integer t_i = 0> (
                    auto & self
                )
                constexpr mutable
                {
                    for (auto const & cell : this->template getOuterNeighbors<1, t_i>())
                    {
                        auto constexpr t_outer_neighbor = ElementTraits<t_element, t_domain>::template getOuterNeighbor<1, t_i>();
                        auto set_faces_unknowns2 = [&] <Integer t_j = 0> (
                            auto & self2
                        )
                        constexpr mutable
                        {
                            for (auto const & face : cell->template getInnerNeighbors<0, t_j>())
                            {
                                auto face_unknown_size = face->degrees_of_freedom_.at(std::string(unknown_label)).template getSize<field, getFaceBasis()>();
                                band_width += face_unknown_size;
                            }
                            if constexpr (t_j < ElementTraits<t_outer_neighbor, t_domain>::template getNumInnerNeighbors<0>() - 1)
                            {
                                self2.template operator ()<t_j + 1>(self2);
                            }
                        };
                        set_faces_unknowns2(set_faces_unknowns2);
                    }
                    if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumOuterNeighbors<1>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                };
                set_faces_unknowns(set_faces_unknowns);
                return band_width;
            }

            void
            hello()
            {

            }

        };

    };
    
} // namespace lolita

#endif /* C17BAFB5_3EEB_4E27_8A93_F8EEC5AAF622 */
