#ifndef C17BAFB5_3EEB_4E27_8A93_F8EEC5AAF622
#define C17BAFB5_3EEB_4E27_8A93_F8EEC5AAF622

#include "lolita_core/lolita.hxx"
#include "lolita_core/lolita_core_n_0000.hxx"
#include "lolita_core/lolita_core_n_1000.hxx"
#include "lolita_core/lolita_core_n_2000.hxx"
#include "lolita_core/lolita_core_n_3000.hxx"

namespace lolita
{

    template<Element t_element, Domain t_domain>
    struct FiniteElementHolder;
        
    template<Basis t_basis>
    struct FiniteElementBasisTraits;

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

            template<Mapping t_mapping, Field t_field>
            RealMatrix<getMappingSize<t_mapping, t_field>(), getNumElementUnknowns<t_field>()>
            getMapping(
                Point const & point
            )
            const
            {
                auto grd = RealMatrix<getMappingSize<t_mapping, t_field>(), getNumElementUnknowns<t_field>()>();
                grd.setZero();
                // getGradientRhs(0, 0);
                auto lhs = getGradientLhs();
                for (auto const & mapping_value : MappingTraits<t_mapping>::template getValues<t_domain, t_field>())
                {
                    auto rhs = getGradientRhs<t_field>(mapping_value.row(), mapping_value.col());
                    auto line = grd.template block<1, getNumElementUnknowns<t_field>()>(mapping_value.rank(), 0);
                    // auto res = lhs * rhs;
                    // auto res2 = this->template getBasisEvaluation<getGradBasis()>(point);
                    // std::cout << "ici rhs : " << rhs.rows() << ", " << rhs.cols() << std::endl;
                    // std::cout << "ici lhs : " << lhs.rows() << ", " << lhs.cols() << std::endl;
                    // std::cout << "ici res2 : " << res2.rows() << ", " << res2.cols() << std::endl;
                    // // auto left_cell_vector = this->template getBasisEvaluation<getCellBasis()>(point);
                    line = mapping_value.value() * this->template getBasisEvaluation<getGradBasis()>(point).transpose() * lhs * rhs;
                }
                return grd;
            }

            // template<Mapping t_mapping, Basis t_gradient_basis>
            // lolita::matrix::Matrix<Real, getMappingSize<t_mapping>(), getNumCols<t_element, t_domain>()>
            // getGradient(
            //     Integer i
            // )
            // const
            // requires(t_mapping.isGradient())
            // {
            //     auto rhs = lolita::matrix::Matrix<Real, getMappingSize<t_mapping>(), getNumCols<t_element, t_domain>()>();
            //     rhs.setZero();
            //     return rhs;
            // }

        };

    };
    
} // namespace lolita

#endif /* C17BAFB5_3EEB_4E27_8A93_F8EEC5AAF622 */
