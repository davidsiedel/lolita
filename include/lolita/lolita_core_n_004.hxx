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
    using RealMatrix = lolita::matrix::Matrix<Real, t_args...>;

    template<auto... t_args>
    using RealVector = lolita::matrix::Vector<Real, t_args...>;

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
        
            lolita::matrix::Vector<Real, getSize()>
            getBasisEvaluation(
                Point const & point
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<Real, getSize()>();
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
            
            lolita::matrix::Vector<Real, getSize()>
            getBasisDerivative(
                Point const & point,
                Integer derivative_direction
            )
            const
            {
                auto basis_vector_values = lolita::matrix::Vector<Real, getSize()>();
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
    
    template<Element t_element, Domain t_domain>
    struct FiniteElementDegreeOfFreedom
    {

        template<Field t_field, Basis t_basis>
        static constexpr
        Integer
        getSize()
        {
            auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
            return t_field_size * t_basis_size;
        }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<Real, getSize<t_field, t_basis>()> const>
        getCoefficients()
        const
        {
            auto const & data = degree_of_freedom_->coefficients_.data() + index_;
            return lolita::matrix::Span<lolita::matrix::Vector<Real, getSize<t_field, t_basis>()> const>(data);
        }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<Real, getSize<t_field, t_basis>()>>
        getCoefficients()
        {
            auto const & data = degree_of_freedom_->coefficients_.data() + index_;
            return lolita::matrix::Span<lolita::matrix::Vector<Real, getSize<t_field, t_basis>()>>(data);
        }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>
        getCoefficients(
            Integer row,
            Integer col
        )
        {
            auto constexpr t_field_shape = FieldTraits<t_field>::template shape<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
            auto const & data = degree_of_freedom_->coefficients_.data() + index_ + (t_field_shape.cols() * row  + col) * t_basis_size;
            return lolita::matrix::Span<lolita::matrix::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>(data);
        }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()> const>
        getCoefficients(
            Integer row,
            Integer col
        )
        const
        {
            auto constexpr t_field_shape = FieldTraits<t_field>::template shape<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
            auto const & data = degree_of_freedom_->coefficients_.data() + index_ + (t_field_shape.cols() * row  + col) * t_basis_size;
            return lolita::matrix::Span<lolita::matrix::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>(data);
        }
        
        std::shared_ptr<DegreeOfFreedom> const &
        getDegreeOfFreedom()
        const
        {
            return degree_of_freedom_;
        }

        Integer index_;

        std::shared_ptr<DegreeOfFreedom> degree_of_freedom_;

    };

    struct Load
    {

        Load(
            Loading const & loading,
            Integer row,
            Integer col
        )
        :
        row_(row),
        col_(col),
        loading_(loading)
        {}

        Load(
            Loading && loading,
            Integer row,
            Integer col
        )
        :
        row_(row),
        col_(col),
        loading_(std::forward<Loading>(loading))
        {}
        
        Real
        getImposedValue(
            Point const & point,
            Real const & time
        )
        const
        {
            return loading_(point, time);
        }

        Integer
        getRow()
        const
        {
            return row_;
        }

        Integer
        getCol()
        const
        {
            return col_;
        }

        Loading loading_;

        Integer row_;

        Integer col_;

    };

    struct QuadratureOperator
    {

        inline
        std::basic_string_view<Character>
        getLabel()
        const
        {
            return label_;
        }

        inline
        lolita::matrix::Matrix<Real>
        getOperator(
            Integer i
        )
        const
        {
            return operators_[i];
        }

        std::basic_string_view<Character> label_;

        std::vector<lolita::matrix::Matrix<Real>> operators_;

    };
    
    template<Domain t_domain>
    struct IntegrationPoint
    {

        IntegrationPoint(
            Point const & coordinates,
            Real const & weight,
            std::shared_ptr<mgis::behaviour::Behaviour> const & behaviour
        )
        :
        coordinates_(coordinates),
        weight_(weight),
        behavior_data_(std::make_unique<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(* behaviour)))
        {}

        std::unique_ptr<mgis::behaviour::BehaviourData> const &
        getMaterialPoint()
        const
        {
            return behavior_data_;
        }

        std::unique_ptr<mgis::behaviour::BehaviourData> &
        getMaterialPoint()
        {
            return behavior_data_;
        }
        
        template<BehaviorConcept auto t_behavior>
        lolita::matrix::Span<lolita::matrix::Vector<Real, BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>()> const>
        getGeneralizedStrain()
        const
        {
            auto constexpr size = BehaviorTraits<t_behavior>::template getGeneralizedStrainSize<t_domain>();
            return lolita::matrix::Span<lolita::matrix::Vector<Real, size> const>(behavior_data_->s1.gradients.data());
        }
        
        template<auto t_finite_element_method>
        lolita::matrix::Span<RealVector<FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()> const>
        getGeneralizedStrain()
        const
        {
            auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
            auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
            return lolita::matrix::Span<lolita::matrix::Vector<Real, size> const>(behavior_data_->s1.gradients.data() + offset);
        }
        
        template<auto t_finite_element_method>
        lolita::matrix::Span<RealVector<FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>()>>
        getGeneralizedStrain()
        {
            auto constexpr size = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainSize<t_domain>();
            auto constexpr offset = FiniteElementMethodTraits<t_finite_element_method>::template getGeneralizedStrainOffset<t_domain>();
            return lolita::matrix::Span<lolita::matrix::Vector<Real, size>>(behavior_data_->s1.gradients.data() + offset);
        }
        
        Point coordinates_;
        
        Real weight_;

        std::vector<QuadratureOperator> quadrature_operators_;
        
        std::unique_ptr<mgis::behaviour::BehaviourData> behavior_data_;

    };

    template<Domain t_domain>
    struct ElementIntegrationPoints
    {

        ElementIntegrationPoints(
            std::shared_ptr<mgis::behaviour::Behaviour> const & behaviour
        )
        :
        behavior_(behaviour)
        {}

        Integer
        getSize()
        const
        {
            return integration_points_.size();
        }

        std::basic_string_view<Character>
        getLabel()
        const
        {
            return behavior_->behaviour;
        }

        void
        integrate()
        {
            for (auto const & integration_point : integration_points_)
            {
                auto behaviour_view = mgis::behaviour::make_view(integration_point->behavior_data_);
                auto result = mgis::behaviour::integrate(behaviour_view, * behavior_);
            }
        }
        
        std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

        std::vector<IntegrationPoint<t_domain>> integration_points_;

    };

    template<Element t_element, Domain t_domain>
    struct FiniteElement
    {

        std::shared_ptr<ElementIntegrationPoints<t_domain>> element_integration_points_;

        std::vector<std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>>> degrees_of_freedom_;

        std::vector<lolita::matrix::Matrix<Real>> element_operators_;

        std::vector<std::shared_ptr<Load>> loads_;

        std::shared_ptr<std::basic_string<Character>> label_;

        FiniteElement(
            std::shared_ptr<std::basic_string<Character>> const & label
        )
        :
        label_(label)
        {}

        std::basic_string_view<Character>
        getLabel()
        const
        {
            return std::basic_string_view<Character>(* label_);
        }

        Boolean
        hasLoad(
            Integer row,
            Integer col
        )
        const
        {
            auto has_load = [&] (
                std::shared_ptr<Load> const & load
            )
            {
                return load->getRow() == row && load->getCol() == col;
            };
            return std::find_if(loads_.begin(), loads_.end(), has_load) != loads_.end();
        }

        Integer
        getLoadIndex(
            Integer row,
            Integer col
        )
        const
        {
            auto has_load = [&] (
                std::shared_ptr<Load> const & load
            )
            {
                return load->getRow() == row && load->getCol() == col;
            };
            if (std::find_if(loads_.begin(), loads_.end(), has_load) != loads_.end())
            {
                return std::distance(loads_.begin(), std::find_if(loads_.begin(), loads_.end(), has_load));
            }
            else
            {
                throw std::runtime_error("NO");
            }
        }
        
        void
        addLoad(
            std::shared_ptr<Load> const & load
        )
        {
            if (!hasLoad(load->getRow(), load->getCol()))
            {
                loads_.push_back(load);
            }
            else
            {
                getLoad(load->getRow(), load->getCol()) = load;
            }
        }

        Boolean
        hasDegreeOfFreedom(
            std::basic_string_view<Character> label
        )
        const
        {
            auto has_label = [&] (
                std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> const & degree_of_freedom
            )
            {
                return degree_of_freedom->getDegreeOfFreedom()->getLabel() == label;
            };
            return std::find_if(degrees_of_freedom_.begin(), degrees_of_freedom_.end(), has_label) != degrees_of_freedom_.end();
        }

        Integer
        getDegreeOfFreedomIndex(
            std::basic_string_view<Character> label
        )
        const
        {
            auto has_label = [&] (
                std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> const & degree_of_freedom
            )
            {
                return degree_of_freedom->getDegreeOfFreedom()->getLabel() == label;
            };
            if (std::find_if(degrees_of_freedom_.begin(), degrees_of_freedom_.end(), has_label) != degrees_of_freedom_.end())
            {
                return std::distance(degrees_of_freedom_.begin(), std::find_if(degrees_of_freedom_.begin(), degrees_of_freedom_.end(), has_label));
            }
            else
            {
                throw std::runtime_error("NO");
            }
        }

        template<Field t_field, Basis t_basis>
        void
        addDegreeOfFreedom(
            std::shared_ptr<DegreeOfFreedom> & degree_of_freedom
        )
        {
            auto constexpr t_size = FiniteElementDegreeOfFreedom<t_element, t_domain>::template getSize<t_field, t_basis>();
            if (!hasDegreeOfFreedom(degree_of_freedom->getLabel()))
            {
                auto element_dof = std::make_unique<FiniteElementDegreeOfFreedom<t_element, t_domain>>();
                element_dof->index_ = degree_of_freedom->coefficients_.size();
                element_dof->degree_of_freedom_ = degree_of_freedom;
                degree_of_freedom->coefficients_.resize(degree_of_freedom->coefficients_.size() + t_size);
                degrees_of_freedom_.push_back(std::move(element_dof));
            }
        }

        std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> const &
        getDegreeOfFreedom(
            std::basic_string_view<Character> label
        )
        const
        {
            return degrees_of_freedom_[getDegreeOfFreedomIndex(label)];
        }

        std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> &
        getDegreeOfFreedom(
            std::basic_string_view<Character> label
        )
        {
            return degrees_of_freedom_[getDegreeOfFreedomIndex(label)];
        }

        std::shared_ptr<Load> const &
        getLoad(
            Integer row,
            Integer col
        )
        const
        {
            return loads_[getLoadIndex(row, col)];
        }

        std::shared_ptr<Load> &
        getLoad(
            Integer row,
            Integer col
        )
        {
            return loads_[getLoadIndex(row, col)];
        }

    };

    template<Element t_element, Domain t_domain>
    struct FiniteElementHolder
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element, t_domain>;
        
        template<Element t__element, Domain t__domain>
        using t_ElementPointer = std::shared_ptr<FiniteElementHolder<t__element, t__domain>>;

    public:
    
        using t_InnerNeighbors = typename t_ElementTraits::template InnerConnectivity<t_ElementPointer>;
        
        using t_OuterNeighbors = typename t_ElementTraits::template OuterConnectivity<t_ElementPointer>;

        template<Basis t_basis>
        using t_Basis = typename FiniteElementBasisTraits<t_basis>::template Implementation<t_element, t_domain>;

        template<auto t_discretization>
        using t_Disc = typename HybridDiscontinuousGalerkinTraits<t_discretization>::template Implementation<t_element, t_domain>;
        
        lolita::boolean
        operator==(
            FiniteElementHolder const & other
        )
        const = default;
        
        lolita::boolean
        operator!=(
            FiniteElementHolder const & other
        )
        const = default;

        lolita::boolean
        isIn(
            std::basic_string_view<lolita::character> domain
        )
        const
        {
            auto has_domain = [&] (
                std::shared_ptr<MeshDomain> const & mesh_domain
            )
            {
                return mesh_domain->tag_ == domain;
            };
            return std::find_if(domains_.begin(), domains_.end(), has_domain) != domains_.end();
        }
        
        std::basic_string<lolita::character>
        getHash()
        const
        requires(t_element.isNode())
        {
            return std::to_string(this->tag_);
        }
        
        std::basic_string<lolita::character>
        getHash()
        const
        {
            std::basic_stringstream<lolita::character> hash;
            auto const & nodes = getInnerNeighbors<t_element.dim_ - 1, 0>();
            for (auto const & node : nodes)
            {
                hash << node->getHash();
            }
            return hash.str();
        }
        
        template<Integer t_i, Integer t_j>
        static constexpr
        Integer
        getInnerNeighborNodeConnection(
            Integer i,
            Integer j
        )
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element, t_domain>::node_connectivity_))[i][j];
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_OuterNeighbors>> &
        getOuterNeighbors()
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_OuterNeighbors>> const &
        getOuterNeighbors()
        const
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_InnerNeighbors>> &
        getInnerNeighbors()
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_InnerNeighbors>> const &
        getInnerNeighbors()
        const
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborIndex(
            Integer i
        )
        const
        requires(!t_element.isNode())
        {
            auto constexpr t_inner_neighbor = t_ElementTraits::template getInnerNeighbor<t_i, t_j>();
            auto constexpr t_coordinates = ElementTraits<t_inner_neighbor, t_domain>::template getOuterNeighborCoordinates<t_element>();
            auto const & items = getInnerNeighbors<t_i, t_j>()[i]->template getOuterNeighbors<t_coordinates.dim_, t_coordinates.tag_>();
            auto is_equal = [&] (t_ElementPointer<t_element, t_domain> const & ptr_element)
            {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborIndex(
            std::shared_ptr<FiniteElementHolder<ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>(), t_domain>> const & ptr_neighbor
        )
        const
        requires(!t_element.isNode())
        {
            auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>();
            auto const & inner_neighbors = getInnerNeighbors<t_i, t_j>();
            auto is_equal = [&] (std::shared_ptr<FiniteElementHolder<t_inner_neighbor, t_domain>> const & neighbor)
            {
                return * neighbor == * ptr_neighbor;
            };
            auto neighbor_index = std::distance(inner_neighbors.begin(), std::find_if(inner_neighbors.begin(), inner_neighbors.end(), is_equal));
            return getInnerNeighborIndex<t_i, t_j>(neighbor_index);
        }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborOrientation(
            Integer i
        )
        const
        requires(!t_element.isNode())
        {
            return getInnerNeighborIndex<t_i, t_j>(i) == 0 ? 1 : -1;
        }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborOrientation(
            std::shared_ptr<FiniteElementHolder<ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>(), t_domain>> const & ptr_neighbor
        )
        const
        requires(!t_element.isNode())
        {
            return getInnerNeighborIndex<t_i, t_j>(ptr_neighbor) == 0 ? 1 : -1;
        }
    
        Point &
        getCurrentCoordinates()
        requires(t_element.isNode())
        {
            return * this->coordinates_;
        }
        
        Point const &
        getCurrentCoordinates()
        const
        requires(t_element.isNode())
        {
            return * this->coordinates_;
        }
        
        lolita::matrix::Matrix<Real, 3, t_element.num_nodes_>
        getCurrentCoordinates()
        const
        requires(!t_element.isNode())
        {
            auto current_nodes_coordinates = lolita::matrix::Matrix<Real, 3, t_element.num_nodes_>();
            auto count = Integer(0);
            for (auto const & node : this->template getInnerNeighbors<t_element.dim_ - 1, 0>())
            {
                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
                count ++;
            }
            return current_nodes_coordinates;
        }
        
        static
        lolita::matrix::Span<lolita::matrix::Matrix<Real, 3, t_element.num_nodes_, lolita::matrix::col_major> const>
        getReferenceCoordinates()
        {
            using t_ReferenceCoordinates = lolita::matrix::Span<lolita::matrix::Matrix<Real, 3, t_element.num_nodes_, lolita::matrix::col_major> const>;
            return t_ReferenceCoordinates(ElementTraits<t_element, t_domain>::reference_nodes_.begin()->begin());
        }
        
        static
        Real
        getShapeMappingEvaluation(
            lolita::matrix::Vector<Real, t_element.num_nodes_> const & nodal_field_values,
            Point const & reference_point
        )
        {
            return t_ElementTraits::getShapeMappingEvaluation(nodal_field_values, reference_point);
        }
        
        static
        Real
        getShapeMappingDerivative(
            lolita::matrix::Vector<Real, t_element.num_nodes_> const & nodal_field_values,
            Point const & reference_point,
            Integer derivative_direction
        )
        {
            return t_ElementTraits::getShapeMappingDerivative(nodal_field_values, reference_point, derivative_direction);
        }
        
        Real
        getShapeMappingDifferential(
            Point const & point
        )
        const
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::matrix::Matrix<Real, 3, 3>();
            auto du = Real(0);
            ru.setZero();
            for (auto i = 0; i < t_domain.dim_; ++i)
            {
                for (auto j = 0; j < t_element.dim_; ++j)
                {
                    ru(i, j) = FiniteElementHolder::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            if constexpr (t_element.dim_ == 3)
            {
                du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
            }
            else if constexpr (t_element.dim_ == 2)
            {
                du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
            }
            else
            {
                du = lolita::numerics::abs(ru.col(0).norm());
            }
            if constexpr (t_domain.frame_ == Domain::Frame::AxiSymmetric)
            {
                Real r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
                if (r0 < 1.e-10)
                {
                    r0 = 1.e-10;
                }
                du *= 2.0 * lolita::numerics::pi * r0;
            }
            return du;
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point,
            Integer direction = -1
        )
        const
        {
            if constexpr (t_ElementTraits::hasDim(0))
            {
                auto const & current_coordinates = this->getCurrentCoordinates();
                auto distance = Real();
                auto mp0 = Point();
                auto mp1 = Point();
                for (auto i = 0; i < t_element.dim_; ++i)
                {
                    mp0(i) = FiniteElementHolder::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                    mp1(i) = FiniteElementHolder::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
                }
                direction == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction);
                return distance;
            }
            else
            {
                using SegmentQuadrature = ElementQuadratureRuleTraits<Element::segment(1), Quadrature::gauss(4)>;
                auto distance = Real(0);
                auto dt = Real();
                auto const current_nodes_coordinates = this->getCurrentCoordinates();
                for (auto q = 0; q < SegmentQuadrature::dim_; ++q)
                {
                    auto pq = SegmentQuadrature::reference_points_[q][0];
                    auto wq = SegmentQuadrature::reference_weights_[q];
                    auto ru = lolita::matrix::Matrix<Real, 3, 3>().setZero();
                    auto difference = second_point - first_point;
                    auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                    for (auto i = 0; i < t_domain.dim_; ++i)
                    {
                        for (auto j = 0; j < t_element.dim_; ++j)
                        {
                            if (direction == -1 || i == static_cast<Integer>(direction))
                            {
                                auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                                auto dx = FiniteElementHolder::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                                ru(i, j) = dx * du;
                            }
                        }
                    }
                    if constexpr (t_element.isFacet())
                    {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        auto Fff = ru.col(0).template dot(ru.col(1));
                        auto Gff = ru.col(1).template dot(ru.col(1));
                        dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                    }
                    else if constexpr (t_element.isCurve())
                    {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        dt = std::sqrt(Eff);
                    }
                    else
                    {
                        dt = 0;
                    }
                    distance += wq * dt;
                }
                return distance;
            }
        }
        
        Point
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = FiniteElementHolder::getReferenceCoordinates();
            auto current_diameters = Point().setZero();
            for (auto i = 0; i < t_element.num_nodes_; ++i)
            {
                for (auto j = i + 1; j < t_element.num_nodes_; ++j)
                {
                    auto const & pt0 = reference_coordinates.col(i);
                    auto const & pt1 = reference_coordinates.col(j);
                    for (auto k = 0; k < 3; ++k)
                    {
                        auto new_value = lolita::numerics::abs(getRiemannianDistance(pt0, pt1, k));
                        auto & current_value = current_diameters(k);
                        if (new_value > current_value)
                        {
                            current_value = new_value;
                        }
                    }
                }
            }
            return current_diameters;
        }
        
        Point
        getCurrentCentroid()
        const
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto barycenter = Point();
            barycenter.setZero();
            for (auto i = 0; i < t_element.num_nodes_; ++i) {
                barycenter += current_nodes_coordinates.col(i);
            }
            barycenter /= Real(t_element.num_nodes_);
            return barycenter;
        }
        
        static
        Point
        getReferenceCentroid()
        {
            auto reference_nodes_coordinates = FiniteElementHolder::getReferenceCoordinates();
            auto barycenter = Point();
            barycenter.setZero();
            for (auto i = 0; i < t_element.num_nodes_; ++i) {
                barycenter += reference_nodes_coordinates.col(i);
            }
            barycenter /= Real(t_element.num_nodes_);
            return barycenter;
        }
        
        static
        Point
        getReferenceDiameters()
        {
            auto dts = Point();
            auto nds = FiniteElementHolder::getReferenceCoordinates();
            dts.setZero();
            for (auto i = 0; i < t_element.num_nodes_; ++i)
            {
                for (auto j = i + 1; j < t_element.num_nodes_; ++j)
                {
                    for (auto k = 0; k < 3; ++k)
                    {
                        auto new_value = lolita::numerics::abs(nds(k, i) - nds(k, j));
                        auto & current_value = dts(k);
                        if (new_value > current_value)
                        {
                            current_value = new_value;
                        }
                    }
                }
            }
            return dts;
        }
        
        Point
        getNormalVector(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1))
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::matrix::Matrix<Real, 3, t_element.dim_>();
            ru.setZero();
            for (auto i = 0; i < 3; ++i)
            {
                for (auto j = 0; j < t_element.dim_; ++j)
                {
                    ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
                }
            }
            if constexpr (t_element.isNode()) {
                return Point{0, 0, 0};
            }
            else if constexpr (t_element.isCurve())
            {
                return Point{ru(1)/ru.norm(), -ru(0)/ru.norm(), 0};
            }
            else
            {
                return (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
            }
        }
        
        Point
        getTangentVector(
            Point const & point,
            Integer direction
        )
        const
        requires(t_element.isCurve())
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = Point();
            for (auto i = 0; i < 3; ++i)
            {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, direction);
            }
            return tangent_vector;
        }
        
        template<Quadrature t_quadrature>
        static
        Real
        getReferenceQuadratureWeight(
            Integer index
        )
        {
            return ElementQuadratureRuleTraits<t_element, t_quadrature>::reference_weights_[index];
        }
        
        template<Quadrature t_quadrature>
        static
        lolita::matrix::Span<Point const>
        getReferenceQuadraturePoint(
            Integer index
        )
        {
            return lolita::matrix::Span<Point const>(ElementQuadratureRuleTraits<t_element, t_quadrature>::reference_points_[index].begin());
        }
        
        template<Quadrature t_quadrature>
        Real
        getCurrentQuadratureWeight(
            Integer index
        )
        const
        {
            auto w = getReferenceQuadratureWeight<t_quadrature>(index);
            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<t_quadrature>(index));
        }
        
        template<Quadrature t_quadrature>
        Point
        getCurrentQuadraturePoint(
            Integer index
        )
        const
        {
            auto p = Point();
            auto const nds = this->getCurrentCoordinates();
            for (auto j = 0; j < 3; ++j)
            {
                p(j) = FiniteElementHolder::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<t_quadrature>(index));
            }
            return p;
        }
        
        template<Quadrature t_quadrature, Integer _i, Integer _j>
        static
        Real
        getInnerNeighborReferenceQuadratureWeight(
            Integer index
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = t_ElementTraits::template getInnerNeighbor<_i, _j>();
            using ComponentGeometry = FiniteElementHolder<_component, t_domain>;
            return ComponentGeometry::template getReferenceQuadratureWeight<t_quadrature>(index);
        }
        
        template<Quadrature t_quadrature, Integer _i, Integer _j>
        static
        Point
        getInnerNeighborReferenceQuadraturePoint(
            Integer component_index,
            Integer index
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = t_ElementTraits ::template getInnerNeighbor<_i, _j>();
            auto p = Point();
            using ComponentGeometry = FiniteElementHolder<_component, t_domain>;
            auto const & elt_reference_nodes = ElementTraits<t_element, t_domain>::reference_nodes_;
            for (auto i = 0; i < 3; ++i)
            {
                auto cpt_coordinates = lolita::matrix::Vector<Real, _component.num_nodes_>();
                for (auto j = 0; j < _component.num_nodes_; ++j)
                {
                    auto const node_tag = getInnerNeighborNodeConnection<_i, _j>(component_index, j);//.get(component_index).get(j);
                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
                }
                auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<t_quadrature>(index);
                p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }
        
        template<Quadrature t_quadrature, Integer t_i, Integer t_j>
        Real
        getInnerNeighborCurrentQuadratureWeight(
            Integer component_index,
            Integer index
        )
        const
        requires(!t_element.isNode())
        {
            auto const & cmp =  this->template getInnerNeighbors<t_i, t_j>()[component_index];//.template get<I>().template get<J>().get(component_index).get();
            return cmp->template getCurrentQuadratureWeight<t_quadrature>(index);
        }
        
        template<Quadrature t_quadrature, Integer t_i, Integer t_j>
        Point
        getInnerNeighborCurrentQuadraturePoint(
            Integer component_index,
            Integer index
        )
        const
        requires(!t_element.isNode())
        {
            auto p = Point();
            auto const cpt_ref_pnt = getInnerNeighborReferenceQuadraturePoint<t_quadrature, t_i, t_j>(component_index, index);
            auto const nds = this->getCurrentCoordinates();
            for (auto j = 0; j < 3; ++j)
            {
                p(j) = FiniteElementHolder::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }

        // -----------------------------------------------------------------------------------------------------------------------------------------------------
        // NEW
        // -----------------------------------------------------------------------------------------------------------------------------------------------------

        Boolean
        hasFiniteElement(
            std::basic_string_view<Character> label
        )
        const
        {
            auto has_label = [&] (
                std::shared_ptr<FiniteElement<t_element, t_domain>> const & finite_element
            )
            {
                return finite_element->getLabel() == label;
            };
            return std::find_if(finite_elements_.begin(), finite_elements_.end(), has_label) != finite_elements_.end();
        }

        Boolean
        hasBehavior(
            std::basic_string_view<Character> label
        )
        const
        {
            auto has_label = [&] (
                std::shared_ptr<ElementIntegrationPoints<t_domain>> const & element_integration_points
            )
            {
                return element_integration_points->getLabel() == label;
            };
            return std::find_if(integration_points_.begin(), integration_points_.end(), has_label) != integration_points_.end();
        }

        Integer
        getFiniteElementIndex(
            std::basic_string_view<Character> label
        )
        const
        {
            auto has_label = [&] (
                std::shared_ptr<FiniteElement<t_element, t_domain>> const & finite_element
            )
            {
                return finite_element->getLabel() == label;
            };
            if (std::find_if(finite_elements_.begin(), finite_elements_.end(), has_label) != finite_elements_.end())
            {
                return std::distance(finite_elements_.begin(), std::find_if(finite_elements_.begin(), finite_elements_.end(), has_label));
            }
            else
            {
                throw std::runtime_error("NO");
            }
        }

        Integer
        getBehaviorIndex(
            std::basic_string_view<Character> label
        )
        const
        {
            auto has_label = [&] (
                std::shared_ptr<ElementIntegrationPoints<t_domain>> const & element_integration_points
            )
            {
                return element_integration_points->getLabel() == label;
            };
            if (std::find_if(integration_points_.begin(), integration_points_.end(), has_label) != integration_points_.end())
            {
                return std::distance(integration_points_.begin(), std::find_if(integration_points_.begin(), integration_points_.end(), has_label));
            }
            else
            {
                throw std::runtime_error("NO");
            }
        }
        
        std::shared_ptr<FiniteElement<t_element, t_domain>> const &
        getFiniteElement(
            std::basic_string_view<Character> label
        )
        const
        {
            return finite_elements_[getFiniteElementIndex(label)];
        }
        
        std::shared_ptr<FiniteElement<t_element, t_domain>> &
        getFiniteElement(
            std::basic_string_view<Character> label
        )
        {
            return finite_elements_[getFiniteElementIndex(label)];
        }
        
        std::shared_ptr<ElementIntegrationPoints<t_domain>> const &
        getIntegrationPoints(
            std::basic_string_view<Character> label
        )
        const
        {
            return integration_points_[getBehaviorIndex(label)];
        }
        
        std::shared_ptr<ElementIntegrationPoints<t_domain>> &
        getIntegrationPoints(
            std::basic_string_view<Character> label
        )
        {
            return integration_points_[getBehaviorIndex(label)];
        }

        void
        addElement(
            std::shared_ptr<std::basic_string<Character>> const & label
        )
        {
            if (!hasFiniteElement(* label))
            {
                finite_elements_.push_back(std::make_shared<FiniteElement<t_element, t_domain>>(label));
            }
        }

        template<Quadrature t_quadrature>
        void
        addBehavior(
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            if (!hasBehavior(behavior->behaviour))
            {
                auto element_integration_points = std::make_shared<ElementIntegrationPoints<t_domain>>(behavior);
                for (auto i = 0; i < ElementQuadratureRuleTraits<t_element, t_quadrature>::getSize(); i++)
                {
                    auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i);
                    auto weight = this->template getCurrentQuadratureWeight<t_quadrature>(i);
                    auto integration_point = IntegrationPoint<t_domain>(point, weight, behavior);
                    element_integration_points->integration_points_.push_back(std::move(integration_point));
                }
                std::cout << "dt : " << element_integration_points->integration_points_[0].behavior_data_->dt << std::endl;
                integration_points_.push_back(element_integration_points);
            }
        }

        void
        addBehavior(
            std::basic_string_view<Character> finite_element_label,
            std::basic_string_view<Character> behavior_label
        )
        {
            getFiniteElement(finite_element_label)->element_integration_points_ = getIntegrationPoints(behavior_label);
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
            Integer derivative_direction
        )
        const
        {
            return static_cast<t_Basis<t_basis> const *>(this)->getBasisDerivative(point, derivative_direction);
        }

        template<Field t_field, Mapping t_mapping, auto t_discretization>
        RealMatrix<MappingTraits<t_mapping>::template getSize<t_domain, t_field>(), t_Disc<t_discretization>::template getNumElementUnknowns<t_field>()>
        getMapping(
            Point const & point
        )
        const
        {
            return static_cast<t_Disc<t_discretization> const *>(this)->template getMapping<t_mapping, t_field>(point);
        }

        template<Field t_field, Quadrature t_quadrature, auto t_discretization>
        void
        makeQuadrature(
            std::basic_string_view<Character> label
        )
        {
            for (auto i = 0; i < getFiniteElement(label)->element_integration_points_->getSize(); i++)
            {
                auto point = this->template getCurrentQuadraturePoint<t_quadrature>(i);
                auto mapp = this->template getMapping<t_field, Mapping::gradient(), t_discretization>(point);
                // element_integration_points_->integration_points_[i].quadrature_operators_[0] = mapp;
            }
            // for (auto const & integration_point : element_integration_points_->integration_points_)
            // {
            //     integration_point.quadrature_operators_[0] = 
            // }
            // std::cout << "mapp : " << std::endl;
            // std::cout << mapp << std::endl;
        }
        
        t_OuterNeighbors outer_neighbors_;
        
        t_InnerNeighbors inner_neighbors_;
        
        lolita::natural tag_;

        std::vector<std::shared_ptr<MeshDomain>> domains_;
        
        std::shared_ptr<Point> coordinates_;

        std::vector<std::shared_ptr<FiniteElement<t_element, t_domain>>> finite_elements_;

        std::vector<std::shared_ptr<ElementIntegrationPoints<t_domain>>> integration_points_;

    };
    
}


#endif /* ECF0EB38_E45D_4E38_ACEB_E40188311824 */
