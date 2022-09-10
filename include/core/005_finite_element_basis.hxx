#ifndef ADE371C2_0D15_4707_B243_25103F8F862B
#define ADE371C2_0D15_4707_B243_25103F8F862B

#include "core/lolita.hxx"
#include "core/000_physics_traits.hxx"
#include "core/001_geometry.hxx"
#include "core/linear_system.hxx"
#include "core/003_quadrature.hxx"
#include "core/004_finite_element.hxx"

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

        public:
        
            Vector<Real, getSize()>
            getBasisEvaluation(
                Point const & point
            )
            const
            {
                auto basis_vector_values = Vector<Real, getSize()>();
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
            
            Vector<Real, getSize()>
            getBasisDerivative(
                Point const & point,
                Integer derivative_direction
            )
            const
            {
                auto basis_vector_values = Vector<Real, getSize()>();
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
    
} // namespace lolita


#endif /* ADE371C2_0D15_4707_B243_25103F8F862B */
