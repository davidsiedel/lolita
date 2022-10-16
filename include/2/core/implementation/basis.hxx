#ifndef CE2653A7_6805_4868_A8D8_8790EF5ABAB1
#define CE2653A7_6805_4868_A8D8_8790EF5ABAB1

#include "2/core/traits/_include.hxx"
#include "2/core/traits/basis_traits.hxx"

namespace lolita::core
{

    template<ShapeConcept auto t_element, MeshConcept auto t_domain, BasisConcept auto t_basis>
    requires(LagrangeBasisConcept<decltype(t_basis)> && t_basis.getOrder() == 1)
    struct BasisImplementation<t_element, t_domain, t_basis> : FiniteElement<t_element, t_domain>
    {

    };

    template<ShapeConcept auto t_element, MeshConcept auto t_domain, BasisConcept auto t_basis>
    requires(MonomialBasisConcept<decltype(t_basis)>)
    struct BasisImplementation<t_element, t_domain, t_basis> : FiniteElement<t_element, t_domain>
    {

    private:

        static constexpr
        Integer
        getSize()
        {
            return BasisTraits<t_basis>::template getSize<t_element>();
        }
    
        static constexpr
        std::array<std::array<Integer, 3>, getSize()>
        getExponents()
        {
            auto exponents = std::array<std::array<Integer, 3>, getSize()>();
            auto row = Integer(0);
            if constexpr (t_element.getDim() == 0)
            {
                exponents[row][0] = 0;
                exponents[row][1] = 0;
                exponents[row][2] = 0;
            }
            else if constexpr (t_element.getDim() == 1)
            {
                for (auto i = 0; i < t_basis.getOrd() + 1; ++i)
                {
                    exponents[row][0] = i;
                    exponents[row][1] = 0;
                    exponents[row][2] = 0;
                    row += 1;
                }
            }
            else if constexpr (t_element.getDim() == 2)
            {
                for (auto i = 0; i < t_basis.getOrd() + 1; ++i)
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
            else if constexpr (t_element.getDim() == 3)
            {
                for (auto i = 0; i < t_basis.getOrd() + 1; ++i)
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
    
        DenseVector<Real, getSize()>
        getBasisEvaluation(
            PointConcept auto const & point
        )
        const
        {
            auto basis_vector_values = DenseVector<Real, getSize()>();
            auto const centroid = this->getReferenceCentroid();
            // auto const centroid = this->getCurrentCentroid();
            // auto const diameters = this->getCurrentDiameters();
            auto const diameters = this->getLocalFrameDiameters();
            for (auto i = 0; i < getSize(); ++i)
            {
                auto value = Real(1);
                for (auto j = 0; j < t_element.getDim(); ++j)
                {
                    // auto dist = this->getRiemannianDistance(centroid, point, j);
                    auto dist = this->getLocalFrameDistance(centroid, point, j);
                    value *= std::pow(2.0 * dist / diameters(j), exponents_[i][j]);
                }
                basis_vector_values(i) = value;
            }
            return basis_vector_values;
        }
        
        DenseVector<Real, getSize()>
        getBasisDerivative(
            PointConcept auto const & point,
            Integer derivative_direction
        )
        const
        {
            auto basis_vector_values = DenseVector<Real, getSize()>();
            auto const centroid = this->getReferenceCentroid();
            // auto const centroid = this->getCurrentCentroid();
            // auto const diameters = this->getCurrentDiameters();
            auto const diameters = this->getLocalFrameDiameters();
            for (auto i = 0; i < getSize(); ++i)
            {
                auto value = Real(1);
                for (auto j = 0; j < t_element.getDim(); ++j)
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
    
} // namespace lolita::core


#endif /* CE2653A7_6805_4868_A8D8_8790EF5ABAB1 */
