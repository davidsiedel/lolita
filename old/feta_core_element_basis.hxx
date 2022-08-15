//
// Created by dsiedel on 02/04/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_BASIS_HXX
#define FETA_FETA_CORE_ELEMENT_BASIS_HXX

#include "new/_feta.hxx"
#include "new/_feta_array.hxx"
#include "new/_feta_matrix.hxx"
#include "new/feta_core_element_basis_description.hxx"
#include "new/feta_core_element_element_description.hxx"
//#include "new/element_geometry.hxx"

namespace feta::core::element::finite_element
{

    namespace detail
    {

        template<ElementDescription E, Basis B, Indx K>
        struct FiniteElementBasisPolicy;

        template<ElementDescription E, Indx K>
        struct FiniteElementBasisPolicy<E, Basis::Monomial, K>
        {

            auto const static constexpr dim_basis = getBinomial(element::ordShape<E>() + K, element::ordShape<E>());

            template<auto M>
            struct Implementation : public Element<E, M>
            {

            private:

                using Self = Implementation;

                static constexpr
                auto
                setExponents()
                {
                    auto exponents_values = Array<Indx, dim_basis, element::dimShape<E>()>();
                    if constexpr (element::ordShape<E>() == 0) {
                        exponents_values(0, 0) = 0;
                    }
                    else if constexpr (element::ordShape<E>() == 1) {
                        auto row = Indx(0);
                        for (Indx i = 0; i < K + 1; ++i) {
                            exponents_values(row, 0) = i;
                            row += 1;
                        }
                    }
                    else if constexpr (element::ordShape<E>() == 2) {
                        auto row = Indx(0);
                        for (Indx i = 0; i < K + 1; ++i) {
                            for (Indx j = 0; j < i + 1; ++j) {
                                exponents_values(row, 0) = i - j;
                                exponents_values(row, 1) = j;
                                row += 1;
                            }
                        }
                    }
                    else if constexpr (element::ordShape<E>() == 3) {
                        auto row = Indx(0);
                        for (Indx i = 0; i < K + 1; ++i) {
                            for (Indx j = 0; j < i + 1; ++j) {
                                for (Indx k = 0; k < i + 1; ++k) {
                                    if (j + k < i + 1) {
                                        exponents_values(row, 0) = i - (j + k);
                                        exponents_values(row, 1) = k;
                                        exponents_values(row, 2) = j;
                                        row += 1;
                                    }
                                }
                            }
                        }
                    }
                    return exponents_values;
                }

                auto const static constexpr exponents = setExponents();

            public:

                auto
                evaluate(
                        auto const &
                        point_arg
                )
                const
                {
                    auto basis_vector_values = Matrix<Real, dim_basis>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (Indx i = 0; i < dim_basis; ++i) {
                        auto value = Real(1);
                        for (Indx j = 0; j < element::dimShape<E>(); ++j) {
                            auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
                            value *= std::pow(2.0 * dist / diameters(j), exponents(i, j));
                        }
                        basis_vector_values(i) = value;
                    }
                    return basis_vector_values;
                }

                auto
                evaluate(
                        auto const &
                        point_arg,
                        auto
                        derivative_direction_arg
                )
                const
                {
                    auto basis_vector_values = Matrix<Real, dim_basis>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (Indx i = 0; i < dim_basis; ++i) {
                        auto value = Real(1);
                        for (Indx j = 0; j < element::dimShape<E>(); ++j) {
                            if (j != derivative_direction_arg) {
                                auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
                                value *= std::pow(2.0 * (dist) / diameters(j), exponents(i, j));
                            } else {
                                if (exponents(i, j) > 0) {
                                    auto c = 2.0 * exponents(i, j) / diameters(j);
                                    auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg,
                                                                                        j);
                                    value *= c * std::pow(2.0 * (dist) / diameters(j), exponents(i, j) - 1);
                                } else {
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

    }

    template<ElementDescription E, BasisDescription B>
    using FiniteElementBasis = detail::FiniteElementBasisPolicy<E, B.basis, B.order>;

    template<ElementDescription E, BasisDescription B>
    constexpr inline
    auto
    dimBasis()
    {
        return FiniteElementBasis<E, B>::dim_basis;
    }

}

#endif //FETA_FETA_CORE_ELEMENT_BASIS_HXX
