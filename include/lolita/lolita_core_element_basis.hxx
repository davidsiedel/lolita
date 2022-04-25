//
// Created by dsiedel on 16/04/22.
//

#ifndef LOLITA_LOLITA_CORE_ELEMENT_BASIS_HXX
#define LOLITA_LOLITA_CORE_ELEMENT_BASIS_HXX

#include "lolita/lolita_matrix.hxx"
#include "lolita/lolita_core.hxx"
#include "lolita/lolita_core_element_geometry.hxx"

namespace lolita::core::element
{

    template<Element E, BasisDescription B>
    struct EB;

    template<Element E, BasisDescription B>
    requires(B.basis == Basis::Monomial)
    struct EB<E, B>
    {

        std::integral auto const static constexpr dim_basis = numerics::binomial(E.dim + B.ord, E.dim);

        template<auto D, auto M>
        struct Implementation : public FiniteElement<E, D, M>
        {

        private:

            using Self = Implementation;

            static constexpr
            auto
            setExponents()
            {
                auto exponents_values = Array<Indx, dim_basis, E.dimPoint()>();
                if constexpr (E.dim == 0) {
                    exponents_values.get(0, 0) = 0;
                }
                else if constexpr (E.dim == 1) {
                    auto row = Indx(0);
                    for (Indx i = 0; i < B.ord + 1; ++i) {
                        exponents_values.get(row, 0) = i;
                        row += 1;
                    }
                }
                else if constexpr (E.dim == 2) {
                    auto row = Indx(0);
                    for (Indx i = 0; i < B.ord + 1; ++i) {
                        for (Indx j = 0; j < i + 1; ++j) {
                            exponents_values.get(row, 0) = i - j;
                            exponents_values.get(row, 1) = j;
                            row += 1;
                        }
                    }
                }
                else if constexpr (E.dim == 3) {
                    auto row = Indx(0);
                    for (Indx i = 0; i < B.ord + 1; ++i) {
                        for (Indx j = 0; j < i + 1; ++j) {
                            for (Indx k = 0; k < i + 1; ++k) {
                                if (j + k < i + 1) {
                                    exponents_values.get(row, 0) = i - (j + k);
                                    exponents_values.get(row, 1) = k;
                                    exponents_values.get(row, 2) = j;
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
                auto basis_vector_values = Vector<Real, dim_basis>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (Indx i = 0; i < dim_basis; ++i) {
                    auto value = Real(1);
                    for (Indx j = 0; j < E.dim; ++j) {
                        auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
                        value *= std::pow(2.0 * dist / diameters(j), exponents.get(i, j));
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
                auto basis_vector_values = Vector<Real, dim_basis>();
                auto const centroid = this->getReferenceCentroid();
                auto const diameters = this->getCurrentDiameters();
                for (Indx i = 0; i < dim_basis; ++i) {
                    auto value = Real(1);
                    for (Indx j = 0; j < E.dim; ++j) {
                        if (j != derivative_direction_arg) {
                            auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
                            value *= std::pow(2.0 * (dist) / diameters(j), exponents.get(i, j));
                        } else {
                            if (exponents.get(i, j) > 0) {
                                auto c = 2.0 * exponents.get(i, j) / diameters(j);
                                auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
//                                value *= c * std::pow(2.0 * (dist) / diameters(j), exponents.get(i, j) - 1);
                                value *= c * numerics::pow(2.0 * (dist) / diameters(j), exponents.get(i, j) - 1);
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

    template<Element E, BasisDescription B>
    static constexpr
    auto
    dimBasis()
    {
        return EB<E, B>::dim_basis;
    }

    namespace detail
    {

        template<Element E, Basis B, Indx K>
        struct FiniteElementBasis;

        template<Element E, Indx K>
        struct FiniteElementBasis<E, Basis::Monomial, K>
        {

            auto const static constexpr dim_basis = numerics::binomial(E.dim + K, E.dim);

            template<auto D, auto M>
            struct Implementation : public FiniteElement<E, D, M>
            {

            private:

                using Self = Implementation;

                static constexpr
                auto
                setExponents()
                {
                    auto exponents_values = Array<Indx, dim_basis, E.dimPoint()>();
                    if constexpr (E.dim == 0) {
                        exponents_values(0, 0) = 0;
                    }
                    else if constexpr (E.dim == 1) {
                        auto row = Indx(0);
                        for (Indx i = 0; i < K + 1; ++i) {
                            exponents_values(row, 0) = i;
                            row += 1;
                        }
                    }
                    else if constexpr (E.dim == 2) {
                        auto row = Indx(0);
                        for (Indx i = 0; i < K + 1; ++i) {
                            for (Indx j = 0; j < i + 1; ++j) {
                                exponents_values(row, 0) = i - j;
                                exponents_values(row, 1) = j;
                                row += 1;
                            }
                        }
                    }
                    else if constexpr (E.dim == 3) {
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
                    auto basis_vector_values = Vector<Real, dim_basis>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (Indx i = 0; i < dim_basis; ++i) {
                        auto value = Real(1);
                        for (Indx j = 0; j < E.dim; ++j) {
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
                    auto basis_vector_values = Vector<Real, dim_basis>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (Indx i = 0; i < dim_basis; ++i) {
                        auto value = Real(1);
                        for (Indx j = 0; j < E.dim; ++j) {
                            if (j != derivative_direction_arg) {
                                auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
                                value *= std::pow(2.0 * (dist) / diameters(j), exponents(i, j));
                            } else {
                                if (exponents(i, j) > 0) {
                                    auto c = 2.0 * exponents(i, j) / diameters(j);
                                    auto dist = this->getDistanceInCurrentConfiguration(centroid, point_arg, j);
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

    template<Element E, BasisDescription B>
    using FiniteElementBasis = detail::FiniteElementBasis<E, B.basis, B.ord>;

    template<Element E, auto D, auto M>
    struct FiniteElementBasis2
    {

        FiniteElementBasis2(
                BasisDescription const &
                basis_description_arg
        )
        :
        basis_description(basis_description_arg)
        {}

        struct Implementation : public FiniteElement<E, D, M>
        {


            auto
            evaluate(
                    auto const &
                    point_arg
            )
            const
            {
                auto v = Vector<Real>();
                if (basis_description == BasisDescription(Basis::Monomial, 1)) {
                    v = detail::FiniteElementBasis<E, Basis::Monomial, 1>::evaluate(point_arg);
                }
                else if (basis_description == BasisDescription(Basis::Monomial, 2)) {
                    v = detail::FiniteElementBasis<E, Basis::Monomial, 2>::evaluate(point_arg);
                }
                else {
                    assert(false);
                }
                return v;
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
                auto v = Vector<Real>();
                if (basis_description == BasisDescription(Basis::Monomial, 1)) {
                    v = detail::FiniteElementBasis<E, Basis::Monomial, 1>::evaluate(point_arg, derivative_direction_arg);
                }
                else if (basis_description == BasisDescription(Basis::Monomial, 2)) {
                    v = detail::FiniteElementBasis<E, Basis::Monomial, 2>::evaluate(point_arg, derivative_direction_arg);
                }
                else {
                    assert(false);
                }
                return v;
            }

        };

        BasisDescription basis_description;

    };

}

#endif //LOLITA_LOLITA_CORE_ELEMENT_BASIS_HXX
