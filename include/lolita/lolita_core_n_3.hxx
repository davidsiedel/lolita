#ifndef E0160644_2FBA_4722_B59E_AA106834F181
#define E0160644_2FBA_4722_B59E_AA106834F181

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_00.hxx"
#include "lolita/lolita_core_n_0.hxx"
#include "lolita/lolita_core_n_1.hxx"
#include "lolita/lolita_core_n_2.hxx"

namespace lolita2::geometry
{

    struct Unknown
    {

        Unknown()
        :
        coefficients_()
        {}

        Unknown(
            lolita::integer size
        )
        :
        coefficients_(size)
        {}

        lolita::matrix::Vector<lolita::real> coefficients_;

    };

    template<Element t_element, Domain t_domain, Field t_field, Basis t_basis>
    struct FiniteElementUnknown
    {

        lolita::integer static constexpr dim_unknown_ = basis::FiniteElementBasisTraits<t_element, t_basis>::dim_;

        lolita::integer static constexpr num_unknown_ = dim_unknown_ * t_field.dim_;

        FiniteElementUnknown(
            std::shared_ptr<FiniteElementGeometry<t_element, t_domain>> const & element,
            std::shared_ptr<Unknown> const & unknown
        )
        :
        element_(element),
        unknown_(unknown)
        {
            unknown->coefficients_.resize(unknown->coefficients_.size() + num_unknown_);
        }

        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, num_unknown_> const>
        getCoefficients()
        const
        {
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, num_unknown_> const>(unknown_->coefficients_.data() + 1);
        }

        std::shared_ptr<FiniteElementGeometry<t_element, t_domain>> element_;

        std::shared_ptr<Unknown> unknown_;

    };

    template<Domain t_domain, Field t_field, Basis t_cell_basis, Basis t_face_basis>
    struct HybridDiscontinuousGalerkinTraits
    {

        template<Element t_element>
        static constexpr
        lolita::integer
        getDimCellBasis()
        {
            return basis::FiniteElementBasisTraits<t_element, t_cell_basis>::dim_;
        }

        template<Element t_element>
        static constexpr
        lolita::integer
        getDimFaceBasis()
        {
            return basis::FiniteElementBasisTraits<t_element, t_face_basis>::dim_;
        }

        // lolita::integer static constexpr dim_cell_unknown_ = basis::FiniteElementBasisTraits<t_element, t_cell_basis>::dim_;

        // lolita::integer static constexpr dim_face_unknown_ = basis::FiniteElementBasisTraits<t_element, t_cell_basis>::dim_;

        // lolita::integer static constexpr dim_unknown_ = basis::FiniteElementBasisTraits<t_element, t_cell_basis>::dim_;

        // lolita::integer static constexpr num_unknown_ = dim_unknown_ * t_field.dim_;

    };

    template<Element t_element, Domain t_domain>
    struct FEBASE : FiniteElementConnectivity<FEBASE, t_element, t_domain>
    {

        FEBASE(
            std::shared_ptr<FiniteElementGeometry<t_element, t_domain>> const & element,
            auto &&... args
        )
        :
        FiniteElementConnectivity<FEBASE, t_element, t_domain>(* element)
        {}

    };

    template<Element t_element, Domain t_domain, Field t_field, Basis t_cell_basis, Basis t_face_basis>
    struct HybridDiscontinuousGalerkin : FEBASE<t_element, t_domain> {};

    template<Element t_element, Domain t_domain, Field t_field, Basis t_cell_basis, Basis t_face_basis>
    requires(t_element.isSub(t_domain, 0))
    struct HybridDiscontinuousGalerkin<t_element, t_domain, t_field, t_cell_basis, t_face_basis> : FEBASE<t_element, t_domain>
    {

        HybridDiscontinuousGalerkin(
            std::shared_ptr<FiniteElementGeometry<t_element, t_domain>> const & element,
            lolita::integer const & a
        )
        :
        FEBASE<t_element, t_domain>(element)
        {}

    };

    template<Element t_element, Domain t_domain, Field t_field, Basis t_cell_basis, Basis t_face_basis>
    requires(t_element.isSub(t_domain, 1))
    struct HybridDiscontinuousGalerkin<t_element, t_domain, t_field, t_cell_basis, t_face_basis> : FEBASE<t_element, t_domain>
    {

        HybridDiscontinuousGalerkin(
            std::shared_ptr<FiniteElementGeometry<t_element, t_domain>> const & element,
            std::shared_ptr<Unknown> const & unknown
        )
        :
        FEBASE<t_element, t_domain>(element)
        {}
        
    };

}


#endif /* E0160644_2FBA_4722_B59E_AA106834F181 */
