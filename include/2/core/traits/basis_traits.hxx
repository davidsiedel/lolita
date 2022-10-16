#ifndef A09243B0_4CD6_437A_8E0A_2484A1307A3C
#define A09243B0_4CD6_437A_8E0A_2484A1307A3C

#include "2/core/traits/_include.hxx"

namespace lolita::core
{
    
    template<ShapeConcept auto t_element, MeshConcept auto t_domain>
    struct FiniteElement;

    template<ShapeConcept auto t_element, MeshConcept auto t_domain, BasisConcept auto t_basis>
    struct BasisImplementation;
    
    template<BasisConcept auto>
    struct BasisTraits;

    template<BasisConcept auto t_basis>
    requires(LagrangeBasisConcept<decltype(t_basis)> && t_basis.getOrder() == 1)
    struct BasisTraits<t_basis>
    {

        template<ShapeConcept auto shape_>
        static constexpr
        Integer
        getSize()
        {
            return 0;
        }

        template<PointShapeConcept auto shape_>
        static constexpr
        Integer
        getSize()
        {
            return 1;
        }

    };

    template<BasisConcept auto t_basis>
    requires(MonomialBasisConcept<decltype(t_basis)>)
    struct BasisTraits<t_basis>
    {
        
        template<ShapeConcept auto shape_>
        static constexpr
        Integer
        getSize()
        {
            return lolita::numerics::binomial(shape_.getDim() + t_basis.getOrd(), shape_.getDim());
        }

    };
    
} // namespace lolita::core



#endif /* A09243B0_4CD6_437A_8E0A_2484A1307A3C */
