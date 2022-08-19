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
    
} // namespace lolita

#endif /* C17BAFB5_3EEB_4E27_8A93_F8EEC5AAF622 */
