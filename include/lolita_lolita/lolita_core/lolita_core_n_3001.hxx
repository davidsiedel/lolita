#ifndef F402E53D_B7DA_4E80_B9B0_6CC57A882459
#define F402E53D_B7DA_4E80_B9B0_6CC57A882459

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_2000.hxx"

namespace lolita
{

    template<Quadrature t_quadrature>
    struct QuadratureTraits;

    template<Quadrature t_quadrature>
    requires(t_quadrature.isGauss())
    struct QuadratureTraits<t_quadrature>
    {

    };
   
} // namespace lolita


#endif /* F402E53D_B7DA_4E80_B9B0_6CC57A882459 */
