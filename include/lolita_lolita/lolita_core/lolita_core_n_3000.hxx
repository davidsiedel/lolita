#ifndef CC79BDC5_49DB_4A81_8A93_18ABD6551AF1
#define CC79BDC5_49DB_4A81_8A93_18ABD6551AF1

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_2000.hxx"

namespace lolita
{
    
    template<Element t_element, Quadrature t_quadrature>
    struct ElementQuadratureRuleTraits;
    
    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.isNode() || !t_element.isNode())
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {
        
        Integer static constexpr dim_ = 1;
        
        Integer static constexpr size_ = 1;

        static constexpr
        Integer
        size()
        {
            return 1;
        }

        static constexpr
        Integer
        getSize()
        {
            return 1;
        }

        std::array<std::array<Real, 3>, dim_> static constexpr reference_points_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };
        
        std::array<Real, dim_> static constexpr reference_weights_ = {
                +1.0000000000000000
        };

    };

} // namespace lolita

#endif /* CC79BDC5_49DB_4A81_8A93_18ABD6551AF1 */
