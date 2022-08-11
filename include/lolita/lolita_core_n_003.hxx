#ifndef CC79BDC5_49DB_4A81_8A93_18ABD6551AF1
#define CC79BDC5_49DB_4A81_8A93_18ABD6551AF1

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_000.hxx"
#include "lolita/lolita_core_n_001.hxx"
#include "lolita/lolita_core_n_002.hxx"

namespace lolita2::geometry
{
    
    template<Element t_element, Quadrature t_quadrature>
    struct ElementQuadratureRuleTraits;
    
    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.isNode() || !t_element.isNode())
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {
        
        lolita::integer static constexpr dim_ = 1;
        
        lolita::integer static constexpr size_ = 1;

        static constexpr
        lolita::integer
        size()
        {
            return 1;
        }

        static constexpr
        lolita::integer
        getSize()
        {
            return 1;
        }

        std::array<std::array<lolita::real, 3>, dim_> static constexpr reference_points_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };
        
        std::array<lolita::real, dim_> static constexpr reference_weights_ = {
                +1.0000000000000000
        };

    };

}

#endif /* CC79BDC5_49DB_4A81_8A93_18ABD6551AF1 */
