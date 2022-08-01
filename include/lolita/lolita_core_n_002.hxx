#ifndef E11520F9_5899_4E9B_8D95_F11263F73062
#define E11520F9_5899_4E9B_8D95_F11263F73062

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_000.hxx"
#include "lolita/lolita_core_n_001.hxx"

namespace lolita2::geometry
{

    struct DegreeOfFreedom
    {

        DegreeOfFreedom(
            std::basic_string_view<lolita::character> label
        )
        :
        label_(label)
        {}
        
        std::basic_string_view<lolita::character>
        getLabel()
        const
        {
            return label_;
        }
        
        lolita::matrix::Vector<lolita::real> const &
        getCoefficients()
        const
        {
            return coefficients_;
        }
        
        lolita::matrix::Vector<lolita::real> &
        getCoefficients()
        {
            return coefficients_;
        }
        
        lolita::matrix::Vector<lolita::real> const &
        getVariations()
        const
        {
            return variations_;
        }
        
        lolita::matrix::Vector<lolita::real> &
        getVariations()
        {
            return variations_;
        }

        std::basic_string_view<lolita::character> label_;

        lolita::matrix::Vector<lolita::real> coefficients_;

        lolita::matrix::Vector<lolita::real> variations_;

    };

}



#endif /* E11520F9_5899_4E9B_8D95_F11263F73062 */
