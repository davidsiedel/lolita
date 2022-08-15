#ifndef E11520F9_5899_4E9B_8D95_F11263F73062
#define E11520F9_5899_4E9B_8D95_F11263F73062

#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_0000.hxx"
#include "lolita/lolita_core_n_1000.hxx"

namespace lolita::core
{

    struct DegreeOfFreedom
    {

        DegreeOfFreedom(
            std::basic_string_view<Character> label
        )
        :
        label_(label)
        {}
        
        std::basic_string_view<Character>
        getLabel()
        const
        {
            return label_;
        }
        
        lolita::matrix::Vector<Real> const &
        getCoefficients()
        const
        {
            return coefficients_;
        }
        
        lolita::matrix::Vector<Real> &
        getCoefficients()
        {
            return coefficients_;
        }
        
        lolita::matrix::Vector<Real> const &
        getVariations()
        const
        {
            return variations_;
        }
        
        lolita::matrix::Vector<Real> &
        getVariations()
        {
            return variations_;
        }

        std::basic_string_view<Character> label_;

        lolita::matrix::Vector<Real> coefficients_;

        lolita::matrix::Vector<Real> variations_;

    };

}



#endif /* E11520F9_5899_4E9B_8D95_F11263F73062 */
