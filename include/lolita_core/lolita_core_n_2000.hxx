#ifndef E11520F9_5899_4E9B_8D95_F11263F73062
#define E11520F9_5899_4E9B_8D95_F11263F73062

#include "lolita_core/lolita.hxx"
#include "lolita_core/lolita_core_n_0000.hxx"
#include "lolita_core/lolita_core_n_1000.hxx"

namespace lolita
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
        
        lolita::algebra::Vector<Real> const &
        getCoefficients()
        const
        {
            return coefficients_;
        }
        
        lolita::algebra::Vector<Real> &
        getCoefficients()
        {
            return coefficients_;
        }
        
        lolita::algebra::Vector<Real> const &
        getVariations()
        const
        {
            return variations_;
        }
        
        lolita::algebra::Vector<Real> &
        getVariations()
        {
            return variations_;
        }

        std::basic_string_view<Character> label_;

        lolita::algebra::Vector<Real> coefficients_;

        lolita::algebra::Vector<Real> variations_;

    };

} // namespace lolita

#endif /* E11520F9_5899_4E9B_8D95_F11263F73062 */
