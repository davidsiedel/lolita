#ifndef CD8635CE_D0A2_4DEA_AA85_D702F7CDE176
#define CD8635CE_D0A2_4DEA_AA85_D702F7CDE176

#include "2/core/traits/_include.hxx"

namespace lolita::core
{

    template<DomainConcept auto domain_>
    struct DomainView
    {
        
        static constexpr
        DomainConcept auto const &
        getDomain()
        {
            return domain_;
        }

    };
    
} // namespace lolita::core

#endif /* CD8635CE_D0A2_4DEA_AA85_D702F7CDE176 */
