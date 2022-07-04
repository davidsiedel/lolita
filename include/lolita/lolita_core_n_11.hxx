#ifndef C31C92CF_B3D3_4D98_9746_7FD9EEFC72FC
#define C31C92CF_B3D3_4D98_9746_7FD9EEFC72FC

#include <execution>
#include "lolita/lolita_core_n_1.hxx"

namespace lolita2::geometry
{
    
    template<Domain _domain>
    struct Mesh;

    template<Element t_element, Domain t_domain>
    struct FiniteElement : FiniteElementConnectivity<t_element, t_domain>
    {

    };

}


#endif /* C31C92CF_B3D3_4D98_9746_7FD9EEFC72FC */
