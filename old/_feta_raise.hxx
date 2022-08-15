//
// Created by dsiedel on 10/02/2022.
//

#ifndef FETA__FETA_RAISE_HXX
#define FETA__FETA_RAISE_HXX

#include "new/_feta.hxx"
//#include <ostream>

namespace feta
{

    inline static
    void
    raise(
            Bool b
    )
    {
        if (b) {
            throw std::runtime_error("");
        }
    }

    template<typename T, typename... Ts>
    inline static
    void
    aassert(
            Bool
            b,
            T &&
            arg,
            Ts &&...
            args
    )
    {
        if (b) {
            StrgStream ss;
            ss << std::forward<T>(arg);
            ((ss << " " << std::forward<Ts>(args)), ...);
            throw std::runtime_error(ss.str());
        }
    }

}

#endif //FETA__FETA_RAISE_HXX
