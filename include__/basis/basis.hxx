/**
 * @file basis.hxx
 * @author David Siedel (davidsiedel@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2023-04-06
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef CAAD12E5_29A8_47C6_8F7D_3B18EC7C1C03
#define CAAD12E5_29A8_47C6_8F7D_3B18EC7C1C03

#include "numerics.hxx"
#include "geometry/domain.hxx"
#include "geometry/point.hxx"
#include "geometry/shape.hxx"
#include "geometry/frame.hxx"

namespace lolita::basis
{
    
    template<Integer k_>
    struct Monomial;

    template<typename T_>
    concept BasisConcept = requires(true);

    template<BasisConcept Basis_>
    struct BasisTraits;
    
} // namespace lolita::basis


#endif /* CAAD12E5_29A8_47C6_8F7D_3B18EC7C1C03 */
