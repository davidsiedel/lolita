#ifndef BDE0C330_8A4B_46B0_8810_63FE2BF2AD37
#define BDE0C330_8A4B_46B0_8810_63FE2BF2AD37

#include "config.hxx"
#include "algebra.hxx"

namespace lolita
{

    using Point = DenseVector<Real, 3>;

    template<typename T>
    concept PointConcept = DenseVectorConcept<T, Real, 3>;

} // namespace lolita

#endif /* BDE0C330_8A4B_46B0_8810_63FE2BF2AD37 */
