#ifndef FBA75E4B_A209_4011_87A7_FF975994C267
#define FBA75E4B_A209_4011_87A7_FF975994C267

namespace lolita::geometry
{

    /**
     * @brief 
     * 
     * @tparam T_ 
     */
    template<typename T_>
    concept FrameConcept = requires
    {

        { T_::getDimEuclidean() } -> std::same_as<lolita::Integer>;
        
    };

    template<Integer dim_>
    struct CartesianFrame
    {

        static constexpr
        Integer
        getDimEuclidean()
        {
            return dim_;
        }

        static constexpr
        Integer dim_euclidean_ = dim_;

    };

    struct AxiSymmetricFrame
    {

        static constexpr
        Integer
        getDimEuclidean()
        {
            return 2;
        }

        static constexpr
        Integer dim_euclidean_ = 2;

    };
    
} // namespace lolita::geometry


#endif /* FBA75E4B_A209_4011_87A7_FF975994C267 */
