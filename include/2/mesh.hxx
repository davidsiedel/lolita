#ifndef FF0E5C46_3F4A_43C2_8E7D_79BAA3FA9102
#define FF0E5C46_3F4A_43C2_8E7D_79BAA3FA9102

#include "config.hxx"

namespace lolita
{

    template<typename Implementation_>
    struct MeshBase;

    namespace detail
    {
        
        struct MeshType
        {

        private:

            constexpr
            MeshType()
            {}

            constexpr
            Boolean
            operator==(
                MeshType const & other
            )
            const
            = default;

            constexpr
            Boolean
            operator!=(
                MeshType const & other
            )
            const
            = default;

            template<typename Implementation_>
            friend class lolita::MeshBase;

        };
        
    } // namespace detail

    template<typename T>
    concept MeshConcept = std::derived_from<std::decay_t<T>, detail::MeshType>;

    template<typename Implementation_>
    struct MeshBase : detail::MeshType
    {

    protected:

        constexpr
        MeshBase(
            Integer dim
        )
        :
        detail::MeshType(),
        dim_(dim)
        {}

    public:

        constexpr
        Boolean
        operator==(
            MeshBase const & other
        )
        const
        = default;

        constexpr
        Boolean
        operator!=(
            MeshBase const & other
        )
        const
        = default;

        template<typename... T>
        constexpr
        Boolean
        operator==(
            MeshBase<T...> const & other
        )
        const
        {
            return utility::areEqual(* this, other);
        }

        template<typename... T>
        constexpr
        Boolean
        operator!=(
            MeshBase<T...> const & other
        )
        const
        {
            return !(* this == other);
        }

        constexpr
        Boolean
        hasDim(
            Integer dim
        )
        const
        {
            return dim_ == dim;
        }

        constexpr
        Integer
        getDim()
        const
        {
            return dim_;
        }

        Integer dim_;

    };

    struct CartesianMesh : MeshBase<CartesianMesh>
    {

    private:

        using Base_ = MeshBase<CartesianMesh>;

    public:

        constexpr explicit
        CartesianMesh(
            Integer dim
        )
        :
        Base_(dim)
        {}

    };

    namespace detail
    {

        template<typename... T>
        struct IsCartesianMeshTraits : std::false_type
        {};

        template<>
        struct IsCartesianMeshTraits<CartesianMesh> : std::true_type
        {};
        
    } // namespace detail

    template<typename T>
    concept CartesianMeshConcept = detail::IsCartesianMeshTraits<T>::value;

    struct AxiSymmetricMesh : MeshBase<AxiSymmetricMesh>
    {

    private:

        using Base_ = MeshBase<AxiSymmetricMesh>;

    public:

        constexpr
        AxiSymmetricMesh()
        :
        Base_(2)
        {}

    };

    namespace detail
    {
        
        template<typename... T>
        struct IsAxiSymmetricMeshTraits : std::false_type
        {};

        template<>
        struct IsAxiSymmetricMeshTraits<AxiSymmetricMesh> : std::true_type
        {};

    } // namespace detail

    template<typename T>
    concept AxiSymmetricMeshConcept = detail::IsAxiSymmetricMeshTraits<T>::value;

} // namespace lolita

#endif /* FF0E5C46_3F4A_43C2_8E7D_79BAA3FA9102 */
