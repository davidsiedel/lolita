//
// Created by dsiedel on 12/02/2022.
//

#ifndef FETA_FETA_CORE_LOAD_HXX
#define FETA_FETA_CORE_LOAD_HXX

//#include "new/load_component.hxx"
#include "new/feta_core_field_description.hxx"

namespace feta::core
{

    enum struct LoadType
    {

        Volumetric,
        Dirichlet,
        Neumann

    };

    template<Indx D>
    struct LoadComponent
    {

        auto const static constexpr dim_euclidean = D;

    private:

        using FunctionT = std::function<Real(Matrix<Real, D> const &, Real const &)>;

    public:

        LoadComponent()
        :
        function([] (Matrix<Real, D> const & point_arg, Real const & time_arg) { return 0.0; }),
        is_null(true),
        direction(0),
        load_type(LoadType::Volumetric)
        {}

        LoadComponent(
                Indx
                direction_arg,
                LoadType
                load_type_arg
        )
        :
        function([] (Matrix<Real, D> const & point_arg, Real const & time_arg) { return 0.0; }),
        is_null(true),
        direction(direction_arg),
        load_type(load_type_arg)
        {}

        LoadComponent(
                FunctionT &&
                function_arg,
                Indx
                direction_arg,
                LoadType
                load_type_arg
        )
        :
        function(function_arg),
        is_null(false),
        direction(direction_arg),
        load_type(load_type_arg)
        {}

        LoadComponent(
                FunctionT const &
                function_arg,
                Indx
                direction_arg,
                LoadType
                load_type_arg
        )
        :
        function(function_arg),
        is_null(false),
        direction(direction_arg),
        load_type(load_type_arg)
        {}

        Real
        getImposedValue(
                Matrix<Real, D> const &
                point_arg,
                Real const &
                time_arg
        )
        const
        {
            return function(point_arg, time_arg);
        }

        FunctionT function;

        Bool is_null;

        Indx direction;

        LoadType load_type;

    };

    /**
     * @brief
     * @tparam Fi
     * @tparam D
     */
    template<Indx D, FieldDescription F>
    struct Load2 : public Array<LoadComponent<D>, F.tensor_description.size>
    {

        using Base = Array<LoadComponent<D>, F.tensor_description.size>;

        Load2()
        {}

        Load2(
                auto ...
                load_component_args
        )
        :
        Base{std::forward<decltype(load_component_args)>(load_component_args)...}
        {
            static_assert(sizeof...(load_component_args) == F.tensor_description.size);
        }

    };

    template<typename... L>
    struct Load : public Array<typename UniqueType<L...>::Type, sizeof...(L)>
    {

        using HHH = typename UniqueType<L...>::Type;

//    LoadComponent<UniqueType<L...>::Type::dim_euclidean>

        static_assert(UniqueType<L...>::value);

        Load()
        {}

        Load(
                L &&...
                load_component_args
        )
        {
            this->data = {std::forward<typename UniqueType<L...>::Type>(load_component_args)...};
//            this->data = {load_component_args...};
        }

    };

    template<auto FiniteElementDescriptionArg>
    using LoadComponent1 = LoadComponent<FiniteElementDescriptionArg.getDimEuclidean()>;


    template<auto FiniteElementDescriptionArg>
    using Load1 = Load2<FiniteElementDescriptionArg.getDimEuclidean(), FiniteElementDescriptionArg.getFieldType()>;


    template<
            auto MixedElementDescriptionArg
    >
    using Load3 = typename RawType<
            decltype(MixedElementDescriptionArg)
    >::template Element<
            Load1
    >;

}

#endif //FETA_FETA_CORE_LOAD_HXX
