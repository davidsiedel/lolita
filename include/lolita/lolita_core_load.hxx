//
// Created by dsiedel on 18/04/22.
//

#ifndef LOLITA_LOLITA_CORE_LOAD_HXX
#define LOLITA_LOLITA_CORE_LOAD_HXX

#include "lolita/lolita_containers.hxx"
#include "lolita/lolita_matrix.hxx"
#include "lolita/lolita_core.hxx"

namespace lolita::core::load
{

    template<Indx D>
    using LoadFuntion = std::function<Real(Vector<Real, D> const &, Real const &)>;

    template<Indx D>
    auto const static constexpr zero = [] (Vector<Real, D> const &, Real const &) { return 0.0; };

    template<Indx D>
    struct LoadComponent
    {

        LoadComponent()
                :
                function([] (Vector<Real, D> const & point_arg, Real const & time_arg) { return 0.0; }),
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
                function([] (Vector<Real, D> const & point_arg, Real const & time_arg) { return 0.0; }),
                is_null(true),
                direction(direction_arg),
                load_type(load_type_arg)
        {}

        LoadComponent(
                LoadFuntion<D> &&
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
                LoadFuntion<D> const &
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

        LoadFuntion<D> function;

        Bool is_null;

        Indx direction;

        LoadType load_type;

    };

    /**
     * @brief
     * @tparam Fi
     * @tparam D
     */
    template<Indx D, FField F>
    struct Load2 : public Array<LoadComponent<D>, F.size()>
    {

        using Base = Array<LoadComponent<D>, F.size()>;

        Load2()
        {}

        Load2(
                auto ...
                load_component_args
        )
                :
                Base{std::forward<decltype(load_component_args)>(load_component_args)...}
        {
            static_assert(sizeof...(load_component_args) == F.size());
        }

    };

}

#endif //LOLITA_LOLITA_CORE_LOAD_HXX
