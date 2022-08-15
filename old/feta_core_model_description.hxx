//
// Created by dsiedel on 19/03/2022.
//

#ifndef FETA_FETA_CORE_MODEL_DESCRIPTION_HXX
#define FETA_FETA_CORE_MODEL_DESCRIPTION_HXX

#include "new/_feta.hxx"
#include "new/feta_core_element_element_reference.hxx"

namespace feta::core::element
{

    enum struct Body
    {

        Cell,
        Face,
        Edge,
        Node

    };

    enum struct Model
    {

        Solid,
        Shell,
        Beam,

    };

    struct ModelDescription
    {

    private:

        static constexpr inline
        auto
        setBody(
                ElementDescription
                element_description_arg,
                Indx
                dim_euclidean_arg
        )
        {
            auto difference = dim_euclidean_arg - element_description_arg.shape_description.ord_shape;
            return static_cast<Body>(difference);
        }

    public:

        constexpr
        ModelDescription()
        :
        body(Body::Cell),
        model(Model::Solid)
        {}

        constexpr
        ModelDescription(
                Body
                body_arg,
                Model
                model_arg
        )
        :
        body(body_arg),
        model(model_arg)
        {}

        constexpr
        ModelDescription(
                ElementDescription
                element_description_arg,
                Indx
                dim_euclidean_arg,
                Model
                model_arg
        )
        :
        body(setBody(element_description_arg, dim_euclidean_arg)),
        model(model_arg)
        {}

        constexpr
        Bool
        operator==(
                ModelDescription const &
                other
        )
        const
        {
            auto eq_0 = body == other.body;
            auto eq_1 = model == other.model;
            return eq_0 && eq_1;
        }

        constexpr
        Bool
        operator!=(
                ModelDescription const &
                other
        )
        const
        {
            return !(other == * this);
        }

        Body body;

        Model model;

    };

    auto const static constexpr cell_solid = ModelDescription(Body::Cell, Model::Solid);
    auto const static constexpr face_solid = ModelDescription(Body::Face, Model::Solid);
    auto const static constexpr edge_solid = ModelDescription(Body::Edge, Model::Solid);
    auto const static constexpr node_solid = ModelDescription(Body::Node, Model::Solid);

}

#endif //FETA_FETA_CORE_MODEL_DESCRIPTION_HXX
