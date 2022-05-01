/*!
* \file   MicromorphicDamageII-generic.cxx
* \brief  This file implements the umat interface for the MicromorphicDamageII behaviour law
* \author Thomas Helfer , Jérémy Bleyer
* \date   21 / 09 / 2021
*/

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif /* NOMINMAX */
#include <windows.h>
#ifdef small
#undef small
#endif /* small */
#endif /* _WIN32 */

#ifndef MFRONT_SHAREDOBJ
#define MFRONT_SHAREDOBJ TFEL_VISIBILITY_EXPORT
#endif /* MFRONT_SHAREDOBJ */

#include<iostream>
#include<cstdlib>
#include"TFEL/Material/OutOfBoundsPolicy.hxx"
#include"TFEL/Math/t2tot2.hxx"
#include"TFEL/Math/t2tost2.hxx"
#include"TFEL/Material/MicromorphicDamageII.hxx"
#include"MFront/GenericBehaviour/Integrate.hxx"

#include"MFront/GenericBehaviour/MicromorphicDamageII-generic.hxx"

static tfel::material::OutOfBoundsPolicy&
MicromorphicDamageII_getOutOfBoundsPolicy(){
using namespace tfel::material;
static OutOfBoundsPolicy policy = None;
return policy;
}

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

MFRONT_SHAREDOBJ const char* 
MicromorphicDamageII_build_id = "";

MFRONT_SHAREDOBJ const char* 
MicromorphicDamageII_mfront_ept = "MicromorphicDamageII";

MFRONT_SHAREDOBJ const char* 
MicromorphicDamageII_tfel_version = "4.0.0-dev";

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_mfront_mkt = 1u;

MFRONT_SHAREDOBJ const char *
MicromorphicDamageII_mfront_interface = "Generic";

MFRONT_SHAREDOBJ const char *
MicromorphicDamageII_src = "bhv_micromorphic_damage.mfront";

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_nModellingHypotheses = 5u;

MFRONT_SHAREDOBJ const char * 
MicromorphicDamageII_ModellingHypotheses[5u] = {"AxisymmetricalGeneralisedPlaneStrain",
"Axisymmetrical",
"PlaneStrain",
"GeneralisedPlaneStrain",
"Tridimensional"};

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_nMainVariables = 2;
MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_nGradients = 2;

MFRONT_SHAREDOBJ int MicromorphicDamageII_GradientsTypes[2] = {2,
0};
MFRONT_SHAREDOBJ const char * MicromorphicDamageII_Gradients[2] = {"MicromorphicDamageGradient",
"MicromorphicDamage"};
MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_nThermodynamicForces = 2;

MFRONT_SHAREDOBJ int MicromorphicDamageII_ThermodynamicForcesTypes[2] = {2,
0};
MFRONT_SHAREDOBJ const char * MicromorphicDamageII_ThermodynamicForces[2] = {"MicromorphicDamageGradientDualForce",
"MicromorphicDamageDualForce"};
MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_nTangentOperatorBlocks = 4;

MFRONT_SHAREDOBJ const char * MicromorphicDamageII_TangentOperatorBlocks[4] = {"MicromorphicDamageGradientDualForce",
"MicromorphicDamageGradient","MicromorphicDamageDualForce","MicromorphicDamage"};
MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_BehaviourType = 0u;

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_BehaviourKinematic = 0u;

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_SymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_ElasticSymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_api_version = 1u;

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_UsableInPurelyImplicitResolution = 0;

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_nMaterialProperties = 3u;

MFRONT_SHAREDOBJ const char *MicromorphicDamageII_MaterialProperties[3u] = {"FractureEnergy",
"CharacteristicLength",
"PenalisationFactor"};

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_nInternalStateVariables = 4;
MFRONT_SHAREDOBJ const char * MicromorphicDamageII_InternalStateVariables[4] = {"Damage",
"DamageDissipatedEnergy","PenalisationEnergy","MicromorphicDamageGradientEnergy"};
MFRONT_SHAREDOBJ int MicromorphicDamageII_InternalStateVariablesTypes [] = {0,0,0,0};

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_nExternalStateVariables = 1;
MFRONT_SHAREDOBJ const char * MicromorphicDamageII_ExternalStateVariables[1] = {"EnergyReleaseRate"};
MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_nParameters = 2;
MFRONT_SHAREDOBJ const char * MicromorphicDamageII_Parameters[2] = {"minimal_time_step_scaling_factor",
"maximal_time_step_scaling_factor"};
MFRONT_SHAREDOBJ int MicromorphicDamageII_ParametersTypes [] = {0,0};

MFRONT_SHAREDOBJ double MicromorphicDamageII_minimal_time_step_scaling_factor_ParameterDefaultValue = 0.1;

MFRONT_SHAREDOBJ double MicromorphicDamageII_maximal_time_step_scaling_factor_ParameterDefaultValue = 1.7976931348623e+308;

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_requiresStiffnessTensor = 0;
MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_requiresThermalExpansionCoefficientTensor = 0;
MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_nInitializeFunctions= 0;

MFRONT_SHAREDOBJ const char * const * MicromorphicDamageII_InitializeFunctions = nullptr;


MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_ComputesInternalEnergy = 0;

MFRONT_SHAREDOBJ unsigned short MicromorphicDamageII_ComputesDissipatedEnergy = 1;

MFRONT_SHAREDOBJ void
MicromorphicDamageII_setOutOfBoundsPolicy(const int p){
if(p==0){
MicromorphicDamageII_getOutOfBoundsPolicy() = tfel::material::None;
} else if(p==1){
MicromorphicDamageII_getOutOfBoundsPolicy() = tfel::material::Warning;
} else if(p==2){
MicromorphicDamageII_getOutOfBoundsPolicy() = tfel::material::Strict;
} else {
std::cerr << "MicromorphicDamageII_setOutOfBoundsPolicy: invalid argument\n";
}
}

MFRONT_SHAREDOBJ int
MicromorphicDamageII_setParameter(const char *const key,const double value){
using tfel::material::MicromorphicDamageIIParametersInitializer;
auto& i = MicromorphicDamageIIParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int MicromorphicDamageII_AxisymmetricalGeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
using Behaviour = MicromorphicDamageII<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, MicromorphicDamageII_getOutOfBoundsPolicy());
return r;
} // end of MicromorphicDamageII_AxisymmetricalGeneralisedPlaneStrain

MFRONT_SHAREDOBJ int MicromorphicDamageII_Axisymmetrical(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICAL;
using Behaviour = MicromorphicDamageII<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, MicromorphicDamageII_getOutOfBoundsPolicy());
return r;
} // end of MicromorphicDamageII_Axisymmetrical

MFRONT_SHAREDOBJ int MicromorphicDamageII_PlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::PLANESTRAIN;
using Behaviour = MicromorphicDamageII<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, MicromorphicDamageII_getOutOfBoundsPolicy());
return r;
} // end of MicromorphicDamageII_PlaneStrain

MFRONT_SHAREDOBJ int MicromorphicDamageII_GeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::GENERALISEDPLANESTRAIN;
using Behaviour = MicromorphicDamageII<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, MicromorphicDamageII_getOutOfBoundsPolicy());
return r;
} // end of MicromorphicDamageII_GeneralisedPlaneStrain

MFRONT_SHAREDOBJ int MicromorphicDamageII_Tridimensional(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::TRIDIMENSIONAL;
using Behaviour = MicromorphicDamageII<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, MicromorphicDamageII_getOutOfBoundsPolicy());
return r;
} // end of MicromorphicDamageII_Tridimensional

#ifdef __cplusplus
}
#endif /* __cplusplus */

