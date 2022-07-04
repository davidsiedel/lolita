/*!
* \file   TFEL/Material/MicromorphicDamageII.hxx
* \brief  this file implements the MicromorphicDamageII Behaviour.
*         File generated by tfel version 4.0.0-dev
* \author Thomas Helfer , Jérémy Bleyer
* \date   21 / 09 / 2021
 */

#ifndef LIB_TFELMATERIAL_MICROMORPHICDAMAGEII_HXX
#define LIB_TFELMATERIAL_MICROMORPHICDAMAGEII_HXX

#include<string>
#include<iostream>
#include<limits>
#include<stdexcept>
#include<algorithm>

#include"TFEL/Raise.hxx"
#include"TFEL/PhysicalConstants.hxx"
#include"TFEL/Config/TFELConfig.hxx"
#include"TFEL/Config/TFELTypes.hxx"
#include"TFEL/TypeTraits/IsFundamentalNumericType.hxx"
#include"TFEL/TypeTraits/IsReal.hxx"
#include"TFEL/Math/General/IEEE754.hxx"
#include"TFEL/Material/MaterialException.hxx"
#include"TFEL/Material/MechanicalBehaviour.hxx"
#include"TFEL/Material/MechanicalBehaviourTraits.hxx"
#include"TFEL/Material/OutOfBoundsPolicy.hxx"
#include"TFEL/Material/BoundsCheck.hxx"
#include"TFEL/Material/IsotropicPlasticity.hxx"
#include"TFEL/Material/Lame.hxx"
#include"TFEL/Material/Hosford1972YieldCriterion.hxx"
#include"TFEL/Material/MicromorphicDamageIIBehaviourData.hxx"
#include"TFEL/Material/MicromorphicDamageIIIntegrationData.hxx"

#include "MFront/GenericBehaviour/State.hxx"
#include "MFront/GenericBehaviour/BehaviourData.hxx"
namespace tfel::material{

struct MicromorphicDamageIIParametersInitializer
{
static MicromorphicDamageIIParametersInitializer&
get();

double minimal_time_step_scaling_factor;
double maximal_time_step_scaling_factor;

void set(const char* const,const double);

/*!
 * \brief convert a string to double
 * \param[in] p : parameter
 * \param[in] v : value
 */
static double getDouble(const std::string&,const std::string&);
private :

MicromorphicDamageIIParametersInitializer();

MicromorphicDamageIIParametersInitializer(const MicromorphicDamageIIParametersInitializer&);

MicromorphicDamageIIParametersInitializer&
operator=(const MicromorphicDamageIIParametersInitializer&);
/*!
 * \brief read the parameters from the given file
 * \param[out] pi : parameters initializer
 * \param[in]  fn : file name
 */
static void readParameters(MicromorphicDamageIIParametersInitializer&,const char* const);
};

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis, typename NumericType, bool use_qt>
class MicromorphicDamageII;

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis hypothesis, typename NumericType>
std::ostream&
 operator <<(std::ostream&,const MicromorphicDamageII<hypothesis, NumericType, false>&);

/*!
* \class MicromorphicDamageII
* \brief This class implements the MicromorphicDamageII behaviour.
* \tparam hypothesis: modelling hypothesis.
* \tparam NumericType: numerical type.
* \author Thomas Helfer , Jérémy Bleyer
* \date   21 / 09 / 2021
*/
template<ModellingHypothesis::Hypothesis hypothesis,typename NumericType>
struct MicromorphicDamageII<hypothesis, NumericType, false> final
: public MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>,
public MicromorphicDamageIIBehaviourData<hypothesis, NumericType, false>,
public MicromorphicDamageIIIntegrationData<hypothesis, NumericType, false>
{

static constexpr unsigned short N = ModellingHypothesisToSpaceDimension<hypothesis>::value;

static_assert(N==1||N==2||N==3);
static_assert(tfel::typetraits::IsFundamentalNumericType<NumericType>::cond);
static_assert(tfel::typetraits::IsReal<NumericType>::cond);

friend std::ostream& operator<< <>(std::ostream&,const MicromorphicDamageII&);

static constexpr unsigned short TVectorSize = N;
typedef tfel::math::StensorDimeToSize<N> StensorDimeToSize;
static constexpr unsigned short StensorSize = StensorDimeToSize::value;
typedef tfel::math::TensorDimeToSize<N> TensorDimeToSize;
static constexpr unsigned short TensorSize = TensorDimeToSize::value;

using ushort =  unsigned short;
using Types = tfel::config::Types<N, NumericType, false>;
using Type = NumericType;
using real = typename Types::real;
using time = typename Types::time;
using length = typename Types::length;
using frequency = typename Types::frequency;
using speed = typename Types::speed;
using stress = typename Types::stress;
using strain = typename Types::strain;
using strainrate = typename Types::strainrate;
using stressrate = typename Types::stressrate;
using temperature = typename Types::temperature;
using thermalexpansion = typename Types::thermalexpansion;
using thermalconductivity = typename Types::thermalconductivity;
using massdensity = typename Types::massdensity;
using energydensity = typename Types::energydensity;
using TVector = typename Types::TVector;
using DisplacementTVector = typename Types::DisplacementTVector;
using ForceTVector = typename Types::ForceTVector;
using HeatFlux = typename Types::HeatFlux;
using TemperatureGradient = typename Types::TemperatureGradient;
using Stensor = typename Types::Stensor;
using StressStensor = typename Types::StressStensor;
using StressRateStensor = typename Types::StressRateStensor;
using StrainStensor = typename Types::StrainStensor;
using StrainRateStensor = typename Types::StrainRateStensor;
using FrequencyStensor = typename Types::FrequencyStensor;
using Tensor = typename Types::Tensor;
using DeformationGradientTensor = typename Types::DeformationGradientTensor;
using StressTensor = typename Types::StressTensor;
using StiffnessTensor = typename Types::StiffnessTensor;
using Stensor4 = typename Types::Stensor4;
using TangentOperator = tfel::math::tvector<(TVectorSize)*(TVectorSize)+(1)*(1),real>;
using PhysicalConstants = tfel::PhysicalConstants<NumericType, false>;

public :

typedef MicromorphicDamageIIBehaviourData<hypothesis, NumericType, false> BehaviourData;
typedef MicromorphicDamageIIIntegrationData<hypothesis, NumericType, false> IntegrationData;
typedef typename MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::SMFlag SMFlag;
typedef typename MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::SMType SMType;
using MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::ELASTIC;
using MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::SECANTOPERATOR;
using MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::TANGENTOPERATOR;
using MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::CONSISTENTTANGENTOPERATOR;
using MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::NOSTIFFNESSREQUESTED;
using IntegrationResult = typename MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::IntegrationResult;

using MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::SUCCESS;
using MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::FAILURE;
using MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::UNRELIABLE_RESULTS;

private :



#line 25 "bhv_micromorphic_damage.mfront"
real dd;
#line 28 "bhv_micromorphic_damage.mfront"
stress dYd;
#line 30 "bhv_micromorphic_damage.mfront"
stress dYtum_03C7__;
#line 32 "bhv_micromorphic_damage.mfront"
stress dYtum_2207__dtum_03C7__;

#line 38 "bhv_micromorphic_damage.mfront"
real Atum_03C7__;
#line 39 "bhv_micromorphic_damage.mfront"
real Htum_03C7__;
#line 40 "bhv_micromorphic_damage.mfront"
real dd_ddtum_03C7__;

time minimal_time_step_scaling_factor;
time maximal_time_step_scaling_factor;

//! Tangent operator;
TangentOperator Dt;
tfel::math::TMatrixView<N,N,real> db_ddtum_2207__dtum_03C7__;
real& da_dddtum_03C7__;
/*!
* \brief Update internal variables at end of integration
*/
void updateIntegrationVariables(){
}

/*!
* \brief Update internal variables at end of integration
*/
void updateStateVariables(){
this->d += this->dd;
this->Yd += this->dYd;
this->Ytum_03C7__ += this->dYtum_03C7__;
this->Ytum_2207__dtum_03C7__ += this->dYtum_2207__dtum_03C7__;
}

/*!
* \brief Update auxiliary state variables at end of integration
*/
void updateAuxiliaryStateVariables()
{}

//! \brief Default constructor (disabled)
MicromorphicDamageII() =delete ;
//! \brief Copy constructor (disabled)
MicromorphicDamageII(const MicromorphicDamageII&) = delete;
//! \brief Assignement operator (disabled)
MicromorphicDamageII& operator = (const MicromorphicDamageII&) = delete;

public:

/*!
* \brief Constructor
*/
MicromorphicDamageII(const MicromorphicDamageIIBehaviourData<hypothesis, NumericType, false>& src1,
const MicromorphicDamageIIIntegrationData<hypothesis, NumericType, false>& src2)
: MicromorphicDamageIIBehaviourData<hypothesis, NumericType, false>(src1),
MicromorphicDamageIIIntegrationData<hypothesis, NumericType, false>(src2),
dd(real(0)),
dYd(stress(0)),
dYtum_03C7__(stress(0)),
dYtum_2207__dtum_03C7__(stress(0)),
db_ddtum_2207__dtum_03C7__(Dt.begin()),
da_dddtum_03C7__(Dt[TVectorSize*TVectorSize])
{
using namespace std;
using namespace tfel::math;
using std::vector;
this->minimal_time_step_scaling_factor = MicromorphicDamageIIParametersInitializer::get().minimal_time_step_scaling_factor;
this->maximal_time_step_scaling_factor = MicromorphicDamageIIParametersInitializer::get().maximal_time_step_scaling_factor;
}

/*
 * \brief constructor for the Generic interface
 * \param[in] mgb_d: behaviour data
 */
MicromorphicDamageII(const mfront::gb::BehaviourData& mgb_d)
: MicromorphicDamageIIBehaviourData<hypothesis, NumericType, false>(mgb_d),
MicromorphicDamageIIIntegrationData<hypothesis, NumericType, false>(mgb_d),
dd(real(0)),
dYd(stress(0)),
dYtum_03C7__(stress(0)),
dYtum_2207__dtum_03C7__(stress(0)),
db_ddtum_2207__dtum_03C7__(Dt.begin()),
da_dddtum_03C7__(Dt[TVectorSize*TVectorSize])
{
using namespace std;
using namespace tfel::math;
using std::vector;
this->minimal_time_step_scaling_factor = MicromorphicDamageIIParametersInitializer::get().minimal_time_step_scaling_factor;
this->maximal_time_step_scaling_factor = MicromorphicDamageIIParametersInitializer::get().maximal_time_step_scaling_factor;
this-> tum_2207__dtum_03C7__ = tfel::math::map<TVector>(mgb_d.s0.gradients);
tfel::fsalgo::transform<TVectorSize>::exe(mgb_d.s1.gradients,mgb_d.s0.gradients,this->dtum_2207__dtum_03C7__.begin(),std::minus<real>());
this-> b = tfel::math::map<TVector>(mgb_d.s0.thermodynamic_forces);
this->dtum_03C7__ = mgb_d.s0.gradients[TVectorSize];
this->ddtum_03C7__ = mgb_d.s1.gradients[TVectorSize] - mgb_d.s0.gradients[TVectorSize];
this->a = mgb_d.s0.thermodynamic_forces[TVectorSize];
}

/*!
 * \ brief initialize the behaviour with user code
 */
void initialize(){
using namespace std;
using namespace tfel::math;
using std::vector;
#line 43 "bhv_micromorphic_damage.mfront"
this->Atum_03C7__ = this->Gc * this->l;
#line 44 "bhv_micromorphic_damage.mfront"
this->Htum_03C7__ = this->beta * this->Gc / this->l;
}

/*!
* \brief set the policy for "out of bounds" conditions
*/
void
setOutOfBoundsPolicy(const OutOfBoundsPolicy policy_value){
this->policy = policy_value;
} // end of setOutOfBoundsPolicy

/*!
* \return the modelling hypothesis
*/
constexpr ModellingHypothesis::Hypothesis
getModellingHypothesis() const{
return hypothesis;
} // end of getModellingHypothesis

/*!
* \brief check bounds
*/
void checkBounds() const{
} // end of checkBounds

IntegrationResult computePredictionOperator(const SMFlag,const SMType) override{
tfel::raise("MicromorphicDamageII::computePredictionOperator: "
"unsupported prediction operator flag");
}

time getMinimalTimeStepScalingFactor() const noexcept override{
  return this->minimal_time_step_scaling_factor;
}

std::pair<bool, time>
computeAPrioriTimeStepScalingFactor(const time current_time_step_scaling_factor) const override{
const auto time_scaling_factor = this->computeAPrioriTimeStepScalingFactorII();
return {time_scaling_factor.first,
        std::min(std::min(std::max(time_scaling_factor.second,
                                   this->minimal_time_step_scaling_factor),
                          this->maximal_time_step_scaling_factor),
                  current_time_step_scaling_factor)};
}

/*!
* \brief Integrate behaviour  over the time step
*/
IntegrationResult
integrate(const SMFlag smflag, const SMType smt) override{
using namespace std;
using namespace tfel::math;
raise_if(smflag!=MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::STANDARDTANGENTOPERATOR,
"invalid tangent operator flag");
bool computeTangentOperator_ = smt!=NOSTIFFNESSREQUESTED;
#line 48 "bhv_micromorphic_damage.mfront"
const auto Y_ets = this->Y + this->dY;
#line 49 "bhv_micromorphic_damage.mfront"
const auto r = 1 / (this->Htum_03C7__ + 2 * Y_ets + (this->Gc / this->l));
#line 51 "bhv_micromorphic_damage.mfront"
const auto d_tr = (2 * Y_ets + this->Htum_03C7__ * (this->dtum_03C7__ + this->ddtum_03C7__)) * r;
#line 52 "bhv_micromorphic_damage.mfront"
if (d_tr > this->d) {
#line 53 "bhv_micromorphic_damage.mfront"
if (d_tr > 1) {
#line 54 "bhv_micromorphic_damage.mfront"
this->d = 1;
#line 55 "bhv_micromorphic_damage.mfront"
this->dd_ddtum_03C7__ = real{};
#line 56 "bhv_micromorphic_damage.mfront"
} else {
#line 57 "bhv_micromorphic_damage.mfront"
this->d = d_tr;
#line 58 "bhv_micromorphic_damage.mfront"
this->dd_ddtum_03C7__ = r * this->Htum_03C7__;
#line 59 "bhv_micromorphic_damage.mfront"
}
#line 60 "bhv_micromorphic_damage.mfront"
} else {
#line 61 "bhv_micromorphic_damage.mfront"
this->dd_ddtum_03C7__ = real{};
#line 62 "bhv_micromorphic_damage.mfront"
}
#line 63 "bhv_micromorphic_damage.mfront"
this->a = -this->Htum_03C7__ * (this->d - this->dtum_03C7__ - this->ddtum_03C7__);
#line 64 "bhv_micromorphic_damage.mfront"
this->b = this->Atum_03C7__ * (this->tum_2207__dtum_03C7__ + this->dtum_2207__dtum_03C7__);
#line 66 "bhv_micromorphic_damage.mfront"
this->Yd = (this->Gc / (2 * this->l)) * power<2>(this->d);
#line 67 "bhv_micromorphic_damage.mfront"
this->Ytum_03C7__ = (this->Htum_03C7__ / 2) * power<2>(this->d - this->dtum_03C7__ - this->ddtum_03C7__);
#line 68 "bhv_micromorphic_damage.mfront"
this->Ytum_2207__dtum_03C7__ = (this->Atum_03C7__ / 2) * ((this->tum_2207__dtum_03C7__ + this->dtum_2207__dtum_03C7__) | (this->tum_2207__dtum_03C7__ + this->dtum_2207__dtum_03C7__));
this->updateIntegrationVariables();
this->updateStateVariables();
this->updateAuxiliaryStateVariables();
if(computeTangentOperator_){
if(!this->computeConsistentTangentOperator(smt)){
return MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::FAILURE;
}
}
return MechanicalBehaviour<MechanicalBehaviourBase::GENERALBEHAVIOUR,hypothesis, NumericType, false>::SUCCESS;
}

std::pair<bool, time>
computeAPosterioriTimeStepScalingFactor(const time current_time_step_scaling_factor) const override{
const auto time_scaling_factor = this->computeAPosterioriTimeStepScalingFactorII();
return {time_scaling_factor.first,
        std::min(std::min(std::max(time_scaling_factor.second,
                                   this->minimal_time_step_scaling_factor),
                          this->maximal_time_step_scaling_factor),
                 current_time_step_scaling_factor)};
}

/*!
* \brief Update the internal energy at end of the time step
* \param[in] Psi_s: internal energy at end of the time step
*/
void computeInternalEnergy(real& Psi_s) const
{
Psi_s=0;
}

/*!
* \brief Update the dissipated energy at end of the time step
* \param[in] Psi_d: dissipated energy at end of the time step
*/
void computeDissipatedEnergy(real& Psi_d) const{
using namespace std;
using namespace tfel::math;
#line 77 "bhv_micromorphic_damage.mfront"
Psi_d = this->Yd + this->Ytum_03C7__ + this->Ytum_2207__dtum_03C7__;
}

bool computeConsistentTangentOperator(const SMType smt){
using namespace std;
using namespace tfel::math;
using std::vector;
#line 72 "bhv_micromorphic_damage.mfront"
this->da_dddtum_03C7__ = this->Htum_03C7__ * (1 - this->dd_ddtum_03C7__);
#line 73 "bhv_micromorphic_damage.mfront"
this->db_ddtum_2207__dtum_03C7__ = this->Atum_03C7__ * tmatrix<N, N, real>::Id();
return true;
}

/*!
* \brief compute the sound velocity
* \param[in] rho_m0: mass density in the reference configuration
*/
speed computeSpeedOfSound(const massdensity&) const {
return speed(0);

}

const TangentOperator& getTangentOperator() const{
return this->Dt;
}

void updateExternalStateVariables(){
this->tum_2207__dtum_03C7__  += this->dtum_2207__dtum_03C7__;
this->dtum_03C7__  += this->ddtum_03C7__;
this->T += this->dT;
this->Y += this->dY;
}

//!
~MicromorphicDamageII()
 override = default;

private:

std::pair<bool, time> computeAPrioriTimeStepScalingFactorII() const{
return {true, this->maximal_time_step_scaling_factor};
}

std::pair<bool, time> computeAPosterioriTimeStepScalingFactorII() const{
return {true,this->maximal_time_step_scaling_factor};
}

//! policy for treating out of bounds conditions
OutOfBoundsPolicy policy = None;
}; // end of MicromorphicDamageII class

template<ModellingHypothesis::Hypothesis hypothesis, typename NumericType>
std::ostream&
operator <<(std::ostream& os,const MicromorphicDamageII<hypothesis, NumericType, false>& b)
{
os << "∇dχ : " << b.tum_2207__dtum_03C7__ << '\n';
os << "Δ∇dχ : " << b.dtum_2207__dtum_03C7__ << '\n';
os << "b : " << b.b << '\n';
os << "dχ : " << b.dtum_03C7__ << '\n';
os << "Δdχ : " << b.ddtum_03C7__ << '\n';
os << "a : " << b.a << '\n';
os << "Δt : " << b.dt << '\n';
os << "Gc : " << b.Gc << '\n';
os << "l : " << b.l << '\n';
os << "beta : " << b.beta << '\n';
os << "d : " << b.d << '\n';
os << "Δd : " << b.dd << '\n';
os << "Yd : " << b.Yd << '\n';
os << "ΔYd : " << b.dYd << '\n';
os << "Yχ : " << b.Ytum_03C7__ << '\n';
os << "ΔYχ : " << b.dYtum_03C7__ << '\n';
os << "Y∇dχ : " << b.Ytum_2207__dtum_03C7__ << '\n';
os << "ΔY∇dχ : " << b.dYtum_2207__dtum_03C7__ << '\n';
os << "T : " << b.T << '\n';
os << "ΔT : " << b.dT << '\n';
os << "Y : " << b.Y << '\n';
os << "ΔY : " << b.dY << '\n';
os << "minimal_time_step_scaling_factor : " << b.minimal_time_step_scaling_factor << '\n';
os << "maximal_time_step_scaling_factor : " << b.maximal_time_step_scaling_factor << '\n';
return os;
}

/*!
* Partial specialisation for MicromorphicDamageII.
*/
template<ModellingHypothesis::Hypothesis hypothesis, typename NumericType>
class MechanicalBehaviourTraits<MicromorphicDamageII<hypothesis, NumericType, false> >
{
static constexpr unsigned short N = ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr unsigned short TVectorSize = N;
typedef tfel::math::StensorDimeToSize<N> StensorDimeToSize;
static constexpr unsigned short StensorSize = StensorDimeToSize::value;
typedef tfel::math::TensorDimeToSize<N> TensorDimeToSize;
static constexpr unsigned short TensorSize = TensorDimeToSize::value;
public:
static constexpr bool is_defined = true;
static constexpr bool use_quantities = false;
static constexpr bool hasStressFreeExpansion = false;
static constexpr bool handlesThermalExpansion = false;
static constexpr unsigned short dimension = N;
static constexpr unsigned short material_properties_nb = 3;
static constexpr unsigned short internal_variables_nb  = 4;
static constexpr unsigned short external_variables_nb  = 2;
static constexpr unsigned short external_variables_nb2 = 1;
static constexpr bool hasConsistentTangentOperator = true;
static constexpr bool isConsistentTangentOperatorSymmetric = false;
static constexpr bool hasPredictionOperator = false;
static constexpr bool hasAPrioriTimeStepScalingFactor = false;
static constexpr bool hasComputeInternalEnergy = false;
static constexpr bool hasComputeDissipatedEnergy = true;
/*!
* \return the name of the class.
*/
static const char* getName(){
return "MicromorphicDamageII";
}

};

/*!
* Partial specialisation for MicromorphicDamageII.
*/
template<typename NumericType>
class MechanicalBehaviourTraits<MicromorphicDamageII<ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS, NumericType, false> >
{
public:
static constexpr bool is_defined = false;
static constexpr bool use_quantities = false;
static constexpr bool hasStressFreeExpansion = false;
static constexpr bool handlesThermalExpansion = false;
static constexpr unsigned short dimension = 0u;
static constexpr unsigned short material_properties_nb = 0;
static constexpr unsigned short internal_variables_nb  = 0;
static constexpr unsigned short external_variables_nb  = 0;
static constexpr unsigned short external_variables_nb2 = 0;
static constexpr bool hasConsistentTangentOperator = false;
static constexpr bool isConsistentTangentOperatorSymmetric = false;
static constexpr bool hasPredictionOperator = false;
static constexpr bool hasAPrioriTimeStepScalingFactor = false;
static constexpr bool hasComputeInternalEnergy = false;
static constexpr bool hasComputeDissipatedEnergy = false;
/*!
* \return the name of the class.
*/
static const char* getName(){
return "MicromorphicDamageII";
}

};

/*!
* Partial specialisation for MicromorphicDamageII.
*/
template<typename NumericType>
class MechanicalBehaviourTraits<MicromorphicDamageII<ModellingHypothesis::PLANESTRESS, NumericType, false> >
{
public:
static constexpr bool is_defined = false;
static constexpr bool use_quantities = false;
static constexpr bool hasStressFreeExpansion = false;
static constexpr bool handlesThermalExpansion = false;
static constexpr unsigned short dimension = 0u;
static constexpr unsigned short material_properties_nb = 0;
static constexpr unsigned short internal_variables_nb  = 0;
static constexpr unsigned short external_variables_nb  = 0;
static constexpr unsigned short external_variables_nb2 = 0;
static constexpr bool hasConsistentTangentOperator = false;
static constexpr bool isConsistentTangentOperatorSymmetric = false;
static constexpr bool hasPredictionOperator = false;
static constexpr bool hasAPrioriTimeStepScalingFactor = false;
static constexpr bool hasComputeInternalEnergy = false;
static constexpr bool hasComputeDissipatedEnergy = false;
/*!
* \return the name of the class.
*/
static const char* getName(){
return "MicromorphicDamageII";
}

};

} // end of namespace tfel::material

#endif /* LIB_TFELMATERIAL_MICROMORPHICDAMAGEII_HXX */