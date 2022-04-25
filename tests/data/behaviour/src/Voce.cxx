/*!
* \file   Voce.cxx
* \brief  this file implements the Voce Behaviour.
*         File generated by tfel version 4.0.0-dev
* \author Ds
* \date   02 / 04 / 2021
 */

#include<string>
#include<cstring>
#include<sstream>
#include<fstream>
#include<stdexcept>

#include"TFEL/Raise.hxx"
#include"TFEL/Material/VoceBehaviourData.hxx"
#include"TFEL/Material/VoceIntegrationData.hxx"
#include"TFEL/Material/Voce.hxx"

namespace tfel::material{

VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer&
VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer::get()
{
static VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer i;
return i;
}

VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer::VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer()
{
// Reading parameters from a file
VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer::readParameters(*this,"Voce-parameters.txt");
VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer::readParameters(*this,"VoceAxisymmetricalGeneralisedPlaneStress-parameters.txt");
}

void
VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer::set(const char* const key,
const double v){
using namespace std;
if(::strcmp("epsilon",key)==0){
VoceParametersInitializer::get().set("epsilon",v);
} else if(::strcmp("theta",key)==0){
VoceParametersInitializer::get().set("theta",v);
} else if(::strcmp("YoungModulus",key)==0){
VoceParametersInitializer::get().set("YoungModulus",v);
} else if(::strcmp("PoissonRatio",key)==0){
VoceParametersInitializer::get().set("PoissonRatio",v);
} else if(::strcmp("RelativeValueForTheEquivalentStressLowerBoundDefinition",key)==0){
VoceParametersInitializer::get().set("RelativeValueForTheEquivalentStressLowerBoundDefinition",v);
} else if(::strcmp("ihr_R0_0",key)==0){
VoceParametersInitializer::get().set("ihr_R0_0",v);
} else if(::strcmp("ihr_Rinf_0",key)==0){
VoceParametersInitializer::get().set("ihr_Rinf_0",v);
} else if(::strcmp("ihr_b_0",key)==0){
VoceParametersInitializer::get().set("ihr_b_0",v);
} else if(::strcmp("ihr_R0_1",key)==0){
VoceParametersInitializer::get().set("ihr_R0_1",v);
} else if(::strcmp("ihr_H_1",key)==0){
VoceParametersInitializer::get().set("ihr_H_1",v);
} else if(::strcmp("minimal_time_step_scaling_factor",key)==0){
VoceParametersInitializer::get().set("minimal_time_step_scaling_factor",v);
} else if(::strcmp("maximal_time_step_scaling_factor",key)==0){
VoceParametersInitializer::get().set("maximal_time_step_scaling_factor",v);
} else if(::strcmp("numerical_jacobian_epsilon",key)==0){
VoceParametersInitializer::get().set("numerical_jacobian_epsilon",v);
} else {
tfel::raise("VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer::set: "
" no parameter named '"+std::string(key)+"'");
}
}

void
VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer::set(const char* const key,
const unsigned short v){
using namespace std;
if(::strcmp("iterMax",key)==0){
VoceParametersInitializer::get().set("iterMax",v);
} else {
tfel::raise("VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer::set: "
"no parameter named '"+std::string(key)+"'");
}
}

void
VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer::readParameters(VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer&,const char* const fn){
auto tokenize = [](const std::string& line){
std::istringstream tokenizer(line);
std::vector<std::string> tokens;
std::copy(std::istream_iterator<std::string>(tokenizer),
std::istream_iterator<std::string>(),
std::back_inserter(tokens));
return tokens;
};
std::ifstream f(fn);
if(!f){
return;
}
size_t ln = 1u;
while(!f.eof()){
auto line = std::string{};
std::getline(f,line);
auto tokens = tokenize(line);
auto throw_if = [ln,line,fn](const bool c,const std::string& m){
tfel::raise_if(c,"VoceAxisymmetricalGeneralisedPlaneStressParametersInitializer::readParameters: "
"error at line '"+std::to_string(ln)+"' "
"while reading parameter file '"+std::string(fn)+"'"
"("+m+")");
};
if(tokens.empty()){
continue;
}
if(tokens[0][0]=='#'){
continue;
}
throw_if(tokens.size()!=2u,"invalid number of tokens");
if("epsilon"==tokens[0]){
VoceParametersInitializer::get().set("epsilon",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("theta"==tokens[0]){
VoceParametersInitializer::get().set("theta",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("YoungModulus"==tokens[0]){
VoceParametersInitializer::get().set("YoungModulus",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("PoissonRatio"==tokens[0]){
VoceParametersInitializer::get().set("PoissonRatio",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("RelativeValueForTheEquivalentStressLowerBoundDefinition"==tokens[0]){
VoceParametersInitializer::get().set("RelativeValueForTheEquivalentStressLowerBoundDefinition",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("ihr_R0_0"==tokens[0]){
VoceParametersInitializer::get().set("ihr_R0_0",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("ihr_Rinf_0"==tokens[0]){
VoceParametersInitializer::get().set("ihr_Rinf_0",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("ihr_b_0"==tokens[0]){
VoceParametersInitializer::get().set("ihr_b_0",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("ihr_R0_1"==tokens[0]){
VoceParametersInitializer::get().set("ihr_R0_1",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("ihr_H_1"==tokens[0]){
VoceParametersInitializer::get().set("ihr_H_1",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("minimal_time_step_scaling_factor"==tokens[0]){
VoceParametersInitializer::get().set("minimal_time_step_scaling_factor",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("maximal_time_step_scaling_factor"==tokens[0]){
VoceParametersInitializer::get().set("maximal_time_step_scaling_factor",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("numerical_jacobian_epsilon"==tokens[0]){
VoceParametersInitializer::get().set("numerical_jacobian_epsilon",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("iterMax"==tokens[0]){
VoceParametersInitializer::get().set("iterMax",
VoceParametersInitializer::getUnsignedShort(tokens[0],tokens[1]));
} else {
throw_if(true,"invalid parameter '"+tokens[0]+"'");
}
}
}

VocePlaneStressParametersInitializer&
VocePlaneStressParametersInitializer::get()
{
static VocePlaneStressParametersInitializer i;
return i;
}

VocePlaneStressParametersInitializer::VocePlaneStressParametersInitializer()
{
// Reading parameters from a file
VocePlaneStressParametersInitializer::readParameters(*this,"Voce-parameters.txt");
VocePlaneStressParametersInitializer::readParameters(*this,"VocePlaneStress-parameters.txt");
}

void
VocePlaneStressParametersInitializer::set(const char* const key,
const double v){
using namespace std;
if(::strcmp("epsilon",key)==0){
VoceParametersInitializer::get().set("epsilon",v);
} else if(::strcmp("theta",key)==0){
VoceParametersInitializer::get().set("theta",v);
} else if(::strcmp("YoungModulus",key)==0){
VoceParametersInitializer::get().set("YoungModulus",v);
} else if(::strcmp("PoissonRatio",key)==0){
VoceParametersInitializer::get().set("PoissonRatio",v);
} else if(::strcmp("RelativeValueForTheEquivalentStressLowerBoundDefinition",key)==0){
VoceParametersInitializer::get().set("RelativeValueForTheEquivalentStressLowerBoundDefinition",v);
} else if(::strcmp("ihr_R0_0",key)==0){
VoceParametersInitializer::get().set("ihr_R0_0",v);
} else if(::strcmp("ihr_Rinf_0",key)==0){
VoceParametersInitializer::get().set("ihr_Rinf_0",v);
} else if(::strcmp("ihr_b_0",key)==0){
VoceParametersInitializer::get().set("ihr_b_0",v);
} else if(::strcmp("ihr_R0_1",key)==0){
VoceParametersInitializer::get().set("ihr_R0_1",v);
} else if(::strcmp("ihr_H_1",key)==0){
VoceParametersInitializer::get().set("ihr_H_1",v);
} else if(::strcmp("minimal_time_step_scaling_factor",key)==0){
VoceParametersInitializer::get().set("minimal_time_step_scaling_factor",v);
} else if(::strcmp("maximal_time_step_scaling_factor",key)==0){
VoceParametersInitializer::get().set("maximal_time_step_scaling_factor",v);
} else if(::strcmp("numerical_jacobian_epsilon",key)==0){
VoceParametersInitializer::get().set("numerical_jacobian_epsilon",v);
} else {
tfel::raise("VocePlaneStressParametersInitializer::set: "
" no parameter named '"+std::string(key)+"'");
}
}

void
VocePlaneStressParametersInitializer::set(const char* const key,
const unsigned short v){
using namespace std;
if(::strcmp("iterMax",key)==0){
VoceParametersInitializer::get().set("iterMax",v);
} else {
tfel::raise("VocePlaneStressParametersInitializer::set: "
"no parameter named '"+std::string(key)+"'");
}
}

void
VocePlaneStressParametersInitializer::readParameters(VocePlaneStressParametersInitializer&,const char* const fn){
auto tokenize = [](const std::string& line){
std::istringstream tokenizer(line);
std::vector<std::string> tokens;
std::copy(std::istream_iterator<std::string>(tokenizer),
std::istream_iterator<std::string>(),
std::back_inserter(tokens));
return tokens;
};
std::ifstream f(fn);
if(!f){
return;
}
size_t ln = 1u;
while(!f.eof()){
auto line = std::string{};
std::getline(f,line);
auto tokens = tokenize(line);
auto throw_if = [ln,line,fn](const bool c,const std::string& m){
tfel::raise_if(c,"VocePlaneStressParametersInitializer::readParameters: "
"error at line '"+std::to_string(ln)+"' "
"while reading parameter file '"+std::string(fn)+"'"
"("+m+")");
};
if(tokens.empty()){
continue;
}
if(tokens[0][0]=='#'){
continue;
}
throw_if(tokens.size()!=2u,"invalid number of tokens");
if("epsilon"==tokens[0]){
VoceParametersInitializer::get().set("epsilon",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("theta"==tokens[0]){
VoceParametersInitializer::get().set("theta",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("YoungModulus"==tokens[0]){
VoceParametersInitializer::get().set("YoungModulus",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("PoissonRatio"==tokens[0]){
VoceParametersInitializer::get().set("PoissonRatio",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("RelativeValueForTheEquivalentStressLowerBoundDefinition"==tokens[0]){
VoceParametersInitializer::get().set("RelativeValueForTheEquivalentStressLowerBoundDefinition",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("ihr_R0_0"==tokens[0]){
VoceParametersInitializer::get().set("ihr_R0_0",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("ihr_Rinf_0"==tokens[0]){
VoceParametersInitializer::get().set("ihr_Rinf_0",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("ihr_b_0"==tokens[0]){
VoceParametersInitializer::get().set("ihr_b_0",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("ihr_R0_1"==tokens[0]){
VoceParametersInitializer::get().set("ihr_R0_1",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("ihr_H_1"==tokens[0]){
VoceParametersInitializer::get().set("ihr_H_1",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("minimal_time_step_scaling_factor"==tokens[0]){
VoceParametersInitializer::get().set("minimal_time_step_scaling_factor",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("maximal_time_step_scaling_factor"==tokens[0]){
VoceParametersInitializer::get().set("maximal_time_step_scaling_factor",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("numerical_jacobian_epsilon"==tokens[0]){
VoceParametersInitializer::get().set("numerical_jacobian_epsilon",
VoceParametersInitializer::getDouble(tokens[0],tokens[1]));
} else if("iterMax"==tokens[0]){
VoceParametersInitializer::get().set("iterMax",
VoceParametersInitializer::getUnsignedShort(tokens[0],tokens[1]));
} else {
throw_if(true,"invalid parameter '"+tokens[0]+"'");
}
}
}

VoceParametersInitializer&
VoceParametersInitializer::get()
{
static VoceParametersInitializer i;
return i;
}

VoceParametersInitializer::VoceParametersInitializer()
{
this->epsilon = 1e-14;
this->theta = 1;
this->young = 206900000000;
this->nu = 0.29;
this->relative_value_for_the_equivalent_stress_lower_bound = 1e-12;
this->ihr_R0_0 = 450000000;
this->ihr_Rinf_0 = 715000000;
this->ihr_b_0 = 16.93;
this->ihr_R0_1 = 0;
this->ihr_H_1 = 129200000;
this->minimal_time_step_scaling_factor = 0.1;
this->maximal_time_step_scaling_factor = 1.7976931348623e+308;
this->numerical_jacobian_epsilon = 1e-15;
this->iterMax = 100;
// Reading parameters from a file
VoceParametersInitializer::readParameters(*this,"Voce-parameters.txt");
}

void
VoceParametersInitializer::set(const char* const key,
const double v){
using namespace std;
if(::strcmp("epsilon",key)==0){
this->epsilon = v;
} else if(::strcmp("theta",key)==0){
this->theta = v;
} else if(::strcmp("YoungModulus",key)==0){
this->young = v;
} else if(::strcmp("PoissonRatio",key)==0){
this->nu = v;
} else if(::strcmp("RelativeValueForTheEquivalentStressLowerBoundDefinition",key)==0){
this->relative_value_for_the_equivalent_stress_lower_bound = v;
} else if(::strcmp("ihr_R0_0",key)==0){
this->ihr_R0_0 = v;
} else if(::strcmp("ihr_Rinf_0",key)==0){
this->ihr_Rinf_0 = v;
} else if(::strcmp("ihr_b_0",key)==0){
this->ihr_b_0 = v;
} else if(::strcmp("ihr_R0_1",key)==0){
this->ihr_R0_1 = v;
} else if(::strcmp("ihr_H_1",key)==0){
this->ihr_H_1 = v;
} else if(::strcmp("minimal_time_step_scaling_factor",key)==0){
this->minimal_time_step_scaling_factor = v;
} else if(::strcmp("maximal_time_step_scaling_factor",key)==0){
this->maximal_time_step_scaling_factor = v;
} else if(::strcmp("numerical_jacobian_epsilon",key)==0){
this->numerical_jacobian_epsilon = v;
} else {
tfel::raise("VoceParametersInitializer::set: "
" no parameter named '"+std::string(key)+"'");
}
}

void
VoceParametersInitializer::set(const char* const key,
const unsigned short v){
using namespace std;
if(::strcmp("iterMax",key)==0){
this->iterMax = v;
} else {
tfel::raise("VoceParametersInitializer::set: "
"no parameter named '"+std::string(key)+"'");
}
}

double
VoceParametersInitializer::getDouble(const std::string& n,
const std::string& v)
{
double value;
std::istringstream converter(v);
converter >> value;
tfel::raise_if(!converter||(!converter.eof()),
"VoceParametersInitializer::getDouble: "
"can't convert '"+v+"' to double for parameter '"+ n+"'");
return value;
}

unsigned short
VoceParametersInitializer::getUnsignedShort(const std::string& n,
const std::string& v)
{
unsigned short value;
std::istringstream converter(v);
converter >> value;
tfel::raise_if(!converter||(!converter.eof()),
"VoceParametersInitializer::getUnsignedShort: "
"can't convert '"+v+"' to unsigned short for parameter '"+ n+"'");
return value;
}

void
VoceParametersInitializer::readParameters(VoceParametersInitializer& pi,const char* const fn){
auto tokenize = [](const std::string& line){
std::istringstream tokenizer(line);
std::vector<std::string> tokens;
std::copy(std::istream_iterator<std::string>(tokenizer),
std::istream_iterator<std::string>(),
std::back_inserter(tokens));
return tokens;
};
std::ifstream f(fn);
if(!f){
return;
}
size_t ln = 1u;
while(!f.eof()){
auto line = std::string{};
std::getline(f,line);
auto tokens = tokenize(line);
auto throw_if = [ln,line,fn](const bool c,const std::string& m){
tfel::raise_if(c,"VoceParametersInitializer::readParameters: "
"error at line '"+std::to_string(ln)+"' "
"while reading parameter file '"+std::string(fn)+"'"
"("+m+")");
};
if(tokens.empty()){
continue;
}
if(tokens[0][0]=='#'){
continue;
}
throw_if(tokens.size()!=2u,"invalid number of tokens");
if("epsilon"==tokens[0]){
pi.epsilon = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("theta"==tokens[0]){
pi.theta = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("YoungModulus"==tokens[0]){
pi.young = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("PoissonRatio"==tokens[0]){
pi.nu = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("RelativeValueForTheEquivalentStressLowerBoundDefinition"==tokens[0]){
pi.relative_value_for_the_equivalent_stress_lower_bound = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("ihr_R0_0"==tokens[0]){
pi.ihr_R0_0 = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("ihr_Rinf_0"==tokens[0]){
pi.ihr_Rinf_0 = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("ihr_b_0"==tokens[0]){
pi.ihr_b_0 = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("ihr_R0_1"==tokens[0]){
pi.ihr_R0_1 = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("ihr_H_1"==tokens[0]){
pi.ihr_H_1 = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("minimal_time_step_scaling_factor"==tokens[0]){
pi.minimal_time_step_scaling_factor = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("maximal_time_step_scaling_factor"==tokens[0]){
pi.maximal_time_step_scaling_factor = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("numerical_jacobian_epsilon"==tokens[0]){
pi.numerical_jacobian_epsilon = VoceParametersInitializer::getDouble(tokens[0],tokens[1]);
} else if("iterMax"==tokens[0]){
pi.iterMax = VoceParametersInitializer::getUnsignedShort(tokens[0],tokens[1]);
} else {
throw_if(true,"invalid parameter '"+tokens[0]+"'");
}
}
}

} // end of namespace tfel::material

