/**
 * \file FEMEnums.h
 * 31/08/2010 KR inital implementation
 *
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef FEMENUMS_H
#define FEMENUMS_H

#include <limits>
#include <string>
#include <list>

namespace FiniteElement
{
enum TemperatureUnit
{
    KELVIN,
    CELSIUS
};

/** \brief Types of physical processes supported by OpenGeoSys.
 * If you change this enum, make sure you apply the changes to
 * the functions convertPorcessType(), convertProcessTypeToString(),
 *  isFlowProcess() and isDeformationProcess()
 */
enum ProcessType
{
    INVALID_PROCESS = 0,  //!< INVALID_PROCESS
    AIR_FLOW,             //!< AIR_FLOW
	BIOLOGICAL,           //Guess what CMCD
    /// M process, single/multi-phase flow
    DEFORMATION,          //!< DEFORMATION
    DEFORMATION_DYNAMIC,  //!< ...
    /// C process, single/multi-phase flow
    DEFORMATION_FLOW,  //!< DEFORMATION_FLOW
    /// H2M monolithic
    DEFORMATION_H2,  //!< DEFORMATION_H2
    FLUID_FLOW,
    FLUID_MOMENTUM,  // BC only
    FLUX,
    /// H process, incompressible flow
    GROUNDWATER_FLOW,  //!< GROUNDWATER_FLOW
    /// T process, single/multi-phase flow
    HEAT_TRANSPORT,  //!< HEAT_TRANSPORT
    /// H process, incompressible flow
    LIQUID_FLOW,       //!< LIQUID_FLOW
    MASS_TRANSPORT,    //!< MASS_TRANSPORT
    MULTI_PHASE_FLOW,  //!< MULTI_PHASE_FLOW
    NO_PCS,            //!< NO_PCS
    /// H process, incompressible flow
    OVERLAND_FLOW,            //!< OVERLAND_FLOW
    PS_GLOBAL,                //!< PS_GLOBAL
    MULTI_COMPONENTIAL_FLOW,  //!< Fluid flow coupled with heat transport
    TNEQ,                     //!< reactive thermal nonequilibrium
    TES,                      //!< reactive thermal equilibrium
    RANDOM_WALK,              //!< RANDOM_WALK
    /// H process, incompressible flow
    RICHARDS_FLOW,  //!< RICHARDS_FLOW
    /// H2 process, compressible flow
    TWO_PHASE_FLOW,  //!< TWO_PHASE_FLOW
    // make sure that this is always the last entry (important for iterating
    // over the enum entries)!
    PROCESS_END
};

/**
 * \brief Convert the given string into the appropriate enum value.
 * @param pcs_type_string string describing a process type
 * @return enum value describing process type
 */
ProcessType convertProcessType(const std::string& pcs_type_string);

/**
 * \brief Convert the given enum value into the appropriate string.
 * @param pcs_type process type described by the enum ProcessType
 * @return string describing the process type
 */
std::string convertProcessTypeToString(ProcessType pcs_type);

/**
 * \brief Checks if the given pcs_type variable corresponds to a flow type of
 * the enum ProcessType.
 * @param pcs_type value of enum ProcessType
 * @return true if pcs_type describes a flow process, else false
 */
bool isFlowProcess(ProcessType pcs_type);

/**
 * \brief Checks if the given pcs_type variable corresponds to a multiphase flow
 * type of the enum ProcessType.
 * @param pcs_type value of enum ProcessType
 * @return true if pcs_type describes a flow process, else false
 */
bool isMultiFlowProcess(ProcessType pcs_type);

/**
 * \brief Checks if the given pcs_type variable corresponds to a deformation
 * type of the enum ProcessType.
 * @param pcs_type value of enum ProcessType
 * @return true if pcs_type describes a deformation process, else false
 */
bool isDeformationProcess(ProcessType pcs_type);

/// Returns a list of strings containing all entries in the ProcessType enum.
const std::list<std::string> getAllProcessNames();

/**
 * \brief Contains all values for primary variables actually handled by OGS.
 */
enum PrimaryVariable
{
    INVALID_PV = 0,   //!< INVALID_PV
    ACCELERATION_X1,  //!< ACCELERATION_X1
    ACCELERATION_Y1,  //!< ACCELERATION_Y1
    ACCELERATION_Z1,  //!< ACCELERATION_Z1
    /// Mass transport
    CONCENTRATION,  //!< CONCENTRATION
    /// Deformation
    DISPLACEMENT_X,  //!< DISPLACEMENT_X
    /// Deformation
    DISPLACEMENT_Y,  //!< DISPLACEMENT_Y
    /// Deformation
    DISPLACEMENT_Z,  //!< DISPLACEMENT_Z
    /// Deformation
    DISPLACEMENT_N,  //!< Normal traction (pressure type load). Only in ST
    EXCAVATION,      // ST
    HEAD,            //!< HEAD
    /// Flow (phase)
    PRESSURE,        //!< PRESSURE
    PRESSURE2,       //!< PRESSURE2
    PRESSURE_RATE1,  // OUT
    SATURATION,      //!< SATURATION
    SATURATION2,     //!< SATURATION2
    STRAIN_XX,       // Output
    STRAIN_XY,       // Output
    STRAIN_XZ,       // Output
    STRAIN_YY,       // Output
    STRAIN_YZ,       // Output
    STRAIN_ZZ,       // Output
    STRAIN_PLS,      // Output
    STRESS_XX,       // IC
    STRESS_XY,       // IC
    STRESS_XZ,       // IC
    STRESS_YY,       // IC
    STRESS_YZ,       // IC
    STRESS_ZZ,       // IC
    /// Heat transport
    TEMPERATURE,    //!< TEMPERATURE
    TEMPERATURE1,   //!< for TNEQ/TES models
    TEMPERATURE2,   //!< for TNEQ model
    VELOCITY_DM_X,  //!< VELOCITY_DM_X
    VELOCITY_DM_Y,  //!< VELOCITY_DM_Y
    VELOCITY_DM_Z,  //!< VELOCITY_DM_Z
    VELOCITY1_X,
    VELOCITY1_Y,
    VELOCITY1_Z,
    CARBON1,
    WATER1,
    METHANE1,
    NITROGEN1,

    // make sure that this is always the last entry (important for iterating
    // over the enum entries)!
    PV_END
};

/**
 * \brief Converts the given string into the appropriate enum value.
 * @param pcs_pv_string string describing the primary variable
 * @return enum value describing the primary variable of the process
 */
//!< PrimaryVariable
PrimaryVariable convertPrimaryVariable(const std::string& pcs_pv_string);

/**
 * \brief Converts the given enum value into the appropriate string.
 * @param pcs_pv primary variable described by the enum ProcessType
 * @return string describing the process type
 */
std::string convertPrimaryVariableToString(PrimaryVariable pcs_pv);

/// Returns a list of strings containing all entries in the PrimaryVariable
/// enum.
const std::list<std::string> getAllPrimaryVariableNames();

enum DistributionType
{
    INVALID_DIS_TYPE = 0,
    ANALYTICAL,  // ST
	BIO, //ST to link to Biological processes
    AVERAGE,
    CONSTANT,  // IC, BC, ST
    CONSTANT_GEO,
    CONSTANT_NEUMANN,  // ST
    SWITCH,
    CRITICALDEPTH,  // ST
    DIRECT,
    RECHARGE_DIRECT,
    FUNCTION,
    GRADIENT,        // IC
    GREEN_AMPT,      // ST
    RESTART,         // IC
    LINEAR,          // BC, ST
    LINEAR_NEUMANN,  // ST
    NORMALDEPTH,     // ST
    POINT,           // BC
	ELEMENT,
    PRECIPITATION,
    SYSTEM_DEPENDENT,  // ST
    TRANSFER_SURROUNDING,
	CAUCHY,
    NODESCONSTANT,  // IC (for DOMAIN keyword)
    CLIMATE,        // climate data (by NB)
    RECHARGE,       // MW
    // Sort of Neumann BC //WW
    // make sure that this is always the last entry (important for iterating
    // over the enum entries)! make sure that this is always the last entry
    // (important for iterating over the enum entries)!
    DIS_END
};

/**
 * \brief Converts the given string into the appropriate enum value.
 * @param pcs_pv_string string describing the primary variable
 * @return enum value describing the primary variable of the process
 */
DistributionType convertDisType(const std::string& dis_type_string);

/**
 * \brief Converts the given enum value into the appropriate string.
 * @param pcs_pv primary variable described by the enum ProcessType
 * @return string describing the process type
 */
std::string convertDisTypeToString(DistributionType dis_type);

/// Returns a list of strings containing all entries in the DistributionType
/// enum.
const std::list<std::string> getAllDistributionNames();

/** \brief Types of error method supported by OpenGeoSys.
 * If you change this enum, make sure you apply the changes to
 * the functions convertErrorMethod(), convertErrorMethodToString()
   Non-Linear and Coupling options (see also
 CRFProcess::CalcIterationNODError()):
   --> LMAX:	max(|x1-x0|)  -- Infinity norm: Local max error (across all
 elements) of solution vector delta (absolute error). Tolerance required for
 each primary variable.
   --> ENORM:	|x1-x0|       -- Euclidian norm: Norm of the solution vector
 delta (absolute error). Norm taken over entire solution vector (all primary
 variables) and checked against a single tolerance.
   --> EVNORM:	|x1-x0|       -- Euclidian varient norm: Norm of the solution
 vector delta (absolute error). Norm taken over solution vector of each primary
 variable, checked againes a tolerence specific to each variable.
   --> ERNORM:	|(x1-x0)/x0)| -- Euclidian Relative norm: Norm of the solution
 vector delta divided by the norm of the solution vector. A single tolerance
 applied to all primary variables.
   --> BNORM:	              -- OGS classic treatment of newton methods. ENORM
 error tolerance plus RHS ("B") control. (note: other error methods (i.e. ENORM)
 will also work well for NEWTON scheme)
 */
enum ErrorMethod
{
    INVALID_ERROR_METHOD = 0,
    LMAX,
    ENORM,
    EVNORM,
    ERNORM,
    BNORM
};

/**
 * \brief Convert the given string into the appropriate enum value.
 * @param pcs_type_string string describing an error method
 * @return enum value describing error method
 */
ErrorMethod convertErrorMethod(const std::string& error_method_string);

// Determines into which energy balance interphase friction is integrated in TES
// models
enum FrictionPhase
{
    INVALID_FRICTION_TYPE = 0,
    SOLID,
    FLUID,
    NONE
};

FrictionPhase convertFrictionPhase(const std::string& friction_string);
std::string convertFrictionPhaseToString(FrictionPhase friction_phase);

// Solid Reactive System
enum SolidReactiveSystem
{
    INVALID_REACTIVE_SYSTEM = 0,
    INERT,
    SINUSOIDAL,
    CaOH2,
    Mn3O4,
    Z13XBF
};

SolidReactiveSystem convertSolidReactiveSystem(
    const std::string& reactive_string);
std::string convertSolidReactiveSystemToString(
    SolidReactiveSystem reactive_system);

enum InitDataReadWriteType
{
    NO_IO = 0,
    READ,
    WRITE,
    READ_WRITE,
};

}  // end namespace FiniteElement

struct TimType
{
    enum type
    {
        INVALID_TIM_TYPE = 0,
        STEADY,
        TRANSIENT,
        PURERWPT
    };
};

TimType::type convertTimType(const std::string& str);

std::string convertTimTypeToString(TimType::type type);

struct IterationType
{
    enum type
    {
        INVALID,
        LINEAR,
        NONLINEAR,
        COUPLED,
        COUPLED_STABLE_ERROR
    };
};

IterationType::type convertIterationType(const std::string& str);

std::string convertIterationTypeToString(IterationType::type type);

struct TimeControlType
{
    enum type
    {
        INVALID,
        FIXED_STEPS,
        PI_AUTO_STEP_SIZE,
        DYNAMIC_VARIABLE,
        DYNAMIC_COURANT,
        DYNAMIC_PRESSURE,
        STEP_SIZE_RESTRICTION,
        NEUMANN,
        ERROR_CONTROL_ADAPTIVE,
        SELF_ADAPTIVE,
        STABLE_ERROR_ADAPTIVE
    };
};

TimeControlType::type convertTimeControlType(const std::string& str);

std::string convertTimeControlTypeToString(TimeControlType::type st_type);

struct ConstrainedType
{
    enum type
    {
        INVALID_CONSTRAINED_TYPE = 0,
        SMALLER,
        GREATER,
        POSITIVE,
        NEGATIVE
    };
};

ConstrainedType::type convertConstrainedType(const std::string& str);
std::string convertConstrainedTypeToString(
    ConstrainedType::type constrained_type);

struct ConstrainedVariable
{
    enum type
    {
        INVALID_CONSTRAINED_VARIABLE = 0,
        VELOCITY
    };
};

ConstrainedVariable::type convertConstrainedVariable(const std::string& str);
std::string convertConstrainedVariableToString(
    ConstrainedVariable::type constrained_variable);

#endif  // FEMENUMS_H
