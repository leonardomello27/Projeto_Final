/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Adaptive.h
 *
 * Code generated for Simulink model 'Adaptive'.
 *
 * Model version                  : 11.0
 * Simulink Coder version         : 9.9 (R2023a) 19-Nov-2022
 * C/C++ source code generated on : Tue Aug  1 17:39:14 2023
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_Adaptive_h_
#define RTW_HEADER_Adaptive_h_
#ifndef Adaptive_COMMON_INCLUDES_
#define Adaptive_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* Adaptive_COMMON_INCLUDES_ */

#include "Adaptive_types.h"
#include "rtGetNaN.h"
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T last_mv_DSTATE;               /* '<S13>/last_mv' */
  real_T last_x_PreviousInput[4];      /* '<S13>/last_x' */
  boolean_T Memory_PreviousInput[96];  /* '<S13>/Memory' */
} DW_Adaptive_T;

/* Invariant block signals (default storage) */
typedef struct {
  const real_T MathFunction[2];        /* '<S13>/Math Function' */
  const real_T MathFunction1;          /* '<S13>/Math Function1' */
  const real_T MathFunction2;          /* '<S13>/Math Function2' */
} ConstB_Adaptive_T;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T set_velocity;                 /* '<Root>/Set velocity' */
  real_T Timegap;                      /* '<Root>/Time gap' */
  real_T v_ego;                        /* '<Root>/Longitudinal velocity' */
  real_T d_rel;                        /* '<Root>/Relative distance' */
  real_T relative_velocity;            /* '<Root>/Relative velocity' */
} ExtU_Adaptive_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T Longitudinalacceleration;     /* '<Root>/Longitudinal acceleration' */
} ExtY_Adaptive_T;

/* Real-time Model Data Structure */
struct tag_RTM_Adaptive_T {
  const char_T * volatile errorStatus;
};

/* Block states (default storage) */
extern DW_Adaptive_T Adaptive_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_Adaptive_T Adaptive_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_Adaptive_T Adaptive_Y;
extern const ConstB_Adaptive_T Adaptive_ConstB;/* constant block i/o */

/* Model entry point functions */
extern void Adaptive_initialize(void);
extern void Adaptive_step(void);
extern void Adaptive_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Adaptive_T *const Adaptive_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S13>/Constant' : Unused code path elimination
 * Block '<S13>/Floor' : Unused code path elimination
 * Block '<S13>/Floor1' : Unused code path elimination
 * Block '<S14>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S15>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S16>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S17>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S18>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S19>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S20>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S21>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S22>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S23>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S24>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S25>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S26>/Vector Dimension Check' : Unused code path elimination
 * Block '<S27>/Vector Dimension Check' : Unused code path elimination
 * Block '<S28>/Vector Dimension Check' : Unused code path elimination
 * Block '<S29>/Vector Dimension Check' : Unused code path elimination
 * Block '<S30>/Vector Dimension Check' : Unused code path elimination
 * Block '<S31>/Vector Dimension Check' : Unused code path elimination
 * Block '<S13>/Min' : Unused code path elimination
 * Block '<S13>/constant' : Unused code path elimination
 * Block '<S32>/Vector Dimension Check' : Unused code path elimination
 * Block '<S13>/umin_scale2' : Unused code path elimination
 * Block '<S13>/umin_scale3' : Unused code path elimination
 * Block '<S13>/umin_scale5' : Unused code path elimination
 * Block '<S13>/ym_zero' : Unused code path elimination
 * Block '<S12>/m_zero' : Unused code path elimination
 * Block '<S12>/p_zero' : Unused code path elimination
 * Block '<S13>/Reshape' : Reshape block reduction
 * Block '<S13>/Reshape1' : Reshape block reduction
 * Block '<S13>/Reshape2' : Reshape block reduction
 * Block '<S13>/Reshape3' : Reshape block reduction
 * Block '<S13>/Reshape4' : Reshape block reduction
 * Block '<S13>/Reshape5' : Reshape block reduction
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Note that this particular code originates from a subsystem build,
 * and has its own system numbers different from the parent model.
 * Refer to the system hierarchy for this subsystem below, and use the
 * MATLAB hilite_system command to trace the generated code back
 * to the parent model.  For example,
 *
 * hilite_system('mpcACCsystem/Adaptive Cruise Control System')    - opens subsystem mpcACCsystem/Adaptive Cruise Control System
 * hilite_system('mpcACCsystem/Adaptive Cruise Control System/Kp') - opens and selects block Kp
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'mpcACCsystem'
 * '<S1>'   : 'mpcACCsystem/Adaptive Cruise Control System'
 * '<S2>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_L0'
 * '<S3>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_amax'
 * '<S4>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_amin'
 * '<S5>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_atrack'
 * '<S6>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_dmin'
 * '<S7>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_optsgn'
 * '<S8>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_reldist'
 * '<S9>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_vego'
 * '<S10>'  : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_vlead'
 * '<S11>'  : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_vset'
 * '<S12>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC'
 * '<S13>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC'
 * '<S14>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check'
 * '<S15>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check1'
 * '<S16>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check2'
 * '<S17>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check'
 * '<S18>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check1'
 * '<S19>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check2'
 * '<S20>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check3'
 * '<S21>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check4'
 * '<S22>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check5'
 * '<S23>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check6'
 * '<S24>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check7'
 * '<S25>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check8'
 * '<S26>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check'
 * '<S27>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check1'
 * '<S28>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check2'
 * '<S29>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check'
 * '<S30>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check1'
 * '<S31>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check6'
 * '<S32>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/moorx'
 * '<S33>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/optimizer'
 * '<S34>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/optimizer/optimizer'
 */
#endif                                 /* RTW_HEADER_Adaptive_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
