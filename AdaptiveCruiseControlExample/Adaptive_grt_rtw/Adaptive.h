/*
 * Adaptive.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "Adaptive".
 *
 * Model version              : 1.1
 * Simulink Coder version : 9.9 (R2023a) 19-Nov-2022
 * C source code generated on : Tue Aug  1 17:41:01 2023
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
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
#include "rt_logging.h"
#endif                                 /* Adaptive_COMMON_INCLUDES_ */

#include "Adaptive_types.h"
#include "rtGetNaN.h"
#include <float.h>
#include <string.h>
#include <stddef.h>
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetRTWLogInfo
#define rtmGetRTWLogInfo(rtm)          ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   ((rtm)->Timing.taskTime0)
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                (&(rtm)->Timing.taskTime0)
#endif

/* Block signals for system '<S1>/DataTypeConversion_L0' */
typedef struct {
  real_T y;                            /* '<S1>/DataTypeConversion_L0' */
} B_DataTypeConversion_L0_Adapt_T;

/* Block signals for system '<S1>/DataTypeConversion_dmin' */
typedef struct {
  real_T y;                            /* '<S1>/DataTypeConversion_dmin' */
} B_DataTypeConversion_dmin_Ada_T;

/* Block signals (default storage) */
typedef struct {
  B_DataTypeConversion_dmin_Ada_T sf_DataTypeConversion_vset;/* '<S1>/DataTypeConversion_vset' */
  B_DataTypeConversion_dmin_Ada_T sf_DataTypeConversion_vlead;/* '<S1>/DataTypeConversion_vlead' */
  B_DataTypeConversion_dmin_Ada_T sf_DataTypeConversion_vego;/* '<S1>/DataTypeConversion_vego' */
  B_DataTypeConversion_dmin_Ada_T sf_DataTypeConversion_reldist;/* '<S1>/DataTypeConversion_reldist' */
  B_DataTypeConversion_dmin_Ada_T sf_DataTypeConversion_dmin;/* '<S1>/DataTypeConversion_dmin' */
  B_DataTypeConversion_L0_Adapt_T sf_DataTypeConversion_atrack;/* '<S1>/DataTypeConversion_atrack' */
  B_DataTypeConversion_L0_Adapt_T sf_DataTypeConversion_amin;/* '<S1>/DataTypeConversion_amin' */
  B_DataTypeConversion_L0_Adapt_T sf_DataTypeConversion_amax;/* '<S1>/DataTypeConversion_amax' */
  B_DataTypeConversion_L0_Adapt_T sf_DataTypeConversion_L0;/* '<S1>/DataTypeConversion_L0' */
} B_Adaptive_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T last_mv_DSTATE;               /* '<S13>/last_mv' */
  real_T last_x_PreviousInput[4];      /* '<S13>/last_x' */
  boolean_T Memory_PreviousInput[96];  /* '<S13>/Memory' */
} DW_Adaptive_T;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T set_velocity;                 /* '<Root>/Set velocity' */
  real_T Timegap;                      /* '<Root>/Time gap' */
  real_T ego_velocity;                 /* '<Root>/Longitudinal velocity' */
  real_T relativedistance;             /* '<Root>/Relative distance' */
  real_T relative_velocity;            /* '<Root>/Relative velocity' */
} ExtU_Adaptive_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T Longitudinalacceleration;     /* '<Root>/Longitudinal acceleration' */
} ExtY_Adaptive_T;

/* Parameters (default storage) */
struct P_Adaptive_T_ {
  real_T AdaptiveCruiseControlSystem_Def;
                              /* Mask Parameter: AdaptiveCruiseControlSystem_Def
                               * Referenced by: '<S1>/Default spacing constant'
                               */
  real_T AdaptiveCruiseControlSystem_Max;
                              /* Mask Parameter: AdaptiveCruiseControlSystem_Max
                               * Referenced by: '<S1>/Maximum longitudinal acceleration constant'
                               */
  real_T AdaptiveCruiseControlSystem_M_j;
                              /* Mask Parameter: AdaptiveCruiseControlSystem_M_j
                               * Referenced by: '<S1>/Maximum velocity constant'
                               */
  real_T AdaptiveCruiseControlSystem_Min;
                              /* Mask Parameter: AdaptiveCruiseControlSystem_Min
                               * Referenced by: '<S1>/Minimum longitudinal acceleration constant'
                               */
  real_T last_x_InitialCondition[4];   /* Expression: lastx+xoff
                                        * Referenced by: '<S13>/last_x'
                                        */
  real_T last_mv_InitialCondition;     /* Expression: lastu+uoff
                                        * Referenced by: '<S13>/last_mv'
                                        */
  real_T Minimumvelocityconstant_Value;/* Expression: MinVelocity
                                        * Referenced by: '<S1>/Minimum velocity constant'
                                        */
  real_T Unconstrained_Value;          /* Expression: 0
                                        * Referenced by: '<S1>/Unconstrained'
                                        */
  real_T E_zero_Value;                 /* Expression: zeros(1,1)
                                        * Referenced by: '<S12>/E_zero'
                                        */
  real_T umin_scale4_Gain;         /* Expression: MVscale(:,ones(1,max(nCC,1)))'
                                    * Referenced by: '<S13>/umin_scale4'
                                    */
  real_T F_zero_Value[2];              /* Expression: zeros(1,2)
                                        * Referenced by: '<S12>/F_zero'
                                        */
  real_T ymin_scale1_Gain[2];       /* Expression: Yscale(:,ones(1,max(nCC,1)))'
                                     * Referenced by: '<S13>/ymin_scale1'
                                     */
  real_T G_zero_Value;                 /* Expression: zeros(1,1)
                                        * Referenced by: '<S12>/G_zero'
                                        */
  real_T S_zero_Value;                 /* Expression: zeros(1,1)
                                        * Referenced by: '<S12>/S_zero'
                                        */
  real_T ymin_scale2_Gain;         /* Expression: MDscale(:,ones(1,max(nCC,1)))'
                                    * Referenced by: '<S13>/ymin_scale2'
                                    */
  real_T Enableoptimizationconstant_Valu;/* Expression: 0
                                          * Referenced by: '<S1>/Enable optimization constant'
                                          */
  real_T extmv_zero_Value;             /* Expression: zeros(1,1)
                                        * Referenced by: '<S12>/ext.mv_zero'
                                        */
  real_T extmv_scale_Gain;             /* Expression: RMVscale
                                        * Referenced by: '<S13>/ext.mv_scale'
                                        */
  real_T mvtarget_zero_Value;          /* Expression: zeros(1,1)
                                        * Referenced by: '<S12>/mv.target_zero'
                                        */
  real_T extmv_scale1_Gain;            /* Expression: RMVscale
                                        * Referenced by: '<S13>/ext.mv_scale1'
                                        */
  real_T ywt_zero_Value[2];            /* Expression: zeros(2,1)
                                        * Referenced by: '<S12>/y.wt_zero'
                                        */
  real_T uwt_zero_Value;               /* Expression: zeros(1,1)
                                        * Referenced by: '<S12>/u.wt_zero'
                                        */
  real_T duwt_zero_Value;              /* Expression: zeros(1,1)
                                        * Referenced by: '<S12>/du.wt_zero'
                                        */
  real_T ecrwt_zero_Value;             /* Expression: zeros(1,1)
                                        * Referenced by: '<S12>/ecr.wt_zero'
                                        */
  real_T umin_scale1_Gain;             /* Expression: MVscale
                                        * Referenced by: '<S13>/umin_scale1'
                                        */
  real_T Externalcontrolsignalconstant_V;/* Expression: 0
                                          * Referenced by: '<S1>/External control signal constant'
                                          */
  boolean_T Memory_InitialCondition[96];/* Expression: iA
                                         * Referenced by: '<S13>/Memory'
                                         */
};

/* Real-time Model Data Structure */
struct tag_RTM_Adaptive_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T taskTime0;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    time_T tFinal;
    boolean_T stopRequestedFlag;
  } Timing;
};

/* Block parameters (default storage) */
extern P_Adaptive_T Adaptive_P;

/* Block signals (default storage) */
extern B_Adaptive_T Adaptive_B;

/* Block states (default storage) */
extern DW_Adaptive_T Adaptive_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_Adaptive_T Adaptive_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_Adaptive_T Adaptive_Y;

/* Model entry point functions */
extern void Adaptive_initialize(void);
extern void Adaptive_step(void);
extern void Adaptive_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Adaptive_T *const Adaptive_M;

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
 * hilite_system('Acc_teste1/Adaptive Cruise Control System')    - opens subsystem Acc_teste1/Adaptive Cruise Control System
 * hilite_system('Acc_teste1/Adaptive Cruise Control System/Kp') - opens and selects block Kp
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Acc_teste1'
 * '<S1>'   : 'Acc_teste1/Adaptive Cruise Control System'
 * '<S2>'   : 'Acc_teste1/Adaptive Cruise Control System/DataTypeConversion_L0'
 * '<S3>'   : 'Acc_teste1/Adaptive Cruise Control System/DataTypeConversion_amax'
 * '<S4>'   : 'Acc_teste1/Adaptive Cruise Control System/DataTypeConversion_amin'
 * '<S5>'   : 'Acc_teste1/Adaptive Cruise Control System/DataTypeConversion_atrack'
 * '<S6>'   : 'Acc_teste1/Adaptive Cruise Control System/DataTypeConversion_dmin'
 * '<S7>'   : 'Acc_teste1/Adaptive Cruise Control System/DataTypeConversion_optsgn'
 * '<S8>'   : 'Acc_teste1/Adaptive Cruise Control System/DataTypeConversion_reldist'
 * '<S9>'   : 'Acc_teste1/Adaptive Cruise Control System/DataTypeConversion_vego'
 * '<S10>'  : 'Acc_teste1/Adaptive Cruise Control System/DataTypeConversion_vlead'
 * '<S11>'  : 'Acc_teste1/Adaptive Cruise Control System/DataTypeConversion_vset'
 * '<S12>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC'
 * '<S13>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC'
 * '<S14>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check'
 * '<S15>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check1'
 * '<S16>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check2'
 * '<S17>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check'
 * '<S18>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check1'
 * '<S19>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check2'
 * '<S20>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check3'
 * '<S21>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check4'
 * '<S22>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check5'
 * '<S23>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check6'
 * '<S24>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check7'
 * '<S25>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check8'
 * '<S26>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check'
 * '<S27>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check1'
 * '<S28>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check2'
 * '<S29>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check'
 * '<S30>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check1'
 * '<S31>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check6'
 * '<S32>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/moorx'
 * '<S33>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/optimizer'
 * '<S34>'  : 'Acc_teste1/Adaptive Cruise Control System/MPC/MPC/optimizer/optimizer'
 */
#endif                                 /* RTW_HEADER_Adaptive_h_ */
