/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: mpcACCsystem.h
 *
 * Code generated for Simulink model 'mpcACCsystem'.
 *
 * Model version                  : 11.0
 * Simulink Coder version         : 9.9 (R2023a) 19-Nov-2022
 * C/C++ source code generated on : Tue Aug  1 17:32:36 2023
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_mpcACCsystem_h_
#define RTW_HEADER_mpcACCsystem_h_
#ifndef mpcACCsystem_COMMON_INCLUDES_
#define mpcACCsystem_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* mpcACCsystem_COMMON_INCLUDES_ */

#include "mpcACCsystem_types.h"
#include "rtGetNaN.h"
#include <string.h>
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
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
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T Sum1;                         /* '<S2>/Sum1' */
  real_T Sum1_m;                       /* '<S3>/Sum1' */
  real_T umin_scale1;                  /* '<S15>/umin_scale1' */
  real_T Integrator;                   /* '<S2>/Integrator' */
  real_T Integrator_c;                 /* '<S3>/Integrator' */
  real_T a_lead;                       /* '<Root>/Sine Wave' */
  real_T xk1[4];                       /* '<S35>/optimizer' */
  real_T u;                            /* '<S35>/optimizer' */
  boolean_T iAout[96];                 /* '<S35>/optimizer' */
} B_mpcACCsystem_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T last_mv_DSTATE;               /* '<S15>/last_mv' */
  real_T last_x_PreviousInput[4];      /* '<S15>/last_x' */
  real_T lastSin;                      /* '<Root>/Sine Wave' */
  real_T lastCos;                      /* '<Root>/Sine Wave' */
  int32_T systemEnable;                /* '<Root>/Sine Wave' */
  boolean_T Memory_PreviousInput[96];  /* '<S15>/Memory' */
} DW_mpcACCsystem_T;

/* Continuous states (default storage) */
typedef struct {
  real_T TransferFcn_CSTATE;           /* '<S2>/Transfer Fcn' */
  real_T Integrator1_CSTATE;           /* '<S3>/Integrator1' */
  real_T Integrator1_CSTATE_p;         /* '<S2>/Integrator1' */
  real_T TransferFcn_CSTATE_i;         /* '<S3>/Transfer Fcn' */
  real_T Integrator_CSTATE;            /* '<S2>/Integrator' */
  real_T Integrator_CSTATE_a;          /* '<S3>/Integrator' */
} X_mpcACCsystem_T;

/* State derivatives (default storage) */
typedef struct {
  real_T TransferFcn_CSTATE;           /* '<S2>/Transfer Fcn' */
  real_T Integrator1_CSTATE;           /* '<S3>/Integrator1' */
  real_T Integrator1_CSTATE_p;         /* '<S2>/Integrator1' */
  real_T TransferFcn_CSTATE_i;         /* '<S3>/Transfer Fcn' */
  real_T Integrator_CSTATE;            /* '<S2>/Integrator' */
  real_T Integrator_CSTATE_a;          /* '<S3>/Integrator' */
} XDot_mpcACCsystem_T;

/* State disabled  */
typedef struct {
  boolean_T TransferFcn_CSTATE;        /* '<S2>/Transfer Fcn' */
  boolean_T Integrator1_CSTATE;        /* '<S3>/Integrator1' */
  boolean_T Integrator1_CSTATE_p;      /* '<S2>/Integrator1' */
  boolean_T TransferFcn_CSTATE_i;      /* '<S3>/Transfer Fcn' */
  boolean_T Integrator_CSTATE;         /* '<S2>/Integrator' */
  boolean_T Integrator_CSTATE_a;       /* '<S3>/Integrator' */
} XDis_mpcACCsystem_T;

/* Invariant block signals (default storage) */
typedef struct {
  const real_T MathFunction[2];        /* '<S15>/Math Function' */
  const real_T MathFunction1;          /* '<S15>/Math Function1' */
  const real_T MathFunction2;          /* '<S15>/Math Function2' */
} ConstB_mpcACCsystem_T;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* Real-time Model Data Structure */
struct tag_RTM_mpcACCsystem_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_mpcACCsystem_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_mpcACCsystem_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[6];
  real_T odeF[3][6];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block signals (default storage) */
extern B_mpcACCsystem_T mpcACCsystem_B;

/* Continuous states (default storage) */
extern X_mpcACCsystem_T mpcACCsystem_X;

/* Block states (default storage) */
extern DW_mpcACCsystem_T mpcACCsystem_DW;
extern const ConstB_mpcACCsystem_T mpcACCsystem_ConstB;/* constant block i/o */

/* Model entry point functions */
extern void mpcACCsystem_initialize(void);
extern void mpcACCsystem_step(void);
extern void mpcACCsystem_terminate(void);

/* Real-time Model object */
extern RT_MODEL_mpcACCsystem_T *const mpcACCsystem_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S15>/Constant' : Unused code path elimination
 * Block '<S15>/Floor' : Unused code path elimination
 * Block '<S15>/Floor1' : Unused code path elimination
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
 * Block '<S26>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S27>/Matrix Dimension Check' : Unused code path elimination
 * Block '<S28>/Vector Dimension Check' : Unused code path elimination
 * Block '<S29>/Vector Dimension Check' : Unused code path elimination
 * Block '<S30>/Vector Dimension Check' : Unused code path elimination
 * Block '<S31>/Vector Dimension Check' : Unused code path elimination
 * Block '<S32>/Vector Dimension Check' : Unused code path elimination
 * Block '<S33>/Vector Dimension Check' : Unused code path elimination
 * Block '<S15>/Min' : Unused code path elimination
 * Block '<S15>/constant' : Unused code path elimination
 * Block '<S34>/Vector Dimension Check' : Unused code path elimination
 * Block '<S15>/umin_scale2' : Unused code path elimination
 * Block '<S15>/umin_scale3' : Unused code path elimination
 * Block '<S15>/umin_scale5' : Unused code path elimination
 * Block '<S15>/ym_zero' : Unused code path elimination
 * Block '<S14>/m_zero' : Unused code path elimination
 * Block '<S14>/p_zero' : Unused code path elimination
 * Block '<S15>/Reshape' : Reshape block reduction
 * Block '<S15>/Reshape1' : Reshape block reduction
 * Block '<S15>/Reshape2' : Reshape block reduction
 * Block '<S15>/Reshape3' : Reshape block reduction
 * Block '<S15>/Reshape4' : Reshape block reduction
 * Block '<S15>/Reshape5' : Reshape block reduction
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'mpcACCsystem'
 * '<S1>'   : 'mpcACCsystem/Adaptive Cruise Control System'
 * '<S2>'   : 'mpcACCsystem/Ego Car'
 * '<S3>'   : 'mpcACCsystem/Lead Car'
 * '<S4>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_L0'
 * '<S5>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_amax'
 * '<S6>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_amin'
 * '<S7>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_atrack'
 * '<S8>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_dmin'
 * '<S9>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_optsgn'
 * '<S10>'  : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_reldist'
 * '<S11>'  : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_vego'
 * '<S12>'  : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_vlead'
 * '<S13>'  : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_vset'
 * '<S14>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC'
 * '<S15>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC'
 * '<S16>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check'
 * '<S17>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check1'
 * '<S18>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check2'
 * '<S19>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check'
 * '<S20>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check1'
 * '<S21>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check2'
 * '<S22>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check3'
 * '<S23>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check4'
 * '<S24>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check5'
 * '<S25>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check6'
 * '<S26>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check7'
 * '<S27>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check8'
 * '<S28>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check'
 * '<S29>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check1'
 * '<S30>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check2'
 * '<S31>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check'
 * '<S32>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check1'
 * '<S33>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check6'
 * '<S34>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/moorx'
 * '<S35>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/optimizer'
 * '<S36>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/optimizer/optimizer'
 */
#endif                                 /* RTW_HEADER_mpcACCsystem_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
