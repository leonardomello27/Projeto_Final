/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: mpcACCsystem_private.h
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

#ifndef RTW_HEADER_mpcACCsystem_private_h_
#define RTW_HEADER_mpcACCsystem_private_h_
#include "rtwtypes.h"
#include "mpcACCsystem_types.h"
#include "mpcACCsystem.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
#define rtmIsMajorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
#define rtmIsMinorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val)           ((rtm)->Timing.t = (val))
#endif

extern real_T rt_roundd_snf(real_T u);
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
extern void mpcACCsys_DataTypeConversion_L0(real_T rtu_u, real_T *rty_y);

/* private model entry point functions */
extern void mpcACCsystem_derivatives(void);

#endif                                 /* RTW_HEADER_mpcACCsystem_private_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
