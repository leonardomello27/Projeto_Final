/*
 * Adaptive_private.h
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

#ifndef RTW_HEADER_Adaptive_private_h_
#define RTW_HEADER_Adaptive_private_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#include "Adaptive.h"
#include "Adaptive_types.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmSetTFinal
#define rtmSetTFinal(rtm, val)         ((rtm)->Timing.tFinal = (val))
#endif

extern real_T rt_roundd_snf(real_T u);
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
extern void Adaptive_DataTypeConversion_L0(real_T rtu_u,
  B_DataTypeConversion_L0_Adapt_T *localB);
extern void Adaptiv_DataTypeConversion_dmin(real_T rtu_u,
  B_DataTypeConversion_dmin_Ada_T *localB);

#endif                                 /* RTW_HEADER_Adaptive_private_h_ */
