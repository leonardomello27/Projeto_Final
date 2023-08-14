/*
 * Adaptive_sf_types.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "Adaptive_sf".
 *
 * Model version              : 11.0
 * Simulink Coder version : 9.9 (R2023a) 19-Nov-2022
 * C source code generated on : Thu Aug  3 20:54:02 2023
 *
 * Target selection: rtwsfcn.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_Adaptive_sf_types_h_
#define RTW_HEADER_Adaptive_sf_types_h_
#include "rtwtypes.h"
#ifndef DEFINED_TYPEDEF_FOR_struct_WTmPWsEMvOzNnnAVv5fQNC_
#define DEFINED_TYPEDEF_FOR_struct_WTmPWsEMvOzNnnAVv5fQNC_

typedef struct {
  int32_T MaxIterations;
  uint8_T sl_padding0[4];
  real_T ConstraintTolerance;
  boolean_T UseWarmStart;
  uint8_T sl_padding1[7];
} struct_WTmPWsEMvOzNnnAVv5fQNC;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_WHjMt45Sk148iktWsfFxl_
#define DEFINED_TYPEDEF_FOR_struct_WHjMt45Sk148iktWsfFxl_

typedef struct {
  int32_T MaxIterations;
  uint8_T sl_padding0[4];
  real_T ConstraintTolerance;
  real_T OptimalityTolerance;
  real_T ComplementarityTolerance;
  real_T StepTolerance;
} struct_WHjMt45Sk148iktWsfFxl;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_lnQ9KXdSZFplhcBp5LBCc_
#define DEFINED_TYPEDEF_FOR_struct_lnQ9KXdSZFplhcBp5LBCc_

typedef struct {
  int32_T MaxIterations;
  uint8_T sl_padding0[4];
  real_T ConstraintTolerance;
  real_T DiscreteConstraintTolerance;
  boolean_T RoundingAtRootNode;
  uint8_T sl_padding1[3];
  int32_T MaxPendingNodes;
} struct_lnQ9KXdSZFplhcBp5LBCc;

#endif
#endif                                 /* RTW_HEADER_Adaptive_sf_types_h_ */
