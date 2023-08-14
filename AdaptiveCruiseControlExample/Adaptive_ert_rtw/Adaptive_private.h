/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Adaptive_private.h
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

#ifndef RTW_HEADER_Adaptive_private_h_
#define RTW_HEADER_Adaptive_private_h_
#include "rtwtypes.h"
#include "Adaptive_types.h"
#include "Adaptive.h"

extern real_T rt_roundd_snf(real_T u);
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
extern void Adaptive_DataTypeConversion_L0(real_T rtu_u, real_T *rty_y);

#endif                                 /* RTW_HEADER_Adaptive_private_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
