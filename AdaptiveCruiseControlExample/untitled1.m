
%{
 % File: Adaptive_Cruise_Control.c
 * Code generated for Simulink model 'Adaptive_Cruise_Control'.
 * Model version                  : 1.2
 * Simulink Coder version         : 9.5 (R2021a) 14-Nov-2020
 * C/C++ source code generated on : Tue Feb 15 18:29:06 2022
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
%}

include "Adaptive_Cruise_Control.h"
include "Adaptive_Cruise_Control_private.h"

% Named constants for Chart: '<S3>/ACC_Logic' */
define Ada_IN_LeadVehicle_Not_Detected ((uint8_T)3U)
define Adaptive_Cr_IN_ACC_STANDBY_MODE ((uint8_T)3U)
define Adaptive_Cru_IN_NO_ACTIVE_CHILD ((uint8_T)0U)
define Adaptive_Cruise_IN_ACC_OFF_MODE ((uint8_T)1U)
define Adaptive_Cruise__IN_ACC_ON_MODE ((uint8_T)2U)
define IN_LeadVehicle_Detected_Follow ((uint8_T)1U)
define IN_LeadVehicle_Detected_Resume ((uint8_T)2U)
define IN_LeadVehicle_Not_Detected_Res ((uint8_T)4U)
define IN_LeadVehicle_Speed_equal_Set_ ((uint8_T)5U)
define IN_LeadVehicle_Speed_lessthan_S ((uint8_T)6U)

% Block states (default storage) */
DW_Adaptive_Cruise_Control_T Adaptive_Cruise_Control_DW;

% External inputs (root inport signals with default storage) */
ExtU_Adaptive_Cruise_Control_T Adaptive_Cruise_Control_U;

% External outputs (root outports fed by signals with default storage) */
ExtY_Adaptive_Cruise_Control_T Adaptive_Cruise_Control_Y;

% Real-time model */
static RT_MODEL_Adaptive_Cruise_Cont_T Adaptive_Cruise_Control_M_;
RT_MODEL_Adaptive_Cruise_Cont_T *const Adaptive_Cruise_Control_M =
  &Adaptive_Cruise_Control_M_;

% Model step function */
void Adaptive_Cruise_Control_step(void)
{
  real_T rtb_Add;
  boolean_T rtb_LogicalOperator;

  %{ Logic: '<S4>/Logical Operator' incorporates:
   %  Inport: '<Root>/CameraInput_LeadVehicle'
   %  Inport: '<Root>/RadarInput_LeadVehicle'
   %}
  rtb_LogicalOperator = (Adaptive_Cruise_Control_U.CameraInput_LeadVehicle &&
    Adaptive_Cruise_Control_U.RadarInput_LeadVehicle);

  % Sum: '<S2>/Add' incorporates:
   %*  Inport: '<Root>/CameraInput_DriveVehicle'
   %*  Inport: '<Root>/RadarInput_DriveVehicle'
   %*  UnitDelay: '<S1>/Unit Delay'
   %
  rtb_Add = (Adaptive_Cruise_Control_U.CameraInput_DriveVehicle +
             Adaptive_Cruise_Control_U.RadarInput_DriveVehicle) +
    Adaptive_Cruise_Control_Y.Acceleration_mode;

  %/* Chart: '<S3>/ACC_Logic' incorporates:
   %*  Inport: '<Root>/CruiseSwitch'
   %*  Inport: '<Root>/RadarInput_DriveVehicle'
   %*  Inport: '<Root>/SetSwitch'
   %*  Inport: '<Root>/Set_Gap'
   %*  Inport: '<Root>/Set_Speed'
   %*  Inport: '<Root>/Time_Gap'
   %*  UnitDelay: '<S1>/Unit Delay'
   %*/
  if (Adaptive_Cruise_Control_DW.is_active_c3_Adaptive_Cruise_Co == 0U) {
    Adaptive_Cruise_Control_DW.is_active_c3_Adaptive_Cruise_Co = 1U;
    Adaptive_Cruise_Control_DW.is_c3_Adaptive_Cruise_Control =
      Adaptive_Cruise_IN_ACC_OFF_MODE;
    Adaptive_Cruise_Control_Y.Acceleration_mode = 0.0;
  } else {
    switch (Adaptive_Cruise_Control_DW.is_c3_Adaptive_Cruise_Control) {
     case Adaptive_Cruise_IN_ACC_OFF_MODE:
      Adaptive_Cruise_Control_Y.Acceleration_mode = 0.0;
      if (Adaptive_Cruise_Control_U.CruiseSwitch == 1.0) {
        Adaptive_Cruise_Control_DW.is_c3_Adaptive_Cruise_Control =
          Adaptive_Cr_IN_ACC_STANDBY_MODE;
        Adaptive_Cruise_Control_Y.Acceleration_mode = 1.0;
      }
      break;

     case Adaptive_Cruise__IN_ACC_ON_MODE:
      if (Adaptive_Cruise_Control_U.SetSwitch == 0.0) {
        Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
          Adaptive_Cru_IN_NO_ACTIVE_CHILD;
        Adaptive_Cruise_Control_DW.is_c3_Adaptive_Cruise_Control =
          Adaptive_Cr_IN_ACC_STANDBY_MODE;
        Adaptive_Cruise_Control_Y.Acceleration_mode = 1.0;
      } else {
        switch (Adaptive_Cruise_Control_DW.is_ACC_ON_MODE) {
         case IN_LeadVehicle_Detected_Follow:
          Adaptive_Cruise_Control_Y.Acceleration_mode = 2.0;
          if (Adaptive_Cruise_Control_U.RadarInput_DriveVehicle == 0.0) {
            Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
              Ada_IN_LeadVehicle_Not_Detected;
            Adaptive_Cruise_Control_Y.Acceleration_mode = 1.0;
          }
          break;

         case IN_LeadVehicle_Detected_Resume:
          Adaptive_Cruise_Control_Y.Acceleration_mode = 3.0;
          if ((rtb_Add == Adaptive_Cruise_Control_U.Set_Speed) && ((real_T)
               rtb_LogicalOperator >= Adaptive_Cruise_Control_U.Set_Speed) &&
              (Adaptive_Cruise_Control_U.Time_Gap >=
               Adaptive_Cruise_Control_U.Set_Gap)) {
            Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
              IN_LeadVehicle_Detected_Follow;
            Adaptive_Cruise_Control_Y.Acceleration_mode = 2.0;
          } else if (Adaptive_Cruise_Control_U.RadarInput_DriveVehicle == 0.0) {
            Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
              IN_LeadVehicle_Not_Detected_Res;
            Adaptive_Cruise_Control_Y.Acceleration_mode = 1.0;
          } else if ((rtb_Add < Adaptive_Cruise_Control_U.Set_Speed) && ((real_T)
                      rtb_LogicalOperator > rtb_Add) &&
                     (Adaptive_Cruise_Control_U.Time_Gap >=
                      Adaptive_Cruise_Control_U.Set_Gap)) {
            Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
              IN_LeadVehicle_Speed_equal_Set_;
            Adaptive_Cruise_Control_Y.Acceleration_mode = 5.0;
          }
          break;

         case Ada_IN_LeadVehicle_Not_Detected:
          Adaptive_Cruise_Control_Y.Acceleration_mode = 1.0;
          if ((Adaptive_Cruise_Control_U.RadarInput_DriveVehicle == 1.0) &&
              (rtb_Add == Adaptive_Cruise_Control_U.Set_Speed) && ((real_T)
               rtb_LogicalOperator >= Adaptive_Cruise_Control_U.Set_Speed) &&
              (Adaptive_Cruise_Control_U.Time_Gap >=
               Adaptive_Cruise_Control_U.Set_Gap)) {
            Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
              IN_LeadVehicle_Detected_Follow;
            Adaptive_Cruise_Control_Y.Acceleration_mode = 2.0;
          } else if (((Adaptive_Cruise_Control_U.RadarInput_DriveVehicle == 1.0)
                      && ((real_T)rtb_LogicalOperator <
                          Adaptive_Cruise_Control_U.Set_Speed)) ||
                     (Adaptive_Cruise_Control_U.Time_Gap <
                      Adaptive_Cruise_Control_U.Set_Gap)) {
            Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
              IN_LeadVehicle_Speed_lessthan_S;
            Adaptive_Cruise_Control_Y.Acceleration_mode = 4.0;
          }
          break;

         case IN_LeadVehicle_Not_Detected_Res:
          Adaptive_Cruise_Control_Y.Acceleration_mode = 1.0;
          break;

         case IN_LeadVehicle_Speed_equal_Set_:
          Adaptive_Cruise_Control_Y.Acceleration_mode = 5.0;
          if ((Adaptive_Cruise_Control_U.RadarInput_DriveVehicle == 0.0) ||
              (rtb_Add <= Adaptive_Cruise_Control_U.Set_Speed)) {
            Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
              IN_LeadVehicle_Not_Detected_Res;
            Adaptive_Cruise_Control_Y.Acceleration_mode = 1.0;
          } else if (((rtb_Add < Adaptive_Cruise_Control_U.Set_Speed) &&
                      ((real_T)rtb_LogicalOperator > rtb_Add)) ||
                     (Adaptive_Cruise_Control_U.Time_Gap >=
                      Adaptive_Cruise_Control_U.Set_Gap)) {
            Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
              IN_LeadVehicle_Detected_Resume;
            Adaptive_Cruise_Control_Y.Acceleration_mode = 3.0;
          } else if ((((real_T)rtb_LogicalOperator <
                       Adaptive_Cruise_Control_U.Set_Speed) && ((real_T)
                       rtb_LogicalOperator < rtb_Add)) || (0.75 *
                      Adaptive_Cruise_Control_U.Set_Gap ==
                      Adaptive_Cruise_Control_U.Time_Gap)) {
            Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
              IN_LeadVehicle_Speed_lessthan_S;
            Adaptive_Cruise_Control_Y.Acceleration_mode = 4.0;
          }
          break;

         default:
          /* case IN_LeadVehicle_Speed_lessthan_Set_Speed: */
          Adaptive_Cruise_Control_Y.Acceleration_mode = 4.0;
          if ((Adaptive_Cruise_Control_U.RadarInput_DriveVehicle == 0.0) &&
              (rtb_Add == Adaptive_Cruise_Control_U.Set_Speed)) {
            Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
              Ada_IN_LeadVehicle_Not_Detected;
            Adaptive_Cruise_Control_Y.Acceleration_mode = 1.0;
          } else if (((real_T)rtb_LogicalOperator * 1.25 >= rtb_Add) && ((real_T)
                      rtb_LogicalOperator * 0.75 <= rtb_Add) && (rtb_Add <
                      Adaptive_Cruise_Control_U.Set_Speed) &&
                     (Adaptive_Cruise_Control_U.Time_Gap <= 1.25 *
                      Adaptive_Cruise_Control_U.Set_Gap) &&
                     (Adaptive_Cruise_Control_U.Time_Gap >= 0.75 *
                      Adaptive_Cruise_Control_U.Set_Gap)) {
            Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
              IN_LeadVehicle_Speed_equal_Set_;
            Adaptive_Cruise_Control_Y.Acceleration_mode = 5.0;
          }
          break;
        }
      }
      break;

     default:
      /* case IN_ACC_STANDBY_MODE: */
      Adaptive_Cruise_Control_Y.Acceleration_mode = 1.0;
      if (Adaptive_Cruise_Control_U.CruiseSwitch == 0.0) {
        Adaptive_Cruise_Control_DW.is_c3_Adaptive_Cruise_Control =
          Adaptive_Cruise_IN_ACC_OFF_MODE;
        Adaptive_Cruise_Control_Y.Acceleration_mode = 0.0;
      } else if (Adaptive_Cruise_Control_U.SetSwitch == 0.0) {
        Adaptive_Cruise_Control_DW.is_c3_Adaptive_Cruise_Control =
          Adaptive_Cruise__IN_ACC_ON_MODE;
        Adaptive_Cruise_Control_DW.is_ACC_ON_MODE =
          IN_LeadVehicle_Detected_Follow;
        Adaptive_Cruise_Control_Y.Acceleration_mode = 2.0;
      }
      break;
    }
  }

  %/* End of Chart: '<S3>/ACC_Logic' */
}

%/* Model initialize function */
void Adaptive_Cruise_Control_initialize(void)
{
  %/* (no initialization code required) */
}

%/* Model terminate function */
void Adaptive_Cruise_Control_terminate(void)
{
  %/* (no terminate code required) */
}

%/*
% * File trailer for generated code.
% *
 %* [EOF]
 %*/