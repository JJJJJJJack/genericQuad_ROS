//
// File: genericQuad_data.cpp
//
// Code generated for Simulink model 'genericQuad'.
//
// Model version                  : 1.51
// Simulink Coder version         : 9.6 (R2021b) 14-May-2021
// C/C++ source code generated on : Sat Dec 24 10:32:34 2022
//
// Target selection: ert.tlc
// Embedded hardware selection: Generic->Unspecified (assume 32-bit Generic)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#include "genericQuad.h"
#include "genericQuad_private.h"

// Block parameters (default storage)
P_genericQuad_T genericQuad_P = {
  // Mask Parameter: PIDVelocityx_D
  //  Referenced by: '<S104>/Derivative Gain'

  0.01,

  // Mask Parameter: PIDVelocityy_D
  //  Referenced by: '<S152>/Derivative Gain'

  0.01,

  // Mask Parameter: PIDVelocityz_D
  //  Referenced by: '<S200>/Derivative Gain'

  2.0,

  // Mask Parameter: PIDattitude_D
  //  Referenced by: '<S392>/Derivative Gain'

  0.0,

  // Mask Parameter: PIDangularroll_D
  //  Referenced by: '<S344>/Derivative Gain'

  0.03,

  // Mask Parameter: PIDangularpitch_D
  //  Referenced by: '<S296>/Derivative Gain'

  0.03,

  // Mask Parameter: PIDangulayaw_D
  //  Referenced by: '<S248>/Derivative Gain'

  0.05,

  // Mask Parameter: PIDattitude_I
  //  Referenced by: '<S395>/Integral Gain'

  0.0,

  // Mask Parameter: PIDangulayaw_I
  //  Referenced by: '<S251>/Integral Gain'

  6.0,

  // Mask Parameter: PIDangularpitch_I
  //  Referenced by: '<S299>/Integral Gain'

  0.5,

  // Mask Parameter: PIDangularroll_I
  //  Referenced by: '<S347>/Integral Gain'

  0.5,

  // Mask Parameter: PIDVelocityx_I
  //  Referenced by: '<S107>/Integral Gain'

  0.0,

  // Mask Parameter: PIDVelocityy_I
  //  Referenced by: '<S155>/Integral Gain'

  0.02,

  // Mask Parameter: PIDVelocityz_I
  //  Referenced by: '<S203>/Integral Gain'

  0.2,

  // Mask Parameter: PIDVelocityx_InitialConditionFo
  //  Referenced by: '<S105>/Filter'

  0.0,

  // Mask Parameter: PIDVelocityy_InitialConditionFo
  //  Referenced by: '<S153>/Filter'

  0.0,

  // Mask Parameter: PIDVelocityz_InitialConditionFo
  //  Referenced by: '<S201>/Filter'

  0.0,

  // Mask Parameter: PIDattitude_InitialConditionFor
  //  Referenced by: '<S393>/Filter'

  0.0,

  // Mask Parameter: PIDangularroll_InitialCondition
  //  Referenced by: '<S345>/Filter'

  0.0,

  // Mask Parameter: PIDangularpitch_InitialConditio
  //  Referenced by: '<S297>/Filter'

  0.0,

  // Mask Parameter: PIDangulayaw_InitialConditionFo
  //  Referenced by: '<S249>/Filter'

  0.0,

  // Mask Parameter: PIDVelocityx_InitialCondition_o
  //  Referenced by: '<S110>/Integrator'

  0.0,

  // Mask Parameter: PIDVelocityy_InitialCondition_i
  //  Referenced by: '<S158>/Integrator'

  0.0,

  // Mask Parameter: PIDVelocityz_InitialCondition_g
  //  Referenced by: '<S206>/Integrator'

  0.0,

  // Mask Parameter: PIDattitude_InitialConditionF_g
  //  Referenced by: '<S398>/Integrator'

  0.0,

  // Mask Parameter: PIDangularroll_InitialConditi_d
  //  Referenced by: '<S350>/Integrator'

  0.0,

  // Mask Parameter: PIDangularpitch_InitialCondit_l
  //  Referenced by: '<S302>/Integrator'

  0.0,

  // Mask Parameter: PIDangulayaw_InitialCondition_k
  //  Referenced by: '<S254>/Integrator'

  0.0,

  // Mask Parameter: PIDVelocityx_N
  //  Referenced by: '<S113>/Filter Coefficient'

  100.0,

  // Mask Parameter: PIDVelocityy_N
  //  Referenced by: '<S161>/Filter Coefficient'

  100.0,

  // Mask Parameter: PIDVelocityz_N
  //  Referenced by: '<S209>/Filter Coefficient'

  100.0,

  // Mask Parameter: PIDattitude_N
  //  Referenced by: '<S401>/Filter Coefficient'

  100.0,

  // Mask Parameter: PIDangularroll_N
  //  Referenced by: '<S353>/Filter Coefficient'

  300.0,

  // Mask Parameter: PIDangularpitch_N
  //  Referenced by: '<S305>/Filter Coefficient'

  300.0,

  // Mask Parameter: PIDangulayaw_N
  //  Referenced by: '<S257>/Filter Coefficient'

  300.0,

  // Mask Parameter: PIDVelocityx_P
  //  Referenced by: '<S115>/Proportional Gain'

  3.0,

  // Mask Parameter: PIDVelocityy_P
  //  Referenced by: '<S163>/Proportional Gain'

  3.0,

  // Mask Parameter: PIDVelocityz_P
  //  Referenced by: '<S211>/Proportional Gain'

  3.0,

  // Mask Parameter: PIDattitude_P
  //  Referenced by: '<S403>/Proportional Gain'

  4.0,

  // Mask Parameter: PIDangularroll_P
  //  Referenced by: '<S348>/Proportional Gain'

  2.0,

  // Mask Parameter: PIDangularpitch_P
  //  Referenced by: '<S300>/Proportional Gain'

  2.0,

  // Mask Parameter: PIDangulayaw_P
  //  Referenced by: '<S252>/Proportional Gain'

  5.0,

  // Mask Parameter: uDOFEulerAngles2_Vm_0
  //  Referenced by: '<S7>/ub,vb,wb'

  { 0.0, 0.0, 0.0 },

  // Mask Parameter: DirectionCosineMatrixtoQuaterni
  //  Referenced by:
  //    '<S43>/Constant'
  //    '<S68>/Constant'
  //    '<S69>/Constant'

  1.0,

  // Mask Parameter: uDOFEulerAngles2_eul_0
  //  Referenced by: '<S23>/phi theta psi'

  { 0.0, 0.0, 0.0 },

  // Mask Parameter: uDOFEulerAngles2_inertia
  //  Referenced by: '<S25>/Constant1'

  { 0.05, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.1 },

  // Mask Parameter: uDOFEulerAngles2_mass_0
  //  Referenced by: '<S25>/Constant'

  2.0,

  // Mask Parameter: uDOFEulerAngles2_pm_0
  //  Referenced by: '<S7>/p,q,r '

  { 0.0, 0.0, 0.0 },

  // Mask Parameter: DirectionCosineMatrixtoQuater_h
  //  Referenced by:
  //    '<S76>/Constant'
  //    '<S78>/Constant'

  4.4408920985006262E-16,

  // Mask Parameter: uDOFEulerAngles2_xme_0
  //  Referenced by: '<S7>/xe,ye,ze'

  { 0.0, 0.0, 0.0 },

  // Computed Parameter: Constant_Value
  //  Referenced by: '<S1>/Constant'

  {
    {
      0U,                              // Seq

      {
        0.0,                           // Sec
        0.0                            // Nsec
      },                               // Stamp

      {
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
        0U, 0U }
      ,                                // FrameId

      {
        0U,                            // CurrentLength
        0U                             // ReceivedLength
      }                                // FrameId_SL_Info
    },                                 // Header

    {
      {
        0.0,                           // X
        0.0,                           // Y
        0.0                            // Z
      },                               // Position

      {
        0.0,                           // X
        0.0,                           // Y
        0.0,                           // Z
        0.0                            // W
      }                                // Orientation
    }                                  // Pose
  },

  // Computed Parameter: Out1_Y0
  //  Referenced by: '<S6>/Out1'

  {
    {
      0.0,                             // X
      0.0,                             // Y
      0.0                              // Z
    },                                 // Linear

    {
      0.0,                             // X
      0.0,                             // Y
      0.0                              // Z
    }                                  // Angular
  },

  // Computed Parameter: Constant_Value_f
  //  Referenced by: '<S3>/Constant'

  {
    {
      0.0,                             // X
      0.0,                             // Y
      0.0                              // Z
    },                                 // Linear

    {
      0.0,                             // X
      0.0,                             // Y
      0.0                              // Z
    }                                  // Angular
  },

  // Expression: 1
  //  Referenced by: '<S42>/Constant'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S42>/Gain'

  0.5,

  // Expression: 2
  //  Referenced by: '<S42>/Gain1'

  2.0,

  // Expression: 1
  //  Referenced by: '<S58>/Constant'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S57>/Constant1'

  0.5,

  // Expression: [0 1]
  //  Referenced by: '<S57>/Constant2'

  { 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S46>/Gain1'

  1.0,

  // Expression: 1
  //  Referenced by: '<S46>/Gain3'

  1.0,

  // Expression: 1
  //  Referenced by: '<S46>/Gain4'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S46>/Gain'

  0.5,

  // Expression: 1
  //  Referenced by: '<S63>/Constant'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S62>/Constant1'

  0.5,

  // Expression: [0 1]
  //  Referenced by: '<S62>/Constant2'

  { 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S47>/Gain1'

  1.0,

  // Expression: 1
  //  Referenced by: '<S47>/Gain2'

  1.0,

  // Expression: 1
  //  Referenced by: '<S47>/Gain3'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S47>/Gain'

  0.5,

  // Expression: 1
  //  Referenced by: '<S53>/Constant'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S52>/Constant1'

  0.5,

  // Expression: [0 1]
  //  Referenced by: '<S52>/Constant2'

  { 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S45>/Gain1'

  1.0,

  // Expression: 1
  //  Referenced by: '<S45>/Gain2'

  1.0,

  // Expression: 1
  //  Referenced by: '<S45>/Gain3'

  1.0,

  // Expression: 0.5
  //  Referenced by: '<S45>/Gain'

  0.5,

  // Expression: 0
  //  Referenced by: '<S69>/Constant1'

  0.0,

  // Expression: 0
  //  Referenced by: '<S68>/Constant1'

  0.0,

  // Expression: -eye(3)
  //  Referenced by: '<S70>/Bias1'

  { -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0 },

  // Expression: -1
  //  Referenced by: '<S71>/Bias'

  -1.0,

  // Expression: zeros(3)
  //  Referenced by: '<S25>/Constant2'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pi
  //  Referenced by: '<S23>/phi theta psi'

  3.1415926535897931,

  // Expression: -pi
  //  Referenced by: '<S23>/phi theta psi'

  -3.1415926535897931,

  // Expression: -1
  //  Referenced by: '<S20>/Gain'

  -1.0,

  // Expression: 0
  //  Referenced by: '<S9>/Constant4'

  0.0,

  // Expression: 0
  //  Referenced by: '<S9>/Constant5'

  0.0,

  // Expression: 2
  //  Referenced by: '<S9>/Mass'

  2.0,

  // Expression: 9.8
  //  Referenced by: '<S9>/g'

  9.8,

  // Computed Parameter: TransferFcn2_A
  //  Referenced by: '<S4>/Transfer Fcn2'

  -24.390243902439025,

  // Computed Parameter: TransferFcn2_C
  //  Referenced by: '<S4>/Transfer Fcn2'

  24.390243902439025,

  // Expression: 10
  //  Referenced by: '<S4>/Saturation'

  10.0,

  // Expression: 0
  //  Referenced by: '<S4>/Saturation'

  0.0,

  // Computed Parameter: TransferFcn1_A
  //  Referenced by: '<S4>/Transfer Fcn1'

  -24.390243902439025,

  // Computed Parameter: TransferFcn1_C
  //  Referenced by: '<S4>/Transfer Fcn1'

  24.390243902439025,

  // Expression: 10
  //  Referenced by: '<S4>/Saturation1'

  10.0,

  // Expression: 0
  //  Referenced by: '<S4>/Saturation1'

  0.0,

  // Computed Parameter: TransferFcn3_A
  //  Referenced by: '<S4>/Transfer Fcn3'

  -24.390243902439025,

  // Computed Parameter: TransferFcn3_C
  //  Referenced by: '<S4>/Transfer Fcn3'

  24.390243902439025,

  // Expression: 10
  //  Referenced by: '<S4>/Saturation2'

  10.0,

  // Expression: 0
  //  Referenced by: '<S4>/Saturation2'

  0.0,

  // Computed Parameter: TransferFcn4_A
  //  Referenced by: '<S4>/Transfer Fcn4'

  -24.390243902439025,

  // Computed Parameter: TransferFcn4_C
  //  Referenced by: '<S4>/Transfer Fcn4'

  24.390243902439025,

  // Expression: 10
  //  Referenced by: '<S4>/Saturation3'

  10.0,

  // Expression: 0
  //  Referenced by: '<S4>/Saturation3'

  0.0,

  // Expression: -1
  //  Referenced by: '<S21>/Gain'

  -1.0,

  // Expression: [1 0 0 0]
  //  Referenced by: '<S8>/Merge'

  { 1.0, 0.0, 0.0, 0.0 },

  // Computed Parameter: Filter_gainval
  //  Referenced by: '<S345>/Filter'

  0.002,

  // Computed Parameter: Integrator_gainval
  //  Referenced by: '<S350>/Integrator'

  0.002,

  // Computed Parameter: Integrator_gainval_n
  //  Referenced by: '<S302>/Integrator'

  0.002,

  // Computed Parameter: Filter_gainval_l
  //  Referenced by: '<S297>/Filter'

  0.002,

  // Computed Parameter: Integrator_gainval_c
  //  Referenced by: '<S254>/Integrator'

  0.002,

  // Computed Parameter: Filter_gainval_b
  //  Referenced by: '<S249>/Filter'

  0.002,

  // Expression: 3.13
  //  Referenced by: '<S4>/Hover throttle'

  3.13,

  // Expression: inf
  //  Referenced by: '<S4>/Saturation4'

  0.0,

  // Expression: 0
  //  Referenced by: '<S4>/Saturation4'

  0.0,

  // Expression: 0.5
  //  Referenced by: '<S4>/Komega'

  0.5,

  // Expression: inf
  //  Referenced by: '<S4>/Saturation5'

  0.0,

  // Expression: 0
  //  Referenced by: '<S4>/Saturation5'

  0.0,

  // Expression: 0.5
  //  Referenced by: '<S4>/Komega1'

  0.5,

  // Expression: inf
  //  Referenced by: '<S4>/Saturation6'

  0.0,

  // Expression: 0
  //  Referenced by: '<S4>/Saturation6'

  0.0,

  // Expression: 0.5
  //  Referenced by: '<S4>/Komega2'

  0.5,

  // Expression: inf
  //  Referenced by: '<S4>/Saturation7'

  0.0,

  // Expression: 0
  //  Referenced by: '<S4>/Saturation7'

  0.0,

  // Expression: 0.5
  //  Referenced by: '<S4>/Komega3'

  0.5
};

//
// File trailer for generated code.
//
// [EOF]
//
