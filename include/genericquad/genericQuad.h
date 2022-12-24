//
// File: genericQuad.h
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
#ifndef RTW_HEADER_genericQuad_h_
#define RTW_HEADER_genericQuad_h_
#include <string.h>
#include <math.h>
#include <stddef.h>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "slros_initialize.h"
#include "genericQuad_types.h"
#include "rt_assert.h"
#include "rtGetInf.h"
#include "rt_nonfinite.h"

// Macros for accessing real-time model data structure
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
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
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

// Block signals (default storage)
struct B_genericQuad_T {
  SL_Bus_genericQuad_geometry_msgs_PoseStamped BusAssignment;// '<Root>/Bus Assignment' 
  real_T VectorConcatenate[18];        // '<S25>/Vector Concatenate'
  real_T Selector1[9];                 // '<S24>/Selector1'
  real_T Selector2[9];                 // '<S24>/Selector2'
  real_T VectorConcatenate_m[9];
  real_T Product_tmp[9];
  real_T Product_tmp_c[9];
  real_T A[9];
  SL_Bus_genericQuad_geometry_msgs_Twist In1;// '<S6>/In1'
  real_T Selector[9];                  // '<S24>/Selector'
  SL_Bus_genericQuad_geometry_msgs_Twist b_varargout_2;
  real_T Product[3];                   // '<S30>/Product'
  real_T Reshape[3];                   // '<S9>/Reshape'
  real_T Sum[3];                       // '<S7>/Sum'
  real_T Product2[3];                  // '<S24>/Product2'
  real_T Sum_d[3];                     // '<Root>/Sum'
  real_T sincos_o2[3];                 // '<S32>/sincos'
  real_T sincos_o1[3];                 // '<S32>/sincos'
  real_T Sum7[3];
  real_T rtb_Saturation2_k[3];
  real_T rtb_Sum_d_c[3];
  boolean_T Compare[9];                // '<S76>/Compare'
  real_T Z;
  real_T FilterCoefficient;            // '<S113>/Filter Coefficient'
  real_T FilterCoefficient_m;          // '<S161>/Filter Coefficient'
  real_T FilterCoefficient_j;          // '<S209>/Filter Coefficient'
  real_T FilterCoefficient_p[3];       // '<S401>/Filter Coefficient'
  real_T IntegralGain[3];              // '<S395>/Integral Gain'
  real_T TmpSignalConversionAtphithetaps[3];// '<S23>/phidot thetadot psidot'
  real_T Merge[4];                     // '<S8>/Merge'
  real_T ProportionalGain;             // '<S348>/Proportional Gain'
  real_T ProportionalGain_a;           // '<S300>/Proportional Gain'
  real_T ProportionalGain_p;           // '<S252>/Proportional Gain'
  real_T Komega;                       // '<S4>/Komega'
  real_T Komega1;                      // '<S4>/Komega1'
  real_T Komega2;                      // '<S4>/Komega2'
  real_T Komega3;                      // '<S4>/Komega3'
  real_T IntegralGain_p;               // '<S107>/Integral Gain'
  real_T IntegralGain_i;               // '<S155>/Integral Gain'
  real_T IntegralGain_g;               // '<S203>/Integral Gain'
  real_T kxj;                          // '<S38>/k x j'
  real_T Saturation;                   // '<S4>/Saturation'
  real_T Saturation1;                  // '<S4>/Saturation1'
  real_T Saturation2;                  // '<S4>/Saturation2'
  real_T Saturation7;                  // '<S4>/Saturation7'
  real_T thrustsetpoint;               // '<S4>/Sum3'
  real_T rtb_Sum2_idx_2;
  real_T rtb_Sum2_idx_0;
  real_T rtb_Sum2_idx_1;
  real_T Selector_tmp;
  real_T rtb_Sum2_idx_1_tmp;
  real_T VectorConcatenate_tmp;
  real_T IntegralGain_g_b;             // '<S347>/Integral Gain'
  real_T IntegralGain_c;               // '<S299>/Integral Gain'
  real_T IntegralGain_pb;              // '<S251>/Integral Gain'
  real_T FilterCoefficient_n;          // '<S257>/Filter Coefficient'
  real_T FilterCoefficient_e;          // '<S305>/Filter Coefficient'
};

// Block states (default storage) for system '<Root>'
struct DW_genericQuad_T {
  ros_slros_internal_block_GetP_T obj; // '<Root>/Get Parameter2'
  ros_slros_internal_block_GetP_T obj_a;// '<Root>/Get Parameter1'
  ros_slros_internal_block_GetP_T obj_p;// '<Root>/Get Parameter'
  ros_slroscpp_internal_block_P_T obj_m;// '<S2>/SinkBlock'
  ros_slroscpp_internal_block_S_T obj_c;// '<S3>/SourceBlock'
  real_T Filter_DSTATE;                // '<S345>/Filter'
  real_T Integrator_DSTATE;            // '<S350>/Integrator'
  real_T Integrator_DSTATE_o;          // '<S302>/Integrator'
  real_T Filter_DSTATE_g;              // '<S297>/Filter'
  real_T Integrator_DSTATE_i;          // '<S254>/Integrator'
  real_T Filter_DSTATE_f;              // '<S249>/Filter'
  real_T Product2_DWORK4[9];           // '<S24>/Product2'
  int8_T If_ActiveSubsystem;           // '<S8>/If'
  int8_T If1_ActiveSubsystem;          // '<S43>/If1'
  int8_T FindMaximumDiagonalValue_Active;// '<S41>/Find Maximum Diagonal Value'
};

// Continuous states (default storage)
struct X_genericQuad_T {
  real_T phithetapsi_CSTATE[3];        // '<S23>/phi theta psi'
  real_T ubvbwb_CSTATE[3];             // '<S7>/ub,vb,wb'
  real_T TransferFcn2_CSTATE;          // '<S4>/Transfer Fcn2'
  real_T TransferFcn1_CSTATE;          // '<S4>/Transfer Fcn1'
  real_T TransferFcn3_CSTATE;          // '<S4>/Transfer Fcn3'
  real_T TransferFcn4_CSTATE;          // '<S4>/Transfer Fcn4'
  real_T pqr_CSTATE[3];                // '<S7>/p,q,r '
  real_T Integrator_CSTATE;            // '<S110>/Integrator'
  real_T Filter_CSTATE;                // '<S105>/Filter'
  real_T Integrator_CSTATE_o;          // '<S158>/Integrator'
  real_T Filter_CSTATE_k;              // '<S153>/Filter'
  real_T Integrator_CSTATE_h;          // '<S206>/Integrator'
  real_T Filter_CSTATE_h;              // '<S201>/Filter'
  real_T Integrator_CSTATE_m[3];       // '<S398>/Integrator'
  real_T Filter_CSTATE_e[3];           // '<S393>/Filter'
  real_T xeyeze_CSTATE[3];             // '<S7>/xe,ye,ze'
};

// Periodic continuous state vector (global)
typedef int_T PeriodicIndX_genericQuad_T[3];
typedef real_T PeriodicRngX_genericQuad_T[6];

// State derivatives (default storage)
struct XDot_genericQuad_T {
  real_T phithetapsi_CSTATE[3];        // '<S23>/phi theta psi'
  real_T ubvbwb_CSTATE[3];             // '<S7>/ub,vb,wb'
  real_T TransferFcn2_CSTATE;          // '<S4>/Transfer Fcn2'
  real_T TransferFcn1_CSTATE;          // '<S4>/Transfer Fcn1'
  real_T TransferFcn3_CSTATE;          // '<S4>/Transfer Fcn3'
  real_T TransferFcn4_CSTATE;          // '<S4>/Transfer Fcn4'
  real_T pqr_CSTATE[3];                // '<S7>/p,q,r '
  real_T Integrator_CSTATE;            // '<S110>/Integrator'
  real_T Filter_CSTATE;                // '<S105>/Filter'
  real_T Integrator_CSTATE_o;          // '<S158>/Integrator'
  real_T Filter_CSTATE_k;              // '<S153>/Filter'
  real_T Integrator_CSTATE_h;          // '<S206>/Integrator'
  real_T Filter_CSTATE_h;              // '<S201>/Filter'
  real_T Integrator_CSTATE_m[3];       // '<S398>/Integrator'
  real_T Filter_CSTATE_e[3];           // '<S393>/Filter'
  real_T xeyeze_CSTATE[3];             // '<S7>/xe,ye,ze'
};

// State disabled
struct XDis_genericQuad_T {
  boolean_T phithetapsi_CSTATE[3];     // '<S23>/phi theta psi'
  boolean_T ubvbwb_CSTATE[3];          // '<S7>/ub,vb,wb'
  boolean_T TransferFcn2_CSTATE;       // '<S4>/Transfer Fcn2'
  boolean_T TransferFcn1_CSTATE;       // '<S4>/Transfer Fcn1'
  boolean_T TransferFcn3_CSTATE;       // '<S4>/Transfer Fcn3'
  boolean_T TransferFcn4_CSTATE;       // '<S4>/Transfer Fcn4'
  boolean_T pqr_CSTATE[3];             // '<S7>/p,q,r '
  boolean_T Integrator_CSTATE;         // '<S110>/Integrator'
  boolean_T Filter_CSTATE;             // '<S105>/Filter'
  boolean_T Integrator_CSTATE_o;       // '<S158>/Integrator'
  boolean_T Filter_CSTATE_k;           // '<S153>/Filter'
  boolean_T Integrator_CSTATE_h;       // '<S206>/Integrator'
  boolean_T Filter_CSTATE_h;           // '<S201>/Filter'
  boolean_T Integrator_CSTATE_m[3];    // '<S398>/Integrator'
  boolean_T Filter_CSTATE_e[3];        // '<S393>/Filter'
  boolean_T xeyeze_CSTATE[3];          // '<S7>/xe,ye,ze'
};

#ifndef ODE3_INTG
#define ODE3_INTG

// ODE3 Integration Data
struct ODE3_IntgData {
  real_T *y;                           // output
  real_T *f[3];                        // derivatives
};

#endif

// Parameters (default storage)
struct P_genericQuad_T_ {
  real_T PIDVelocityx_D;               // Mask Parameter: PIDVelocityx_D
                                          //  Referenced by: '<S104>/Derivative Gain'

  real_T PIDVelocityy_D;               // Mask Parameter: PIDVelocityy_D
                                          //  Referenced by: '<S152>/Derivative Gain'

  real_T PIDVelocityz_D;               // Mask Parameter: PIDVelocityz_D
                                          //  Referenced by: '<S200>/Derivative Gain'

  real_T PIDattitude_D;                // Mask Parameter: PIDattitude_D
                                          //  Referenced by: '<S392>/Derivative Gain'

  real_T PIDangularroll_D;             // Mask Parameter: PIDangularroll_D
                                          //  Referenced by: '<S344>/Derivative Gain'

  real_T PIDangularpitch_D;            // Mask Parameter: PIDangularpitch_D
                                          //  Referenced by: '<S296>/Derivative Gain'

  real_T PIDangulayaw_D;               // Mask Parameter: PIDangulayaw_D
                                          //  Referenced by: '<S248>/Derivative Gain'

  real_T PIDattitude_I;                // Mask Parameter: PIDattitude_I
                                          //  Referenced by: '<S395>/Integral Gain'

  real_T PIDangulayaw_I;               // Mask Parameter: PIDangulayaw_I
                                          //  Referenced by: '<S251>/Integral Gain'

  real_T PIDangularpitch_I;            // Mask Parameter: PIDangularpitch_I
                                          //  Referenced by: '<S299>/Integral Gain'

  real_T PIDangularroll_I;             // Mask Parameter: PIDangularroll_I
                                          //  Referenced by: '<S347>/Integral Gain'

  real_T PIDVelocityx_I;               // Mask Parameter: PIDVelocityx_I
                                          //  Referenced by: '<S107>/Integral Gain'

  real_T PIDVelocityy_I;               // Mask Parameter: PIDVelocityy_I
                                          //  Referenced by: '<S155>/Integral Gain'

  real_T PIDVelocityz_I;               // Mask Parameter: PIDVelocityz_I
                                          //  Referenced by: '<S203>/Integral Gain'

  real_T PIDVelocityx_InitialConditionFo;
                              // Mask Parameter: PIDVelocityx_InitialConditionFo
                                 //  Referenced by: '<S105>/Filter'

  real_T PIDVelocityy_InitialConditionFo;
                              // Mask Parameter: PIDVelocityy_InitialConditionFo
                                 //  Referenced by: '<S153>/Filter'

  real_T PIDVelocityz_InitialConditionFo;
                              // Mask Parameter: PIDVelocityz_InitialConditionFo
                                 //  Referenced by: '<S201>/Filter'

  real_T PIDattitude_InitialConditionFor;
                              // Mask Parameter: PIDattitude_InitialConditionFor
                                 //  Referenced by: '<S393>/Filter'

  real_T PIDangularroll_InitialCondition;
                              // Mask Parameter: PIDangularroll_InitialCondition
                                 //  Referenced by: '<S345>/Filter'

  real_T PIDangularpitch_InitialConditio;
                              // Mask Parameter: PIDangularpitch_InitialConditio
                                 //  Referenced by: '<S297>/Filter'

  real_T PIDangulayaw_InitialConditionFo;
                              // Mask Parameter: PIDangulayaw_InitialConditionFo
                                 //  Referenced by: '<S249>/Filter'

  real_T PIDVelocityx_InitialCondition_o;
                              // Mask Parameter: PIDVelocityx_InitialCondition_o
                                 //  Referenced by: '<S110>/Integrator'

  real_T PIDVelocityy_InitialCondition_i;
                              // Mask Parameter: PIDVelocityy_InitialCondition_i
                                 //  Referenced by: '<S158>/Integrator'

  real_T PIDVelocityz_InitialCondition_g;
                              // Mask Parameter: PIDVelocityz_InitialCondition_g
                                 //  Referenced by: '<S206>/Integrator'

  real_T PIDattitude_InitialConditionF_g;
                              // Mask Parameter: PIDattitude_InitialConditionF_g
                                 //  Referenced by: '<S398>/Integrator'

  real_T PIDangularroll_InitialConditi_d;
                              // Mask Parameter: PIDangularroll_InitialConditi_d
                                 //  Referenced by: '<S350>/Integrator'

  real_T PIDangularpitch_InitialCondit_l;
                              // Mask Parameter: PIDangularpitch_InitialCondit_l
                                 //  Referenced by: '<S302>/Integrator'

  real_T PIDangulayaw_InitialCondition_k;
                              // Mask Parameter: PIDangulayaw_InitialCondition_k
                                 //  Referenced by: '<S254>/Integrator'

  real_T PIDVelocityx_N;               // Mask Parameter: PIDVelocityx_N
                                          //  Referenced by: '<S113>/Filter Coefficient'

  real_T PIDVelocityy_N;               // Mask Parameter: PIDVelocityy_N
                                          //  Referenced by: '<S161>/Filter Coefficient'

  real_T PIDVelocityz_N;               // Mask Parameter: PIDVelocityz_N
                                          //  Referenced by: '<S209>/Filter Coefficient'

  real_T PIDattitude_N;                // Mask Parameter: PIDattitude_N
                                          //  Referenced by: '<S401>/Filter Coefficient'

  real_T PIDangularroll_N;             // Mask Parameter: PIDangularroll_N
                                          //  Referenced by: '<S353>/Filter Coefficient'

  real_T PIDangularpitch_N;            // Mask Parameter: PIDangularpitch_N
                                          //  Referenced by: '<S305>/Filter Coefficient'

  real_T PIDangulayaw_N;               // Mask Parameter: PIDangulayaw_N
                                          //  Referenced by: '<S257>/Filter Coefficient'

  real_T PIDVelocityx_P;               // Mask Parameter: PIDVelocityx_P
                                          //  Referenced by: '<S115>/Proportional Gain'

  real_T PIDVelocityy_P;               // Mask Parameter: PIDVelocityy_P
                                          //  Referenced by: '<S163>/Proportional Gain'

  real_T PIDVelocityz_P;               // Mask Parameter: PIDVelocityz_P
                                          //  Referenced by: '<S211>/Proportional Gain'

  real_T PIDattitude_P;                // Mask Parameter: PIDattitude_P
                                          //  Referenced by: '<S403>/Proportional Gain'

  real_T PIDangularroll_P;             // Mask Parameter: PIDangularroll_P
                                          //  Referenced by: '<S348>/Proportional Gain'

  real_T PIDangularpitch_P;            // Mask Parameter: PIDangularpitch_P
                                          //  Referenced by: '<S300>/Proportional Gain'

  real_T PIDangulayaw_P;               // Mask Parameter: PIDangulayaw_P
                                          //  Referenced by: '<S252>/Proportional Gain'

  real_T uDOFEulerAngles2_Vm_0[3];     // Mask Parameter: uDOFEulerAngles2_Vm_0
                                          //  Referenced by: '<S7>/ub,vb,wb'

  real_T DirectionCosineMatrixtoQuaterni;
                              // Mask Parameter: DirectionCosineMatrixtoQuaterni
                                 //  Referenced by:
                                 //    '<S43>/Constant'
                                 //    '<S68>/Constant'
                                 //    '<S69>/Constant'

  real_T uDOFEulerAngles2_eul_0[3];    // Mask Parameter: uDOFEulerAngles2_eul_0
                                          //  Referenced by: '<S23>/phi theta psi'

  real_T uDOFEulerAngles2_inertia[9];// Mask Parameter: uDOFEulerAngles2_inertia
                                        //  Referenced by: '<S25>/Constant1'

  real_T uDOFEulerAngles2_mass_0;     // Mask Parameter: uDOFEulerAngles2_mass_0
                                         //  Referenced by: '<S25>/Constant'

  real_T uDOFEulerAngles2_pm_0[3];     // Mask Parameter: uDOFEulerAngles2_pm_0
                                          //  Referenced by: '<S7>/p,q,r '

  real_T DirectionCosineMatrixtoQuater_h;
                              // Mask Parameter: DirectionCosineMatrixtoQuater_h
                                 //  Referenced by:
                                 //    '<S76>/Constant'
                                 //    '<S78>/Constant'

  real_T uDOFEulerAngles2_xme_0[3];    // Mask Parameter: uDOFEulerAngles2_xme_0
                                          //  Referenced by: '<S7>/xe,ye,ze'

  SL_Bus_genericQuad_geometry_msgs_PoseStamped Constant_Value;// Computed Parameter: Constant_Value
                                                                 //  Referenced by: '<S1>/Constant'

  SL_Bus_genericQuad_geometry_msgs_Twist Out1_Y0;// Computed Parameter: Out1_Y0
                                                    //  Referenced by: '<S6>/Out1'

  SL_Bus_genericQuad_geometry_msgs_Twist Constant_Value_f;// Computed Parameter: Constant_Value_f
                                                             //  Referenced by: '<S3>/Constant'

  real_T Constant_Value_o;             // Expression: 1
                                          //  Referenced by: '<S42>/Constant'

  real_T Gain_Gain;                    // Expression: 0.5
                                          //  Referenced by: '<S42>/Gain'

  real_T Gain1_Gain;                   // Expression: 2
                                          //  Referenced by: '<S42>/Gain1'

  real_T Constant_Value_n;             // Expression: 1
                                          //  Referenced by: '<S58>/Constant'

  real_T Constant1_Value;              // Expression: 0.5
                                          //  Referenced by: '<S57>/Constant1'

  real_T Constant2_Value[2];           // Expression: [0 1]
                                          //  Referenced by: '<S57>/Constant2'

  real_T Gain1_Gain_k;                 // Expression: 1
                                          //  Referenced by: '<S46>/Gain1'

  real_T Gain3_Gain;                   // Expression: 1
                                          //  Referenced by: '<S46>/Gain3'

  real_T Gain4_Gain;                   // Expression: 1
                                          //  Referenced by: '<S46>/Gain4'

  real_T Gain_Gain_g;                  // Expression: 0.5
                                          //  Referenced by: '<S46>/Gain'

  real_T Constant_Value_p;             // Expression: 1
                                          //  Referenced by: '<S63>/Constant'

  real_T Constant1_Value_h;            // Expression: 0.5
                                          //  Referenced by: '<S62>/Constant1'

  real_T Constant2_Value_k[2];         // Expression: [0 1]
                                          //  Referenced by: '<S62>/Constant2'

  real_T Gain1_Gain_l;                 // Expression: 1
                                          //  Referenced by: '<S47>/Gain1'

  real_T Gain2_Gain;                   // Expression: 1
                                          //  Referenced by: '<S47>/Gain2'

  real_T Gain3_Gain_a;                 // Expression: 1
                                          //  Referenced by: '<S47>/Gain3'

  real_T Gain_Gain_d;                  // Expression: 0.5
                                          //  Referenced by: '<S47>/Gain'

  real_T Constant_Value_fv;            // Expression: 1
                                          //  Referenced by: '<S53>/Constant'

  real_T Constant1_Value_b;            // Expression: 0.5
                                          //  Referenced by: '<S52>/Constant1'

  real_T Constant2_Value_o[2];         // Expression: [0 1]
                                          //  Referenced by: '<S52>/Constant2'

  real_T Gain1_Gain_b;                 // Expression: 1
                                          //  Referenced by: '<S45>/Gain1'

  real_T Gain2_Gain_d;                 // Expression: 1
                                          //  Referenced by: '<S45>/Gain2'

  real_T Gain3_Gain_az;                // Expression: 1
                                          //  Referenced by: '<S45>/Gain3'

  real_T Gain_Gain_i;                  // Expression: 0.5
                                          //  Referenced by: '<S45>/Gain'

  real_T Constant1_Value_d;            // Expression: 0
                                          //  Referenced by: '<S69>/Constant1'

  real_T Constant1_Value_j;            // Expression: 0
                                          //  Referenced by: '<S68>/Constant1'

  real_T Bias1_Bias[9];                // Expression: -eye(3)
                                          //  Referenced by: '<S70>/Bias1'

  real_T Bias_Bias;                    // Expression: -1
                                          //  Referenced by: '<S71>/Bias'

  real_T Constant2_Value_n[9];         // Expression: zeros(3)
                                          //  Referenced by: '<S25>/Constant2'

  real_T phithetapsi_WrappedStateUpperVa;// Expression: pi
                                            //  Referenced by: '<S23>/phi theta psi'

  real_T phithetapsi_WrappedStateLowerVa;// Expression: -pi
                                            //  Referenced by: '<S23>/phi theta psi'

  real_T Gain_Gain_ip;                 // Expression: -1
                                          //  Referenced by: '<S20>/Gain'

  real_T Constant4_Value;              // Expression: 0
                                          //  Referenced by: '<S9>/Constant4'

  real_T Constant5_Value;              // Expression: 0
                                          //  Referenced by: '<S9>/Constant5'

  real_T Mass_Value;                   // Expression: 2
                                          //  Referenced by: '<S9>/Mass'

  real_T g_Gain;                       // Expression: 9.8
                                          //  Referenced by: '<S9>/g'

  real_T TransferFcn2_A;               // Computed Parameter: TransferFcn2_A
                                          //  Referenced by: '<S4>/Transfer Fcn2'

  real_T TransferFcn2_C;               // Computed Parameter: TransferFcn2_C
                                          //  Referenced by: '<S4>/Transfer Fcn2'

  real_T Saturation_UpperSat;          // Expression: 10
                                          //  Referenced by: '<S4>/Saturation'

  real_T Saturation_LowerSat;          // Expression: 0
                                          //  Referenced by: '<S4>/Saturation'

  real_T TransferFcn1_A;               // Computed Parameter: TransferFcn1_A
                                          //  Referenced by: '<S4>/Transfer Fcn1'

  real_T TransferFcn1_C;               // Computed Parameter: TransferFcn1_C
                                          //  Referenced by: '<S4>/Transfer Fcn1'

  real_T Saturation1_UpperSat;         // Expression: 10
                                          //  Referenced by: '<S4>/Saturation1'

  real_T Saturation1_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S4>/Saturation1'

  real_T TransferFcn3_A;               // Computed Parameter: TransferFcn3_A
                                          //  Referenced by: '<S4>/Transfer Fcn3'

  real_T TransferFcn3_C;               // Computed Parameter: TransferFcn3_C
                                          //  Referenced by: '<S4>/Transfer Fcn3'

  real_T Saturation2_UpperSat;         // Expression: 10
                                          //  Referenced by: '<S4>/Saturation2'

  real_T Saturation2_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S4>/Saturation2'

  real_T TransferFcn4_A;               // Computed Parameter: TransferFcn4_A
                                          //  Referenced by: '<S4>/Transfer Fcn4'

  real_T TransferFcn4_C;               // Computed Parameter: TransferFcn4_C
                                          //  Referenced by: '<S4>/Transfer Fcn4'

  real_T Saturation3_UpperSat;         // Expression: 10
                                          //  Referenced by: '<S4>/Saturation3'

  real_T Saturation3_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S4>/Saturation3'

  real_T Gain_Gain_m;                  // Expression: -1
                                          //  Referenced by: '<S21>/Gain'

  real_T Merge_InitialOutput[4];       // Expression: [1 0 0 0]
                                          //  Referenced by: '<S8>/Merge'

  real_T Filter_gainval;               // Computed Parameter: Filter_gainval
                                          //  Referenced by: '<S345>/Filter'

  real_T Integrator_gainval;           // Computed Parameter: Integrator_gainval
                                          //  Referenced by: '<S350>/Integrator'

  real_T Integrator_gainval_n;       // Computed Parameter: Integrator_gainval_n
                                        //  Referenced by: '<S302>/Integrator'

  real_T Filter_gainval_l;             // Computed Parameter: Filter_gainval_l
                                          //  Referenced by: '<S297>/Filter'

  real_T Integrator_gainval_c;       // Computed Parameter: Integrator_gainval_c
                                        //  Referenced by: '<S254>/Integrator'

  real_T Filter_gainval_b;             // Computed Parameter: Filter_gainval_b
                                          //  Referenced by: '<S249>/Filter'

  real_T Hoverthrottle_Value;          // Expression: 3.13
                                          //  Referenced by: '<S4>/Hover throttle'

  real_T Saturation4_UpperSat;         // Expression: inf
                                          //  Referenced by: '<S4>/Saturation4'

  real_T Saturation4_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S4>/Saturation4'

  real_T Komega_Gain;                  // Expression: 0.5
                                          //  Referenced by: '<S4>/Komega'

  real_T Saturation5_UpperSat;         // Expression: inf
                                          //  Referenced by: '<S4>/Saturation5'

  real_T Saturation5_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S4>/Saturation5'

  real_T Komega1_Gain;                 // Expression: 0.5
                                          //  Referenced by: '<S4>/Komega1'

  real_T Saturation6_UpperSat;         // Expression: inf
                                          //  Referenced by: '<S4>/Saturation6'

  real_T Saturation6_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S4>/Saturation6'

  real_T Komega2_Gain;                 // Expression: 0.5
                                          //  Referenced by: '<S4>/Komega2'

  real_T Saturation7_UpperSat;         // Expression: inf
                                          //  Referenced by: '<S4>/Saturation7'

  real_T Saturation7_LowerSat;         // Expression: 0
                                          //  Referenced by: '<S4>/Saturation7'

  real_T Komega3_Gain;                 // Expression: 0.5
                                          //  Referenced by: '<S4>/Komega3'

};

// Real-time Model Data Structure
struct tag_RTM_genericQuad_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_genericQuad_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[28];
  real_T odeF[3][28];
  ODE3_IntgData intgData;

  //
  //  Sizes:
  //  The following substructure contains sizes information
  //  for many of the model attributes such as inputs, outputs,
  //  dwork, sample times, etc.

  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  //
  //  Timing:
  //  The following substructure contains information regarding
  //  the timing information for the model.

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

// Block parameters (default storage)
#ifdef __cplusplus

extern "C" {

#endif

  extern P_genericQuad_T genericQuad_P;

#ifdef __cplusplus

}
#endif

// Block signals (default storage)
#ifdef __cplusplus

extern "C" {

#endif

  extern struct B_genericQuad_T genericQuad_B;

#ifdef __cplusplus

}
#endif

// Continuous states (default storage)
extern X_genericQuad_T genericQuad_X;

// Block states (default storage)
extern struct DW_genericQuad_T genericQuad_DW;

#ifdef __cplusplus

extern "C" {

#endif

  // Model entry point functions
  extern void genericQuad_initialize(void);
  extern void genericQuad_step(void);
  extern void genericQuad_terminate(void);

#ifdef __cplusplus

}
#endif

// Real-time Model object
#ifdef __cplusplus

extern "C" {

#endif

  extern RT_MODEL_genericQuad_T *const genericQuad_M;

#ifdef __cplusplus

}
#endif

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S27>/Unit Conversion' : Unused code path elimination
//  Block '<S29>/Unit Conversion' : Unused code path elimination
//  Block '<S4>/Scope2' : Unused code path elimination
//  Block '<S4>/attitude' : Unused code path elimination
//  Block '<S4>/omega' : Unused code path elimination
//  Block '<S4>/pose' : Unused code path elimination
//  Block '<S4>/velocity' : Unused code path elimination
//  Block '<S4>/velocity1' : Unused code path elimination
//  Block '<S33>/Reshape (9) to [3x3] column-major' : Reshape block reduction
//  Block '<S35>/Reshape1' : Reshape block reduction
//  Block '<S35>/Reshape2' : Reshape block reduction
//  Block '<S36>/Reshape1' : Reshape block reduction
//  Block '<S36>/Reshape2' : Reshape block reduction
//  Block '<S24>/Reshape' : Reshape block reduction
//  Block '<S24>/Reshape1' : Reshape block reduction
//  Block '<S28>/Unit Conversion' : Eliminated nontunable gain of 1
//  Block '<S30>/Reshape1' : Reshape block reduction
//  Block '<S30>/Reshape2' : Reshape block reduction
//  Block '<S8>/Reshape 3x3 -> 9' : Reshape block reduction
//  Block '<S70>/Reshape' : Reshape block reduction


//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Use the MATLAB hilite_system command to trace the generated code back
//  to the model.  For example,
//
//  hilite_system('<S3>')    - opens system 3
//  hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'genericQuad'
//  '<S1>'   : 'genericQuad/Blank Message'
//  '<S2>'   : 'genericQuad/Publish'
//  '<S3>'   : 'genericQuad/Subscribe1'
//  '<S4>'   : 'genericQuad/Subsystem1'
//  '<S5>'   : 'genericQuad/time to sec & nsec'
//  '<S6>'   : 'genericQuad/Subscribe1/Enabled Subsystem'
//  '<S7>'   : 'genericQuad/Subsystem1/6DOF (Euler Angles)2'
//  '<S8>'   : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions'
//  '<S9>'   : 'genericQuad/Subsystem1/Gravity1'
//  '<S10>'  : 'genericQuad/Subsystem1/Mapping'
//  '<S11>'  : 'genericQuad/Subsystem1/Mixer'
//  '<S12>'  : 'genericQuad/Subsystem1/PID Velocity x'
//  '<S13>'  : 'genericQuad/Subsystem1/PID Velocity y'
//  '<S14>'  : 'genericQuad/Subsystem1/PID Velocity z'
//  '<S15>'  : 'genericQuad/Subsystem1/PID angula yaw'
//  '<S16>'  : 'genericQuad/Subsystem1/PID angular pitch'
//  '<S17>'  : 'genericQuad/Subsystem1/PID angular roll'
//  '<S18>'  : 'genericQuad/Subsystem1/PID attitude'
//  '<S19>'  : 'genericQuad/Subsystem1/R_EW'
//  '<S20>'  : 'genericQuad/Subsystem1/Subsystem'
//  '<S21>'  : 'genericQuad/Subsystem1/Subsystem1'
//  '<S22>'  : 'genericQuad/Subsystem1/euler to body'
//  '<S23>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Calculate DCM & Euler Angles'
//  '<S24>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot'
//  '<S25>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Determine Force,  Mass & Inertia'
//  '<S26>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Vbxw'
//  '<S27>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Velocity Conversion'
//  '<S28>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Velocity Conversion1'
//  '<S29>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Velocity Conversion2'
//  '<S30>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/transform to Inertial axes '
//  '<S31>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix'
//  '<S32>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Calculate DCM & Euler Angles/phidot thetadot psidot'
//  '<S33>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
//  '<S34>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot/3x3 Cross Product'
//  '<S35>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot/I x w'
//  '<S36>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot/I x w1'
//  '<S37>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot/3x3 Cross Product/Subsystem'
//  '<S38>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Calculate omega_dot/3x3 Cross Product/Subsystem1'
//  '<S39>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Vbxw/Subsystem'
//  '<S40>'  : 'genericQuad/Subsystem1/6DOF (Euler Angles)2/Vbxw/Subsystem1'
//  '<S41>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace'
//  '<S42>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Positive Trace'
//  '<S43>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM'
//  '<S44>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/trace(DCM)'
//  '<S45>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)'
//  '<S46>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)'
//  '<S47>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)'
//  '<S48>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/diag(DCM)'
//  '<S49>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) -sin(theta)'
//  '<S50>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(theta)sin(phi) - (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
//  '<S51>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(theta)sin(psi) + (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
//  '<S52>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/if s~=0; s=0.5//s'
//  '<S53>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/u(1) -(u(5)+u(9)) +1'
//  '<S54>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) +sin(theta)'
//  '<S55>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(theta)sin(phi) + (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
//  '<S56>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(theta)sin(psi) + (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
//  '<S57>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/if s~=0; s=0.5//s'
//  '<S58>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/u(5) -(u(1)+u(9)) +1'
//  '<S59>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) -sin(theta)'
//  '<S60>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(theta)sin(phi) + (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
//  '<S61>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(theta)sin(psi) - (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
//  '<S62>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/if s~=0; s=0.5//s'
//  '<S63>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/u(9) -(u(1)+u(5)) +1'
//  '<S64>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) +sin(theta)'
//  '<S65>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(theta)sin(phi) - (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
//  '<S66>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(theta)sin(psi) - (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
//  '<S67>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error'
//  '<S68>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal'
//  '<S69>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper'
//  '<S70>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotOrthogonal'
//  '<S71>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper'
//  '<S72>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal/Error'
//  '<S73>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal/Warning'
//  '<S74>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper/Error'
//  '<S75>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper/Warning'
//  '<S76>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotOrthogonal/transpose*dcm ~= eye(3)'
//  '<S77>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper/Determinant of 3x3 Matrix'
//  '<S78>'  : 'genericQuad/Subsystem1/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper/determinant does not equal 1'
//  '<S79>'  : 'genericQuad/Subsystem1/PID Velocity x/Anti-windup'
//  '<S80>'  : 'genericQuad/Subsystem1/PID Velocity x/D Gain'
//  '<S81>'  : 'genericQuad/Subsystem1/PID Velocity x/Filter'
//  '<S82>'  : 'genericQuad/Subsystem1/PID Velocity x/Filter ICs'
//  '<S83>'  : 'genericQuad/Subsystem1/PID Velocity x/I Gain'
//  '<S84>'  : 'genericQuad/Subsystem1/PID Velocity x/Ideal P Gain'
//  '<S85>'  : 'genericQuad/Subsystem1/PID Velocity x/Ideal P Gain Fdbk'
//  '<S86>'  : 'genericQuad/Subsystem1/PID Velocity x/Integrator'
//  '<S87>'  : 'genericQuad/Subsystem1/PID Velocity x/Integrator ICs'
//  '<S88>'  : 'genericQuad/Subsystem1/PID Velocity x/N Copy'
//  '<S89>'  : 'genericQuad/Subsystem1/PID Velocity x/N Gain'
//  '<S90>'  : 'genericQuad/Subsystem1/PID Velocity x/P Copy'
//  '<S91>'  : 'genericQuad/Subsystem1/PID Velocity x/Parallel P Gain'
//  '<S92>'  : 'genericQuad/Subsystem1/PID Velocity x/Reset Signal'
//  '<S93>'  : 'genericQuad/Subsystem1/PID Velocity x/Saturation'
//  '<S94>'  : 'genericQuad/Subsystem1/PID Velocity x/Saturation Fdbk'
//  '<S95>'  : 'genericQuad/Subsystem1/PID Velocity x/Sum'
//  '<S96>'  : 'genericQuad/Subsystem1/PID Velocity x/Sum Fdbk'
//  '<S97>'  : 'genericQuad/Subsystem1/PID Velocity x/Tracking Mode'
//  '<S98>'  : 'genericQuad/Subsystem1/PID Velocity x/Tracking Mode Sum'
//  '<S99>'  : 'genericQuad/Subsystem1/PID Velocity x/Tsamp - Integral'
//  '<S100>' : 'genericQuad/Subsystem1/PID Velocity x/Tsamp - Ngain'
//  '<S101>' : 'genericQuad/Subsystem1/PID Velocity x/postSat Signal'
//  '<S102>' : 'genericQuad/Subsystem1/PID Velocity x/preSat Signal'
//  '<S103>' : 'genericQuad/Subsystem1/PID Velocity x/Anti-windup/Passthrough'
//  '<S104>' : 'genericQuad/Subsystem1/PID Velocity x/D Gain/Internal Parameters'
//  '<S105>' : 'genericQuad/Subsystem1/PID Velocity x/Filter/Cont. Filter'
//  '<S106>' : 'genericQuad/Subsystem1/PID Velocity x/Filter ICs/Internal IC - Filter'
//  '<S107>' : 'genericQuad/Subsystem1/PID Velocity x/I Gain/Internal Parameters'
//  '<S108>' : 'genericQuad/Subsystem1/PID Velocity x/Ideal P Gain/Passthrough'
//  '<S109>' : 'genericQuad/Subsystem1/PID Velocity x/Ideal P Gain Fdbk/Disabled'
//  '<S110>' : 'genericQuad/Subsystem1/PID Velocity x/Integrator/Continuous'
//  '<S111>' : 'genericQuad/Subsystem1/PID Velocity x/Integrator ICs/Internal IC'
//  '<S112>' : 'genericQuad/Subsystem1/PID Velocity x/N Copy/Disabled'
//  '<S113>' : 'genericQuad/Subsystem1/PID Velocity x/N Gain/Internal Parameters'
//  '<S114>' : 'genericQuad/Subsystem1/PID Velocity x/P Copy/Disabled'
//  '<S115>' : 'genericQuad/Subsystem1/PID Velocity x/Parallel P Gain/Internal Parameters'
//  '<S116>' : 'genericQuad/Subsystem1/PID Velocity x/Reset Signal/Disabled'
//  '<S117>' : 'genericQuad/Subsystem1/PID Velocity x/Saturation/Passthrough'
//  '<S118>' : 'genericQuad/Subsystem1/PID Velocity x/Saturation Fdbk/Disabled'
//  '<S119>' : 'genericQuad/Subsystem1/PID Velocity x/Sum/Sum_PID'
//  '<S120>' : 'genericQuad/Subsystem1/PID Velocity x/Sum Fdbk/Disabled'
//  '<S121>' : 'genericQuad/Subsystem1/PID Velocity x/Tracking Mode/Disabled'
//  '<S122>' : 'genericQuad/Subsystem1/PID Velocity x/Tracking Mode Sum/Passthrough'
//  '<S123>' : 'genericQuad/Subsystem1/PID Velocity x/Tsamp - Integral/Passthrough'
//  '<S124>' : 'genericQuad/Subsystem1/PID Velocity x/Tsamp - Ngain/Passthrough'
//  '<S125>' : 'genericQuad/Subsystem1/PID Velocity x/postSat Signal/Forward_Path'
//  '<S126>' : 'genericQuad/Subsystem1/PID Velocity x/preSat Signal/Forward_Path'
//  '<S127>' : 'genericQuad/Subsystem1/PID Velocity y/Anti-windup'
//  '<S128>' : 'genericQuad/Subsystem1/PID Velocity y/D Gain'
//  '<S129>' : 'genericQuad/Subsystem1/PID Velocity y/Filter'
//  '<S130>' : 'genericQuad/Subsystem1/PID Velocity y/Filter ICs'
//  '<S131>' : 'genericQuad/Subsystem1/PID Velocity y/I Gain'
//  '<S132>' : 'genericQuad/Subsystem1/PID Velocity y/Ideal P Gain'
//  '<S133>' : 'genericQuad/Subsystem1/PID Velocity y/Ideal P Gain Fdbk'
//  '<S134>' : 'genericQuad/Subsystem1/PID Velocity y/Integrator'
//  '<S135>' : 'genericQuad/Subsystem1/PID Velocity y/Integrator ICs'
//  '<S136>' : 'genericQuad/Subsystem1/PID Velocity y/N Copy'
//  '<S137>' : 'genericQuad/Subsystem1/PID Velocity y/N Gain'
//  '<S138>' : 'genericQuad/Subsystem1/PID Velocity y/P Copy'
//  '<S139>' : 'genericQuad/Subsystem1/PID Velocity y/Parallel P Gain'
//  '<S140>' : 'genericQuad/Subsystem1/PID Velocity y/Reset Signal'
//  '<S141>' : 'genericQuad/Subsystem1/PID Velocity y/Saturation'
//  '<S142>' : 'genericQuad/Subsystem1/PID Velocity y/Saturation Fdbk'
//  '<S143>' : 'genericQuad/Subsystem1/PID Velocity y/Sum'
//  '<S144>' : 'genericQuad/Subsystem1/PID Velocity y/Sum Fdbk'
//  '<S145>' : 'genericQuad/Subsystem1/PID Velocity y/Tracking Mode'
//  '<S146>' : 'genericQuad/Subsystem1/PID Velocity y/Tracking Mode Sum'
//  '<S147>' : 'genericQuad/Subsystem1/PID Velocity y/Tsamp - Integral'
//  '<S148>' : 'genericQuad/Subsystem1/PID Velocity y/Tsamp - Ngain'
//  '<S149>' : 'genericQuad/Subsystem1/PID Velocity y/postSat Signal'
//  '<S150>' : 'genericQuad/Subsystem1/PID Velocity y/preSat Signal'
//  '<S151>' : 'genericQuad/Subsystem1/PID Velocity y/Anti-windup/Passthrough'
//  '<S152>' : 'genericQuad/Subsystem1/PID Velocity y/D Gain/Internal Parameters'
//  '<S153>' : 'genericQuad/Subsystem1/PID Velocity y/Filter/Cont. Filter'
//  '<S154>' : 'genericQuad/Subsystem1/PID Velocity y/Filter ICs/Internal IC - Filter'
//  '<S155>' : 'genericQuad/Subsystem1/PID Velocity y/I Gain/Internal Parameters'
//  '<S156>' : 'genericQuad/Subsystem1/PID Velocity y/Ideal P Gain/Passthrough'
//  '<S157>' : 'genericQuad/Subsystem1/PID Velocity y/Ideal P Gain Fdbk/Disabled'
//  '<S158>' : 'genericQuad/Subsystem1/PID Velocity y/Integrator/Continuous'
//  '<S159>' : 'genericQuad/Subsystem1/PID Velocity y/Integrator ICs/Internal IC'
//  '<S160>' : 'genericQuad/Subsystem1/PID Velocity y/N Copy/Disabled'
//  '<S161>' : 'genericQuad/Subsystem1/PID Velocity y/N Gain/Internal Parameters'
//  '<S162>' : 'genericQuad/Subsystem1/PID Velocity y/P Copy/Disabled'
//  '<S163>' : 'genericQuad/Subsystem1/PID Velocity y/Parallel P Gain/Internal Parameters'
//  '<S164>' : 'genericQuad/Subsystem1/PID Velocity y/Reset Signal/Disabled'
//  '<S165>' : 'genericQuad/Subsystem1/PID Velocity y/Saturation/Passthrough'
//  '<S166>' : 'genericQuad/Subsystem1/PID Velocity y/Saturation Fdbk/Disabled'
//  '<S167>' : 'genericQuad/Subsystem1/PID Velocity y/Sum/Sum_PID'
//  '<S168>' : 'genericQuad/Subsystem1/PID Velocity y/Sum Fdbk/Disabled'
//  '<S169>' : 'genericQuad/Subsystem1/PID Velocity y/Tracking Mode/Disabled'
//  '<S170>' : 'genericQuad/Subsystem1/PID Velocity y/Tracking Mode Sum/Passthrough'
//  '<S171>' : 'genericQuad/Subsystem1/PID Velocity y/Tsamp - Integral/Passthrough'
//  '<S172>' : 'genericQuad/Subsystem1/PID Velocity y/Tsamp - Ngain/Passthrough'
//  '<S173>' : 'genericQuad/Subsystem1/PID Velocity y/postSat Signal/Forward_Path'
//  '<S174>' : 'genericQuad/Subsystem1/PID Velocity y/preSat Signal/Forward_Path'
//  '<S175>' : 'genericQuad/Subsystem1/PID Velocity z/Anti-windup'
//  '<S176>' : 'genericQuad/Subsystem1/PID Velocity z/D Gain'
//  '<S177>' : 'genericQuad/Subsystem1/PID Velocity z/Filter'
//  '<S178>' : 'genericQuad/Subsystem1/PID Velocity z/Filter ICs'
//  '<S179>' : 'genericQuad/Subsystem1/PID Velocity z/I Gain'
//  '<S180>' : 'genericQuad/Subsystem1/PID Velocity z/Ideal P Gain'
//  '<S181>' : 'genericQuad/Subsystem1/PID Velocity z/Ideal P Gain Fdbk'
//  '<S182>' : 'genericQuad/Subsystem1/PID Velocity z/Integrator'
//  '<S183>' : 'genericQuad/Subsystem1/PID Velocity z/Integrator ICs'
//  '<S184>' : 'genericQuad/Subsystem1/PID Velocity z/N Copy'
//  '<S185>' : 'genericQuad/Subsystem1/PID Velocity z/N Gain'
//  '<S186>' : 'genericQuad/Subsystem1/PID Velocity z/P Copy'
//  '<S187>' : 'genericQuad/Subsystem1/PID Velocity z/Parallel P Gain'
//  '<S188>' : 'genericQuad/Subsystem1/PID Velocity z/Reset Signal'
//  '<S189>' : 'genericQuad/Subsystem1/PID Velocity z/Saturation'
//  '<S190>' : 'genericQuad/Subsystem1/PID Velocity z/Saturation Fdbk'
//  '<S191>' : 'genericQuad/Subsystem1/PID Velocity z/Sum'
//  '<S192>' : 'genericQuad/Subsystem1/PID Velocity z/Sum Fdbk'
//  '<S193>' : 'genericQuad/Subsystem1/PID Velocity z/Tracking Mode'
//  '<S194>' : 'genericQuad/Subsystem1/PID Velocity z/Tracking Mode Sum'
//  '<S195>' : 'genericQuad/Subsystem1/PID Velocity z/Tsamp - Integral'
//  '<S196>' : 'genericQuad/Subsystem1/PID Velocity z/Tsamp - Ngain'
//  '<S197>' : 'genericQuad/Subsystem1/PID Velocity z/postSat Signal'
//  '<S198>' : 'genericQuad/Subsystem1/PID Velocity z/preSat Signal'
//  '<S199>' : 'genericQuad/Subsystem1/PID Velocity z/Anti-windup/Passthrough'
//  '<S200>' : 'genericQuad/Subsystem1/PID Velocity z/D Gain/Internal Parameters'
//  '<S201>' : 'genericQuad/Subsystem1/PID Velocity z/Filter/Cont. Filter'
//  '<S202>' : 'genericQuad/Subsystem1/PID Velocity z/Filter ICs/Internal IC - Filter'
//  '<S203>' : 'genericQuad/Subsystem1/PID Velocity z/I Gain/Internal Parameters'
//  '<S204>' : 'genericQuad/Subsystem1/PID Velocity z/Ideal P Gain/Passthrough'
//  '<S205>' : 'genericQuad/Subsystem1/PID Velocity z/Ideal P Gain Fdbk/Disabled'
//  '<S206>' : 'genericQuad/Subsystem1/PID Velocity z/Integrator/Continuous'
//  '<S207>' : 'genericQuad/Subsystem1/PID Velocity z/Integrator ICs/Internal IC'
//  '<S208>' : 'genericQuad/Subsystem1/PID Velocity z/N Copy/Disabled'
//  '<S209>' : 'genericQuad/Subsystem1/PID Velocity z/N Gain/Internal Parameters'
//  '<S210>' : 'genericQuad/Subsystem1/PID Velocity z/P Copy/Disabled'
//  '<S211>' : 'genericQuad/Subsystem1/PID Velocity z/Parallel P Gain/Internal Parameters'
//  '<S212>' : 'genericQuad/Subsystem1/PID Velocity z/Reset Signal/Disabled'
//  '<S213>' : 'genericQuad/Subsystem1/PID Velocity z/Saturation/Passthrough'
//  '<S214>' : 'genericQuad/Subsystem1/PID Velocity z/Saturation Fdbk/Disabled'
//  '<S215>' : 'genericQuad/Subsystem1/PID Velocity z/Sum/Sum_PID'
//  '<S216>' : 'genericQuad/Subsystem1/PID Velocity z/Sum Fdbk/Disabled'
//  '<S217>' : 'genericQuad/Subsystem1/PID Velocity z/Tracking Mode/Disabled'
//  '<S218>' : 'genericQuad/Subsystem1/PID Velocity z/Tracking Mode Sum/Passthrough'
//  '<S219>' : 'genericQuad/Subsystem1/PID Velocity z/Tsamp - Integral/Passthrough'
//  '<S220>' : 'genericQuad/Subsystem1/PID Velocity z/Tsamp - Ngain/Passthrough'
//  '<S221>' : 'genericQuad/Subsystem1/PID Velocity z/postSat Signal/Forward_Path'
//  '<S222>' : 'genericQuad/Subsystem1/PID Velocity z/preSat Signal/Forward_Path'
//  '<S223>' : 'genericQuad/Subsystem1/PID angula yaw/Anti-windup'
//  '<S224>' : 'genericQuad/Subsystem1/PID angula yaw/D Gain'
//  '<S225>' : 'genericQuad/Subsystem1/PID angula yaw/Filter'
//  '<S226>' : 'genericQuad/Subsystem1/PID angula yaw/Filter ICs'
//  '<S227>' : 'genericQuad/Subsystem1/PID angula yaw/I Gain'
//  '<S228>' : 'genericQuad/Subsystem1/PID angula yaw/Ideal P Gain'
//  '<S229>' : 'genericQuad/Subsystem1/PID angula yaw/Ideal P Gain Fdbk'
//  '<S230>' : 'genericQuad/Subsystem1/PID angula yaw/Integrator'
//  '<S231>' : 'genericQuad/Subsystem1/PID angula yaw/Integrator ICs'
//  '<S232>' : 'genericQuad/Subsystem1/PID angula yaw/N Copy'
//  '<S233>' : 'genericQuad/Subsystem1/PID angula yaw/N Gain'
//  '<S234>' : 'genericQuad/Subsystem1/PID angula yaw/P Copy'
//  '<S235>' : 'genericQuad/Subsystem1/PID angula yaw/Parallel P Gain'
//  '<S236>' : 'genericQuad/Subsystem1/PID angula yaw/Reset Signal'
//  '<S237>' : 'genericQuad/Subsystem1/PID angula yaw/Saturation'
//  '<S238>' : 'genericQuad/Subsystem1/PID angula yaw/Saturation Fdbk'
//  '<S239>' : 'genericQuad/Subsystem1/PID angula yaw/Sum'
//  '<S240>' : 'genericQuad/Subsystem1/PID angula yaw/Sum Fdbk'
//  '<S241>' : 'genericQuad/Subsystem1/PID angula yaw/Tracking Mode'
//  '<S242>' : 'genericQuad/Subsystem1/PID angula yaw/Tracking Mode Sum'
//  '<S243>' : 'genericQuad/Subsystem1/PID angula yaw/Tsamp - Integral'
//  '<S244>' : 'genericQuad/Subsystem1/PID angula yaw/Tsamp - Ngain'
//  '<S245>' : 'genericQuad/Subsystem1/PID angula yaw/postSat Signal'
//  '<S246>' : 'genericQuad/Subsystem1/PID angula yaw/preSat Signal'
//  '<S247>' : 'genericQuad/Subsystem1/PID angula yaw/Anti-windup/Passthrough'
//  '<S248>' : 'genericQuad/Subsystem1/PID angula yaw/D Gain/Internal Parameters'
//  '<S249>' : 'genericQuad/Subsystem1/PID angula yaw/Filter/Disc. Forward Euler Filter'
//  '<S250>' : 'genericQuad/Subsystem1/PID angula yaw/Filter ICs/Internal IC - Filter'
//  '<S251>' : 'genericQuad/Subsystem1/PID angula yaw/I Gain/Internal Parameters'
//  '<S252>' : 'genericQuad/Subsystem1/PID angula yaw/Ideal P Gain/Internal Parameters'
//  '<S253>' : 'genericQuad/Subsystem1/PID angula yaw/Ideal P Gain Fdbk/Disabled'
//  '<S254>' : 'genericQuad/Subsystem1/PID angula yaw/Integrator/Discrete'
//  '<S255>' : 'genericQuad/Subsystem1/PID angula yaw/Integrator ICs/Internal IC'
//  '<S256>' : 'genericQuad/Subsystem1/PID angula yaw/N Copy/Disabled'
//  '<S257>' : 'genericQuad/Subsystem1/PID angula yaw/N Gain/Internal Parameters'
//  '<S258>' : 'genericQuad/Subsystem1/PID angula yaw/P Copy/Disabled'
//  '<S259>' : 'genericQuad/Subsystem1/PID angula yaw/Parallel P Gain/Passthrough'
//  '<S260>' : 'genericQuad/Subsystem1/PID angula yaw/Reset Signal/Disabled'
//  '<S261>' : 'genericQuad/Subsystem1/PID angula yaw/Saturation/Passthrough'
//  '<S262>' : 'genericQuad/Subsystem1/PID angula yaw/Saturation Fdbk/Disabled'
//  '<S263>' : 'genericQuad/Subsystem1/PID angula yaw/Sum/Sum_PID'
//  '<S264>' : 'genericQuad/Subsystem1/PID angula yaw/Sum Fdbk/Disabled'
//  '<S265>' : 'genericQuad/Subsystem1/PID angula yaw/Tracking Mode/Disabled'
//  '<S266>' : 'genericQuad/Subsystem1/PID angula yaw/Tracking Mode Sum/Passthrough'
//  '<S267>' : 'genericQuad/Subsystem1/PID angula yaw/Tsamp - Integral/Passthrough'
//  '<S268>' : 'genericQuad/Subsystem1/PID angula yaw/Tsamp - Ngain/Passthrough'
//  '<S269>' : 'genericQuad/Subsystem1/PID angula yaw/postSat Signal/Forward_Path'
//  '<S270>' : 'genericQuad/Subsystem1/PID angula yaw/preSat Signal/Forward_Path'
//  '<S271>' : 'genericQuad/Subsystem1/PID angular pitch/Anti-windup'
//  '<S272>' : 'genericQuad/Subsystem1/PID angular pitch/D Gain'
//  '<S273>' : 'genericQuad/Subsystem1/PID angular pitch/Filter'
//  '<S274>' : 'genericQuad/Subsystem1/PID angular pitch/Filter ICs'
//  '<S275>' : 'genericQuad/Subsystem1/PID angular pitch/I Gain'
//  '<S276>' : 'genericQuad/Subsystem1/PID angular pitch/Ideal P Gain'
//  '<S277>' : 'genericQuad/Subsystem1/PID angular pitch/Ideal P Gain Fdbk'
//  '<S278>' : 'genericQuad/Subsystem1/PID angular pitch/Integrator'
//  '<S279>' : 'genericQuad/Subsystem1/PID angular pitch/Integrator ICs'
//  '<S280>' : 'genericQuad/Subsystem1/PID angular pitch/N Copy'
//  '<S281>' : 'genericQuad/Subsystem1/PID angular pitch/N Gain'
//  '<S282>' : 'genericQuad/Subsystem1/PID angular pitch/P Copy'
//  '<S283>' : 'genericQuad/Subsystem1/PID angular pitch/Parallel P Gain'
//  '<S284>' : 'genericQuad/Subsystem1/PID angular pitch/Reset Signal'
//  '<S285>' : 'genericQuad/Subsystem1/PID angular pitch/Saturation'
//  '<S286>' : 'genericQuad/Subsystem1/PID angular pitch/Saturation Fdbk'
//  '<S287>' : 'genericQuad/Subsystem1/PID angular pitch/Sum'
//  '<S288>' : 'genericQuad/Subsystem1/PID angular pitch/Sum Fdbk'
//  '<S289>' : 'genericQuad/Subsystem1/PID angular pitch/Tracking Mode'
//  '<S290>' : 'genericQuad/Subsystem1/PID angular pitch/Tracking Mode Sum'
//  '<S291>' : 'genericQuad/Subsystem1/PID angular pitch/Tsamp - Integral'
//  '<S292>' : 'genericQuad/Subsystem1/PID angular pitch/Tsamp - Ngain'
//  '<S293>' : 'genericQuad/Subsystem1/PID angular pitch/postSat Signal'
//  '<S294>' : 'genericQuad/Subsystem1/PID angular pitch/preSat Signal'
//  '<S295>' : 'genericQuad/Subsystem1/PID angular pitch/Anti-windup/Passthrough'
//  '<S296>' : 'genericQuad/Subsystem1/PID angular pitch/D Gain/Internal Parameters'
//  '<S297>' : 'genericQuad/Subsystem1/PID angular pitch/Filter/Disc. Forward Euler Filter'
//  '<S298>' : 'genericQuad/Subsystem1/PID angular pitch/Filter ICs/Internal IC - Filter'
//  '<S299>' : 'genericQuad/Subsystem1/PID angular pitch/I Gain/Internal Parameters'
//  '<S300>' : 'genericQuad/Subsystem1/PID angular pitch/Ideal P Gain/Internal Parameters'
//  '<S301>' : 'genericQuad/Subsystem1/PID angular pitch/Ideal P Gain Fdbk/Disabled'
//  '<S302>' : 'genericQuad/Subsystem1/PID angular pitch/Integrator/Discrete'
//  '<S303>' : 'genericQuad/Subsystem1/PID angular pitch/Integrator ICs/Internal IC'
//  '<S304>' : 'genericQuad/Subsystem1/PID angular pitch/N Copy/Disabled'
//  '<S305>' : 'genericQuad/Subsystem1/PID angular pitch/N Gain/Internal Parameters'
//  '<S306>' : 'genericQuad/Subsystem1/PID angular pitch/P Copy/Disabled'
//  '<S307>' : 'genericQuad/Subsystem1/PID angular pitch/Parallel P Gain/Passthrough'
//  '<S308>' : 'genericQuad/Subsystem1/PID angular pitch/Reset Signal/Disabled'
//  '<S309>' : 'genericQuad/Subsystem1/PID angular pitch/Saturation/Passthrough'
//  '<S310>' : 'genericQuad/Subsystem1/PID angular pitch/Saturation Fdbk/Disabled'
//  '<S311>' : 'genericQuad/Subsystem1/PID angular pitch/Sum/Sum_PID'
//  '<S312>' : 'genericQuad/Subsystem1/PID angular pitch/Sum Fdbk/Disabled'
//  '<S313>' : 'genericQuad/Subsystem1/PID angular pitch/Tracking Mode/Disabled'
//  '<S314>' : 'genericQuad/Subsystem1/PID angular pitch/Tracking Mode Sum/Passthrough'
//  '<S315>' : 'genericQuad/Subsystem1/PID angular pitch/Tsamp - Integral/Passthrough'
//  '<S316>' : 'genericQuad/Subsystem1/PID angular pitch/Tsamp - Ngain/Passthrough'
//  '<S317>' : 'genericQuad/Subsystem1/PID angular pitch/postSat Signal/Forward_Path'
//  '<S318>' : 'genericQuad/Subsystem1/PID angular pitch/preSat Signal/Forward_Path'
//  '<S319>' : 'genericQuad/Subsystem1/PID angular roll/Anti-windup'
//  '<S320>' : 'genericQuad/Subsystem1/PID angular roll/D Gain'
//  '<S321>' : 'genericQuad/Subsystem1/PID angular roll/Filter'
//  '<S322>' : 'genericQuad/Subsystem1/PID angular roll/Filter ICs'
//  '<S323>' : 'genericQuad/Subsystem1/PID angular roll/I Gain'
//  '<S324>' : 'genericQuad/Subsystem1/PID angular roll/Ideal P Gain'
//  '<S325>' : 'genericQuad/Subsystem1/PID angular roll/Ideal P Gain Fdbk'
//  '<S326>' : 'genericQuad/Subsystem1/PID angular roll/Integrator'
//  '<S327>' : 'genericQuad/Subsystem1/PID angular roll/Integrator ICs'
//  '<S328>' : 'genericQuad/Subsystem1/PID angular roll/N Copy'
//  '<S329>' : 'genericQuad/Subsystem1/PID angular roll/N Gain'
//  '<S330>' : 'genericQuad/Subsystem1/PID angular roll/P Copy'
//  '<S331>' : 'genericQuad/Subsystem1/PID angular roll/Parallel P Gain'
//  '<S332>' : 'genericQuad/Subsystem1/PID angular roll/Reset Signal'
//  '<S333>' : 'genericQuad/Subsystem1/PID angular roll/Saturation'
//  '<S334>' : 'genericQuad/Subsystem1/PID angular roll/Saturation Fdbk'
//  '<S335>' : 'genericQuad/Subsystem1/PID angular roll/Sum'
//  '<S336>' : 'genericQuad/Subsystem1/PID angular roll/Sum Fdbk'
//  '<S337>' : 'genericQuad/Subsystem1/PID angular roll/Tracking Mode'
//  '<S338>' : 'genericQuad/Subsystem1/PID angular roll/Tracking Mode Sum'
//  '<S339>' : 'genericQuad/Subsystem1/PID angular roll/Tsamp - Integral'
//  '<S340>' : 'genericQuad/Subsystem1/PID angular roll/Tsamp - Ngain'
//  '<S341>' : 'genericQuad/Subsystem1/PID angular roll/postSat Signal'
//  '<S342>' : 'genericQuad/Subsystem1/PID angular roll/preSat Signal'
//  '<S343>' : 'genericQuad/Subsystem1/PID angular roll/Anti-windup/Passthrough'
//  '<S344>' : 'genericQuad/Subsystem1/PID angular roll/D Gain/Internal Parameters'
//  '<S345>' : 'genericQuad/Subsystem1/PID angular roll/Filter/Disc. Forward Euler Filter'
//  '<S346>' : 'genericQuad/Subsystem1/PID angular roll/Filter ICs/Internal IC - Filter'
//  '<S347>' : 'genericQuad/Subsystem1/PID angular roll/I Gain/Internal Parameters'
//  '<S348>' : 'genericQuad/Subsystem1/PID angular roll/Ideal P Gain/Internal Parameters'
//  '<S349>' : 'genericQuad/Subsystem1/PID angular roll/Ideal P Gain Fdbk/Disabled'
//  '<S350>' : 'genericQuad/Subsystem1/PID angular roll/Integrator/Discrete'
//  '<S351>' : 'genericQuad/Subsystem1/PID angular roll/Integrator ICs/Internal IC'
//  '<S352>' : 'genericQuad/Subsystem1/PID angular roll/N Copy/Disabled'
//  '<S353>' : 'genericQuad/Subsystem1/PID angular roll/N Gain/Internal Parameters'
//  '<S354>' : 'genericQuad/Subsystem1/PID angular roll/P Copy/Disabled'
//  '<S355>' : 'genericQuad/Subsystem1/PID angular roll/Parallel P Gain/Passthrough'
//  '<S356>' : 'genericQuad/Subsystem1/PID angular roll/Reset Signal/Disabled'
//  '<S357>' : 'genericQuad/Subsystem1/PID angular roll/Saturation/Passthrough'
//  '<S358>' : 'genericQuad/Subsystem1/PID angular roll/Saturation Fdbk/Disabled'
//  '<S359>' : 'genericQuad/Subsystem1/PID angular roll/Sum/Sum_PID'
//  '<S360>' : 'genericQuad/Subsystem1/PID angular roll/Sum Fdbk/Disabled'
//  '<S361>' : 'genericQuad/Subsystem1/PID angular roll/Tracking Mode/Disabled'
//  '<S362>' : 'genericQuad/Subsystem1/PID angular roll/Tracking Mode Sum/Passthrough'
//  '<S363>' : 'genericQuad/Subsystem1/PID angular roll/Tsamp - Integral/Passthrough'
//  '<S364>' : 'genericQuad/Subsystem1/PID angular roll/Tsamp - Ngain/Passthrough'
//  '<S365>' : 'genericQuad/Subsystem1/PID angular roll/postSat Signal/Forward_Path'
//  '<S366>' : 'genericQuad/Subsystem1/PID angular roll/preSat Signal/Forward_Path'
//  '<S367>' : 'genericQuad/Subsystem1/PID attitude/Anti-windup'
//  '<S368>' : 'genericQuad/Subsystem1/PID attitude/D Gain'
//  '<S369>' : 'genericQuad/Subsystem1/PID attitude/Filter'
//  '<S370>' : 'genericQuad/Subsystem1/PID attitude/Filter ICs'
//  '<S371>' : 'genericQuad/Subsystem1/PID attitude/I Gain'
//  '<S372>' : 'genericQuad/Subsystem1/PID attitude/Ideal P Gain'
//  '<S373>' : 'genericQuad/Subsystem1/PID attitude/Ideal P Gain Fdbk'
//  '<S374>' : 'genericQuad/Subsystem1/PID attitude/Integrator'
//  '<S375>' : 'genericQuad/Subsystem1/PID attitude/Integrator ICs'
//  '<S376>' : 'genericQuad/Subsystem1/PID attitude/N Copy'
//  '<S377>' : 'genericQuad/Subsystem1/PID attitude/N Gain'
//  '<S378>' : 'genericQuad/Subsystem1/PID attitude/P Copy'
//  '<S379>' : 'genericQuad/Subsystem1/PID attitude/Parallel P Gain'
//  '<S380>' : 'genericQuad/Subsystem1/PID attitude/Reset Signal'
//  '<S381>' : 'genericQuad/Subsystem1/PID attitude/Saturation'
//  '<S382>' : 'genericQuad/Subsystem1/PID attitude/Saturation Fdbk'
//  '<S383>' : 'genericQuad/Subsystem1/PID attitude/Sum'
//  '<S384>' : 'genericQuad/Subsystem1/PID attitude/Sum Fdbk'
//  '<S385>' : 'genericQuad/Subsystem1/PID attitude/Tracking Mode'
//  '<S386>' : 'genericQuad/Subsystem1/PID attitude/Tracking Mode Sum'
//  '<S387>' : 'genericQuad/Subsystem1/PID attitude/Tsamp - Integral'
//  '<S388>' : 'genericQuad/Subsystem1/PID attitude/Tsamp - Ngain'
//  '<S389>' : 'genericQuad/Subsystem1/PID attitude/postSat Signal'
//  '<S390>' : 'genericQuad/Subsystem1/PID attitude/preSat Signal'
//  '<S391>' : 'genericQuad/Subsystem1/PID attitude/Anti-windup/Passthrough'
//  '<S392>' : 'genericQuad/Subsystem1/PID attitude/D Gain/Internal Parameters'
//  '<S393>' : 'genericQuad/Subsystem1/PID attitude/Filter/Cont. Filter'
//  '<S394>' : 'genericQuad/Subsystem1/PID attitude/Filter ICs/Internal IC - Filter'
//  '<S395>' : 'genericQuad/Subsystem1/PID attitude/I Gain/Internal Parameters'
//  '<S396>' : 'genericQuad/Subsystem1/PID attitude/Ideal P Gain/Passthrough'
//  '<S397>' : 'genericQuad/Subsystem1/PID attitude/Ideal P Gain Fdbk/Disabled'
//  '<S398>' : 'genericQuad/Subsystem1/PID attitude/Integrator/Continuous'
//  '<S399>' : 'genericQuad/Subsystem1/PID attitude/Integrator ICs/Internal IC'
//  '<S400>' : 'genericQuad/Subsystem1/PID attitude/N Copy/Disabled'
//  '<S401>' : 'genericQuad/Subsystem1/PID attitude/N Gain/Internal Parameters'
//  '<S402>' : 'genericQuad/Subsystem1/PID attitude/P Copy/Disabled'
//  '<S403>' : 'genericQuad/Subsystem1/PID attitude/Parallel P Gain/Internal Parameters'
//  '<S404>' : 'genericQuad/Subsystem1/PID attitude/Reset Signal/Disabled'
//  '<S405>' : 'genericQuad/Subsystem1/PID attitude/Saturation/Passthrough'
//  '<S406>' : 'genericQuad/Subsystem1/PID attitude/Saturation Fdbk/Disabled'
//  '<S407>' : 'genericQuad/Subsystem1/PID attitude/Sum/Sum_PID'
//  '<S408>' : 'genericQuad/Subsystem1/PID attitude/Sum Fdbk/Disabled'
//  '<S409>' : 'genericQuad/Subsystem1/PID attitude/Tracking Mode/Disabled'
//  '<S410>' : 'genericQuad/Subsystem1/PID attitude/Tracking Mode Sum/Passthrough'
//  '<S411>' : 'genericQuad/Subsystem1/PID attitude/Tsamp - Integral/Passthrough'
//  '<S412>' : 'genericQuad/Subsystem1/PID attitude/Tsamp - Ngain/Passthrough'
//  '<S413>' : 'genericQuad/Subsystem1/PID attitude/postSat Signal/Forward_Path'
//  '<S414>' : 'genericQuad/Subsystem1/PID attitude/preSat Signal/Forward_Path'

#endif                                 // RTW_HEADER_genericQuad_h_

//
// File trailer for generated code.
//
// [EOF]
//
