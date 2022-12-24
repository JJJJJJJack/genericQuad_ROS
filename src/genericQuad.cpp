//
// File: genericQuad.cpp
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

// Block signals (default storage)
B_genericQuad_T genericQuad_B;

// Continuous states
X_genericQuad_T genericQuad_X;

// Block states (default storage)
DW_genericQuad_T genericQuad_DW;

// Periodic continuous states
PeriodicIndX_genericQuad_T genericQuad_PeriodicIndX;
PeriodicRngX_genericQuad_T genericQuad_PeriodicRngX;

// Real-time model
RT_MODEL_genericQuad_T genericQuad_M_ = RT_MODEL_genericQuad_T();
RT_MODEL_genericQuad_T *const genericQuad_M = &genericQuad_M_;

// Forward declaration for local functions
static void rt_mrdivide_U1d1x3_U2d_9vOrDY_k(const real_T u0[3], const real_T u1
  [9], real_T y[3]);

// State reduction function
void local_stateReduction(real_T* x, int_T* p, int_T n, real_T* r)
{
  int_T i, j;
  for (i = 0, j = 0; i < n; ++i, ++j) {
    int_T k = p[i];
    real_T lb = r[j++];
    real_T xk = x[k]-lb;
    real_T rk = r[j]-lb;
    int_T q = (int_T) floor(xk/rk);
    if (q) {
      x[k] = xk-q*rk+lb;
    }
  }
}

//
// This function updates continuous states using the ODE3 fixed-step
// solver algorithm
//
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  // Solver Matrices
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = static_cast<ODE3_IntgData *>(rtsiGetSolverData(si));
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 28;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  // Save the state values at time t in y, we'll use x as ynew.
  (void) memcpy(y, x,
                static_cast<uint_T>(nXc)*sizeof(real_T));

  // Assumes that rtsiSetT and ModelOutputs are up-to-date
  // f0 = f(t,y)
  rtsiSetdX(si, f0);
  genericQuad_derivatives();

  // f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*));
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  genericQuad_step();
  genericQuad_derivatives();

  // f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*));
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  genericQuad_step();
  genericQuad_derivatives();

  // tnew = t + hA(3);
  // ynew = y + f*hB(:,3);
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  local_stateReduction(rtsiGetContStates(si), rtsiGetPeriodicContStateIndices(si),
                       3,
                       rtsiGetPeriodicContStateRanges(si));
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

static void rt_mrdivide_U1d1x3_U2d_9vOrDY_k(const real_T u0[3], const real_T u1
  [9], real_T y[3])
{
  real_T a21;
  real_T maxval;
  int32_T r1;
  int32_T r2;
  int32_T r3;
  memcpy(&genericQuad_B.A[0], &u1[0], 9U * sizeof(real_T));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = fabs(u1[0]);
  a21 = fabs(u1[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (fabs(u1[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  genericQuad_B.A[r2] = u1[r2] / u1[r1];
  genericQuad_B.A[r3] /= genericQuad_B.A[r1];
  genericQuad_B.A[r2 + 3] -= genericQuad_B.A[r1 + 3] * genericQuad_B.A[r2];
  genericQuad_B.A[r3 + 3] -= genericQuad_B.A[r1 + 3] * genericQuad_B.A[r3];
  genericQuad_B.A[r2 + 6] -= genericQuad_B.A[r1 + 6] * genericQuad_B.A[r2];
  genericQuad_B.A[r3 + 6] -= genericQuad_B.A[r1 + 6] * genericQuad_B.A[r3];
  if (fabs(genericQuad_B.A[r3 + 3]) > fabs(genericQuad_B.A[r2 + 3])) {
    int32_T rtemp;
    rtemp = r2 + 1;
    r2 = r3;
    r3 = rtemp - 1;
  }

  genericQuad_B.A[r3 + 3] /= genericQuad_B.A[r2 + 3];
  genericQuad_B.A[r3 + 6] -= genericQuad_B.A[r3 + 3] * genericQuad_B.A[r2 + 6];
  y[r1] = u0[0] / genericQuad_B.A[r1];
  y[r2] = u0[1] - genericQuad_B.A[r1 + 3] * y[r1];
  y[r3] = u0[2] - genericQuad_B.A[r1 + 6] * y[r1];
  y[r2] /= genericQuad_B.A[r2 + 3];
  y[r3] -= genericQuad_B.A[r2 + 6] * y[r2];
  y[r3] /= genericQuad_B.A[r3 + 6];
  y[r2] -= genericQuad_B.A[r3 + 3] * y[r3];
  y[r1] -= y[r3] * genericQuad_B.A[r3];
  y[r1] -= y[r2] * genericQuad_B.A[r2];
}

// Model step function
void genericQuad_step(void)
{
  // local block i/o variables
  real_T rtb_FilterCoefficient;
  int32_T i;
  int32_T rtb_VectorConcatenate_tmp;
  int32_T rtb_VectorConcatenate_tmp_0;
  boolean_T tmp;
  if (rtmIsMajorTimeStep(genericQuad_M)) {
    // set solver stop time
    rtsiSetSolverStopTime(&genericQuad_M->solverInfo,
                          ((genericQuad_M->Timing.clockTick0+1)*
      genericQuad_M->Timing.stepSize0));
  }                                    // end MajorTimeStep

  // Update absolute time of base rate at minor time step
  if (rtmIsMinorTimeStep(genericQuad_M)) {
    genericQuad_M->Timing.t[0] = rtsiGetT(&genericQuad_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(genericQuad_M)) {
    for (i = 0; i < 3; i++) {
      int32_T rtb_VectorConcatenate_tmp_1;

      // Concatenate: '<S25>/Vector Concatenate' incorporates:
      //   Constant: '<S25>/Constant1'
      //   Constant: '<S25>/Constant2'
      //   Selector: '<S24>/Selector1'

      genericQuad_B.VectorConcatenate[6 * i] =
        genericQuad_P.uDOFEulerAngles2_inertia[3 * i];
      rtb_VectorConcatenate_tmp = 6 * i + 3;
      genericQuad_B.VectorConcatenate[rtb_VectorConcatenate_tmp] =
        genericQuad_P.Constant2_Value_n[3 * i];

      // Selector: '<S24>/Selector' incorporates:
      //   Concatenate: '<S25>/Vector Concatenate'
      //   Selector: '<S24>/Selector2'

      genericQuad_B.Selector_tmp = genericQuad_B.VectorConcatenate[6 * i];

      // Selector: '<S24>/Selector'
      genericQuad_B.Selector[3 * i] = genericQuad_B.Selector_tmp;

      // Selector: '<S24>/Selector1' incorporates:
      //   Concatenate: '<S25>/Vector Concatenate'

      genericQuad_B.Selector1[3 * i] =
        genericQuad_B.VectorConcatenate[rtb_VectorConcatenate_tmp];

      // Selector: '<S24>/Selector2'
      genericQuad_B.Selector2[3 * i] = genericQuad_B.Selector_tmp;

      // Concatenate: '<S25>/Vector Concatenate' incorporates:
      //   Constant: '<S25>/Constant1'
      //   Constant: '<S25>/Constant2'
      //   Selector: '<S24>/Selector'
      //   Selector: '<S24>/Selector1'
      //   Selector: '<S24>/Selector2'

      rtb_VectorConcatenate_tmp = 3 * i + 1;
      rtb_VectorConcatenate_tmp_0 = 6 * i + 1;
      genericQuad_B.VectorConcatenate[rtb_VectorConcatenate_tmp_0] =
        genericQuad_P.uDOFEulerAngles2_inertia[rtb_VectorConcatenate_tmp];
      rtb_VectorConcatenate_tmp_1 = 6 * i + 4;
      genericQuad_B.VectorConcatenate[rtb_VectorConcatenate_tmp_1] =
        genericQuad_P.Constant2_Value_n[rtb_VectorConcatenate_tmp];

      // Selector: '<S24>/Selector' incorporates:
      //   Concatenate: '<S25>/Vector Concatenate'
      //   Selector: '<S24>/Selector2'

      genericQuad_B.Selector_tmp =
        genericQuad_B.VectorConcatenate[rtb_VectorConcatenate_tmp_0];

      // Selector: '<S24>/Selector'
      genericQuad_B.Selector[rtb_VectorConcatenate_tmp] =
        genericQuad_B.Selector_tmp;

      // Selector: '<S24>/Selector1' incorporates:
      //   Concatenate: '<S25>/Vector Concatenate'

      genericQuad_B.Selector1[rtb_VectorConcatenate_tmp] =
        genericQuad_B.VectorConcatenate[rtb_VectorConcatenate_tmp_1];

      // Selector: '<S24>/Selector2'
      genericQuad_B.Selector2[rtb_VectorConcatenate_tmp] =
        genericQuad_B.Selector_tmp;

      // Concatenate: '<S25>/Vector Concatenate' incorporates:
      //   Constant: '<S25>/Constant1'
      //   Constant: '<S25>/Constant2'
      //   Selector: '<S24>/Selector'
      //   Selector: '<S24>/Selector1'
      //   Selector: '<S24>/Selector2'

      rtb_VectorConcatenate_tmp = 3 * i + 2;
      rtb_VectorConcatenate_tmp_0 = 6 * i + 2;
      genericQuad_B.VectorConcatenate[rtb_VectorConcatenate_tmp_0] =
        genericQuad_P.uDOFEulerAngles2_inertia[rtb_VectorConcatenate_tmp];
      rtb_VectorConcatenate_tmp_1 = 6 * i + 5;
      genericQuad_B.VectorConcatenate[rtb_VectorConcatenate_tmp_1] =
        genericQuad_P.Constant2_Value_n[rtb_VectorConcatenate_tmp];

      // Selector: '<S24>/Selector' incorporates:
      //   Concatenate: '<S25>/Vector Concatenate'
      //   Selector: '<S24>/Selector2'

      genericQuad_B.Selector_tmp =
        genericQuad_B.VectorConcatenate[rtb_VectorConcatenate_tmp_0];

      // Selector: '<S24>/Selector'
      genericQuad_B.Selector[rtb_VectorConcatenate_tmp] =
        genericQuad_B.Selector_tmp;

      // Selector: '<S24>/Selector1' incorporates:
      //   Concatenate: '<S25>/Vector Concatenate'

      genericQuad_B.Selector1[rtb_VectorConcatenate_tmp] =
        genericQuad_B.VectorConcatenate[rtb_VectorConcatenate_tmp_1];

      // Selector: '<S24>/Selector2'
      genericQuad_B.Selector2[rtb_VectorConcatenate_tmp] =
        genericQuad_B.Selector_tmp;
    }
  }

  // Trigonometry: '<S31>/sincos' incorporates:
  //   Integrator: '<S23>/phi theta psi'
  //   MATLAB Function: '<S4>/euler to body'
  //   SignalConversion generated from: '<S31>/sincos'

  genericQuad_B.rtb_Sum2_idx_0 = cos(genericQuad_X.phithetapsi_CSTATE[2]);
  genericQuad_B.rtb_Sum2_idx_2 = sin(genericQuad_X.phithetapsi_CSTATE[2]);
  genericQuad_B.rtb_Sum2_idx_1_tmp = cos(genericQuad_X.phithetapsi_CSTATE[1]);
  genericQuad_B.Selector_tmp = sin(genericQuad_X.phithetapsi_CSTATE[1]);
  genericQuad_B.Saturation7 = cos(genericQuad_X.phithetapsi_CSTATE[0]);
  genericQuad_B.thrustsetpoint = sin(genericQuad_X.phithetapsi_CSTATE[0]);

  // Fcn: '<S31>/Fcn11' incorporates:
  //   Concatenate: '<S33>/Vector Concatenate'
  //   Trigonometry: '<S31>/sincos'

  genericQuad_B.VectorConcatenate_m[0] = genericQuad_B.rtb_Sum2_idx_0 *
    genericQuad_B.rtb_Sum2_idx_1_tmp;

  // Fcn: '<S31>/Fcn21' incorporates:
  //   Concatenate: '<S33>/Vector Concatenate'
  //   Fcn: '<S31>/Fcn22'
  //   Trigonometry: '<S31>/sincos'

  genericQuad_B.rtb_Sum2_idx_1 = genericQuad_B.Selector_tmp *
    genericQuad_B.thrustsetpoint;
  genericQuad_B.VectorConcatenate_m[1] = genericQuad_B.rtb_Sum2_idx_1 *
    genericQuad_B.rtb_Sum2_idx_0 - genericQuad_B.rtb_Sum2_idx_2 *
    genericQuad_B.Saturation7;

  // Fcn: '<S31>/Fcn31' incorporates:
  //   Concatenate: '<S33>/Vector Concatenate'
  //   Fcn: '<S31>/Fcn32'
  //   Trigonometry: '<S31>/sincos'

  genericQuad_B.kxj = genericQuad_B.Selector_tmp * genericQuad_B.Saturation7;
  genericQuad_B.VectorConcatenate_m[2] = genericQuad_B.kxj *
    genericQuad_B.rtb_Sum2_idx_0 + genericQuad_B.rtb_Sum2_idx_2 *
    genericQuad_B.thrustsetpoint;

  // Fcn: '<S31>/Fcn12' incorporates:
  //   Concatenate: '<S33>/Vector Concatenate'
  //   Trigonometry: '<S31>/sincos'

  genericQuad_B.VectorConcatenate_m[3] = genericQuad_B.rtb_Sum2_idx_2 *
    genericQuad_B.rtb_Sum2_idx_1_tmp;

  // Fcn: '<S31>/Fcn22' incorporates:
  //   Concatenate: '<S33>/Vector Concatenate'
  //   Trigonometry: '<S31>/sincos'

  genericQuad_B.VectorConcatenate_m[4] = genericQuad_B.rtb_Sum2_idx_1 *
    genericQuad_B.rtb_Sum2_idx_2 + genericQuad_B.rtb_Sum2_idx_0 *
    genericQuad_B.Saturation7;

  // Fcn: '<S31>/Fcn32' incorporates:
  //   Concatenate: '<S33>/Vector Concatenate'
  //   Trigonometry: '<S31>/sincos'

  genericQuad_B.VectorConcatenate_m[5] = genericQuad_B.kxj *
    genericQuad_B.rtb_Sum2_idx_2 - genericQuad_B.rtb_Sum2_idx_0 *
    genericQuad_B.thrustsetpoint;

  // Fcn: '<S31>/Fcn13' incorporates:
  //   Concatenate: '<S33>/Vector Concatenate'
  //   Trigonometry: '<S31>/sincos'

  genericQuad_B.VectorConcatenate_m[6] = -genericQuad_B.Selector_tmp;

  // Fcn: '<S31>/Fcn23' incorporates:
  //   Concatenate: '<S33>/Vector Concatenate'
  //   MATLAB Function: '<S4>/euler to body'
  //   Trigonometry: '<S31>/sincos'

  genericQuad_B.VectorConcatenate_tmp = genericQuad_B.rtb_Sum2_idx_1_tmp *
    genericQuad_B.thrustsetpoint;
  genericQuad_B.VectorConcatenate_m[7] = genericQuad_B.VectorConcatenate_tmp;

  // Fcn: '<S31>/Fcn33' incorporates:
  //   Concatenate: '<S33>/Vector Concatenate'
  //   Trigonometry: '<S31>/sincos'

  genericQuad_B.VectorConcatenate_m[8] = genericQuad_B.rtb_Sum2_idx_1_tmp *
    genericQuad_B.Saturation7;

  // Outputs for IfAction SubSystem: '<S43>/If Warning//Error' incorporates:
  //   ActionPort: '<S67>/if'

  for (i = 0; i < 3; i++) {
    // If: '<S43>/If1' incorporates:
    //   Concatenate: '<S33>/Vector Concatenate'
    //   Math: '<S70>/Math Function'
    //   Math: '<S7>/Transpose'

    genericQuad_B.Product_tmp[3 * i] = genericQuad_B.VectorConcatenate_m[i];
    genericQuad_B.Product_tmp[3 * i + 1] = genericQuad_B.VectorConcatenate_m[i +
      3];
    genericQuad_B.Product_tmp[3 * i + 2] = genericQuad_B.VectorConcatenate_m[i +
      6];
  }

  // End of Outputs for SubSystem: '<S43>/If Warning//Error'
  for (i = 0; i < 3; i++) {
    // Product: '<S30>/Product' incorporates:
    //   Integrator: '<S7>/ub,vb,wb'
    //   Math: '<S7>/Transpose'

    genericQuad_B.Product[i] = 0.0;
    genericQuad_B.Product[i] += genericQuad_B.Product_tmp[i] *
      genericQuad_X.ubvbwb_CSTATE[0];
    genericQuad_B.Product[i] += genericQuad_B.Product_tmp[i + 3] *
      genericQuad_X.ubvbwb_CSTATE[1];
    genericQuad_B.Product[i] += genericQuad_B.Product_tmp[i + 6] *
      genericQuad_X.ubvbwb_CSTATE[2];
  }

  if (rtmIsMajorTimeStep(genericQuad_M)) {
    // Outputs for Atomic SubSystem: '<Root>/Subscribe1'
    // MATLABSystem: '<S3>/SourceBlock' incorporates:
    //   Inport: '<S6>/In1'

    tmp = Sub_genericQuad_426.getLatestMessage(&genericQuad_B.b_varargout_2);

    // Outputs for Enabled SubSystem: '<S3>/Enabled Subsystem' incorporates:
    //   EnablePort: '<S6>/Enable'

    if (tmp) {
      genericQuad_B.In1 = genericQuad_B.b_varargout_2;
    }

    // End of MATLABSystem: '<S3>/SourceBlock'
    // End of Outputs for SubSystem: '<S3>/Enabled Subsystem'
    // End of Outputs for SubSystem: '<Root>/Subscribe1'

    // Reshape: '<S9>/Reshape' incorporates:
    //   Constant: '<S9>/Constant4'
    //   Constant: '<S9>/Constant5'
    //   Constant: '<S9>/Mass'
    //   Gain: '<S9>/g'

    genericQuad_B.Reshape[0] = genericQuad_P.Constant4_Value;
    genericQuad_B.Reshape[1] = genericQuad_P.Constant5_Value;
    genericQuad_B.Reshape[2] = genericQuad_P.g_Gain * genericQuad_P.Mass_Value;
  }

  // Sum: '<S4>/Sum2' incorporates:
  //   Gain: '<S20>/Gain'

  genericQuad_B.rtb_Sum2_idx_0 = genericQuad_B.In1.Linear.X -
    genericQuad_B.Product[0];
  genericQuad_B.rtb_Sum2_idx_1 = genericQuad_B.In1.Linear.Y -
    genericQuad_B.Product[1];
  genericQuad_B.rtb_Sum2_idx_2 = genericQuad_B.In1.Linear.Z -
    genericQuad_P.Gain_Gain_ip * genericQuad_B.Product[2];

  // TransferFcn: '<S4>/Transfer Fcn2'
  genericQuad_B.kxj = genericQuad_P.TransferFcn2_C *
    genericQuad_X.TransferFcn2_CSTATE;

  // Saturate: '<S4>/Saturation'
  if (genericQuad_B.kxj > genericQuad_P.Saturation_UpperSat) {
    genericQuad_B.Saturation = genericQuad_P.Saturation_UpperSat;
  } else if (genericQuad_B.kxj < genericQuad_P.Saturation_LowerSat) {
    genericQuad_B.Saturation = genericQuad_P.Saturation_LowerSat;
  } else {
    genericQuad_B.Saturation = genericQuad_B.kxj;
  }

  // End of Saturate: '<S4>/Saturation'

  // TransferFcn: '<S4>/Transfer Fcn1'
  genericQuad_B.kxj = genericQuad_P.TransferFcn1_C *
    genericQuad_X.TransferFcn1_CSTATE;

  // Saturate: '<S4>/Saturation1'
  if (genericQuad_B.kxj > genericQuad_P.Saturation1_UpperSat) {
    genericQuad_B.Saturation1 = genericQuad_P.Saturation1_UpperSat;
  } else if (genericQuad_B.kxj < genericQuad_P.Saturation1_LowerSat) {
    genericQuad_B.Saturation1 = genericQuad_P.Saturation1_LowerSat;
  } else {
    genericQuad_B.Saturation1 = genericQuad_B.kxj;
  }

  // End of Saturate: '<S4>/Saturation1'

  // TransferFcn: '<S4>/Transfer Fcn3'
  genericQuad_B.kxj = genericQuad_P.TransferFcn3_C *
    genericQuad_X.TransferFcn3_CSTATE;

  // Saturate: '<S4>/Saturation2'
  if (genericQuad_B.kxj > genericQuad_P.Saturation2_UpperSat) {
    genericQuad_B.Saturation2 = genericQuad_P.Saturation2_UpperSat;
  } else if (genericQuad_B.kxj < genericQuad_P.Saturation2_LowerSat) {
    genericQuad_B.Saturation2 = genericQuad_P.Saturation2_LowerSat;
  } else {
    genericQuad_B.Saturation2 = genericQuad_B.kxj;
  }

  // End of Saturate: '<S4>/Saturation2'

  // TransferFcn: '<S4>/Transfer Fcn4'
  genericQuad_B.kxj = genericQuad_P.TransferFcn4_C *
    genericQuad_X.TransferFcn4_CSTATE;

  // Saturate: '<S4>/Saturation3'
  if (genericQuad_B.kxj > genericQuad_P.Saturation3_UpperSat) {
    genericQuad_B.kxj = genericQuad_P.Saturation3_UpperSat;
  } else if (genericQuad_B.kxj < genericQuad_P.Saturation3_LowerSat) {
    genericQuad_B.kxj = genericQuad_P.Saturation3_LowerSat;
  }

  // End of Saturate: '<S4>/Saturation3'

  // Sum: '<S4>/Sum15' incorporates:
  //   MATLAB Function: '<S4>/Mapping'

  genericQuad_B.Sum7[0] = 0.0;
  genericQuad_B.Sum7[1] = 0.0;
  genericQuad_B.Sum7[2] = -(((genericQuad_B.Saturation +
    genericQuad_B.Saturation1) + genericQuad_B.Saturation2) + genericQuad_B.kxj);

  // Sum: '<S26>/Sum' incorporates:
  //   Integrator: '<S7>/p,q,r '
  //   Integrator: '<S7>/ub,vb,wb'
  //   Product: '<S39>/i x j'
  //   Product: '<S39>/j x k'
  //   Product: '<S39>/k x i'
  //   Product: '<S40>/i x k'
  //   Product: '<S40>/j x i'
  //   Product: '<S40>/k x j'

  genericQuad_B.rtb_Saturation2_k[0] = genericQuad_X.ubvbwb_CSTATE[1] *
    genericQuad_X.pqr_CSTATE[2];
  genericQuad_B.rtb_Saturation2_k[1] = genericQuad_X.pqr_CSTATE[0] *
    genericQuad_X.ubvbwb_CSTATE[2];
  genericQuad_B.rtb_Saturation2_k[2] = genericQuad_X.ubvbwb_CSTATE[0] *
    genericQuad_X.pqr_CSTATE[1];
  genericQuad_B.rtb_Sum_d_c[0] = genericQuad_X.pqr_CSTATE[1] *
    genericQuad_X.ubvbwb_CSTATE[2];
  genericQuad_B.rtb_Sum_d_c[1] = genericQuad_X.ubvbwb_CSTATE[0] *
    genericQuad_X.pqr_CSTATE[2];
  genericQuad_B.rtb_Sum_d_c[2] = genericQuad_X.pqr_CSTATE[0] *
    genericQuad_X.ubvbwb_CSTATE[1];
  for (i = 0; i < 3; i++) {
    // Sum: '<S7>/Sum' incorporates:
    //   Concatenate: '<S33>/Vector Concatenate'
    //   Constant: '<S25>/Constant'
    //   Product: '<S7>/Product'
    //   Product: '<S9>/Product'
    //   Reshape: '<S9>/Reshape'
    //   Sum: '<S26>/Sum'
    //   Sum: '<S4>/Sum15'

    genericQuad_B.Sum[i] = (((genericQuad_B.VectorConcatenate_m[i + 3] *
      genericQuad_B.Reshape[1] + genericQuad_B.VectorConcatenate_m[i] *
      genericQuad_B.Reshape[0]) + genericQuad_B.VectorConcatenate_m[i + 6] *
      genericQuad_B.Reshape[2]) + genericQuad_B.Sum7[i]) /
      genericQuad_P.uDOFEulerAngles2_mass_0 + (genericQuad_B.rtb_Saturation2_k[i]
      - genericQuad_B.rtb_Sum_d_c[i]);

    // Product: '<S35>/Product' incorporates:
    //   Integrator: '<S7>/p,q,r '
    //   Selector: '<S24>/Selector'
    //   Sum: '<Root>/Sum'

    genericQuad_B.Sum_d[i] = (genericQuad_B.Selector[i + 3] *
      genericQuad_X.pqr_CSTATE[1] + genericQuad_B.Selector[i] *
      genericQuad_X.pqr_CSTATE[0]) + genericQuad_B.Selector[i + 6] *
      genericQuad_X.pqr_CSTATE[2];
  }

  // MATLAB Function: '<S4>/Mapping'
  genericQuad_B.rtb_Saturation2_k[0] = (((genericQuad_B.Saturation2 +
    genericQuad_B.kxj) - genericQuad_B.Saturation) - genericQuad_B.Saturation1) *
    0.15;
  genericQuad_B.rtb_Saturation2_k[1] = (((genericQuad_B.Saturation1 +
    genericQuad_B.kxj) - genericQuad_B.Saturation) - genericQuad_B.Saturation2) *
    0.15;
  genericQuad_B.rtb_Saturation2_k[2] = (((genericQuad_B.Saturation1 +
    genericQuad_B.Saturation2) - genericQuad_B.Saturation) - genericQuad_B.kxj) *
    0.2 * 0.2;

  // Sum: '<S34>/Sum' incorporates:
  //   Integrator: '<S7>/p,q,r '
  //   Product: '<S37>/i x j'
  //   Product: '<S37>/j x k'
  //   Product: '<S37>/k x i'
  //   Product: '<S38>/i x k'
  //   Product: '<S38>/j x i'
  //   Product: '<S38>/k x j'

  genericQuad_B.Sum7[0] = genericQuad_X.pqr_CSTATE[1] * genericQuad_B.Sum_d[2];
  genericQuad_B.Sum7[1] = genericQuad_B.Sum_d[0] * genericQuad_X.pqr_CSTATE[2];
  genericQuad_B.Sum7[2] = genericQuad_X.pqr_CSTATE[0] * genericQuad_B.Sum_d[1];
  genericQuad_B.rtb_Sum_d_c[0] = genericQuad_B.Sum_d[1] *
    genericQuad_X.pqr_CSTATE[2];
  genericQuad_B.rtb_Sum_d_c[1] = genericQuad_X.pqr_CSTATE[0] *
    genericQuad_B.Sum_d[2];
  genericQuad_B.rtb_Sum_d_c[2] = genericQuad_B.Sum_d[0] *
    genericQuad_X.pqr_CSTATE[1];
  for (i = 0; i < 3; i++) {
    // Sum: '<S24>/Sum2' incorporates:
    //   Integrator: '<S7>/p,q,r '
    //   Product: '<S36>/Product'
    //   Selector: '<S24>/Selector1'
    //   Sum: '<S34>/Sum'

    genericQuad_B.Sum_d[i] = (genericQuad_B.rtb_Saturation2_k[i] -
      ((genericQuad_B.Selector1[i + 3] * genericQuad_X.pqr_CSTATE[1] +
        genericQuad_B.Selector1[i] * genericQuad_X.pqr_CSTATE[0]) +
       genericQuad_B.Selector1[i + 6] * genericQuad_X.pqr_CSTATE[2])) -
      (genericQuad_B.Sum7[i] - genericQuad_B.rtb_Sum_d_c[i]);
  }

  // Product: '<S24>/Product2' incorporates:
  //   Selector: '<S24>/Selector2'

  rt_mrdivide_U1d1x3_U2d_9vOrDY_k(genericQuad_B.Sum_d, genericQuad_B.Selector2,
    genericQuad_B.Product2);
  if (rtmIsMajorTimeStep(genericQuad_M)) {
    // SignalConversion generated from: '<Root>/Bus Selector3'
    genericQuad_B.Z = genericQuad_B.In1.Angular.Z;
  }

  // Gain: '<S113>/Filter Coefficient' incorporates:
  //   Gain: '<S104>/Derivative Gain'
  //   Integrator: '<S105>/Filter'
  //   Sum: '<S105>/SumD'

  genericQuad_B.FilterCoefficient = (genericQuad_P.PIDVelocityx_D *
    genericQuad_B.rtb_Sum2_idx_0 - genericQuad_X.Filter_CSTATE) *
    genericQuad_P.PIDVelocityx_N;

  // Gain: '<S161>/Filter Coefficient' incorporates:
  //   Gain: '<S152>/Derivative Gain'
  //   Integrator: '<S153>/Filter'
  //   Sum: '<S153>/SumD'

  genericQuad_B.FilterCoefficient_m = (genericQuad_P.PIDVelocityy_D *
    genericQuad_B.rtb_Sum2_idx_1 - genericQuad_X.Filter_CSTATE_k) *
    genericQuad_P.PIDVelocityy_N;

  // Gain: '<S209>/Filter Coefficient' incorporates:
  //   Gain: '<S200>/Derivative Gain'
  //   Integrator: '<S201>/Filter'
  //   Sum: '<S201>/SumD'

  genericQuad_B.FilterCoefficient_j = (genericQuad_P.PIDVelocityz_D *
    genericQuad_B.rtb_Sum2_idx_2 - genericQuad_X.Filter_CSTATE_h) *
    genericQuad_P.PIDVelocityz_N;

  // SignalConversion generated from: '<S19>/ SFunction ' incorporates:
  //   Gain: '<S115>/Proportional Gain'
  //   Gain: '<S163>/Proportional Gain'
  //   Integrator: '<S110>/Integrator'
  //   Integrator: '<S158>/Integrator'
  //   MATLAB Function: '<S4>/R_EW'
  //   Sum: '<S119>/Sum'
  //   Sum: '<S167>/Sum'

  genericQuad_B.sincos_o1[0] = (genericQuad_P.PIDVelocityx_P *
    genericQuad_B.rtb_Sum2_idx_0 + genericQuad_X.Integrator_CSTATE) +
    genericQuad_B.FilterCoefficient;
  genericQuad_B.sincos_o1[1] = (genericQuad_P.PIDVelocityy_P *
    genericQuad_B.rtb_Sum2_idx_1 + genericQuad_X.Integrator_CSTATE_o) +
    genericQuad_B.FilterCoefficient_m;

  // MATLAB Function: '<S4>/R_EW' incorporates:
  //   Concatenate: '<S33>/Vector Concatenate'

  for (i = 0; i < 3; i++) {
    genericQuad_B.Sum_d[i] = ((genericQuad_B.VectorConcatenate_m[i + 3] *
      genericQuad_B.sincos_o1[1] + genericQuad_B.VectorConcatenate_m[i] *
      genericQuad_B.sincos_o1[0]) + genericQuad_B.VectorConcatenate_m[i + 6] *
      0.0) * 0.017453292519943295;
  }

  // Sum: '<S4>/Sum1' incorporates:
  //   Integrator: '<S23>/phi theta psi'
  //   MATLAB Function: '<S4>/R_EW'

  genericQuad_B.sincos_o2[0] = genericQuad_B.Sum_d[1] -
    genericQuad_X.phithetapsi_CSTATE[0];
  genericQuad_B.sincos_o2[1] = -genericQuad_B.Sum_d[0] -
    genericQuad_X.phithetapsi_CSTATE[1];
  genericQuad_B.sincos_o2[2] = genericQuad_B.Z -
    genericQuad_X.phithetapsi_CSTATE[2];

  // MATLAB Function: '<S4>/euler to body' incorporates:
  //   Integrator: '<S23>/phi theta psi'

  genericQuad_B.Product_tmp_c[0] = 1.0;
  genericQuad_B.Product_tmp_c[3] = 0.0;
  genericQuad_B.Product_tmp_c[6] = -genericQuad_B.Selector_tmp;
  genericQuad_B.Product_tmp_c[1] = 0.0;
  genericQuad_B.Product_tmp_c[4] = genericQuad_B.Saturation7;
  genericQuad_B.Product_tmp_c[7] = genericQuad_B.VectorConcatenate_tmp;
  genericQuad_B.Product_tmp_c[2] = 0.0;
  genericQuad_B.Product_tmp_c[5] = -genericQuad_B.thrustsetpoint;
  genericQuad_B.Product_tmp_c[8] = cos(genericQuad_X.phithetapsi_CSTATE[0] *
    genericQuad_B.rtb_Sum2_idx_1_tmp);

  // Gain: '<S401>/Filter Coefficient' incorporates:
  //   Gain: '<S392>/Derivative Gain'
  //   Integrator: '<S393>/Filter'
  //   Sum: '<S393>/SumD'

  genericQuad_B.FilterCoefficient_p[0] = (genericQuad_P.PIDattitude_D *
    genericQuad_B.sincos_o2[0] - genericQuad_X.Filter_CSTATE_e[0]) *
    genericQuad_P.PIDattitude_N;

  // Sum: '<S407>/Sum' incorporates:
  //   Gain: '<S392>/Derivative Gain'
  //   Gain: '<S403>/Proportional Gain'
  //   Integrator: '<S398>/Integrator'

  genericQuad_B.rtb_Sum2_idx_1_tmp = (genericQuad_P.PIDattitude_P *
    genericQuad_B.sincos_o2[0] + genericQuad_X.Integrator_CSTATE_m[0]) +
    genericQuad_B.FilterCoefficient_p[0];

  // Gain: '<S401>/Filter Coefficient' incorporates:
  //   Gain: '<S392>/Derivative Gain'
  //   Integrator: '<S393>/Filter'
  //   Sum: '<S393>/SumD'

  genericQuad_B.FilterCoefficient_p[1] = (genericQuad_P.PIDattitude_D *
    genericQuad_B.sincos_o2[1] - genericQuad_X.Filter_CSTATE_e[1]) *
    genericQuad_P.PIDattitude_N;

  // Sum: '<S407>/Sum' incorporates:
  //   Gain: '<S392>/Derivative Gain'
  //   Gain: '<S403>/Proportional Gain'
  //   Integrator: '<S398>/Integrator'

  genericQuad_B.Saturation7 = (genericQuad_P.PIDattitude_P *
    genericQuad_B.sincos_o2[1] + genericQuad_X.Integrator_CSTATE_m[1]) +
    genericQuad_B.FilterCoefficient_p[1];

  // Gain: '<S401>/Filter Coefficient' incorporates:
  //   Gain: '<S392>/Derivative Gain'
  //   Integrator: '<S393>/Filter'
  //   Sum: '<S393>/SumD'

  genericQuad_B.FilterCoefficient_p[2] = (genericQuad_P.PIDattitude_D *
    genericQuad_B.sincos_o2[2] - genericQuad_X.Filter_CSTATE_e[2]) *
    genericQuad_P.PIDattitude_N;

  // Sum: '<S407>/Sum' incorporates:
  //   Gain: '<S392>/Derivative Gain'
  //   Gain: '<S403>/Proportional Gain'
  //   Integrator: '<S398>/Integrator'

  genericQuad_B.thrustsetpoint = (genericQuad_P.PIDattitude_P *
    genericQuad_B.sincos_o2[2] + genericQuad_X.Integrator_CSTATE_m[2]) +
    genericQuad_B.FilterCoefficient_p[2];
  for (i = 0; i < 3; i++) {
    // Sum: '<S4>/Sum7' incorporates:
    //   Integrator: '<S7>/p,q,r '
    //   MATLAB Function: '<S4>/euler to body'

    genericQuad_B.Sum7[i] = ((genericQuad_B.Product_tmp_c[i + 3] *
      genericQuad_B.Saturation7 + genericQuad_B.Product_tmp_c[i] *
      genericQuad_B.rtb_Sum2_idx_1_tmp) + genericQuad_B.Product_tmp_c[i + 6] *
      genericQuad_B.thrustsetpoint) - genericQuad_X.pqr_CSTATE[i];

    // Gain: '<S395>/Integral Gain'
    genericQuad_B.IntegralGain[i] = genericQuad_P.PIDattitude_I *
      genericQuad_B.sincos_o2[i];

    // Trigonometry: '<S32>/sincos' incorporates:
    //   Integrator: '<S23>/phi theta psi'

    genericQuad_B.sincos_o1[i] = sin(genericQuad_X.phithetapsi_CSTATE[i]);
    genericQuad_B.sincos_o2[i] = cos(genericQuad_X.phithetapsi_CSTATE[i]);
  }

  // MATLABSystem: '<Root>/Get Parameter'
  ParamGet_genericQuad_459.get_parameter(&genericQuad_B.rtb_Sum2_idx_1_tmp);

  // MATLABSystem: '<Root>/Get Parameter1'
  ParamGet_genericQuad_460.get_parameter(&genericQuad_B.Saturation7);

  // MATLABSystem: '<Root>/Get Parameter2'
  ParamGet_genericQuad_461.get_parameter(&genericQuad_B.thrustsetpoint);

  // Fcn: '<S32>/phidot' incorporates:
  //   Fcn: '<S32>/psidot'
  //   Integrator: '<S7>/p,q,r '

  genericQuad_B.kxj = genericQuad_B.sincos_o1[0] * genericQuad_X.pqr_CSTATE[1] +
    genericQuad_B.sincos_o2[0] * genericQuad_X.pqr_CSTATE[2];

  // SignalConversion generated from: '<S23>/phi theta psi' incorporates:
  //   Fcn: '<S32>/phidot'
  //   Fcn: '<S32>/psidot'
  //   Fcn: '<S32>/thetadot'
  //   Integrator: '<S7>/p,q,r '

  genericQuad_B.TmpSignalConversionAtphithetaps[0] = genericQuad_B.sincos_o1[1] /
    genericQuad_B.sincos_o2[1] * genericQuad_B.kxj + genericQuad_X.pqr_CSTATE[0];
  genericQuad_B.TmpSignalConversionAtphithetaps[1] = genericQuad_B.sincos_o2[0] *
    genericQuad_X.pqr_CSTATE[1] - genericQuad_B.sincos_o1[0] *
    genericQuad_X.pqr_CSTATE[2];
  genericQuad_B.TmpSignalConversionAtphithetaps[2] = genericQuad_B.kxj /
    genericQuad_B.sincos_o2[1];

  // Sum: '<S44>/Add'
  genericQuad_B.kxj = (genericQuad_B.VectorConcatenate_m[0] +
                       genericQuad_B.VectorConcatenate_m[4]) +
    genericQuad_B.VectorConcatenate_m[8];

  // If: '<S8>/If' incorporates:
  //   Sum: '<S44>/Add'

  if (rtmIsMajorTimeStep(genericQuad_M)) {
    genericQuad_DW.If_ActiveSubsystem = static_cast<int8_T>(!(genericQuad_B.kxj >
      0.0));
  }

  switch (genericQuad_DW.If_ActiveSubsystem) {
   case 0:
    // Outputs for IfAction SubSystem: '<S8>/Positive Trace' incorporates:
    //   ActionPort: '<S42>/Action Port'

    // Sqrt: '<S42>/sqrt' incorporates:
    //   Constant: '<S42>/Constant'
    //   Sum: '<S42>/Sum'
    //   Sum: '<S44>/Add'

    genericQuad_B.kxj = sqrt(genericQuad_B.kxj + genericQuad_P.Constant_Value_o);

    // Gain: '<S42>/Gain' incorporates:
    //   Merge: '<S8>/Merge'

    genericQuad_B.Merge[0] = genericQuad_P.Gain_Gain * genericQuad_B.kxj;

    // Gain: '<S42>/Gain1'
    genericQuad_B.kxj *= genericQuad_P.Gain1_Gain;

    // Product: '<S42>/Product' incorporates:
    //   Fcn: '<S31>/Fcn13'
    //   Fcn: '<S31>/Fcn23'
    //   Merge: '<S8>/Merge'
    //   Sum: '<S64>/Add'
    //   Sum: '<S65>/Add'
    //   Sum: '<S66>/Add'
    //   Trigonometry: '<S31>/sincos'

    genericQuad_B.Merge[1] = (genericQuad_B.VectorConcatenate_tmp -
      genericQuad_B.VectorConcatenate_m[5]) / genericQuad_B.kxj;
    genericQuad_B.Merge[2] = (genericQuad_B.VectorConcatenate_m[2] -
      (-genericQuad_B.Selector_tmp)) / genericQuad_B.kxj;
    genericQuad_B.Merge[3] = (genericQuad_B.VectorConcatenate_m[3] -
      genericQuad_B.VectorConcatenate_m[1]) / genericQuad_B.kxj;

    // End of Outputs for SubSystem: '<S8>/Positive Trace'
    break;

   case 1:
    // Outputs for IfAction SubSystem: '<S8>/Negative Trace' incorporates:
    //   ActionPort: '<S41>/Action Port'

    // If: '<S41>/Find Maximum Diagonal Value'
    if (rtmIsMajorTimeStep(genericQuad_M)) {
      if ((genericQuad_B.VectorConcatenate_m[4] >
           genericQuad_B.VectorConcatenate_m[0]) &&
          (genericQuad_B.VectorConcatenate_m[4] >
           genericQuad_B.VectorConcatenate_m[8])) {
        genericQuad_DW.FindMaximumDiagonalValue_Active = 0;
      } else if (genericQuad_B.VectorConcatenate_m[8] >
                 genericQuad_B.VectorConcatenate_m[0]) {
        genericQuad_DW.FindMaximumDiagonalValue_Active = 1;
      } else {
        genericQuad_DW.FindMaximumDiagonalValue_Active = 2;
      }
    }

    switch (genericQuad_DW.FindMaximumDiagonalValue_Active) {
     case 0:
      // Outputs for IfAction SubSystem: '<S41>/Maximum Value at DCM(2,2)' incorporates:
      //   ActionPort: '<S46>/Action Port'

      // Sqrt: '<S46>/sqrt' incorporates:
      //   Constant: '<S58>/Constant'
      //   Sum: '<S58>/Add'

      genericQuad_B.Saturation1 = sqrt(((genericQuad_B.VectorConcatenate_m[4] -
        genericQuad_B.VectorConcatenate_m[0]) -
        genericQuad_B.VectorConcatenate_m[8]) + genericQuad_P.Constant_Value_n);

      // Switch: '<S57>/Switch' incorporates:
      //   Constant: '<S57>/Constant1'
      //   Constant: '<S57>/Constant2'

      if (genericQuad_B.Saturation1 != 0.0) {
        genericQuad_B.kxj = genericQuad_P.Constant1_Value;
        genericQuad_B.Saturation = genericQuad_B.Saturation1;
      } else {
        genericQuad_B.kxj = genericQuad_P.Constant2_Value[0];
        genericQuad_B.Saturation = genericQuad_P.Constant2_Value[1];
      }

      // End of Switch: '<S57>/Switch'

      // Product: '<S57>/Product'
      genericQuad_B.kxj /= genericQuad_B.Saturation;

      // Gain: '<S46>/Gain1' incorporates:
      //   Merge: '<S8>/Merge'
      //   Product: '<S46>/Product'
      //   Sum: '<S56>/Add'

      genericQuad_B.Merge[1] = (genericQuad_B.VectorConcatenate_m[1] +
        genericQuad_B.VectorConcatenate_m[3]) * genericQuad_B.kxj *
        genericQuad_P.Gain1_Gain_k;

      // Gain: '<S46>/Gain3' incorporates:
      //   Fcn: '<S31>/Fcn23'
      //   Merge: '<S8>/Merge'
      //   Product: '<S46>/Product'
      //   Sum: '<S55>/Add'

      genericQuad_B.Merge[3] = (genericQuad_B.VectorConcatenate_m[5] +
        genericQuad_B.VectorConcatenate_tmp) * genericQuad_B.kxj *
        genericQuad_P.Gain3_Gain;

      // Gain: '<S46>/Gain4' incorporates:
      //   Fcn: '<S31>/Fcn13'
      //   Merge: '<S8>/Merge'
      //   Product: '<S46>/Product'
      //   Sum: '<S54>/Add'
      //   Trigonometry: '<S31>/sincos'

      genericQuad_B.Merge[0] = (genericQuad_B.VectorConcatenate_m[2] -
        (-genericQuad_B.Selector_tmp)) * genericQuad_B.kxj *
        genericQuad_P.Gain4_Gain;

      // Gain: '<S46>/Gain' incorporates:
      //   Merge: '<S8>/Merge'

      genericQuad_B.Merge[2] = genericQuad_P.Gain_Gain_g *
        genericQuad_B.Saturation1;

      // End of Outputs for SubSystem: '<S41>/Maximum Value at DCM(2,2)'
      break;

     case 1:
      // Outputs for IfAction SubSystem: '<S41>/Maximum Value at DCM(3,3)' incorporates:
      //   ActionPort: '<S47>/Action Port'

      // Sqrt: '<S47>/sqrt' incorporates:
      //   Constant: '<S63>/Constant'
      //   Sum: '<S63>/Add'

      genericQuad_B.Saturation1 = sqrt(((genericQuad_B.VectorConcatenate_m[8] -
        genericQuad_B.VectorConcatenate_m[0]) -
        genericQuad_B.VectorConcatenate_m[4]) + genericQuad_P.Constant_Value_p);

      // Switch: '<S62>/Switch' incorporates:
      //   Constant: '<S62>/Constant1'
      //   Constant: '<S62>/Constant2'

      if (genericQuad_B.Saturation1 != 0.0) {
        genericQuad_B.kxj = genericQuad_P.Constant1_Value_h;
        genericQuad_B.Saturation = genericQuad_B.Saturation1;
      } else {
        genericQuad_B.kxj = genericQuad_P.Constant2_Value_k[0];
        genericQuad_B.Saturation = genericQuad_P.Constant2_Value_k[1];
      }

      // End of Switch: '<S62>/Switch'

      // Product: '<S62>/Product'
      genericQuad_B.kxj /= genericQuad_B.Saturation;

      // Gain: '<S47>/Gain1' incorporates:
      //   Fcn: '<S31>/Fcn13'
      //   Merge: '<S8>/Merge'
      //   Product: '<S47>/Product'
      //   Sum: '<S59>/Add'
      //   Trigonometry: '<S31>/sincos'

      genericQuad_B.Merge[1] = (genericQuad_B.VectorConcatenate_m[2] +
        -genericQuad_B.Selector_tmp) * genericQuad_B.kxj *
        genericQuad_P.Gain1_Gain_l;

      // Gain: '<S47>/Gain2' incorporates:
      //   Fcn: '<S31>/Fcn23'
      //   Merge: '<S8>/Merge'
      //   Product: '<S47>/Product'
      //   Sum: '<S60>/Add'

      genericQuad_B.Merge[2] = (genericQuad_B.VectorConcatenate_m[5] +
        genericQuad_B.VectorConcatenate_tmp) * genericQuad_B.kxj *
        genericQuad_P.Gain2_Gain;

      // Gain: '<S47>/Gain3' incorporates:
      //   Merge: '<S8>/Merge'
      //   Product: '<S47>/Product'
      //   Sum: '<S61>/Add'

      genericQuad_B.Merge[0] = (genericQuad_B.VectorConcatenate_m[3] -
        genericQuad_B.VectorConcatenate_m[1]) * genericQuad_B.kxj *
        genericQuad_P.Gain3_Gain_a;

      // Gain: '<S47>/Gain' incorporates:
      //   Merge: '<S8>/Merge'

      genericQuad_B.Merge[3] = genericQuad_P.Gain_Gain_d *
        genericQuad_B.Saturation1;

      // End of Outputs for SubSystem: '<S41>/Maximum Value at DCM(3,3)'
      break;

     case 2:
      // Outputs for IfAction SubSystem: '<S41>/Maximum Value at DCM(1,1)' incorporates:
      //   ActionPort: '<S45>/Action Port'

      // Sqrt: '<S45>/sqrt' incorporates:
      //   Constant: '<S53>/Constant'
      //   Sum: '<S53>/Add'

      genericQuad_B.Saturation1 = sqrt(((genericQuad_B.VectorConcatenate_m[0] -
        genericQuad_B.VectorConcatenate_m[4]) -
        genericQuad_B.VectorConcatenate_m[8]) + genericQuad_P.Constant_Value_fv);

      // Switch: '<S52>/Switch' incorporates:
      //   Constant: '<S52>/Constant1'
      //   Constant: '<S52>/Constant2'

      if (genericQuad_B.Saturation1 != 0.0) {
        genericQuad_B.kxj = genericQuad_P.Constant1_Value_b;
        genericQuad_B.Saturation = genericQuad_B.Saturation1;
      } else {
        genericQuad_B.kxj = genericQuad_P.Constant2_Value_o[0];
        genericQuad_B.Saturation = genericQuad_P.Constant2_Value_o[1];
      }

      // End of Switch: '<S52>/Switch'

      // Product: '<S52>/Product'
      genericQuad_B.kxj /= genericQuad_B.Saturation;

      // Gain: '<S45>/Gain1' incorporates:
      //   Merge: '<S8>/Merge'
      //   Product: '<S45>/Product'
      //   Sum: '<S51>/Add'

      genericQuad_B.Merge[2] = (genericQuad_B.VectorConcatenate_m[1] +
        genericQuad_B.VectorConcatenate_m[3]) * genericQuad_B.kxj *
        genericQuad_P.Gain1_Gain_b;

      // Gain: '<S45>/Gain2' incorporates:
      //   Fcn: '<S31>/Fcn13'
      //   Merge: '<S8>/Merge'
      //   Product: '<S45>/Product'
      //   Sum: '<S49>/Add'
      //   Trigonometry: '<S31>/sincos'

      genericQuad_B.Merge[3] = (genericQuad_B.VectorConcatenate_m[2] +
        -genericQuad_B.Selector_tmp) * genericQuad_B.kxj *
        genericQuad_P.Gain2_Gain_d;

      // Gain: '<S45>/Gain3' incorporates:
      //   Fcn: '<S31>/Fcn23'
      //   Merge: '<S8>/Merge'
      //   Product: '<S45>/Product'
      //   Sum: '<S50>/Add'

      genericQuad_B.Merge[0] = (genericQuad_B.VectorConcatenate_tmp -
        genericQuad_B.VectorConcatenate_m[5]) * genericQuad_B.kxj *
        genericQuad_P.Gain3_Gain_az;

      // Gain: '<S45>/Gain' incorporates:
      //   Merge: '<S8>/Merge'

      genericQuad_B.Merge[1] = genericQuad_P.Gain_Gain_i *
        genericQuad_B.Saturation1;

      // End of Outputs for SubSystem: '<S41>/Maximum Value at DCM(1,1)'
      break;
    }

    // End of If: '<S41>/Find Maximum Diagonal Value'
    // End of Outputs for SubSystem: '<S8>/Negative Trace'
    break;
  }

  // End of If: '<S8>/If'

  // Clock: '<Root>/Clock'
  genericQuad_B.kxj = genericQuad_M->Timing.t[0];

  // MATLAB Function: '<Root>/time to sec & nsec'
  genericQuad_B.Saturation = floor(genericQuad_B.kxj);

  // BusAssignment: '<Root>/Bus Assignment' incorporates:
  //   Constant: '<S1>/Constant'
  //   Gain: '<S21>/Gain'
  //   Integrator: '<S7>/xe,ye,ze'
  //   MATLAB Function: '<Root>/time to sec & nsec'
  //   MATLABSystem: '<Root>/Get Parameter'
  //   MATLABSystem: '<Root>/Get Parameter1'
  //   MATLABSystem: '<Root>/Get Parameter2'
  //   Sum: '<Root>/Sum'

  genericQuad_B.BusAssignment = genericQuad_P.Constant_Value;
  genericQuad_B.BusAssignment.Header.Stamp.Sec = genericQuad_B.Saturation;
  genericQuad_B.BusAssignment.Header.Stamp.Nsec = (genericQuad_B.kxj -
    genericQuad_B.Saturation) * 1.0E+9;
  genericQuad_B.BusAssignment.Pose.Position.X = genericQuad_B.rtb_Sum2_idx_1_tmp
    + genericQuad_X.xeyeze_CSTATE[0];
  genericQuad_B.BusAssignment.Pose.Position.Y = genericQuad_B.Saturation7 +
    genericQuad_X.xeyeze_CSTATE[1];
  genericQuad_B.BusAssignment.Pose.Position.Z = genericQuad_P.Gain_Gain_m *
    genericQuad_X.xeyeze_CSTATE[2] + genericQuad_B.thrustsetpoint;
  genericQuad_B.BusAssignment.Pose.Orientation.W = genericQuad_B.Merge[0];
  genericQuad_B.BusAssignment.Pose.Orientation.X = genericQuad_B.Merge[1];
  genericQuad_B.BusAssignment.Pose.Orientation.Y = genericQuad_B.Merge[2];
  genericQuad_B.BusAssignment.Pose.Orientation.Z = genericQuad_B.Merge[3];

  // Outputs for Atomic SubSystem: '<Root>/Publish'
  // MATLABSystem: '<S2>/SinkBlock'
  Pub_genericQuad_438.publish(&genericQuad_B.BusAssignment);

  // End of Outputs for SubSystem: '<Root>/Publish'
  if (rtmIsMajorTimeStep(genericQuad_M)) {
    // Gain: '<S353>/Filter Coefficient' incorporates:
    //   DiscreteIntegrator: '<S345>/Filter'
    //   Gain: '<S344>/Derivative Gain'
    //   Sum: '<S345>/SumD'

    rtb_FilterCoefficient = (genericQuad_P.PIDangularroll_D *
      genericQuad_B.Sum7[0] - genericQuad_DW.Filter_DSTATE) *
      genericQuad_P.PIDangularroll_N;

    // Gain: '<S348>/Proportional Gain' incorporates:
    //   DiscreteIntegrator: '<S350>/Integrator'
    //   Sum: '<S359>/Sum'

    genericQuad_B.ProportionalGain = ((genericQuad_B.Sum7[0] +
      genericQuad_DW.Integrator_DSTATE) + rtb_FilterCoefficient) *
      genericQuad_P.PIDangularroll_P;

    // Gain: '<S305>/Filter Coefficient' incorporates:
    //   DiscreteIntegrator: '<S297>/Filter'
    //   Gain: '<S296>/Derivative Gain'
    //   Sum: '<S297>/SumD'

    genericQuad_B.FilterCoefficient_e = (genericQuad_P.PIDangularpitch_D *
      genericQuad_B.Sum7[1] - genericQuad_DW.Filter_DSTATE_g) *
      genericQuad_P.PIDangularpitch_N;

    // Gain: '<S300>/Proportional Gain' incorporates:
    //   DiscreteIntegrator: '<S302>/Integrator'
    //   Sum: '<S311>/Sum'

    genericQuad_B.ProportionalGain_a = ((genericQuad_B.Sum7[1] +
      genericQuad_DW.Integrator_DSTATE_o) + genericQuad_B.FilterCoefficient_e) *
      genericQuad_P.PIDangularpitch_P;

    // Gain: '<S257>/Filter Coefficient' incorporates:
    //   DiscreteIntegrator: '<S249>/Filter'
    //   Gain: '<S248>/Derivative Gain'
    //   Sum: '<S249>/SumD'

    genericQuad_B.FilterCoefficient_n = (genericQuad_P.PIDangulayaw_D *
      genericQuad_B.Sum7[2] - genericQuad_DW.Filter_DSTATE_f) *
      genericQuad_P.PIDangulayaw_N;

    // Gain: '<S252>/Proportional Gain' incorporates:
    //   DiscreteIntegrator: '<S254>/Integrator'
    //   Sum: '<S263>/Sum'

    genericQuad_B.ProportionalGain_p = ((genericQuad_B.Sum7[2] +
      genericQuad_DW.Integrator_DSTATE_i) + genericQuad_B.FilterCoefficient_n) *
      genericQuad_P.PIDangulayaw_P;
  }

  // Sum: '<S4>/Sum3' incorporates:
  //   Constant: '<S4>/Hover throttle'
  //   Gain: '<S211>/Proportional Gain'
  //   Integrator: '<S206>/Integrator'
  //   Sum: '<S215>/Sum'

  genericQuad_B.thrustsetpoint = ((genericQuad_P.PIDVelocityz_P *
    genericQuad_B.rtb_Sum2_idx_2 + genericQuad_X.Integrator_CSTATE_h) +
    genericQuad_B.FilterCoefficient_j) + genericQuad_P.Hoverthrottle_Value;

  // MATLAB Function: '<S4>/Mixer'
  genericQuad_B.rtb_Sum2_idx_1_tmp = genericQuad_B.thrustsetpoint
    - genericQuad_B.ProportionalGain;
  genericQuad_B.Saturation7 = (genericQuad_B.rtb_Sum2_idx_1_tmp -
    genericQuad_B.ProportionalGain_a) - genericQuad_B.ProportionalGain_p;

  // Saturate: '<S4>/Saturation4'
  if (genericQuad_B.Saturation7 > genericQuad_P.Saturation4_UpperSat) {
    genericQuad_B.Saturation7 = genericQuad_P.Saturation4_UpperSat;
  } else if (genericQuad_B.Saturation7 < genericQuad_P.Saturation4_LowerSat) {
    genericQuad_B.Saturation7 = genericQuad_P.Saturation4_LowerSat;
  }

  // End of Saturate: '<S4>/Saturation4'

  // Gain: '<S4>/Komega' incorporates:
  //   Math: '<S4>/Square'

  genericQuad_B.Komega = genericQuad_B.Saturation7 * genericQuad_B.Saturation7 *
    genericQuad_P.Komega_Gain;

  // MATLAB Function: '<S4>/Mixer'
  genericQuad_B.Saturation7 = (genericQuad_B.rtb_Sum2_idx_1_tmp +
    genericQuad_B.ProportionalGain_a) + genericQuad_B.ProportionalGain_p;

  // Saturate: '<S4>/Saturation5'
  if (genericQuad_B.Saturation7 > genericQuad_P.Saturation5_UpperSat) {
    genericQuad_B.Saturation7 = genericQuad_P.Saturation5_UpperSat;
  } else if (genericQuad_B.Saturation7 < genericQuad_P.Saturation5_LowerSat) {
    genericQuad_B.Saturation7 = genericQuad_P.Saturation5_LowerSat;
  }

  // End of Saturate: '<S4>/Saturation5'

  // Gain: '<S4>/Komega1' incorporates:
  //   Math: '<S4>/Square1'

  genericQuad_B.Komega1 = genericQuad_B.Saturation7 * genericQuad_B.Saturation7 *
    genericQuad_P.Komega1_Gain;

  // MATLAB Function: '<S4>/Mixer'
  genericQuad_B.rtb_Sum2_idx_1_tmp = genericQuad_B.thrustsetpoint
    + genericQuad_B.ProportionalGain;
  genericQuad_B.Saturation7 = (genericQuad_B.rtb_Sum2_idx_1_tmp -
    genericQuad_B.ProportionalGain_a) + genericQuad_B.ProportionalGain_p;

  // Saturate: '<S4>/Saturation6'
  if (genericQuad_B.Saturation7 > genericQuad_P.Saturation6_UpperSat) {
    genericQuad_B.Saturation7 = genericQuad_P.Saturation6_UpperSat;
  } else if (genericQuad_B.Saturation7 < genericQuad_P.Saturation6_LowerSat) {
    genericQuad_B.Saturation7 = genericQuad_P.Saturation6_LowerSat;
  }

  // End of Saturate: '<S4>/Saturation6'

  // Gain: '<S4>/Komega2' incorporates:
  //   Math: '<S4>/Square2'

  genericQuad_B.Komega2 = genericQuad_B.Saturation7 * genericQuad_B.Saturation7 *
    genericQuad_P.Komega2_Gain;

  // MATLAB Function: '<S4>/Mixer'
  genericQuad_B.Saturation7 = (genericQuad_B.rtb_Sum2_idx_1_tmp +
    genericQuad_B.ProportionalGain_a) - genericQuad_B.ProportionalGain_p;

  // Saturate: '<S4>/Saturation7'
  if (genericQuad_B.Saturation7 > genericQuad_P.Saturation7_UpperSat) {
    genericQuad_B.Saturation7 = genericQuad_P.Saturation7_UpperSat;
  } else if (genericQuad_B.Saturation7 < genericQuad_P.Saturation7_LowerSat) {
    genericQuad_B.Saturation7 = genericQuad_P.Saturation7_LowerSat;
  }

  // End of Saturate: '<S4>/Saturation7'

  // Gain: '<S4>/Komega3' incorporates:
  //   Math: '<S4>/Square3'

  genericQuad_B.Komega3 = genericQuad_B.Saturation7 * genericQuad_B.Saturation7 *
    genericQuad_P.Komega3_Gain;

  // Gain: '<S107>/Integral Gain'
  genericQuad_B.IntegralGain_p = genericQuad_P.PIDVelocityx_I *
    genericQuad_B.rtb_Sum2_idx_0;

  // Gain: '<S155>/Integral Gain'
  genericQuad_B.IntegralGain_i = genericQuad_P.PIDVelocityy_I *
    genericQuad_B.rtb_Sum2_idx_1;

  // Gain: '<S203>/Integral Gain'
  genericQuad_B.IntegralGain_g = genericQuad_P.PIDVelocityz_I *
    genericQuad_B.rtb_Sum2_idx_2;
  if (rtmIsMajorTimeStep(genericQuad_M)) {
    int8_T rtAction;

    // Gain: '<S251>/Integral Gain'
    genericQuad_B.IntegralGain_pb = genericQuad_P.PIDangulayaw_I *
      genericQuad_B.Sum7[2];

    // Gain: '<S299>/Integral Gain'
    genericQuad_B.IntegralGain_c = genericQuad_P.PIDangularpitch_I *
      genericQuad_B.Sum7[1];

    // Gain: '<S347>/Integral Gain'
    genericQuad_B.IntegralGain_g_b = genericQuad_P.PIDangularroll_I *
      genericQuad_B.Sum7[0];

    // If: '<S43>/If1' incorporates:
    //   Constant: '<S43>/Constant'

    rtAction = -1;
    if (rtmIsMajorTimeStep(genericQuad_M)) {
      if (genericQuad_P.DirectionCosineMatrixtoQuaterni != 1.0) {
        rtAction = 0;
      }

      genericQuad_DW.If1_ActiveSubsystem = rtAction;
    } else {
      rtAction = genericQuad_DW.If1_ActiveSubsystem;
    }

    if (rtAction == 0) {
      // Outputs for IfAction SubSystem: '<S43>/If Warning//Error' incorporates:
      //   ActionPort: '<S67>/if'

      // Bias: '<S70>/Bias1' incorporates:
      //   Concatenate: '<S33>/Vector Concatenate'
      //   Math: '<S70>/Math Function'
      //   Product: '<S70>/Product'

      for (i = 0; i < 3; i++) {
        for (rtb_VectorConcatenate_tmp = 0; rtb_VectorConcatenate_tmp < 3;
             rtb_VectorConcatenate_tmp++) {
          rtb_VectorConcatenate_tmp_0 = 3 * rtb_VectorConcatenate_tmp + i;
          genericQuad_B.Product_tmp_c[rtb_VectorConcatenate_tmp_0] =
            ((genericQuad_B.VectorConcatenate_m[3 * rtb_VectorConcatenate_tmp +
              1] * genericQuad_B.Product_tmp[i + 3] +
              genericQuad_B.VectorConcatenate_m[3 * rtb_VectorConcatenate_tmp] *
              genericQuad_B.Product_tmp[i]) + genericQuad_B.VectorConcatenate_m
             [3 * rtb_VectorConcatenate_tmp + 2] * genericQuad_B.Product_tmp[i +
             6]) + genericQuad_P.Bias1_Bias[rtb_VectorConcatenate_tmp_0];
        }
      }

      // End of Bias: '<S70>/Bias1'

      // RelationalOperator: '<S76>/Compare' incorporates:
      //   Abs: '<S70>/Abs2'
      //   Constant: '<S76>/Constant'

      for (i = 0; i < 9; i++) {
        genericQuad_B.Compare[i] = (fabs(genericQuad_B.Product_tmp_c[i]) >
          genericQuad_P.DirectionCosineMatrixtoQuater_h);
      }

      // End of RelationalOperator: '<S76>/Compare'

      // Logic: '<S70>/Logical Operator1' incorporates:
      //   RelationalOperator: '<S76>/Compare'

      tmp = genericQuad_B.Compare[0];
      for (i = 0; i < 8; i++) {
        tmp = (tmp || genericQuad_B.Compare[i + 1]);
      }

      // If: '<S67>/If' incorporates:
      //   Abs: '<S71>/Abs1'
      //   Bias: '<S71>/Bias'
      //   Concatenate: '<S33>/Vector Concatenate'
      //   Constant: '<S78>/Constant'
      //   Fcn: '<S31>/Fcn13'
      //   Fcn: '<S31>/Fcn23'
      //   Logic: '<S70>/Logical Operator1'
      //   Product: '<S77>/Product'
      //   Product: '<S77>/Product1'
      //   Product: '<S77>/Product2'
      //   Product: '<S77>/Product3'
      //   Product: '<S77>/Product4'
      //   Product: '<S77>/Product5'
      //   RelationalOperator: '<S78>/Compare'
      //   Reshape: '<S77>/Reshape'
      //   Sum: '<S77>/Sum'
      //   Trigonometry: '<S31>/sincos'

      if (fabs((((((genericQuad_B.VectorConcatenate_m[0] *
                    genericQuad_B.VectorConcatenate_m[4] *
                    genericQuad_B.VectorConcatenate_m[8] -
                    genericQuad_B.VectorConcatenate_m[0] *
                    genericQuad_B.VectorConcatenate_m[5] *
                    genericQuad_B.VectorConcatenate_tmp) -
                   genericQuad_B.VectorConcatenate_m[1] *
                   genericQuad_B.VectorConcatenate_m[3] *
                   genericQuad_B.VectorConcatenate_m[8]) +
                  genericQuad_B.VectorConcatenate_m[2] *
                  genericQuad_B.VectorConcatenate_m[3] *
                  genericQuad_B.VectorConcatenate_tmp) +
                 genericQuad_B.VectorConcatenate_m[1] *
                 genericQuad_B.VectorConcatenate_m[5] *
                 -genericQuad_B.Selector_tmp) -
                genericQuad_B.VectorConcatenate_m[2] *
                genericQuad_B.VectorConcatenate_m[4] *
                -genericQuad_B.Selector_tmp) + genericQuad_P.Bias_Bias) >
          genericQuad_P.DirectionCosineMatrixtoQuater_h) {
        // Outputs for IfAction SubSystem: '<S67>/If Not Proper' incorporates:
        //   ActionPort: '<S69>/Action Port'

        // If: '<S69>/If' incorporates:
        //   Constant: '<S69>/Constant'

        if (genericQuad_P.DirectionCosineMatrixtoQuaterni == 2.0) {
          // Outputs for IfAction SubSystem: '<S69>/Warning' incorporates:
          //   ActionPort: '<S75>/Action Port'

          // Assertion: '<S75>/Assertion' incorporates:
          //   Constant: '<S69>/Constant1'

          utAssert(genericQuad_P.Constant1_Value_d != 0.0);

          // End of Outputs for SubSystem: '<S69>/Warning'
        } else if (genericQuad_P.DirectionCosineMatrixtoQuaterni == 3.0) {
          // Outputs for IfAction SubSystem: '<S69>/Error' incorporates:
          //   ActionPort: '<S74>/Action Port'

          // Assertion: '<S74>/Assertion' incorporates:
          //   Constant: '<S69>/Constant1'

          utAssert(genericQuad_P.Constant1_Value_d != 0.0);

          // End of Outputs for SubSystem: '<S69>/Error'
        }

        // End of If: '<S69>/If'
        // End of Outputs for SubSystem: '<S67>/If Not Proper'
      } else if (tmp) {
        // Outputs for IfAction SubSystem: '<S67>/Else If Not Orthogonal' incorporates:
        //   ActionPort: '<S68>/Action Port'

        // If: '<S68>/If' incorporates:
        //   Constant: '<S68>/Constant'

        if (genericQuad_P.DirectionCosineMatrixtoQuaterni == 2.0) {
          // Outputs for IfAction SubSystem: '<S68>/Warning' incorporates:
          //   ActionPort: '<S73>/Action Port'

          // Assertion: '<S73>/Assertion' incorporates:
          //   Constant: '<S68>/Constant1'

          utAssert(genericQuad_P.Constant1_Value_j != 0.0);

          // End of Outputs for SubSystem: '<S68>/Warning'
        } else if (genericQuad_P.DirectionCosineMatrixtoQuaterni == 3.0) {
          // Outputs for IfAction SubSystem: '<S68>/Error' incorporates:
          //   ActionPort: '<S72>/Action Port'

          // Assertion: '<S72>/Assertion' incorporates:
          //   Constant: '<S68>/Constant1'

          utAssert(genericQuad_P.Constant1_Value_j != 0.0);

          // End of Outputs for SubSystem: '<S68>/Error'
        }

        // End of If: '<S68>/If'
        // End of Outputs for SubSystem: '<S67>/Else If Not Orthogonal'
      }

      // End of If: '<S67>/If'
      // End of Outputs for SubSystem: '<S43>/If Warning//Error'
    }
  }

  if (rtmIsMajorTimeStep(genericQuad_M)) {
    if (rtmIsMajorTimeStep(genericQuad_M)) {
      // Update for DiscreteIntegrator: '<S345>/Filter'
      genericQuad_DW.Filter_DSTATE += genericQuad_P.Filter_gainval *
        rtb_FilterCoefficient;

      // Update for DiscreteIntegrator: '<S350>/Integrator'
      genericQuad_DW.Integrator_DSTATE += genericQuad_P.Integrator_gainval *
        genericQuad_B.IntegralGain_g_b;

      // Update for DiscreteIntegrator: '<S302>/Integrator'
      genericQuad_DW.Integrator_DSTATE_o += genericQuad_P.Integrator_gainval_n *
        genericQuad_B.IntegralGain_c;

      // Update for DiscreteIntegrator: '<S297>/Filter'
      genericQuad_DW.Filter_DSTATE_g += genericQuad_P.Filter_gainval_l *
        genericQuad_B.FilterCoefficient_e;

      // Update for DiscreteIntegrator: '<S254>/Integrator'
      genericQuad_DW.Integrator_DSTATE_i += genericQuad_P.Integrator_gainval_c *
        genericQuad_B.IntegralGain_pb;

      // Update for DiscreteIntegrator: '<S249>/Filter'
      genericQuad_DW.Filter_DSTATE_f += genericQuad_P.Filter_gainval_b *
        genericQuad_B.FilterCoefficient_n;
    }
  }                                    // end MajorTimeStep

  if (rtmIsMajorTimeStep(genericQuad_M)) {
    rt_ertODEUpdateContinuousStates(&genericQuad_M->solverInfo);

    // Update absolute time for base rate
    // The "clockTick0" counts the number of times the code of this task has
    //  been executed. The absolute time is the multiplication of "clockTick0"
    //  and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
    //  overflow during the application lifespan selected.

    ++genericQuad_M->Timing.clockTick0;
    genericQuad_M->Timing.t[0] = rtsiGetSolverStopTime
      (&genericQuad_M->solverInfo);

    {
      // Update absolute timer for sample time: [0.002s, 0.0s]
      // The "clockTick1" counts the number of times the code of this task has
      //  been executed. The resolution of this integer timer is 0.002, which is the step size
      //  of the task. Size of "clockTick1" ensures timer will not overflow during the
      //  application lifespan selected.

      genericQuad_M->Timing.clockTick1++;
    }
  }                                    // end MajorTimeStep
}

// Derivatives for root system: '<Root>'
void genericQuad_derivatives(void)
{
  XDot_genericQuad_T *_rtXdot;
  _rtXdot = ((XDot_genericQuad_T *) genericQuad_M->derivs);

  // Derivatives for TransferFcn: '<S4>/Transfer Fcn2'
  _rtXdot->TransferFcn2_CSTATE = 0.0;
  _rtXdot->TransferFcn2_CSTATE += genericQuad_P.TransferFcn2_A *
    genericQuad_X.TransferFcn2_CSTATE;
  _rtXdot->TransferFcn2_CSTATE += genericQuad_B.Komega;

  // Derivatives for TransferFcn: '<S4>/Transfer Fcn1'
  _rtXdot->TransferFcn1_CSTATE = 0.0;
  _rtXdot->TransferFcn1_CSTATE += genericQuad_P.TransferFcn1_A *
    genericQuad_X.TransferFcn1_CSTATE;
  _rtXdot->TransferFcn1_CSTATE += genericQuad_B.Komega1;

  // Derivatives for TransferFcn: '<S4>/Transfer Fcn3'
  _rtXdot->TransferFcn3_CSTATE = 0.0;
  _rtXdot->TransferFcn3_CSTATE += genericQuad_P.TransferFcn3_A *
    genericQuad_X.TransferFcn3_CSTATE;
  _rtXdot->TransferFcn3_CSTATE += genericQuad_B.Komega2;

  // Derivatives for TransferFcn: '<S4>/Transfer Fcn4'
  _rtXdot->TransferFcn4_CSTATE = 0.0;
  _rtXdot->TransferFcn4_CSTATE += genericQuad_P.TransferFcn4_A *
    genericQuad_X.TransferFcn4_CSTATE;
  _rtXdot->TransferFcn4_CSTATE += genericQuad_B.Komega3;

  // Derivatives for Integrator: '<S110>/Integrator'
  _rtXdot->Integrator_CSTATE = genericQuad_B.IntegralGain_p;

  // Derivatives for Integrator: '<S105>/Filter'
  _rtXdot->Filter_CSTATE = genericQuad_B.FilterCoefficient;

  // Derivatives for Integrator: '<S158>/Integrator'
  _rtXdot->Integrator_CSTATE_o = genericQuad_B.IntegralGain_i;

  // Derivatives for Integrator: '<S153>/Filter'
  _rtXdot->Filter_CSTATE_k = genericQuad_B.FilterCoefficient_m;

  // Derivatives for Integrator: '<S206>/Integrator'
  _rtXdot->Integrator_CSTATE_h = genericQuad_B.IntegralGain_g;

  // Derivatives for Integrator: '<S201>/Filter'
  _rtXdot->Filter_CSTATE_h = genericQuad_B.FilterCoefficient_j;

  // Derivatives for Integrator: '<S23>/phi theta psi'
  _rtXdot->phithetapsi_CSTATE[0] =
    genericQuad_B.TmpSignalConversionAtphithetaps[0];

  // Derivatives for Integrator: '<S7>/ub,vb,wb'
  _rtXdot->ubvbwb_CSTATE[0] = genericQuad_B.Sum[0];

  // Derivatives for Integrator: '<S7>/p,q,r '
  _rtXdot->pqr_CSTATE[0] = genericQuad_B.Product2[0];

  // Derivatives for Integrator: '<S398>/Integrator'
  _rtXdot->Integrator_CSTATE_m[0] = genericQuad_B.IntegralGain[0];

  // Derivatives for Integrator: '<S393>/Filter'
  _rtXdot->Filter_CSTATE_e[0] = genericQuad_B.FilterCoefficient_p[0];

  // Derivatives for Integrator: '<S7>/xe,ye,ze'
  _rtXdot->xeyeze_CSTATE[0] = genericQuad_B.Product[0];

  // Derivatives for Integrator: '<S23>/phi theta psi'
  _rtXdot->phithetapsi_CSTATE[1] =
    genericQuad_B.TmpSignalConversionAtphithetaps[1];

  // Derivatives for Integrator: '<S7>/ub,vb,wb'
  _rtXdot->ubvbwb_CSTATE[1] = genericQuad_B.Sum[1];

  // Derivatives for Integrator: '<S7>/p,q,r '
  _rtXdot->pqr_CSTATE[1] = genericQuad_B.Product2[1];

  // Derivatives for Integrator: '<S398>/Integrator'
  _rtXdot->Integrator_CSTATE_m[1] = genericQuad_B.IntegralGain[1];

  // Derivatives for Integrator: '<S393>/Filter'
  _rtXdot->Filter_CSTATE_e[1] = genericQuad_B.FilterCoefficient_p[1];

  // Derivatives for Integrator: '<S7>/xe,ye,ze'
  _rtXdot->xeyeze_CSTATE[1] = genericQuad_B.Product[1];

  // Derivatives for Integrator: '<S23>/phi theta psi'
  _rtXdot->phithetapsi_CSTATE[2] =
    genericQuad_B.TmpSignalConversionAtphithetaps[2];

  // Derivatives for Integrator: '<S7>/ub,vb,wb'
  _rtXdot->ubvbwb_CSTATE[2] = genericQuad_B.Sum[2];

  // Derivatives for Integrator: '<S7>/p,q,r '
  _rtXdot->pqr_CSTATE[2] = genericQuad_B.Product2[2];

  // Derivatives for Integrator: '<S398>/Integrator'
  _rtXdot->Integrator_CSTATE_m[2] = genericQuad_B.IntegralGain[2];

  // Derivatives for Integrator: '<S393>/Filter'
  _rtXdot->Filter_CSTATE_e[2] = genericQuad_B.FilterCoefficient_p[2];

  // Derivatives for Integrator: '<S7>/xe,ye,ze'
  _rtXdot->xeyeze_CSTATE[2] = genericQuad_B.Product[2];
}

// Model initialize function
void genericQuad_initialize(void)
{
  // Registration code

  // initialize non-finites
  rt_InitInfAndNaN(sizeof(real_T));

  // non-finite (run-time) assignments
  genericQuad_P.Saturation4_UpperSat = rtInf;
  genericQuad_P.Saturation5_UpperSat = rtInf;
  genericQuad_P.Saturation6_UpperSat = rtInf;
  genericQuad_P.Saturation7_UpperSat = rtInf;

  {
    // Setup solver object
    rtsiSetSimTimeStepPtr(&genericQuad_M->solverInfo,
                          &genericQuad_M->Timing.simTimeStep);
    rtsiSetTPtr(&genericQuad_M->solverInfo, &rtmGetTPtr(genericQuad_M));
    rtsiSetStepSizePtr(&genericQuad_M->solverInfo,
                       &genericQuad_M->Timing.stepSize0);
    rtsiSetdXPtr(&genericQuad_M->solverInfo, &genericQuad_M->derivs);
    rtsiSetContStatesPtr(&genericQuad_M->solverInfo, (real_T **)
                         &genericQuad_M->contStates);
    rtsiSetNumContStatesPtr(&genericQuad_M->solverInfo,
      &genericQuad_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&genericQuad_M->solverInfo,
      &genericQuad_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&genericQuad_M->solverInfo,
      &genericQuad_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&genericQuad_M->solverInfo,
      &genericQuad_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&genericQuad_M->solverInfo, (&rtmGetErrorStatus
      (genericQuad_M)));
    rtsiSetRTModelPtr(&genericQuad_M->solverInfo, genericQuad_M);
  }

  rtsiSetSimTimeStep(&genericQuad_M->solverInfo, MAJOR_TIME_STEP);
  genericQuad_M->intgData.y = genericQuad_M->odeY;
  genericQuad_M->intgData.f[0] = genericQuad_M->odeF[0];
  genericQuad_M->intgData.f[1] = genericQuad_M->odeF[1];
  genericQuad_M->intgData.f[2] = genericQuad_M->odeF[2];
  genericQuad_M->contStates = ((X_genericQuad_T *) &genericQuad_X);
  genericQuad_M->periodicContStateIndices = ((int_T*) genericQuad_PeriodicIndX);
  genericQuad_M->periodicContStateRanges = ((real_T*) genericQuad_PeriodicRngX);
  rtsiSetSolverData(&genericQuad_M->solverInfo, static_cast<void *>
                    (&genericQuad_M->intgData));
  rtsiSetSolverName(&genericQuad_M->solverInfo,"ode3");
  rtmSetTPtr(genericQuad_M, &genericQuad_M->Timing.tArray[0]);
  genericQuad_M->Timing.stepSize0 = 0.002;

  {
    int32_T i;
    char_T b_zeroDelimName[12];
    char_T b_zeroDelimTopic[8];
    char_T b_zeroDelimTopic_0[5];
    static const char_T tmp[7] = { 'c', 'm', 'd', '_', 'v', 'e', 'l' };

    static const char_T tmp_0[11] = { 'i', 'n', 'i', 't', '_', 'p', 'o', 's',
      'e', '_', 'x' };

    static const char_T tmp_1[11] = { 'i', 'n', 'i', 't', '_', 'p', 'o', 's',
      'e', '_', 'y' };

    static const char_T tmp_2[11] = { 'i', 'n', 'i', 't', '_', 'p', 'o', 's',
      'e', '_', 'z' };

    // Start for If: '<S8>/If'
    genericQuad_DW.If_ActiveSubsystem = -1;

    // Start for If: '<S43>/If1'
    genericQuad_DW.If1_ActiveSubsystem = -1;

    // InitializeConditions for TransferFcn: '<S4>/Transfer Fcn2'
    genericQuad_X.TransferFcn2_CSTATE = 0.0;

    // InitializeConditions for TransferFcn: '<S4>/Transfer Fcn1'
    genericQuad_X.TransferFcn1_CSTATE = 0.0;

    // InitializeConditions for TransferFcn: '<S4>/Transfer Fcn3'
    genericQuad_X.TransferFcn3_CSTATE = 0.0;

    // InitializeConditions for TransferFcn: '<S4>/Transfer Fcn4'
    genericQuad_X.TransferFcn4_CSTATE = 0.0;

    // InitializeConditions for Integrator: '<S110>/Integrator'
    genericQuad_X.Integrator_CSTATE =
      genericQuad_P.PIDVelocityx_InitialCondition_o;

    // InitializeConditions for Integrator: '<S105>/Filter'
    genericQuad_X.Filter_CSTATE = genericQuad_P.PIDVelocityx_InitialConditionFo;

    // InitializeConditions for Integrator: '<S158>/Integrator'
    genericQuad_X.Integrator_CSTATE_o =
      genericQuad_P.PIDVelocityy_InitialCondition_i;

    // InitializeConditions for Integrator: '<S153>/Filter'
    genericQuad_X.Filter_CSTATE_k =
      genericQuad_P.PIDVelocityy_InitialConditionFo;

    // InitializeConditions for Integrator: '<S206>/Integrator'
    genericQuad_X.Integrator_CSTATE_h =
      genericQuad_P.PIDVelocityz_InitialCondition_g;

    // InitializeConditions for Integrator: '<S201>/Filter'
    genericQuad_X.Filter_CSTATE_h =
      genericQuad_P.PIDVelocityz_InitialConditionFo;

    // InitializeConditions for Integrator: '<S23>/phi theta psi'
    genericQuad_X.phithetapsi_CSTATE[0] = genericQuad_P.uDOFEulerAngles2_eul_0[0];

    // InitializeConditions for Integrator: '<S7>/ub,vb,wb'
    genericQuad_X.ubvbwb_CSTATE[0] = genericQuad_P.uDOFEulerAngles2_Vm_0[0];

    // InitializeConditions for Integrator: '<S7>/p,q,r '
    genericQuad_X.pqr_CSTATE[0] = genericQuad_P.uDOFEulerAngles2_pm_0[0];

    // InitializeConditions for Integrator: '<S398>/Integrator'
    genericQuad_X.Integrator_CSTATE_m[0] =
      genericQuad_P.PIDattitude_InitialConditionF_g;

    // InitializeConditions for Integrator: '<S393>/Filter'
    genericQuad_X.Filter_CSTATE_e[0] =
      genericQuad_P.PIDattitude_InitialConditionFor;

    // InitializeConditions for Integrator: '<S7>/xe,ye,ze'
    genericQuad_X.xeyeze_CSTATE[0] = genericQuad_P.uDOFEulerAngles2_xme_0[0];

    // InitializeConditions for Integrator: '<S23>/phi theta psi'
    genericQuad_X.phithetapsi_CSTATE[1] = genericQuad_P.uDOFEulerAngles2_eul_0[1];

    // InitializeConditions for Integrator: '<S7>/ub,vb,wb'
    genericQuad_X.ubvbwb_CSTATE[1] = genericQuad_P.uDOFEulerAngles2_Vm_0[1];

    // InitializeConditions for Integrator: '<S7>/p,q,r '
    genericQuad_X.pqr_CSTATE[1] = genericQuad_P.uDOFEulerAngles2_pm_0[1];

    // InitializeConditions for Integrator: '<S398>/Integrator'
    genericQuad_X.Integrator_CSTATE_m[1] =
      genericQuad_P.PIDattitude_InitialConditionF_g;

    // InitializeConditions for Integrator: '<S393>/Filter'
    genericQuad_X.Filter_CSTATE_e[1] =
      genericQuad_P.PIDattitude_InitialConditionFor;

    // InitializeConditions for Integrator: '<S7>/xe,ye,ze'
    genericQuad_X.xeyeze_CSTATE[1] = genericQuad_P.uDOFEulerAngles2_xme_0[1];

    // InitializeConditions for Integrator: '<S23>/phi theta psi'
    genericQuad_X.phithetapsi_CSTATE[2] = genericQuad_P.uDOFEulerAngles2_eul_0[2];

    // InitializeConditions for Integrator: '<S7>/ub,vb,wb'
    genericQuad_X.ubvbwb_CSTATE[2] = genericQuad_P.uDOFEulerAngles2_Vm_0[2];

    // InitializeConditions for Integrator: '<S7>/p,q,r '
    genericQuad_X.pqr_CSTATE[2] = genericQuad_P.uDOFEulerAngles2_pm_0[2];

    // InitializeConditions for Integrator: '<S398>/Integrator'
    genericQuad_X.Integrator_CSTATE_m[2] =
      genericQuad_P.PIDattitude_InitialConditionF_g;

    // InitializeConditions for Integrator: '<S393>/Filter'
    genericQuad_X.Filter_CSTATE_e[2] =
      genericQuad_P.PIDattitude_InitialConditionFor;

    // InitializeConditions for Integrator: '<S7>/xe,ye,ze'
    genericQuad_X.xeyeze_CSTATE[2] = genericQuad_P.uDOFEulerAngles2_xme_0[2];

    // InitializeConditions for DiscreteIntegrator: '<S345>/Filter'
    genericQuad_DW.Filter_DSTATE = genericQuad_P.PIDangularroll_InitialCondition;

    // InitializeConditions for DiscreteIntegrator: '<S350>/Integrator'
    genericQuad_DW.Integrator_DSTATE =
      genericQuad_P.PIDangularroll_InitialConditi_d;

    // InitializeConditions for DiscreteIntegrator: '<S302>/Integrator'
    genericQuad_DW.Integrator_DSTATE_o =
      genericQuad_P.PIDangularpitch_InitialCondit_l;

    // InitializeConditions for DiscreteIntegrator: '<S297>/Filter'
    genericQuad_DW.Filter_DSTATE_g =
      genericQuad_P.PIDangularpitch_InitialConditio;

    // InitializeConditions for DiscreteIntegrator: '<S254>/Integrator'
    genericQuad_DW.Integrator_DSTATE_i =
      genericQuad_P.PIDangulayaw_InitialCondition_k;

    // InitializeConditions for DiscreteIntegrator: '<S249>/Filter'
    genericQuad_DW.Filter_DSTATE_f =
      genericQuad_P.PIDangulayaw_InitialConditionFo;

    // SystemInitialize for Atomic SubSystem: '<Root>/Subscribe1'
    // SystemInitialize for Enabled SubSystem: '<S3>/Enabled Subsystem'
    // SystemInitialize for Outport: '<S6>/Out1' incorporates:
    //   Inport: '<S6>/In1'

    genericQuad_B.In1 = genericQuad_P.Out1_Y0;

    // End of SystemInitialize for SubSystem: '<S3>/Enabled Subsystem'

    // Start for MATLABSystem: '<S3>/SourceBlock'
    genericQuad_DW.obj_c.matlabCodegenIsDeleted = false;
    genericQuad_DW.obj_c.isInitialized = 1;
    for (i = 0; i < 7; i++) {
      b_zeroDelimTopic[i] = tmp[i];
    }

    b_zeroDelimTopic[7] = '\x00';
    Sub_genericQuad_426.createSubscriber(&b_zeroDelimTopic[0], 51);
    genericQuad_DW.obj_c.isSetupComplete = true;

    // End of Start for MATLABSystem: '<S3>/SourceBlock'
    // End of SystemInitialize for SubSystem: '<Root>/Subscribe1'

    // SystemInitialize for IfAction SubSystem: '<S8>/Negative Trace'
    // Start for If: '<S41>/Find Maximum Diagonal Value'
    genericQuad_DW.FindMaximumDiagonalValue_Active = -1;

    // End of SystemInitialize for SubSystem: '<S8>/Negative Trace'

    // SystemInitialize for Atomic SubSystem: '<Root>/Publish'
    // Start for MATLABSystem: '<S2>/SinkBlock'
    genericQuad_DW.obj_m.matlabCodegenIsDeleted = false;
    genericQuad_DW.obj_m.isInitialized = 1;

    // End of SystemInitialize for SubSystem: '<Root>/Publish'

    // SystemInitialize for Merge: '<S8>/Merge'
    genericQuad_B.Merge[0] = genericQuad_P.Merge_InitialOutput[0];

    // SystemInitialize for Atomic SubSystem: '<Root>/Publish'
    // Start for MATLABSystem: '<S2>/SinkBlock'
    b_zeroDelimTopic_0[0] = 'p';

    // End of SystemInitialize for SubSystem: '<Root>/Publish'

    // SystemInitialize for Merge: '<S8>/Merge'
    genericQuad_B.Merge[1] = genericQuad_P.Merge_InitialOutput[1];

    // SystemInitialize for Atomic SubSystem: '<Root>/Publish'
    // Start for MATLABSystem: '<S2>/SinkBlock'
    b_zeroDelimTopic_0[1] = 'o';

    // End of SystemInitialize for SubSystem: '<Root>/Publish'

    // SystemInitialize for Merge: '<S8>/Merge'
    genericQuad_B.Merge[2] = genericQuad_P.Merge_InitialOutput[2];

    // SystemInitialize for Atomic SubSystem: '<Root>/Publish'
    // Start for MATLABSystem: '<S2>/SinkBlock'
    b_zeroDelimTopic_0[2] = 's';

    // End of SystemInitialize for SubSystem: '<Root>/Publish'

    // SystemInitialize for Merge: '<S8>/Merge'
    genericQuad_B.Merge[3] = genericQuad_P.Merge_InitialOutput[3];

    // SystemInitialize for Atomic SubSystem: '<Root>/Publish'
    // Start for MATLABSystem: '<S2>/SinkBlock'
    b_zeroDelimTopic_0[3] = 'e';
    b_zeroDelimTopic_0[4] = '\x00';
    Pub_genericQuad_438.createPublisher(&b_zeroDelimTopic_0[0], 1);
    genericQuad_DW.obj_m.isSetupComplete = true;

    // End of SystemInitialize for SubSystem: '<Root>/Publish'

    // Start for MATLABSystem: '<Root>/Get Parameter'
    genericQuad_DW.obj_p.matlabCodegenIsDeleted = false;
    genericQuad_DW.obj_p.isInitialized = 1;
    for (i = 0; i < 11; i++) {
      b_zeroDelimName[i] = tmp_0[i];
    }

    b_zeroDelimName[11] = '\x00';
    ParamGet_genericQuad_459.initialize(&b_zeroDelimName[0]);
    ParamGet_genericQuad_459.initialize_error_codes(0, 1, 2, 3);
    ParamGet_genericQuad_459.set_initial_value(0.0);
    genericQuad_DW.obj_p.isSetupComplete = true;

    // End of Start for MATLABSystem: '<Root>/Get Parameter'

    // Start for MATLABSystem: '<Root>/Get Parameter1'
    genericQuad_DW.obj_a.matlabCodegenIsDeleted = false;
    genericQuad_DW.obj_a.isInitialized = 1;
    for (i = 0; i < 11; i++) {
      b_zeroDelimName[i] = tmp_1[i];
    }

    b_zeroDelimName[11] = '\x00';
    ParamGet_genericQuad_460.initialize(&b_zeroDelimName[0]);
    ParamGet_genericQuad_460.initialize_error_codes(0, 1, 2, 3);
    ParamGet_genericQuad_460.set_initial_value(0.0);
    genericQuad_DW.obj_a.isSetupComplete = true;

    // End of Start for MATLABSystem: '<Root>/Get Parameter1'

    // Start for MATLABSystem: '<Root>/Get Parameter2'
    genericQuad_DW.obj.matlabCodegenIsDeleted = false;
    genericQuad_DW.obj.isInitialized = 1;
    for (i = 0; i < 11; i++) {
      b_zeroDelimName[i] = tmp_2[i];
    }

    b_zeroDelimName[11] = '\x00';
    ParamGet_genericQuad_461.initialize(&b_zeroDelimName[0]);
    ParamGet_genericQuad_461.initialize_error_codes(0, 1, 2, 3);
    ParamGet_genericQuad_461.set_initial_value(0.0);
    genericQuad_DW.obj.isSetupComplete = true;

    // End of Start for MATLABSystem: '<Root>/Get Parameter2'

    // InitializeConditions for root-level periodic continuous states
    {
      int_T rootPeriodicContStateIndices[3] = { 0, 1, 2 };

      real_T rootPeriodicContStateRanges[6] = { -3.1415926535897931,
        3.1415926535897931, -3.1415926535897931, 3.1415926535897931,
        -3.1415926535897931, 3.1415926535897931 };

      (void) memcpy((void*)genericQuad_PeriodicIndX,
                    rootPeriodicContStateIndices,
                    3*sizeof(int_T));
      (void) memcpy((void*)genericQuad_PeriodicRngX, rootPeriodicContStateRanges,
                    6*sizeof(real_T));
    }
  }
}

// Model terminate function
void genericQuad_terminate(void)
{
  // Terminate for Atomic SubSystem: '<Root>/Subscribe1'
  // Terminate for MATLABSystem: '<S3>/SourceBlock'
  if (!genericQuad_DW.obj_c.matlabCodegenIsDeleted) {
    genericQuad_DW.obj_c.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<S3>/SourceBlock'
  // End of Terminate for SubSystem: '<Root>/Subscribe1'

  // Terminate for MATLABSystem: '<Root>/Get Parameter'
  if (!genericQuad_DW.obj_p.matlabCodegenIsDeleted) {
    genericQuad_DW.obj_p.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<Root>/Get Parameter'

  // Terminate for MATLABSystem: '<Root>/Get Parameter1'
  if (!genericQuad_DW.obj_a.matlabCodegenIsDeleted) {
    genericQuad_DW.obj_a.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<Root>/Get Parameter1'

  // Terminate for MATLABSystem: '<Root>/Get Parameter2'
  if (!genericQuad_DW.obj.matlabCodegenIsDeleted) {
    genericQuad_DW.obj.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<Root>/Get Parameter2'

  // Terminate for Atomic SubSystem: '<Root>/Publish'
  // Terminate for MATLABSystem: '<S2>/SinkBlock'
  if (!genericQuad_DW.obj_m.matlabCodegenIsDeleted) {
    genericQuad_DW.obj_m.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<S2>/SinkBlock'
  // End of Terminate for SubSystem: '<Root>/Publish'
}

//
// File trailer for generated code.
//
// [EOF]
//
