//
// File: genericQuad_types.h
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
#ifndef RTW_HEADER_genericQuad_types_h_
#define RTW_HEADER_genericQuad_types_h_
#include "rtwtypes.h"

// Model Code Variants
#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_ros_time_Time_
#define DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_ros_time_Time_

// MsgType=ros_time/Time
struct SL_Bus_genericQuad_ros_time_Time
{
  real_T Sec;
  real_T Nsec;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_ROSVariableLengthArrayInfo_
#define DEFINED_TYPEDEF_FOR_SL_Bus_ROSVariableLengthArrayInfo_

struct SL_Bus_ROSVariableLengthArrayInfo
{
  uint32_T CurrentLength;
  uint32_T ReceivedLength;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_std_msgs_Header_
#define DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_std_msgs_Header_

// MsgType=std_msgs/Header
struct SL_Bus_genericQuad_std_msgs_Header
{
  uint32_T Seq;

  // MsgType=ros_time/Time
  SL_Bus_genericQuad_ros_time_Time Stamp;

  // PrimitiveROSType=string:IsVarLen=1:VarLenCategory=data:VarLenElem=FrameId_SL_Info:TruncateAction=warn 
  uint8_T FrameId[128];

  // IsVarLen=1:VarLenCategory=length:VarLenElem=FrameId
  SL_Bus_ROSVariableLengthArrayInfo FrameId_SL_Info;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_Point_
#define DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_Point_

// MsgType=geometry_msgs/Point
struct SL_Bus_genericQuad_geometry_msgs_Point
{
  real_T X;
  real_T Y;
  real_T Z;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_Quaternion_
#define DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_Quaternion_

// MsgType=geometry_msgs/Quaternion
struct SL_Bus_genericQuad_geometry_msgs_Quaternion
{
  real_T X;
  real_T Y;
  real_T Z;
  real_T W;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_Pose_
#define DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_Pose_

// MsgType=geometry_msgs/Pose
struct SL_Bus_genericQuad_geometry_msgs_Pose
{
  // MsgType=geometry_msgs/Point
  SL_Bus_genericQuad_geometry_msgs_Point Position;

  // MsgType=geometry_msgs/Quaternion
  SL_Bus_genericQuad_geometry_msgs_Quaternion Orientation;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_PoseStamped_
#define DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_PoseStamped_

// MsgType=geometry_msgs/PoseStamped
struct SL_Bus_genericQuad_geometry_msgs_PoseStamped
{
  // MsgType=std_msgs/Header
  SL_Bus_genericQuad_std_msgs_Header Header;

  // MsgType=geometry_msgs/Pose
  SL_Bus_genericQuad_geometry_msgs_Pose Pose;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_Vector3_
#define DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_Vector3_

// MsgType=geometry_msgs/Vector3
struct SL_Bus_genericQuad_geometry_msgs_Vector3
{
  real_T X;
  real_T Y;
  real_T Z;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_Twist_
#define DEFINED_TYPEDEF_FOR_SL_Bus_genericQuad_geometry_msgs_Twist_

// MsgType=geometry_msgs/Twist
struct SL_Bus_genericQuad_geometry_msgs_Twist
{
  // MsgType=geometry_msgs/Vector3
  SL_Bus_genericQuad_geometry_msgs_Vector3 Linear;

  // MsgType=geometry_msgs/Vector3
  SL_Bus_genericQuad_geometry_msgs_Vector3 Angular;
};

#endif

#ifndef struct_f_robotics_slcore_internal_bl_T
#define struct_f_robotics_slcore_internal_bl_T

struct f_robotics_slcore_internal_bl_T
{
  int32_T __dummy;
};

#endif                                // struct_f_robotics_slcore_internal_bl_T

#ifndef struct_ros_slros_internal_block_GetP_T
#define struct_ros_slros_internal_block_GetP_T

struct ros_slros_internal_block_GetP_T
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  f_robotics_slcore_internal_bl_T SampleTimeHandler;
};

#endif                                // struct_ros_slros_internal_block_GetP_T

#ifndef struct_ros_slroscpp_internal_block_P_T
#define struct_ros_slroscpp_internal_block_P_T

struct ros_slroscpp_internal_block_P_T
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
};

#endif                                // struct_ros_slroscpp_internal_block_P_T

#ifndef struct_ros_slroscpp_internal_block_S_T
#define struct_ros_slroscpp_internal_block_S_T

struct ros_slroscpp_internal_block_S_T
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
};

#endif                                // struct_ros_slroscpp_internal_block_S_T

// Parameters (default storage)
typedef struct P_genericQuad_T_ P_genericQuad_T;

// Forward declaration for rtModel
typedef struct tag_RTM_genericQuad_T RT_MODEL_genericQuad_T;

#endif                                 // RTW_HEADER_genericQuad_types_h_

//
// File trailer for generated code.
//
// [EOF]
//
