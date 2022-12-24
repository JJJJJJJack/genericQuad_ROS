#ifndef _SLROS_INITIALIZE_H_
#define _SLROS_INITIALIZE_H_

#include "slros_busmsg_conversion.h"
#include "slros_generic.h"
#include "genericQuad_types.h"

extern ros::NodeHandle * SLROSNodePtr;
extern const std::string SLROSNodeName;

// For Block genericQuad/Subscribe1
extern SimulinkSubscriber<geometry_msgs::Twist, SL_Bus_genericQuad_geometry_msgs_Twist> Sub_genericQuad_426;

// For Block genericQuad/Publish
extern SimulinkPublisher<geometry_msgs::PoseStamped, SL_Bus_genericQuad_geometry_msgs_PoseStamped> Pub_genericQuad_438;

// For Block genericQuad/Get Parameter
extern SimulinkParameterGetter<real64_T, double> ParamGet_genericQuad_459;

// For Block genericQuad/Get Parameter1
extern SimulinkParameterGetter<real64_T, double> ParamGet_genericQuad_460;

// For Block genericQuad/Get Parameter2
extern SimulinkParameterGetter<real64_T, double> ParamGet_genericQuad_461;

void slros_node_init(int argc, char** argv);

#endif
