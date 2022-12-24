#include "slros_initialize.h"

ros::NodeHandle * SLROSNodePtr;
const std::string SLROSNodeName = "genericQuad";

// For Block genericQuad/Subscribe1
SimulinkSubscriber<geometry_msgs::Twist, SL_Bus_genericQuad_geometry_msgs_Twist> Sub_genericQuad_426;

// For Block genericQuad/Publish
SimulinkPublisher<geometry_msgs::PoseStamped, SL_Bus_genericQuad_geometry_msgs_PoseStamped> Pub_genericQuad_438;

// For Block genericQuad/Get Parameter
SimulinkParameterGetter<real64_T, double> ParamGet_genericQuad_459;

// For Block genericQuad/Get Parameter1
SimulinkParameterGetter<real64_T, double> ParamGet_genericQuad_460;

// For Block genericQuad/Get Parameter2
SimulinkParameterGetter<real64_T, double> ParamGet_genericQuad_461;

void slros_node_init(int argc, char** argv)
{
  ros::init(argc, argv, SLROSNodeName);
  SLROSNodePtr = new ros::NodeHandle("~");
}

