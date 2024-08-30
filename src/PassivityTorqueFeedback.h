/*
 * Copyright 2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_control/GlobalPlugin.h>
#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/VirtualTorqueSensor.h>
#include <mc_rbdyn/ExternalTorqueSensor.h>
#include <RBDyn/Coriolis.h>
#include <RBDyn/FD.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/FA.h>
#include <RBDyn/ID.h>
#include <Eigen/src/Core/IO.h>
#include <jrl-qp/experimental/BoxAndSingleConstraintSolver.h>

#include "utils/ROSSubscriber.h"
namespace mc_plugin
{

struct PassivityTorqueFeedback : public mc_control::GlobalPlugin
{
  void init(mc_control::MCGlobalController & controller, const mc_rtc::Configuration & config) override;

  void reset(mc_control::MCGlobalController & controller) override;

  void before(mc_control::MCGlobalController & controller) override;

  void after(mc_control::MCGlobalController & controller) override;

  mc_control::GlobalPlugin::GlobalPluginConfiguration configuration() override;

  ~PassivityTorqueFeedback() override;

private:
  void addGUI(mc_control::MCGlobalController & controller);
  void addLOG(mc_control::MCGlobalController & controller);
  void removeLOG(mc_control::MCGlobalController & controller);

  std::string getIntegralTermType(void);
  void setIntegralTermType(std::string type);
  void torque_continuity(double lambda_massmatrix,double lambda_id,double lambda_diag_massmatrix);
  void torque_activation(mc_control::MCGlobalController & controller);

private:
  // /** instance of integral term anti windup */
  // torque_control::IntegralTermAntiWindup integralTermAntiWindup_;
  /** Deactivate the output */
  bool is_active_;
  /** Deactivate the output */
  bool estimation_;
  /** Disable unecessary console outputs */
  bool verbose_;
  /** Track if the value of perc_ is changing */
  bool is_changing_;
  /** simple or filtered integration */
  bool filtered_activated_;
  /** Count to track time */
  int count_;
  /** Bool for taking the coriolis effects in consideration or not*/
  bool coriolis_indicator_;
  /** Corresponding value of the booleen for taking the coriolis effects in consideration or not*/
  double coriolis_indicator_value_;
  /** Control timestep used for integration*/
  double dt_;
  /** Coriolis object for computation of coriolis matrix from robot configuration */
  rbd::Coriolis * coriolis_;
  /** Forward Dynamics */
  rbd::ForwardDynamics * fd_;
  /** Inverse Dynamics */
  rbd::InverseDynamics * id_;
  /** Virtual torque sensor used to hold value of feedback term */
  mc_rbdyn::VirtualTorqueSensor * virtual_torque_sensor_;
  /** External torque sensor used to hold value of feedback term */
  mc_rbdyn::ExternalTorqueSensor * external_torque_sensor_;
  /** Positive gain value for mass matrix*/
  double lambda_massmatrix_;
  /** Positive gain value for identity matrix */
  double lambda_id_;
  /** Positive gain value for the massmatrix diagonale */
  double lambda_diag_massmatrix_; 
  /** Slow filter time constant */
  double phi_slow_;
  /** Fast filter time constant */
  double phi_fast_;
  /** Exponential of slow filter time constant */
  double exp_phi_slow_;
  /** Exponential of fast filter time constant */
  double exp_phi_fast_;
  /** Weight for fast filtered s (mu) used only for filtered integral type */
  double fast_filter_weight_;
  /** Percentage for the torque limit in the antiwindup regarding the physical limits of the robot */
  double perc_;
  /** Percentage for the torque limit in the antiwindup regarding the physical limits of the robot, used for the exponential growing due to change in the GUI */
  double perc_target_;
  /** Power */
  double power_;
  /** Configuration velocity reference from the integral of configuration acceleration reference*/
  Eigen::VectorXd alpha_r_;
  /** Velocity gain matrix */
  Eigen::MatrixXd K_;
  /** Gain Matrix for the comparison with External Torque Estimation */
  Eigen::MatrixXd K_test_;
  /** Velocity gain matrix diagonal*/
  Eigen::MatrixXd D_;
  /** M dot test*/
  // Eigen::MatrixXd M_dot_;
  // /** M old*/
  // Eigen::MatrixXd M_old_;
  // /** M new*/
  // Eigen::MatrixXd M_new_;
  /** Coriolis matrix */
  Eigen::MatrixXd C_;
  /** Vector of velocity error using integration of acceleration reference for velocity reference */
  Eigen::VectorXd s_;
  /** Newest value of vector of velocity error used only for filtered integral type */
  Eigen::VectorXd new_s_;
  /** Previous value of vector of velocity error used only for filtered integral type */
  Eigen::VectorXd prev_s_;
  /** Slow filtered value of vector of velocity error used only for filtered integral type */
  Eigen::VectorXd slow_filtered_s_;
  /** Fast filtered value of vector of velocity error used only for filtered integral type */
  Eigen::VectorXd fast_filtered_s_;
  /** Current from the motor of each joint */
  Eigen::VectorXd motor_current_;
  /** Torque integral term */
  Eigen::VectorXd tau_;
  /** Torque QP term */
  Eigen::VectorXd tau_qp_;
  /** Torque send*/
  Eigen::VectorXd tau_Out_;
  /** Torque QP and Integral and Coriolis */
  Eigen::VectorXd tau_sum_;
  /** Torque coriolis term */
  Eigen::VectorXd tau_coriolis_;
  /** Torque current motor */
  Eigen::VectorXd tau_current_;
  /** Torque test for comparison*/
  Eigen::VectorXd tau_test_;
  /** Torque lower bound */
  Eigen::VectorXd torqueL_;
  /** Torque lower bound */
  Eigen::VectorXd torqueU_;
  /** maximal angular acceleration */
  Eigen::Vector3d maxAngAcc_ ; 
  /** maximal linear acceleration */
  Eigen::Vector3d maxLinAcc_ ; 
  /** solver for antiwindup */
  std::shared_ptr<jrl::qp::experimental::BoxAndSingleConstraintSolver> solver_;


  Eigen::VectorXd integralTerm;
  Eigen::VectorXd residual;
  rbd::Jacobian jac;
  double residualGains;
  std::string referenceFrame;
  std::string force_sensor_topic_;
  ROSWrenchStampedSubscriber wrench_sub_;
  bool ros_force_sensor_;
  Eigen::VectorXd FTSensorTorques;
  Eigen::VectorXd prevFTSensorTorques;
  Eigen::VectorXd filteredFTSensorTorques;
  Eigen::VectorXd newExternalTorques;
  Eigen::VectorXd externalTorques;
  Eigen::VectorXd externalTorquesMinus;
  Eigen::VectorXd filteredExternalTorques;
  sva::ForceVecd externalForces;
  sva::ForceVecd externalForcesResidual;
  sva::ForceVecd newExternalForces;
  sva::ForceVecd filteredFTSensorForces;
  Eigen::Vector6d externalForcesFT;


Eigen::IOFormat format;
};

}

