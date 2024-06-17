/*
 * Copyright 2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_control/GlobalPlugin.h>
#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/VirtualTorqueSensor.h>
#include <RBDyn/Coriolis.h>
#include <RBDyn/FD.h>
#include <jrl-qp/experimental/BoxAndSingleConstraintSolver.h>


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

public:
  enum IntegralTermType
  {
    Simple,
    Filtered
  };

private:
  // /** instance of integral term anti windup */
  // torque_control::IntegralTermAntiWindup integralTermAntiWindup_; // est-ce que je mets une etoile ?
  /** Deactivate the output */
  bool is_active_;
  /** Disable unecessary console outputs */
  bool verbose_;
  /** Count to track time */
  int count_;
  /** Control timestep used for integration*/
  double dt_;
  /** Coriolis object for computation of coriolis matrix from robot configuration */
  rbd::Coriolis * coriolis_;
  /** Coriolis object for computation of coriolis matrix from robot configuration */
  rbd::ForwardDynamics * fd_;
  /** Virtual torque sensor used to hold value of feedback term */
  mc_rbdyn::VirtualTorqueSensor * virtual_torque_sensor_;
  /** Type of integral used */
  IntegralTermType integralType_;
  /** Integral term lower limit*/
  double integral_boundl_;
  /** Integral term upper limit*/
  double integral_boundu_;
  /** Positive velocity gain value */
  double lambda_massmatrix_;
  /** Positive velocity gain value */
  double lambda_id_;
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
  /** Configuration velocity reference from the integral of configuration acceleration reference*/
  Eigen::VectorXd alpha_r_;
  /** Velocity gain matrix */
  Eigen::MatrixXd K_;
  /** Integral gain matrix */
  Eigen::MatrixXd L_;
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
  /** Torque feedback term */
  Eigen::VectorXd tau_;
  /** Torque reference term */
  Eigen::VectorXd tau_ref_;
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
};

}

