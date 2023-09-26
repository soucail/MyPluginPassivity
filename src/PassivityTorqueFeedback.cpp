#include "PassivityTorqueFeedback.h"

#include <mc_control/GlobalPluginMacros.h>
#include <mc_tvm/Robot.h>

namespace mc_plugin
{

PassivityTorqueFeedback::~PassivityTorqueFeedback() = default;

void PassivityTorqueFeedback::init(mc_control::MCGlobalController & controller, const mc_rtc::Configuration & config)
{
  mc_rtc::log::info("[PassivityTorqueFeedback][init] Init called with configuration:\n{}", config.dump(true, true));

  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);
  auto & robot = ctl.robot(ctl.robots()[0].name());
  auto nrDof = robot.mb().nrDof();

  // ========== Check for Virtual Torque Sensor to hold feedback term value ========== //
  if (!robot.hasDevice<mc_rbdyn::VirtualTorqueSensor>("virtualTorqueSensor"))
  {
    mc_rtc::log::error_and_throw<std::runtime_error>("[PassivityTorqueFeedback][Init] No \"VirtualTorqueSensor\" with the name \"virtualTorqueSensor\" found in the robot module, please add one to the robot's RobotModule.");
  }
  virtual_torque_sensor_ = &robot.device<mc_rbdyn::VirtualTorqueSensor>("virtualTorqueSensor");

  dt_ = ctl.timestep();
  coriolis_ = new rbd::Coriolis(robot.mb());
  fd_ = new rbd::ForwardDynamics(robot.mb());

  // ====================  Load  config  ==================== //
  verbose_ = config("verbose", false);
  auto integrationTypeStr = config("integration_type",(std::string)"Simple");
  auto velGainTypeStr = config("velocity_gain_type",(std::string)"Diagonal");
  lambda_ = config("lambda", 1.0);
  fast_filter_weight_ = config("fast_filter_weight",0.9);
  phi_slow_ = config("phi_slow",0.1);
  phi_fast_ = config("phi_fast",100.0);
  integral_boundl_ = config("integral_bound")("low",-20.0);
  integral_boundu_ = config("integral_bound")("up",-20.0);
  // ==================== Config loaded ==================== //

  // Integral term type
  if (integrationTypeStr.compare("Simple") == 0) integralType_ = IntegralTermType::Simple;
  else if (integrationTypeStr.compare("Filtered") == 0) integralType_ = IntegralTermType::Filtered;
  else mc_rtc::log::error_and_throw("[PassivityTorqueFeedback][init] Integration type loaded from config does not match any implemented type.\nPlease use one of the following type : {Simple, Filtered}");

  // Velocity gain type
  if (velGainTypeStr.compare("Diagonal") == 0) velGainType_ = VelocityGainType::Diagonal;
  else if (velGainTypeStr.compare("MassDiagonal") == 0) velGainType_ = VelocityGainType::MassDiagonal;
  else if (velGainTypeStr.compare("MassMatrix") == 0) velGainType_ = VelocityGainType::MassMatrix;
  else mc_rtc::log::error_and_throw("[PassivityTorqueFeedback][init] Velocity gain type loaded from config does not match any implemented type.\nPlease use one of the following type : {Diagonal, MassDiagonal, MassMatrix}");

  // Initializing all variables
  exp_phi_slow_ = exp(-dt_*phi_slow_);
  exp_phi_fast_ = exp(-dt_*phi_fast_);
  alpha_r_ = Eigen::VectorXd::Zero(nrDof);
  K_ = Eigen::MatrixXd::Zero(nrDof,nrDof);
  L_ = Eigen::MatrixXd::Zero(nrDof,nrDof);
  s_ = Eigen::VectorXd::Zero(nrDof);
  new_s_ = Eigen::VectorXd::Zero(nrDof);
  prev_s_ = Eigen::VectorXd::Zero(nrDof);
  slow_filtered_s_ = Eigen::VectorXd::Zero(nrDof);
  fast_filtered_s_ = Eigen::VectorXd::Zero(nrDof);
  tau_ = Eigen::VectorXd::Zero(nrDof);

  addGUI(controller);
  addLOG(controller);

  count_ = 0;

  mc_rtc::log::success("[PassivityTorqueFeedback][init] Initialization completed");
}

void PassivityTorqueFeedback::reset(mc_control::MCGlobalController & controller)
{
  removeLOG(controller);
  mc_rtc::log::info("[PassivityTorqueFeedback] Reset called");
}

void PassivityTorqueFeedback::before(mc_control::MCGlobalController & controller)
{

  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);
  auto & robot = ctl.robot("kinova");
  auto & realRobot = ctl.realRobot("kinova");

  if (robot.encoderVelocities().empty())
  {
    return;
  }


  Eigen::VectorXd alpha_d(robot.mb().nrDof());
  Eigen::VectorXd alpha(robot.mb().nrDof());

  rbd::paramToVector(robot.alphaD(),alpha_d);
  rbd::paramToVector(robot.alpha(),alpha);
  fd_->forwardDynamics(robot.mb(),robot.mbc());

  alpha_r_ +=  alpha_d*dt_;

  if (velGainType_ == VelocityGainType::Diagonal) K_ = lambda_*Eigen::MatrixXd::Identity(robot.mb().nrDof(),robot.mb().nrDof());
  else if (velGainType_ == VelocityGainType::MassDiagonal) K_ = lambda_* fd_->H().diagonal();
  else if (velGainType_ == VelocityGainType::MassMatrix) K_ = lambda_* fd_->H();

  L_ = coriolis_->coriolis(robot.mb(), robot.mbc()) + K_;

  if (integralType_ == IntegralTermType::Simple)
  {
    s_ = alpha_r_ - alpha;
  }
  else if (integralType_ == IntegralTermType::Filtered)
  {
    new_s_ = alpha_r_ - alpha;

    slow_filtered_s_ = exp_phi_slow_*slow_filtered_s_ + new_s_ - prev_s_;
    fast_filtered_s_ = exp_phi_fast_*fast_filtered_s_ + new_s_ - prev_s_;

    prev_s_ = new_s_;

    s_ = fast_filter_weight_ * fast_filtered_s_ + (1 - fast_filter_weight_) * slow_filtered_s_;
    // mc_rtc::log::info("Filtered integration s: {}", s_.transpose());
  }

  if(is_active_)
  {
    tau_ = L_*s_;
    tau_ = tau_.cwiseMax(integral_boundl_).cwiseMin(integral_boundu_);
  }
  else
  {
    tau_.setZero();
  }

  if((count_%1000) == 0 && verbose_)
  {
    mc_rtc::log::info("qdd : {}",alpha_d.transpose());
    mc_rtc::log::info("qd : {}",alpha.transpose());
    mc_rtc::log::info("Filtered integration s: {}", s_.transpose());
    mc_rtc::log::info("Tau : {}", tau_.transpose());
    mc_rtc::log::info("Tau clamped: {}", tau_.transpose());
  }

  virtual_torque_sensor_->torques(tau_);

  count_++;
  // mc_rtc::log::info("PassivityTorqueFeedback::before");
}

void PassivityTorqueFeedback::after(mc_control::MCGlobalController & controller)
{
  mc_rtc::log::info("PassivityTorqueFeedback::after");
} 

void PassivityTorqueFeedback::addGUI(mc_control::MCGlobalController & controller)
{
  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);
  auto gui = ctl.controller().gui();

  gui->addElement({"Plugins","Integral term feedback","Configure"},
    mc_rtc::gui::Checkbox("Is active", this->is_active_)
  );

  gui->addElement({"Plugins","Integral term feedback","Configure"},
    mc_rtc::gui::ComboInput("Integral term type", {"Simple","Filtered"},
      [this]() { return getIntegralTermType(); },
      [this] (const std::string & t) { setIntegralTermType(t); }
    )
  );
  gui->addElement({"Plugins","Integral term feedback","Configure"},
    mc_rtc::gui::ComboInput("Velocity gain type", {"Diagonal","MassDiagonal","MassMatrix"},
      [this]() { return getVelocityGainType(); },
      [this] (const std::string & t) { setVelocityGainType(t); }
    )
  );
  gui->addElement({"Plugins","Integral term feedback","Configure"},
    mc_rtc::gui::NumberInput("Gain", this->lambda_)
  );
  gui->addElement({"Plugins","Integral term feedback","Configure"},
    mc_rtc::gui::NumberInput("Slow filter phi",
      [this]() { return phi_slow_; },
      [this](double phi) {
        phi_slow_ = phi;
        exp_phi_slow_ = exp(-phi*dt_);
      }
    )
  );
  gui->addElement({"Plugins","Integral term feedback","Configure"},
    mc_rtc::gui::NumberInput("Fast filter phi",
      [this]() { return phi_fast_; },
      [this](double phi) {
        phi_fast_ = phi;
        exp_phi_fast_ = exp(-phi*dt_);
      }
    )
  );
  gui->addElement({"Plugins","Integral term feedback","Configure"},
    mc_rtc::gui::NumberInput("Fast filter weight", this->fast_filter_weight_)
  );
  gui->addElement({"Plugins","Integral term feedback","Configure","Bounds"},
    mc_rtc::gui::NumberInput("Lower", this->integral_boundl_)
  );
  gui->addElement({"Plugins","Integral term feedback","Configure","Bounds"},
    mc_rtc::gui::NumberInput("Upper", this->integral_boundu_)
  );

  // Logs
  gui->addElement({"Plugins","Integral term feedback","Logs"},
    mc_rtc::gui::ArrayLabel("Ref velocity",
      ctl.controller().robot().refJointOrder(),
      [this]() { return alpha_r_;}
    )
  );

  gui->addElement({"Plugins","Integral term feedback","Logs"},
    mc_rtc::gui::ArrayLabel("New non filtered integral term",
      ctl.controller().robot().refJointOrder(),
      [this]() { return new_s_;}
    )
  );
  
  gui->addElement({"Plugins","Integral term feedback","Logs"},
    mc_rtc::gui::ArrayLabel("Past non filtered integral term",
      ctl.controller().robot().refJointOrder(),
      [this]() { return prev_s_;}
    )
  );

  gui->addElement({"Plugins","Integral term feedback","Logs"},
    mc_rtc::gui::ArrayLabel("Slow filtered integral term",
      ctl.controller().robot().refJointOrder(),
      [this]() { return slow_filtered_s_;}
    )
  );

  gui->addElement({"Plugins","Integral term feedback","Logs"},
    mc_rtc::gui::ArrayLabel("Fast filtered integral term",
      ctl.controller().robot().refJointOrder(),
      [this]() { return fast_filtered_s_;}
    )
  );

  gui->addElement({"Plugins","Integral term feedback","Logs"},
    mc_rtc::gui::ArrayLabel("Integral term",
      ctl.controller().robot().refJointOrder(),
      [this]() { return s_;}
    )
  );

  gui->addElement({"Plugins","Integral term feedback","Logs"},
    mc_rtc::gui::ArrayLabel("Torque sent",
      ctl.controller().robot().refJointOrder(),
      [this]() { return tau_;}
    )
  );

  // Plots
  gui->addElement({"Plugins","Integral term feedback","Plots"},
    mc_rtc::gui::ElementsStacking::Horizontal,
    mc_rtc::gui::Button("Plot integral term",
      [this, gui]() {
        gui->addPlot("Integral term feedback",
          mc_rtc::gui::plot::X("t", [this]() { return static_cast<double>(count_) * dt_; }),
          mc_rtc::gui::plot::Y("Integral term q(0)",[this]() { return s_[0]; }, mc_rtc::gui::Color::Blue),
          mc_rtc::gui::plot::Y("Integral term q(1)",[this]() { return s_[1]; }, mc_rtc::gui::Color::Red),
          mc_rtc::gui::plot::Y("Integral term q(2)",[this]() { return s_[2]; }, mc_rtc::gui::Color::Green),
          mc_rtc::gui::plot::Y("Integral term q(3)",[this]() { return s_[3]; }, mc_rtc::gui::Color::Yellow),
          mc_rtc::gui::plot::Y("Integral term q(4)",[this]() { return s_[4]; }, mc_rtc::gui::Color::Magenta),
          mc_rtc::gui::plot::Y("Integral term q(5)",[this]() { return s_[5]; }, mc_rtc::gui::Color::Cyan),
          mc_rtc::gui::plot::Y("Integral term q(6)",[this]() { return s_[6]; }, mc_rtc::gui::Color::Black)
        );
      }
    ),
    mc_rtc::gui::Button("Remove integral term plot", [this,gui]() { gui->removePlot("Integral term feedback"); })
  );
}

void PassivityTorqueFeedback::addLOG(mc_control::MCGlobalController & controller)
{
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_active", [&, this]() { return this->is_active_; });

  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_type", [&, this]() { return static_cast<int>(this->integralType_); });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_gain_type", [&, this]() { return static_cast<int>(this->velGainType_); });

  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_bounds_lower", [&, this]() { return this->integral_boundl_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_bounds_up", [&, this]() { return this->integral_boundu_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_gain_scalar", [&, this]() { return this->lambda_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_filter_slow_phi", [&, this]() { return this->phi_slow_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_filter_fast_phi", [&, this]() { return this->phi_fast_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_filter_slow_exp_phi", [&, this]() { return this->exp_phi_slow_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_filter_fast_exp_phi", [&, this]() { return this->exp_phi_fast_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_filter_weight", [&, this]() { return this->fast_filter_weight_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_integral_of_reference_acceleration", [&, this]() { return this->alpha_r_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_K", [&, this]() { return static_cast<Eigen::VectorXd>(this->K_.diagonal()); });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_L", [&, this]() { return static_cast<Eigen::VectorXd>(this->L_.diagonal()); });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_s", [&, this]() { return this->s_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_filter_new_s", [&, this]() { return this->new_s_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_filter_previous_s", [&, this]() { return this->prev_s_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_filter_slow_s", [&, this]() { return this->slow_filtered_s_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_filter_fast_s", [&, this]() { return this->fast_filtered_s_; });
  controller.controller().logger().addLogEntry("PassivityTorqueFeedback_torque", [&, this]() { return this->tau_; });
}

void PassivityTorqueFeedback::removeLOG(mc_control::MCGlobalController & controller)
{
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_active");

  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_type");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_gain_type");

  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_bounds_lower");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_bounds_up");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_gain_scalar");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_filter_slow_phi");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_filter_fast_phi");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_filter_slow_exp_phi");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_filter_fast_exp_phi");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_filter_weight");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_integral_of_reference_acceleration");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_K");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_L");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_s");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_filter_new_s");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_filter_previous_s");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_filter_slow_s");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_filter_fast_s");
  controller.controller().logger().removeLogEntry("PassivityTorqueFeedback_torque");
}

std::string PassivityTorqueFeedback::getIntegralTermType(void)
{
  switch(integralType_)
  {
    case IntegralTermType::Simple:
      return "Simple";
    case IntegralTermType::Filtered:
      return "Filtered";
    default:
      return "";
  }
}

std::string PassivityTorqueFeedback::getVelocityGainType(void)
{
  switch(integralType_)
  {
    case VelocityGainType::Diagonal:
      return "Diagonal";
    case VelocityGainType::MassDiagonal:
      return "MassDiagonal";
    case VelocityGainType::MassMatrix:
      return "MassMatrix";
    default:
      return "";
  }
}

void PassivityTorqueFeedback::setIntegralTermType(std::string type)
{
  if (type.compare("Simple") == 0) integralType_ = IntegralTermType::Simple;
  else if (type.compare("Filtered") == 0) integralType_ = IntegralTermType::Filtered;
}

void PassivityTorqueFeedback::setVelocityGainType(std::string type)
{
  if (type.compare("Diagonal") == 0) velGainType_ = VelocityGainType::Diagonal;
  else if (type.compare("MassDiagonal") == 0) velGainType_ = VelocityGainType::MassDiagonal;
  else if (type.compare("MassMatrix") == 0) velGainType_ = VelocityGainType::MassMatrix;
}

mc_control::GlobalPlugin::GlobalPluginConfiguration PassivityTorqueFeedback::configuration()
{
  mc_control::GlobalPlugin::GlobalPluginConfiguration out;
  out.should_run_before = true;
  out.should_run_after = false;
  out.should_always_run = false;
  return out;
}

} // namespace mc_plugin

EXPORT_MC_RTC_PLUGIN("PassivityTorqueFeedback", mc_plugin::PassivityTorqueFeedback)
