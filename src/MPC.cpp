#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"


using CppAD::AD;

// Timestep length and duration
size_t N = 8;
double dt = 0.3;



// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// Starting indices of the variables in the combined optimization vector
  size_t x_start=0;
  size_t y_start=N;
  size_t psi_start=2*N;
  size_t v_start=3*N; 
  size_t cte_start=4*N;
  size_t epsi_start=5*N;
  size_t delta_start=6*N;
  size_t a_start=delta_start+(N-1);
size_t ref_v=20; // this is a reasonably fast velocity

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // Implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
  
    // Cost function

   fg[0]=0;
	//error from the waypoints
    for (unsigned int t = 0; t < N; t++) {
      fg[0] += 1*CppAD::pow(vars[cte_start + t], 2);
      fg[0] += 1*CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += 0.005*CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (unsigned int t = 0; t < N - 1; t++) {
      fg[0] += 100*CppAD::pow(vars[delta_start + t], 2);
      fg[0] += 0.01*CppAD::pow(vars[a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (unsigned int t = 0; t < N - 2; t++) {
      fg[0] += 10000*CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += 1*CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }
    // The starting values of the simulation need to stay constant
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // Set the kinematic equations as constraints
   for (unsigned int t = 1; t < N; t++) {
  // First calculate the state at time t+1 .
  AD<double> x1 = vars[x_start + t];
  AD<double> y1 = vars[y_start + t];
  AD<double> psi1 = vars[psi_start + t];
  AD<double> v1 = vars[v_start + t];
  AD<double> cte1 = vars[cte_start + t];
  AD<double> epsi1 = vars[epsi_start + t];

  // And the state at time t.
  AD<double> x0 = vars[x_start + t - 1];
  AD<double> y0 = vars[y_start + t - 1];
  AD<double> psi0 = vars[psi_start + t - 1];
  AD<double> v0 = vars[v_start + t - 1];
  AD<double> cte0 = vars[cte_start + t - 1];
  AD<double> epsi0 = vars[epsi_start + t - 1];

  // Only consider the actuation at time t.
  AD<double> delta0 = vars[delta_start + t - 1];
  AD<double> a0 = vars[a_start + t - 1];

  AD<double> f0 = coeffs[0] + coeffs[1] * x0;
  AD<double> psides0 = CppAD::atan(coeffs[1]);

  // Set the constraints 
  fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
  fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
  
  fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
  fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));

  fg[1 + psi_start + t] = psi1 - (psi0 - v0 * delta0 / Lf * dt);
  fg[1 + epsi_start + t] =   epsi1 - ((psi0 - psides0) - v0 * delta0 / Lf * dt);
  }

}
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

  vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  
  // Number of model variables (includes both states and inputs).
   size_t n_vars = 6*N+2*(N-1);
  // The number of constraints
  size_t n_constraints = 6*N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (unsigned int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // Lower and upper limits for variables.
  double degrees_25_torad=0.436332313;
  double inf_bound=1000;
  for (unsigned int t = 1; t < N; t++) {
    vars_lowerbound[x_start+t]=-inf_bound;
    vars_lowerbound[y_start+t]=-inf_bound;
    vars_lowerbound[psi_start+t]=-3.14;
    vars_lowerbound[v_start+t]=-inf_bound;
    vars_lowerbound[cte_start+t]=-inf_bound;
    vars_lowerbound[epsi_start+t]=-inf_bound;

    vars_upperbound[x_start+t]=inf_bound;
    vars_upperbound[y_start+t]=inf_bound;
    vars_upperbound[psi_start+t]=3.14;
    vars_upperbound[v_start+t]=inf_bound;
    vars_upperbound[cte_start+t]=inf_bound;
    vars_upperbound[epsi_start+t]=inf_bound;
  }

  for (unsigned int t = 0; t < N-1; t++) {
    vars_lowerbound[delta_start+t]=-degrees_25_torad;
    vars_lowerbound[a_start+t]=-1;
    vars_upperbound[delta_start+t]=degrees_25_torad;
    vars_upperbound[a_start+t]=1;
  }
  // Fix the initial state values
  vars_lowerbound[x_start]=state[0];
  vars_lowerbound[y_start]=state[1];
  vars_lowerbound[psi_start]=state[2];
  vars_lowerbound[v_start]=state[3];
  vars_lowerbound[cte_start]=state[4];
  vars_lowerbound[epsi_start]=state[5];

  vars_upperbound[x_start]=state[0];
  vars_upperbound[y_start]=state[1];
  vars_upperbound[psi_start]=state[2];
  vars_upperbound[v_start]=state[3];
  vars_upperbound[cte_start]=state[4];
  vars_upperbound[epsi_start]=state[5];

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (unsigned int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start]=state[0];
  constraints_lowerbound[y_start]=state[1];
  constraints_lowerbound[psi_start]=state[2];
  constraints_lowerbound[v_start]=state[3];
  constraints_lowerbound[cte_start]=state[4];
  constraints_lowerbound[epsi_start]=state[5];

  constraints_upperbound[x_start]=state[0];
  constraints_upperbound[y_start]=state[1];
  constraints_upperbound[psi_start]=state[2];
  constraints_upperbound[v_start]=state[3];
  constraints_upperbound[cte_start]=state[4];
  constraints_upperbound[epsi_start]=state[5];

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem

  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  //auto cost = solution.obj_value;
  //std::cout << "Cost " << cost << std::endl;

  // Return the first actuator values. 
  auto x1={solution.x[delta_start],solution.x[a_start]};
  // Return the predicted path to visualize it in the simulator
  this->plottable_trajectory_x = {};
  this->plottable_trajectory_y = {};
  for (unsigned int i = 0; i < N; i++) {
    this->plottable_trajectory_x.push_back(solution.x[x_start + i]);
    this->plottable_trajectory_y.push_back(solution.x[y_start + i]);
  }
  return {x1};
}









