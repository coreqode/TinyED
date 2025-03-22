#include <iostream>
#include "common.h"
#include "spring.h"
#include "util.h"

using namespace tinyed;

int main()
{

  std::vector<VectorXd> rest_vertices = {
    VectorXd::Zero(2),
    (VectorXd(2) << 1.0, 0.0).finished()
};

  std::vector<std::pair<int, int>> constraints = {{0, 1}};
  double stiffness = 10.0;
  // SpringModel model(rest_vertices, constraints, stiffness);
  tinyed::SpringModel model(rest_vertices, constraints, stiffness);

  std::vector<VectorXd> x = {
    VectorXd::Zero(2),
    (VectorXd(2) << 1.2, 0.0).finished()
};

  // Analytical computations
  double energy = model.getEnergy(0, x);
  // VectorXd force;  force.setZero();
  auto force = model.getForce(0, x);
  MatrixXd jacobian = model.getForceJacobian(0, x);

  // Numerical tests
  VectorXd numerical_force = -numericalGradient([&](const std::vector<VectorXd>& x) { return model.getEnergy(0, x); }, x);
  MatrixXd numerical_jacobian = numericalJacobian([&](const std::vector<VectorXd>& x) { return model.getForce(0, x); }, x);

  std::cout << "\n================ Derivative Test Results ================\n";
  std::cout << "Analytical Force: \n" << force.transpose() << "\n";
  std::cout << "Numerical Force: \n" << numerical_force.transpose() << "\n\n";
  std::cout << "Analytical Jacobian: \n" << jacobian << "\n";
  std::cout << "Numerical Jacobian: \n" << numerical_jacobian << "\n";

  return 0;
}
