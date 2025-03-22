#include "spring.h"

namespace tinyed {
  
        Real SpringModel::getEnergy(int idx, const std::vector<VectorXd>& x) const {
            const auto& [i, j] = constraints[idx];
            Real cur_length = (x[i] - x[j]).norm();
            Real rest_length = data[idx].rest_length;
            return 0.5 * data[idx].stiffness * std::pow(cur_length - rest_length, 2);
        }

        // void SpringModel::getForce(int idx, const std::vector<VectorXd>& x, VectorXd& force) const {
        VectorXd SpringModel::getForce(int idx, const std::vector<VectorXd>& x) const {
            const auto& [i, j] = constraints[idx];
            VectorXd x21 = x[j] - x[i];
            Real cur_length = x21.norm();
            Real force_magnitude = data[idx].stiffness * (cur_length - data[idx].rest_length);
            VectorXd force = force_magnitude * (x21 / cur_length);
            VectorXd forces(2 * dim);
            forces << force, -force;
            return forces;
        }

        MatrixXd SpringModel::getForceJacobian(int idx, const std::vector<VectorXd>& x) const {
          const auto& [i, j] = constraints[idx];
          VectorXd x21 = x[j] - x[i];
          Real cur_length = x21.norm();
          Real rest_length = data[idx].rest_length;

          // Normalize x21 before computing the dyadic product
          VectorXd x21_normalized = x21 / cur_length;
          MatrixXd x21x21Dyadic = x21_normalized * x21_normalized.transpose();

          // Compute stiffness matrix following Python’s logic
          MatrixXd mat = -data[idx].stiffness * (
              (1 - rest_length / cur_length) * (MatrixXd::Identity(dim, dim) - x21x21Dyadic)
              + x21x21Dyadic
          );

          // Construct the full Jacobian matrix
          MatrixXd jacobian(2 * dim, 2 * dim);
          jacobian.topLeftCorner(dim, dim) = mat;
          jacobian.topRightCorner(dim, dim) = -mat;
          jacobian.bottomLeftCorner(dim, dim) = -mat;
          jacobian.bottomRightCorner(dim, dim) = mat;

          return jacobian;
      }
}
