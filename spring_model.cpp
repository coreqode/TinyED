#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include <iomanip>
// using namespace Eigen;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Map;
class EnergyModel {
public:
    virtual ~EnergyModel() = default;

    virtual double getEnergy(int idx, const std::vector<VectorXd>& x) const = 0;
    virtual VectorXd getForce(int idx, const std::vector<VectorXd>& x) const = 0;
    virtual MatrixXd getForceJacobian(int idx, const std::vector<VectorXd>& x) const = 0;
    virtual VectorXd compute_dS(int idx, const std::vector<VectorXd>& x) const = 0;
    virtual int nonRigidDofs() const = 0;
    virtual int totalConstraints() const = 0;

    double getTotalEnergy(const std::vector<VectorXd>& x) const {
        double totalEnergy = 0.0;
        for (int idx = 0; idx < totalConstraints(); ++idx) {
            totalEnergy += getEnergy(idx, x);
        }
        return totalEnergy;
    }

    VectorXd getTotalForce(const std::vector<std::vector<int>>& constraints, const std::vector<VectorXd>& x) const {
        assert(constraints.size() == static_cast<size_t>(totalConstraints()));
        int dim = x[0].size();
        VectorXd totalForce = VectorXd::Zero(x.size() * dim);

        for (int idx = 0; idx < totalConstraints(); ++idx) {
            VectorXd force = getForce(idx, x);
            for (size_t n1 = 0; n1 < constraints[idx].size(); ++n1) {
                int v1 = constraints[idx][n1];
                totalForce.segment(dim * v1, dim) += force.segment(dim * n1, dim);
            }
        }
        return totalForce;
    }

    std::pair<VectorXd, MatrixXd> getTotalForceAndJacobian(const std::vector<std::vector<int>>& constraints, const std::vector<VectorXd>& x) const {
        assert(constraints.size() == static_cast<size_t>(totalConstraints()));
        int dim = x[0].size();
        VectorXd totalForce = VectorXd::Zero(x.size() * dim);
        MatrixXd totalForceJacobian = MatrixXd::Zero(totalForce.size(), totalForce.size());

        for (int idx = 0; idx < totalConstraints(); ++idx) {
            std::vector<VectorXd> pos;
            for (int v : constraints[idx]) {
                if (v >= 0) {
                    pos.push_back(x[v]);
                } else {
                    pos.push_back(VectorXd::Zero(dim));
                }
            }

            VectorXd force = getForce(idx, pos);
            for (size_t n1 = 0; n1 < constraints[idx].size(); ++n1) {
                int v1 = constraints[idx][n1];
                if (v1 < 0) continue;
                totalForce.segment(dim * v1, dim) += force.segment(dim * n1, dim);
            }

            MatrixXd forceJacobian = getForceJacobian(idx, pos);
            for (size_t n1 = 0; n1 < constraints[idx].size(); ++n1) {
                int v1 = constraints[idx][n1];
                if (v1 < 0) continue;
                for (size_t n2 = 0; n2 < constraints[idx].size(); ++n2) {
                    int v2 = constraints[idx][n2];
                    if (v2 < 0) continue;
                    totalForceJacobian.block(dim * v1, dim * v2, dim, dim) += forceJacobian.block(dim * n1, dim * n2, dim, dim);
                }
            }
        }
        return {totalForce, totalForceJacobian};
    }
};
struct SpringData {
  double rest_length;
  double stiffness;
};
class SpringModel : public EnergyModel {
  private:
      std::vector<SpringData> data;
      std::vector<std::pair<int, int>> constraints;
      int dim;
  
  public:
      SpringModel(const std::vector<VectorXd>& rest_vertices, const std::vector<std::pair<int, int>>& constraints, double k)
          : constraints(constraints), dim(rest_vertices[0].size()) {
          for (const auto& edge : constraints) {
              double rest_length = (rest_vertices[edge.first] - rest_vertices[edge.second]).norm();
              data.push_back({rest_length, k});
          }
      }
  
      double getEnergy(int idx, const std::vector<VectorXd>& x) const override {
          const auto& [i, j] = constraints[idx];
          double cur_length = (x[i] - x[j]).norm();
          double rest_length = data[idx].rest_length;
          return 0.5 * data[idx].stiffness * std::pow(cur_length - rest_length, 2);
      }
  
      VectorXd getForce(int idx, const std::vector<VectorXd>& x) const override {
          const auto& [i, j] = constraints[idx];
          VectorXd x21 = x[j] - x[i];
          double cur_length = x21.norm();
          double force_magnitude = data[idx].stiffness * (cur_length - data[idx].rest_length);
          VectorXd force = force_magnitude * (x21 / cur_length);
          VectorXd forces(2 * dim);
          forces << force, -force;
          return forces;
      }
  
      MatrixXd getForceJacobian(int idx, const std::vector<VectorXd>& x) const override {
        const auto& [i, j] = constraints[idx];
        VectorXd x21 = x[j] - x[i];
        double cur_length = x21.norm();
        double rest_length = data[idx].rest_length;
        
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
    
      
    VectorXd compute_dS(int idx, const std::vector<VectorXd>& x) const override {
        const auto& [i, j] = constraints[idx];
        VectorXd x21 = x[i] - x[j];
        double cur_length = x21.norm();
        VectorXd x21_cap = x21 / cur_length;
        VectorXd dS(2 * dim);
        dS << x21_cap, -x21_cap;
        return dS;
    }

    // Return the number of non-rigid degrees of freedom
    int nonRigidDofs() const override {
 
        return 1;
    }
  
    int totalConstraints() const override {
          return static_cast<int>(constraints.size());
      }
  };

// Function to print a 2D Eigen Matrix
void printMatrix(const MatrixXd& mat, const std::string& label) {
  std::cout << label << ":\n" << mat << "\n\n";
}
// Function for numerical differentiation
VectorXd numericalGradient(std::function<double(const std::vector<VectorXd>&)> func, std::vector<VectorXd>& x, double eps = 1e-5) {
    VectorXd grad(x.size() * x[0].size());
    int dim = x[0].size();
 
    for (size_t i = 0; i < x.size(); ++i) {
        for (int j = 0; j < dim; ++j) {
            x[i][j] += eps;
            double E_plus = func(x);
            x[i][j] -= 2 * eps;
            double E_minus = func(x);
            x[i][j] += eps;
            grad[i * dim + j] = (E_plus - E_minus) / (2 * eps);
        }
    }
    return grad;
}

// Function for numerical Jacobian
MatrixXd numericalJacobian(std::function<VectorXd(const std::vector<VectorXd>&)> func, std::vector<VectorXd>& x, double eps = 1e-5) {
    int dim = x[0].size();
    VectorXd f0 = func(x);
    MatrixXd J(f0.size(), f0.size());
    
    for (size_t i = 0; i < x.size(); ++i) {
        for (int j = 0; j < dim; ++j) {
            x[i][j] += eps;
            VectorXd f_plus = func(x);
            
            x[i][j] -= 2 * eps;
            VectorXd f_minus = func(x);
            x[i][j] += eps;
            J.col(i * dim + j) = (f_plus - f_minus) / (2 * eps);
            
            std::cout<<"fplus"<<f_plus<<"\n";
            std::cout<<"f_minus"<<f_minus<<"\n";

            std::cout<<"huh"<<(f_plus - f_minus)<<"\n";
            std::cout << "Jacobian column " << (i * dim + j) << ":\n" << J.col(i * dim + j) << "\n\n";
        }
    }
    return J;
}
// int main() {
//   std::vector<VectorXd> rest_vertices = {
//       VectorXd::Zero(2),
//       (VectorXd(2) << 1.0, 0.0).finished()
//   };

//   std::vector<std::pair<int, int>> constraints = {{0, 1}};
//   double stiffness = 10.0;
//   SpringModel model(rest_vertices, constraints, stiffness);

//   std::vector<VectorXd> x = {
//       VectorXd::Zero(2),
//       (VectorXd(2) << 1.2, 0.0).finished()
//   };

//   double energy = model.getEnergy(0, x);
//   VectorXd force = model.getForce(0, x);
//   MatrixXd jacobian = model.getForceJacobian(0, x);

//   // Reshape the force vector into a 2x2 matrix
//   MatrixXd forceMatrix = Map<MatrixXd>(force.data(), 2, 2);

//   // Display results in matrix-like format
//   std::cout << "\n================ Spring Model Results ================\n";
//   std::cout << "Spring 0 Energy: " << std::fixed << std::setprecision(4) << energy << " J\n\n";
  
//   printMatrix(forceMatrix, "Force exerted by spring 0");
//   printMatrix(jacobian, "Force Jacobian (stiffness matrix) of spring 0");

//   return 0;
// }
int main() {
    std::vector<VectorXd> rest_vertices = {
        VectorXd::Zero(2),
        (VectorXd(2) << 1.0, 0.0).finished()
    };

    std::vector<std::pair<int, int>> constraints = {{0, 1}};
    double stiffness = 10.0;
    SpringModel model(rest_vertices, constraints, stiffness);

    std::vector<VectorXd> x = {
        VectorXd::Zero(2),
        (VectorXd(2) << 1.2, 0.0).finished()
    };

    // Analytical computations
    double energy = model.getEnergy(0, x);
    VectorXd force = model.getForce(0, x);
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