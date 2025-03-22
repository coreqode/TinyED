#include "common.h"

namespace tinyed {

class EnergyModel {
    public:
      virtual ~EnergyModel() = default;
      virtual Real getEnergy(int idx, const std::vector<VectorXd>& x) const = 0;
      // virtual void getForce(int idx, const std::vector<VectorXd>& x, VectorXd& force) const = 0;
      virtual VectorXd getForce(int idx, const std::vector<VectorXd>& x) const = 0;
      virtual MatrixXd getForceJacobian(int idx, const std::vector<VectorXd>& x) const = 0;
      int totalConstraints() const;

      Real getTotalEnergy(const std::vector<VectorXd>& x) const {
          Real totalEnergy = 0.0;
          for (int idx = 0; idx < totalConstraints(); ++idx) {
              totalEnergy += getEnergy(idx, x);
          }
          return totalEnergy;
      }


      // TODO: pass by reference for force and jacobians
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
  
}
