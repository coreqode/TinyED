#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include <iomanip>
// using namespace Eigen;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Map;
#include "functions.h"
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
