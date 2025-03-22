#include "common.h"

namespace tinyed
    {
       VectorXd numericalGradient(std::function<Real(const std::vector<VectorXd>&)> func, std::vector<VectorXd>& x, Real eps = 1e-5);
       MatrixXd numericalJacobian(std::function<VectorXd(const std::vector<VectorXd>&)> func, std::vector<VectorXd>& x, Real eps = 1e-5);
}