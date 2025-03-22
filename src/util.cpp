#include "util.h"

namespace tinyed
    {

// Function for numerical differentiation
VectorXd numericalGradient(std::function<Real(const std::vector<VectorXd>&)> func, std::vector<VectorXd>& x, Real eps) {
    VectorXd grad(x.size() * x[0].size());
    int dim = x[0].size();
 
    for (size_t i = 0; i < x.size(); ++i) {
        for (int j = 0; j < dim; ++j) {
            x[i][j] += eps;
            Real E_plus = func(x);
            x[i][j] -= 2 * eps;
            Real E_minus = func(x);
            x[i][j] += eps;
            grad[i * dim + j] = (E_plus - E_minus) / (2 * eps);
        }
    }
    return grad;
}

// Function for numerical Jacobian
MatrixXd numericalJacobian(std::function<VectorXd(const std::vector<VectorXd>&)> func, std::vector<VectorXd>& x, Real eps) {
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
            
//            std::cout<<"fplus"<<f_plus<<"\n";
//            std::cout<<"f_minus"<<f_minus<<"\n";
//
//            std::cout<<"huh"<<(f_plus - f_minus)<<"\n";
//            std::cout << "Jacobian column " << (i * dim + j) << ":\n" << J.col(i * dim + j) << "\n\n";
        }
    }
    return J;
}

}
