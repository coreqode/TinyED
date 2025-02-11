#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include <iomanip>
using namespace Eigen;
// Abstract base class for energy models
class EnergyModel {
public:
    virtual double getEnergy(int idx, const std::vector<VectorXd>& x) = 0;
    virtual VectorXd getForce(int idx, const std::vector<VectorXd>& x) = 0;
    virtual MatrixXd getForceJacobian(int idx, const std::vector<VectorXd>& x) = 0;
    virtual ~EnergyModel() {}
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

    double getEnergy(int idx, const std::vector<VectorXd>& x) override {
        const auto& [i, j] = constraints[idx];
        double cur_length = (x[i] - x[j]).norm();
        double rest_length = data[idx].rest_length;
        return 0.5 * data[idx].stiffness * std::pow(cur_length - rest_length, 2);
    }

    VectorXd getForce(int idx, const std::vector<VectorXd>& x) override {
        const auto& [i, j] = constraints[idx];
        VectorXd x21 = x[j] - x[i];
        double cur_length = x21.norm();
        double force_magnitude = data[idx].stiffness * (cur_length - data[idx].rest_length);
        VectorXd force = force_magnitude * (x21 / cur_length);
        VectorXd forces(2 * dim);
        forces << force, -force;
        return forces;
    }

    MatrixXd getForceJacobian(int idx, const std::vector<VectorXd>& x) override {
        const auto& [i, j] = constraints[idx];
        VectorXd x21 = x[j] - x[i];
        double cur_length = x21.norm();
        double rest_length = data[idx].rest_length;
        MatrixXd x21x21Dyadic = x21 * x21.transpose();
        MatrixXd mat = -data[idx].stiffness * (MatrixXd::Identity(dim, dim) - x21x21Dyadic / (cur_length * cur_length)) * (1 - rest_length / cur_length);
        MatrixXd jacobian = MatrixXd::Zero(2 * dim, 2 * dim);
        jacobian.topLeftCorner(dim, dim) = mat;
        jacobian.topRightCorner(dim, dim) = -mat;
        jacobian.bottomLeftCorner(dim, dim) = -mat;
        jacobian.bottomRightCorner(dim, dim) = mat;
        return jacobian;
    }
};



// Function to print a 2D Eigen Matrix
void printMatrix(const MatrixXd& mat, const std::string& label) {
    std::cout << label << ":\n" << mat << "\n\n";
}

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

    double energy = model.getEnergy(0, x);
    VectorXd force = model.getForce(0, x);
    MatrixXd jacobian = model.getForceJacobian(0, x);

    // Reshape the force vector into a 2x2 matrix
    MatrixXd forceMatrix = Map<MatrixXd>(force.data(), 2, 2);

    // Display results in matrix-like format
    std::cout << "\n================ Spring Model Results ================\n";
    std::cout << "Spring 0 Energy: " << std::fixed << std::setprecision(4) << energy << " J\n\n";
    
    printMatrix(forceMatrix, "Force exerted by spring 0");
    printMatrix(jacobian, "Force Jacobian (stiffness matrix) of spring 0");

    return 0;
}

