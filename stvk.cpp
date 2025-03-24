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

    // std::pair<VectorXd, MatrixXd> getTotalForceAndJacobian(const std::vector<std::vector<int>>& constraints, const std::vector<VectorXd>& x) const {
    //     assert(constraints.size() == static_cast<size_t>(totalConstraints()));
    //     int dim = x[0].size();
    //     VectorXd totalForce = VectorXd::Zero(x.size() * dim);
    //     MatrixXd totalForceJacobian = MatrixXd::Zero(totalForce.size(), totalForce.size());

    //     for (int idx = 0; idx < totalConstraints(); ++idx) {
    //         std::vector<VectorXd> pos;
    //         for (int v : constraints[idx]) {
    //             if (v >= 0) {
    //                 pos.push_back(x[v]);
    //             } else {
    //                 pos.push_back(VectorXd::Zero(dim));
    //             }
    //         }

    //         VectorXd force = getForce(idx, pos);
    //         for (size_t n1 = 0; n1 < constraints[idx].size(); ++n1) {
    //             int v1 = constraints[idx][n1];
    //             if (v1 < 0) continue;
    //             totalForce.segment(dim * v1, dim) += force.segment(dim * n1, dim);
    //         }

    //         MatrixXd forceJacobian = getForceJacobian(idx, pos);
    //         for (size_t n1 = 0; n1 < constraints[idx].size(); ++n1) {
    //             int v1 = constraints[idx][n1];
    //             if (v1 < 0) continue;
    //             for (size_t n2 = 0; n2 < constraints[idx].size(); ++n2) {
    //                 int v2 = constraints[idx][n2];
    //                 if (v2 < 0) continue;
    //                 totalForceJacobian.block(dim * v1, dim * v2, dim, dim) += forceJacobian.block(dim * n1, dim * n2, dim, dim);
    //             }
    //         }
    //     }
    //     return {totalForce, totalForceJacobian};
    // }
};


//fem_model.cpp  
struct FemData {
    double mu_;
    double lambda_;
    MatrixXd invRestMat;
    double restVolume;
};

class FemModel : public EnergyModel {
public:
    FemModel(const std::vector<VectorXd>& rest_vertices,
             const std::vector<std::vector<int>>& constraints,
             double poisson_ratio,
             double youngs_modulus)
        : dim_(rest_vertices[0].size()) {
        
        mu_ = youngs_modulus / (2 * (1 + poisson_ratio));
        lambda_ = (youngs_modulus * poisson_ratio) /
                  ((1 + poisson_ratio) * (1 - 2 * poisson_ratio));
        
        for (const auto& elem : constraints) {
           
            FemData d;
            d.mu_ = mu_;
            d.lambda_ = lambda_;
            d.invRestMat = compute_inv_rest_mat(rest_vertices, elem);
            d.restVolume = nd_volume(rest_vertices, elem);
            data_.push_back(d);
        }
    }

    int totalConstraints() const override {
        return static_cast<int>(data_.size());
    }

    VectorXd compute_dS(int idx, const std::vector<VectorXd>& x) const override {
        return compute_dG(x, data_[idx].invRestMat);
    }

    int nonRigidDofs() const override {
        if (dim_ == 3) {
            return 3;
        } else if (dim_ == 4) {
            return 6;
        } else {
            throw std::invalid_argument("Unsupported element type");
        }
    }

    // Implement other pure virtual functions from EnergyModel
    double getEnergy(int idx, const std::vector<VectorXd>& x) const override {
        // Implementation here
        return 0;
    }

    VectorXd getForce(int idx, const std::vector<VectorXd>& x) const override {
        // Implementation here
        return VectorXd::Zero(x[0].size()); 
    }

    MatrixXd getForceJacobian(int idx, const std::vector<VectorXd>& x) const override {
        // Implementation here
        return MatrixXd::Zero(x[0].size(), x[0].size());
    }

protected:
    int dim_;
    double mu_;
    double lambda_;
    std::vector<FemData> data_;

    
    
    MatrixXd compute_inv_rest_mat(const std::vector<VectorXd>& rest_vertices, const std::vector<int>& elem) const {
        // MatrixXd x_rest=rest_vertices[elem[0]];
        if (elem.size() < 3) {
            std::cerr << "Error: elem size is too small!\n";
            return MatrixXd::Zero(3, 3);  // Or handle appropriately
        }
       
        if (elem.size() == 3) {
            MatrixXd H(3, 2);
            H << -1, -1,
                  1,  0,
                  0,  1;
            MatrixXd X(2, 3);
            X.row(0) = (rest_vertices[elem[1]] - rest_vertices[elem[0]]).transpose();
            X.row(1) = (rest_vertices[elem[2]] - rest_vertices[elem[0]]).transpose();
            return (X * H).inverse();
        } else if (elem.size() == 4) {
            MatrixXd H(4, 3);
            H <<  1,  0,  0,
                  0,  1,  0,
                  0,  0,  1,
                 -1, -1, -1;
            MatrixXd X(3, 4);
            for (int i = 0; i < 3; ++i) {
                X.row(i) = rest_vertices[elem[i]].transpose();
            }
            return (X * H).inverse();
        } else {
            throw std::invalid_argument("Unsupported element size");
        }
    }

    // Compute the volume based on rest_vertices and element

    double nd_volume(const std::vector<VectorXd>& rest_vertices,
                         const std::vector<int>& elem) const {
        
        

            if (elem.size() == 2) {
                return (rest_vertices[elem[0]] - rest_vertices[elem[1]]).norm();
            } else if (elem.size() == 3) {
               
                VectorXd vec1 = rest_vertices[elem[1]] - rest_vertices[elem[0]];
                VectorXd vec2 = rest_vertices[elem[2]] - rest_vertices[elem[0]];


                Eigen::VectorXd cross_product(3);
                cross_product << vec1(1) * vec2(2) - vec1(2) * vec2(1),
                                vec1(2) * vec2(0) - vec1(0) * vec2(2),
                                vec1(0) * vec2(1) - vec1(1) * vec2(0);

                double result = 0.5 * cross_product.norm();
                return result;
            } else if (elem.size() == 4) {
                MatrixXd mat(4, 4);
                for (int i = 0; i < 4; ++i) {
                    mat(i, 0) = rest_vertices[elem[i]][0];
                    mat(i, 1) = rest_vertices[elem[i]][1];
                    mat(i, 2) = rest_vertices[elem[i]][2];
                    mat(i, 3) = 1.0;
                }
                return std::abs(mat.determinant()) / 6.0;
            } else {
                throw std::invalid_argument("Unsupported element type");
            }
}
    

    // Compute dG based on x and invRestMat
    VectorXd compute_dG(const std::vector<VectorXd>& x,
                        const MatrixXd& invRestMat) const {
        int num_nodes = x.size();
        MatrixXd Ds(dim_, num_nodes - 1);

        // Construct the Ds matrix
        for (int i = 1; i < num_nodes; ++i) {
            Ds.col(i - 1) = x[i] - x[0];
        }

        // Compute the deformation gradient F
        MatrixXd F = Ds * invRestMat;

        // Flatten F into a vector dG
        VectorXd dG = Map<VectorXd>(F.data(), F.size());

        return dG;
    }
};



class StVKModel : public FemModel {
    private:
        mutable VectorXd jacobian_x;
        mutable MatrixXd jacobian;
    
    public:
        StVKModel(const std::vector<VectorXd>& rest_vertices, const std::vector<std::vector<int>>& constraints, double poisson_ratio, double youngs_modulus)
            : FemModel(rest_vertices, constraints, poisson_ratio, youngs_modulus) {}
    
        double getEnergy(int idx, const std::vector<Eigen::VectorXd>& x) const override {
                // Determine the total size needed for the flattened vector
                size_t total_size = 0;
                for (const auto& vec : x) {
                    total_size += vec.size();
                }
                // std::cout << "Total size: " << total_size << std::endl;
                // Create a flattened vector and copy data
                Eigen::VectorXd flattened_x(total_size);
                size_t pos = 0;
                for (const auto& vec : x) {
                    flattened_x.segment(pos, vec.size()) = vec;
                    pos += vec.size();
                }
                // std::cout << "Flattened x: " << flattened_x << std::endl;
                // Now pass flattened_x to the energy computation function
                 return FEM::stvkEnergy(flattened_x, data_[idx].invRestMat, data_[idx].restVolume, data_[idx].mu_, data_[idx].lambda_);
        }
            
    
        VectorXd getForce(int idx, const std::vector<VectorXd>& x) const override {
            jacobian_x = Eigen::Map<const VectorXd>(x[0].data(), x.size() * x[0].size());
            // VectorXd force;
            Vectorr<9> force; force.setZero(); 
            Matrixr<9,9> forceJacobians; forceJacobians.setZero(); 
          
            FEM::stvkForceJacobians(jacobian_x, data_[idx].invRestMat, data_[idx].restVolume, data_[idx].mu_, data_[idx].lambda_,force, forceJacobians); 
        
            std::make_tuple(force, forceJacobians); 
            // std::cout << "Force: " << force.transpose() << std::endl;
            return VectorXd(force);
            // Convert force (Vectorr<9>) to VectorXd
            
           
        }
    
        MatrixXd getForceJacobian(int idx, const std::vector<VectorXd>& x) const override {
            jacobian_x = Eigen::Map<const VectorXd>(x[0].data(), x.size() * x[0].size());
            // VectorXd force;
            Vectorr<9> force; force.setZero(); 
            Matrixr<9,9> forceJacobians; forceJacobians.setZero(); 
          
            FEM::stvkForceJacobians(jacobian_x, data_[idx].invRestMat, data_[idx].restVolume, data_[idx].mu_, data_[idx].lambda_,force, forceJacobians); 
        
            std::make_tuple(force, forceJacobians); 
            // std::cout << "Force: " << force.transpose() << std::endl;
            // return VectorXd(force);
            // Convert force (Vectorr<9>) to VectorXd
            assert((Eigen::Map<const VectorXd>(x[0].data(), x.size() * x[0].size()) - jacobian_x).norm() < 1e-9);
            return forceJacobians;
        }
    
       
};

int main() {
    // Define initial vertex positions using VectorXd (dynamic size)
    Eigen::VectorXd v1(3), v2(3), v3(3);
    v1 << 0.0, 0.0, 0.0;
    v2 << 1.0, 0.0, 0.0;
    v3 << 0.5, 1.0, 0.5;

    std::vector<Eigen::VectorXd> vertices = {v1, v2, v3};

    // Constraints for the FEM model (simple triangular element)
    std::vector<std::vector<int>> femConstraints = {{0, 1, 2}};

    // Material properties
    double youngsModulus = 200.0;
    double poissonRatio = 0.3;

    // Initialize StVK Model
    StVKModel stvkModel(vertices, femConstraints, poissonRatio, youngsModulus);

    // Compute energy, force, and Jacobian before deformation
    std::cout << "=== Before Deformation ===" << std::endl;
    std::cout << "Energy: " << stvkModel.getEnergy(0, vertices) << std::endl;
    std::cout << "Force:\n" << stvkModel.getForce(0, vertices) << "\n";
    std::cout << "Jacobian:\n" << stvkModel.getForceJacobian(0, vertices) << "\n";

    // Apply deformation (Modify the existing vertices)
    vertices[1] += Eigen::VectorXd::Map((Eigen::Vector3d(0.3, 0.2, 0.1)).data(), 3);
    vertices[2] += Eigen::VectorXd::Map((Eigen::Vector3d(-0.1, 0.4, -0.2)).data(), 3);

    // Compute energy, force, and Jacobian after deformation
    std::cout << "\n=== After Deformation ===" << std::endl;
    std::cout << "Energy: " << stvkModel.getEnergy(0, vertices) << std::endl;
    std::cout << "Force:\n" << stvkModel.getForce(0, vertices) << "\n";
    std::cout << "Jacobian:\n" << stvkModel.getForceJacobian(0, vertices) << "\n";

    return 0;
}