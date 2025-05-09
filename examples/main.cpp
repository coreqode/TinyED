#include <iostream>
#include "materials.h"

using namespace TinyED;

// int main()
// {
//     // Define the positions of two points in 3D space (flattened into VecS<6>)
//     VecS<6> x;
//     x << 2, 0,0,   // First point (x1, y1, z1)
//          0, 0, 0;   // Second point (x2, y2, z2)

//     Real l0 = 1.0;  // Rest length
//     Real k = 100.0; // Stiffness constant

//     // Compute energy
//     const Real energy = SpringModel::getEnergy(x, l0, k);
//     std::cout << "Computed energy: " << energy << std::endl;

//     // Compute force
//     VecS<6> force;
//     SpringModel::getForce(x, l0, k, force);
//     std::cout << "Computed force: " << force.transpose() << std::endl;

//     // Compute force Jacobian
//     MatS<6, 6> jacobian;
//     SpringModel::getForceJacobian(x, l0, k, jacobian);
//     std::cout << "Computed force Jacobian:\n" << jacobian << std::endl;

//     // Compute stress
//     VecS<3> stress;
//     SpringModel::getStress(x, l0, k, stress);
//     std::cout << "Computed stress: " << stress.transpose() << std::endl;
//     //compute DF/dx
//     MatS<3, 3> dFdx;
//     SpringModel::get_dFdx(x, l0, k, dFdx);
//     std::cout << "Computed dF/dx:\n" << dFdx << std::endl;

  
//     return 0;
// }

int main() {
    // Define rest state (undeformed tetrahedron)
    VecS<12> rest_vertices;
    rest_vertices << 0, 0, 0,   // Vertex 0
                     1, 0, 0,   // Vertex 1
                     0, 1, 0,   // Vertex 2
                     0, 0, 1;   // Vertex 3

    // Define deformed state (stretched tetrahedron)
    VecS<12> deformed_vertices;
    deformed_vertices  << 0, 0, 0,   // Vertex 0
                     1.1, 0, 0,   // Vertex 1
                     0, 1.0, 0,   // Vertex 2
                     0, 0, 1.1; 

    // Material properties
    Real lambda = 1.0;
    Real mu = 0.5;
    Real volume = 1.0 / 6.0;  // Tetrahedron rest volume (simple unit tetrahedron)

    // Compute inverse rest shape matrix
    MatS<3, 3> invRestMat;
    computeInvRestMat<4>(rest_vertices, invRestMat);
    std::cout << "Inverse Rest Shape Matrix:\n" << invRestMat << std::endl;

    // Compute energy
    Real energy = NeoHookeanModel::getEnergy(deformed_vertices, invRestMat, volume, mu, lambda);
    std::cout << "Neo-Hookean Energy: " << energy << std::endl;

    // // Compute deformation gradient
    MatS<3, 3> F;
    // computeDeformationGradient<4>(deformed_vertices, invRestMat, dXInv, F);

    // // Compute stress
    MatS<3, 3> sigma;
    NeoHookeanModel::getStress(deformed_vertices, invRestMat, mu, lambda,sigma);
    std::cout << "Neo-Hookean Stress:\n" << sigma << std::endl;
   
    // Compute force
    VecS<12> force;
    NeoHookeanModel::getForce<4>(deformed_vertices, invRestMat, volume, mu, lambda, force);
    std::cout << "Neo-Hookean Force Vector:\n" << force.transpose() << std::endl;

    // Compute force Jacobian
    MatS<12, 12> J;
    NeoHookeanModel::getForceJacobian(deformed_vertices, invRestMat, volume, mu, lambda, force, J);
    std::cout << "Neo-Hookean Force Jacobian:\n" << J << std::endl;

    return 0;
}
