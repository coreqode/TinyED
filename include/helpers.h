#ifndef HELPERS_H
#define HELPERS_H
#include "common.h"
#include <iostream>
namespace TinyED
{
    // void computeInvRestMat(const VecS<9>& rest_vertices , MatS<2, 2>& invRestMat) {
    // //2d triangle
    //     // Define transformation matrix H
        
    //     MatS<3, 2> H;
    //     H << -1, -1,
    //           1,  0,
    //           0,  1;
    
    //     // Compute rest shape matrix X
    //     MatS<3, 2> X;
    //     if (rest_vertices[1] == rest_vertices[4] && rest_vertices[4] == rest_vertices[7]) {
    //         // All y-coordinates are the same → XZ plane
    //         X << rest_vertices[0], rest_vertices[2],  // Vertex 0 (x, z)
    //              rest_vertices[3], rest_vertices[5],  // Vertex 1 (x, z)
    //              rest_vertices[6], rest_vertices[8];  // Vertex 2 (x, z)
    //     } else if (rest_vertices[2] == rest_vertices[5] && rest_vertices[5] == rest_vertices[8]) {
    //         // All z-coordinates are the same → XY plane
    //         X << rest_vertices[0], rest_vertices[1],  // Vertex 0 (x, y)
    //              rest_vertices[3], rest_vertices[4],  // Vertex 1 (x, y)
    //              rest_vertices[6], rest_vertices[7];  // Vertex 2 (x, y)
    //     } else {
    //         throw std::runtime_error("Only XY and XZ rest alignment are supported");
    //     }
    //     // Compute inverse of X * H
    //     invRestMat =  (X.transpose() * H).inverse();

    //     //print X and H too
    //     std::cout << "Rest shape matrix X:\n" << X << std::endl;
    //     std::cout << "Transformation matrix H:\n" << H << std::endl;
       
    // }
template <unsigned int N>
void computeInvRestMat(const VecS< N * 3>& rest_vertices, MatS<N - 1, N - 1>& invRestMat) {
    static_assert(N == 3 || N == 4, "Only triangles (N=3) and tetrahedra (N=4) are supported");

    MatS<N, N - 1> H;
    MatS<N, N - 1> X;

    if constexpr (N == 3) {
        // H for triangle
        H << -1, -1,
              1,  0,
              0,  1;

        // Determine the 2D plane (XY or XZ)
        if (rest_vertices[1] == rest_vertices[4] && rest_vertices[4] == rest_vertices[7]) {
            // All y are same → XZ plane
            X << rest_vertices[0], rest_vertices[2],
                 rest_vertices[3], rest_vertices[5],
                 rest_vertices[6], rest_vertices[8];
        } else if (rest_vertices[2] == rest_vertices[5] && rest_vertices[5] == rest_vertices[8]) {
            // All z are same → XY plane
            X << rest_vertices[0], rest_vertices[1],
                 rest_vertices[3], rest_vertices[4],
                 rest_vertices[6], rest_vertices[7];
        } else {
            throw std::runtime_error("Only XY or XZ aligned triangle rest shapes are supported.");
        invRestMat = (X.transpose() * H).inverse();
        }
    }
     else if constexpr (N == 4) {
        // H for tetrahedron
        H << 1, 0, 0,
             0, 1, 0,
             0, 0, 1,
            -1, -1, -1;

        // X is transpose of rest_vertices
        MatS<3, 4> X;
        X << rest_vertices[0],  rest_vertices[3],  rest_vertices[6],  rest_vertices[9],
        rest_vertices[1],  rest_vertices[4],  rest_vertices[7],  rest_vertices[10],
        rest_vertices[2],  rest_vertices[5],  rest_vertices[8],  rest_vertices[11];
        
       
        invRestMat = (X*H).inverse();}
}


template <unsigned int N>
void computeDeformationGradient(
    const VecS<N * 3>& positions,
    const MatS<N - 1, N - 1>& invRestMat,
    MatS<N, N - 1>& deformedCoordsToDeformationGradient,
    MatS<3, N - 1>& F)
{
    MatS<N, N - 1> B;

    if constexpr (N == 3) {
        // Triangle: 3 nodes, 2D simplex
        B << -1, -1,
              1,  0,
              0,  1;
    }
    else if constexpr (N == 4) {
        // Tetrahedron: 4 nodes, 3D simplex
        B << 1, 0, 0,
				0, 1, 0,
				0, 0, 1,
				-1, -1, -1;
    }
    else {
        static_assert(N == 3 || N == 4, "Only triangles (N=3) and tetrahedra (N=4) are supported.");
    }

    deformedCoordsToDeformationGradient = B * invRestMat;

    MatS<3, N> positionsMat;
    // for (unsigned int p = 0; p < N; ++p) {
    //     positionsMat.col(p) = positions.template segment<3>(p * 3);
    // }
    
	for(unsigned int p = 0; p < N; p++)
	{
		positionsMat.col(p) = positions.segment(p*3, 3);
	}
    std::cout<<"positionsMat: " << positionsMat << std::endl;
    F = positionsMat * deformedCoordsToDeformationGradient;
  
}

    void computeAllForces()
    {

    }

    void computeAllForceJacobians()
    {

    }


}

#endif //HELPERS_H
