#ifndef HELPERS_H
#define HELPERS_H
#include <iostream>
#include "common.h"

namespace TinyED
{
    void computeInvRestMat(const VecS<9>& rest_vertices , MatS<2, 2>& invRestMat) {
    //2d triangle
        // Define transformation matrix H
        MatS<3, 2> H;
        H << -1, -1,
              1,  0,
              0,  1;
    
        // Compute rest shape matrix X
        MatS<3, 2> X;
        if (rest_vertices[1] == rest_vertices[4] && rest_vertices[4] == rest_vertices[7]) {
            // All y-coordinates are the same → XZ plane
            X << rest_vertices[0], rest_vertices[2],  // Vertex 0 (x, z)
                 rest_vertices[3], rest_vertices[5],  // Vertex 1 (x, z)
                 rest_vertices[6], rest_vertices[8];  // Vertex 2 (x, z)
        } else if (rest_vertices[2] == rest_vertices[5] && rest_vertices[5] == rest_vertices[8]) {
            // All z-coordinates are the same → XY plane
            X << rest_vertices[0], rest_vertices[1],  // Vertex 0 (x, y)
                 rest_vertices[3], rest_vertices[4],  // Vertex 1 (x, y)
                 rest_vertices[6], rest_vertices[7];  // Vertex 2 (x, y)
        } else {
            throw std::runtime_error("Only XY and XZ rest alignment are supported");
        }
        // Compute inverse of X * H
        invRestMat =  (X.transpose() * H).inverse();

        //print X and H too
        // std::cout << "Rest shape matrix X:\n" << X << std::endl;
        // std::cout << "Transformation matrix H:\n" << H << std::endl;
       
    }
    

    void computeDeformationGradient(
        const VecS<9>& positions, 
        const MatS<2,2>& invRestMat,
        MatS<3,2>& B,
        MatS<3,2>& F)
    {
        const Mat32 H{
            {-1, -1},
            {1, 0},
            {0, 1}
        };

        B = H * invRestMat;

        MatS<3,3> positionsMat;
    
        for(unsigned int p = 0; p < 3; p++) {
            positionsMat.col(p) = positions.segment(p * 3, 3);
        }
    
        F = positionsMat * B;
    }


}

#endif //HELPERS_H
