#include "common.h"
// #include "hylc/MathematicaDefinitions.h"
// #include "hylc/strain/strain.hpp"
// #include "hylc/vectors.h"
// #include "hylc/hylc_conf.hpp"
#include "Eigen/Dense"
class FEM
{
  static void stvkPk1Stress( const Real mu, const Real lambda, const Matrixr<3, 2>& F, Matrixr<3, 2>& sigma);
  static void stableNeohookeanPk1Stress( const Real mu, const Real lambda, const Matrixr<3, 3>& F, Matrixr<3, 3>& sigma);
  static void neohookeanPk1Stress( const Real mu, const Real lambda, const Matrixr<3, 3>& F, Matrixr<3, 3>& sigma);

    static void femStretchForceTriangles(
        const Vectorr<9>& positions,
        const Matrixr<2,2>& invRestMat,
        const Real restVolume,
        const Real mu, const Real lambda,
        Vectorr<9>& force);

    static void femStretchForceTets(
        const Vectorr<12>& positions,
        const Matrixr<3,3>& invRestMat,
        const Real restVolume,
        const Real mu, const Real lambda,
        Vectorr<12>& force); 

  public:
    static void CalculateDeformationGradientTriangles(const Vectorr<9>& positions, const Matrixr<2,2>& invRestMat,  Matrixr<3,2>& deformedCoordsToDeformationGradient, Matrixr<3, 2>& F);
    static void CalculateDeformationGradientTets(const Vectorr<12>& positions, const Matrixr<3,3>& invRestMat,  Matrixr<4,3>& deformedCoordsToDeformationGradient, Matrixr<3, 3>& F);

    static Real stvkEnergy(
        const Vectorr<9>& positions,
        const Matrixr<2,2>& invRestMat,
        const Real restVolume,
        const Real mu, const Real lambda);

    static void stvkForceJacobians(const Vectorr<9>& positions,
        const Matrixr<2,2>& invRestMat,
        const Real volume,
        const Real mu, const Real lambda,
        Vectorr<9> &force,
        Matrixr<9,9>& forcesPositionDeriv);


    static Real stableNeohookeanEnergy( const Vectorr<12>& positions, const Matrixr<3, 3>& invRestMat, const Real restVolume, const Real mu, const Real lambda); 
    static void stableNeohookeanForceJacobians(const Vectorr<12>& positions, const Matrixr<3, 3>& invRestMat, const Real volume, const Real mu, const Real lambda, const bool get_force_jacobian, Vector12r &force, Matrix1212r& forcesPositionDeriv); 

    static Real neohookeanEnergy( const Vectorr<12>& positions, const Matrixr<3, 3>& invRestMat, const Real restVolume, const Real mu, const Real lambda); 
    static void neohookeanForceJacobians(const Vectorr<12>& positions, const Matrixr<3, 3>& invRestMat, const Real volume, const Real mu, const Real lambda, const bool get_force_jacobian, Vector12r& force, Matrix1212r& forcesPositionDeriv); 

};

// class HYLCMaterial
// 	{
// 		public:
// 			static Real calculateEnergy(
// 				const Vectorr<18> &positions,
// 				const bool nn0_exists, const bool nn1_exists, const bool nn2_exists,
// 				const Matrixr<2,2>& invRestMat,
// 				const Real restVolume,
// 				const hylc::Config *hylc_config
// 				);

// 			static void calculateForceJacobian(
// 				const Vectorr<18> &positions,
// 				const bool nn0_exists, const bool nn1_exists, const bool nn2_exists,
// 				const Matrix2r& invRestMat,
// 				const Real restVolume,
// 				hylc::Config *hylc_config,
// 				Vector18r& phi,
// 				Matrixr<18,18>& jacobianFull
// 				);
// 	};
