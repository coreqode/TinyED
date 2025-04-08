#ifndef MATERIALS_H
#define MATERIALS_H

#include "common.h"
#include "helpers.h"
namespace TinyED
{
    struct SpringModel {
        static Real getEnergy(const VecS<6> &x, const Real rest_length, const Real k) {
            const Real l = (x.segment<3>(0) - x.segment<3>(3)).norm();
            return 0.5 * k * (l - rest_length) * (l - rest_length);
        }
    
        static void getForce(const VecS<6> &x, const Real l0, const Real k, VecS<6> &force) {
            VecS<3> x1 = x.segment<3>(0);
            VecS<3> x2 = x.segment<3>(3);
            VecS<3> x21 = x2 - x1;
            Real l = x21.norm();
    
            if (l < 1e-8) {
                force.setZero();
              
            }
    
            Real force_magnitude = k * (l - l0);
            VecS<3> forceVec = force_magnitude * (x21 / l);
    
            force.segment<3>(0) = -forceVec; // Force on first point
            force.segment<3>(3) = forceVec;  // Force on second point
        }
    
        static void get_dFdx(const VecS<6> &x, const Real l0, const Real k, MatS<3, 3> &dFdx) {
            VecS<3> x1 = x.segment<3>(0);
            VecS<3> x2 = x.segment<3>(3);
            VecS<3> x21 = x2 - x1;
            Real l = x21.norm();
    
            if (l < 1e-8) {
                dFdx.setZero();
                return;
            }
    
            MatS<3, 3> I = MatS<3, 3>::Identity();
            MatS<3, 3> x_hat_x_hatT = (x21 * x21.transpose()) / (l * l);
            dFdx = -k * ((1 - (l0 / l)) * (I - x_hat_x_hatT) + x_hat_x_hatT);
        }
    
        static void getForceJacobian(const VecS<6> &x, const Real l0, const Real k, MatS<6, 6> &jacobian) {
            MatS<3, 3> dFdx;
            get_dFdx(x, l0, k, dFdx);
    
            jacobian.block<3, 3>(0, 0) = dFdx;
            jacobian.block<3, 3>(3, 3) = dFdx;
            jacobian.block<3, 3>(0, 3) = -dFdx;
            jacobian.block<3, 3>(3, 0) = -dFdx;
        }
    
        static void getStress(const VecS<6> &x, const Real l0, const Real k, VecS<3> &stress) {
            VecS<3> x1 = x.segment<3>(0);
            VecS<3> x2 = x.segment<3>(3);
            VecS<3> x21 = x2 - x1;
            Real l = x21.norm();
    
            if (l < 1e-8) {
                stress.setZero();
               
            }
    
            Real strain = (l - l0) / l0;
            Real stress_magnitude = k * strain;
            stress = stress_magnitude * (x21 / l);
        }
    };
    

    struct StVKModel
    {
        static Real getEnergy(
            const VecS<9>& positions,
            const MatS<2,2>& invRestMat,
            const Real restVolume,
            const Real mu, const Real lambda) 
        {
            std::cout << "Calculating Energy...\n";
        
            // Compute Deformation Gradient
            MatS<3,2> F;
            MatS<3,2> deformedCoordsToDeformationGradient;
            computeDeformationGradient(positions, invRestMat, deformedCoordsToDeformationGradient, F);
         
        
            // Compute Green Strain
            MatS<2,2> greenStrain = (F.transpose() * F - MatS<2,2>::Identity()) * 0.5;
            std::cout << "Green Strain:\n" << greenStrain << std::endl;
        
            // Compute Energy
            Real trace = greenStrain.trace();
            Real energy = restVolume * (mu * greenStrain.squaredNorm() + lambda * 0.5 * trace * trace);
        
            return energy;
        }
      
        static void get_stvk_stress(
                const Real mu, const Real lambda, const MatS<3, 2>& F, MatS<3, 2>& sigma)
        {
            MatS<2,2> greenStrain = (F.transpose() * F - MatS<2, 2>::Identity()) * 0.5;
            // St. Venant-Kirchhoff model
            sigma = F * (2 * mu * greenStrain + lambda * greenStrain.trace() * MatS<2,2>::Identity());
        }
        static void get_stvk_Force(
            const VecS<9>& x,
                const MatS<2,2>& invRestMat,
                const Real restVolume,
                const Real mu, const Real lambda,
                
                VecS<9>& force)
        {
            MatS<3, 2> sigma;

            MatS<3,2> deformedCoordsToDeformationGradient;
            MatS<3,2> F;

            computeDeformationGradient(x, invRestMat, deformedCoordsToDeformationGradient, F);

            get_stvk_stress(mu, lambda, F, sigma);
            
            auto force_map = Eigen::Map<MatS<3, 3>>(reinterpret_cast<Real*>(force.data())); 
            force_map = sigma * deformedCoordsToDeformationGradient.transpose() * restVolume;
            force = -force;
        }


        // static void getstvkForceJacobian(const VecS<9>& x, const Real area, const Real lambda, const Real mu, MatS<9, 9>&
        // force)
        // {

        // }
        static void getstvkForceJacobians(
            const VecS<9>& positions, 
            const MatS<2, 2>& invRestMat, 
            const Real volume, 
            const Real mu, 
            const Real lambda, 
            VecS<9>& force, 
            MatS<9, 9>& forcesPositionDeriv)
        {
            get_stvk_Force(positions, invRestMat, volume, mu, lambda,  force);
        
            MatS<3, 2> deformedCoordsToDeformationGradient;
            MatS<3, 2> F;
            computeDeformationGradient(positions, invRestMat, deformedCoordsToDeformationGradient, F);
            const MatS<2, 3> BT = deformedCoordsToDeformationGradient.transpose();
            const MatS<2, 2> greenStrain = (F.transpose() * F - MatS<2, 2>::Identity()) * 0.5;
        
            const Real greenStrainTrace = greenStrain.trace();
            const Real lambdaGreenTrace = lambda * greenStrainTrace;
        
            for (int s = 0; s < 2; s++) {
                const Eigen::Block<MatS<3, 2>, 3, 1, true> delF_rsdelx_j = deformedCoordsToDeformationGradient.col(s);
        
                for (int l = 0; l < 2; l++) {
                    Real delP_ildelF_rs_p1_cached = 2 * mu * greenStrain(s, l);
                    if (s == l) {
                        delP_ildelF_rs_p1_cached += lambdaGreenTrace;
                    }
        
                    for (int r = 0; r < 3; r++) {
                        Real deltrEdelF_rs = F(r, s);
                        MatS<2, 1> delE_kldelF_rs_p1 = MatS<2, 1>::Zero();
        
                        if (l == s) {
                            for (int k = 0; k < 2; k++) {
                                delE_kldelF_rs_p1(k) = F(r, k);
                            }
                        }
        
                        for (int i = 0; i < 3; i++) {
                            Real delP_ildelF_rs_p1 = (i == r) ? delP_ildelF_rs_p1_cached : 0.0;
                            Real delP_ildelF_rs_p2 = 0.0;
        
                            if (l == s) {
                                delP_ildelF_rs_p2 += mu * F.row(i).dot(delE_kldelF_rs_p1);
                            }
        
                            delP_ildelF_rs_p2 += F(i, s) * mu * F(r, l);
                            delP_ildelF_rs_p2 += F(i, l) * lambda * deltrEdelF_rs;
        
                            Real delP_ildelF_rs = delP_ildelF_rs_p1 + delP_ildelF_rs_p2;
        
                            for (int x = 0; x < 3; x++) {
                                for (int y = 0; y < 3; y++) {
                                    forcesPositionDeriv(i + 3 * x, r + 3 * y) += 
                                        -volume * delP_ildelF_rs * delF_rsdelx_j(y) * deformedCoordsToDeformationGradient(x, l);
                                }
                            }
                        }
                    }
                }
            }
        }
        
      
};       
     struct NeoHookeanModel
    {
        static Real getEnergy(const VecS<12>& x, const Real area, const Real lambda, const Real mu)
        {
            return 0;
        }

        static void getForce(const VecS<12>& x, const Real area, const Real lambda, const Real mu, VecS<12>& force)
        {

        }

        static void getForceJacobian(const VecS<12>& x, const Real area, const Real lambda, const Real mu, MatS<12, 12>&
        force)
        {

        }
    };

}

#endif //MATERIALS_H






