#ifndef MATERIALS_H
#define MATERIALS_H

#include "common.h"

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
        static Real getEnergy(const VecS<9>& x, const Real area, const Real lambda, const Real mu)
        {
            return 0;
        }

        static void getForce(const VecS<9>& x, const Real area, const Real lambda, const Real mu, VecS<9>& force)
        {

        }

        static void getForceJacobian(const VecS<9>& x, const Real area, const Real lambda, const Real mu, MatS<9, 9>&
        force)
        {

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
