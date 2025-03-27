#ifndef MATERIALS_H
#define MATERIALS_H

#include "common.h"

namespace TinyED
{
    struct SpringModel{

        static Real getEnergy(const VecS<6> &x, const Real l0, const Real k)
        {
            const Real l = (x.segment<3>(0) - x.segment<3>(3)).norm();
            return 0.5 * k * (l - l0) * (l - l0);
        }

        static void getForce(const VecS<6> &x, const Real l0, const Real k, VecS<6> &force)
        {
        }

        static void getForceJacobian(const VecS<6> &x, const Real l0, const Real k, MatS<6, 6> &jacobian)
        {
        }

        static void getStress(const VecS<6> &x, const Real l0, const Real k, VecS<3> &stress)
        {
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
