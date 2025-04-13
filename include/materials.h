#ifndef MATERIALS_H
#define MATERIALS_H

#include "common.h"
#include "helpers.h"
#include <variant>

namespace TinyED
{
    struct SpringData
    {
        Real rest_length;
        Real k;
    };

    struct FEMTriData
    {
        Real mu;
        Real lambda;
        Real restVolume;
        Mat2 invRestMat;
    };

    struct SpringModel {
        static inline Real getEnergy(const VecS<6> &x, const SpringData& data) {
            const Real l = (x.segment<3>(0) - x.segment<3>(3)).norm();
            return 0.5 * data.k * (l - data.rest_length) * (l - data.rest_length);
        }
    
        static inline void getForce(const VecS<6> &x, const SpringData& data, VecS<6> &force) {
            VecS<3> x1 = x.segment<3>(0);
            VecS<3> x2 = x.segment<3>(3);
            VecS<3> x21 = x2 - x1;
            Real l = x21.norm();
    
            if (l < 1e-8) {
                force.setZero();
            }
    
            Real force_magnitude = data.k * (l - data.rest_length);
            VecS<3> forceVec = force_magnitude * (x21 / l);
    
            force.segment<3>(0) = -forceVec; // Force on first point
            force.segment<3>(3) = forceVec;  // Force on second point
        }
    
        static inline void get_dFdx(const VecS<6> &x, const SpringData& data, MatS<3, 3> &dFdx) {
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
            dFdx = -data.k * ((1 - (data.rest_length / l)) * (I - x_hat_x_hatT) + x_hat_x_hatT);
        }
    
        static inline void getForceJacobian(const VecS<6> &x, const SpringData& data, MatS<6, 6> &jacobian) {
            MatS<3, 3> dFdx;
            get_dFdx(x, data, dFdx);
    
            jacobian.block<3, 3>(0, 0) = dFdx;
            jacobian.block<3, 3>(3, 3) = dFdx;
            jacobian.block<3, 3>(0, 3) = -dFdx;
            jacobian.block<3, 3>(3, 0) = -dFdx;
        }
    
        static inline void getStress(const VecS<6> &x, const SpringData& data, VecS<3> &stress) {
            VecS<3> x1 = x.segment<3>(0);
            VecS<3> x2 = x.segment<3>(3);
            VecS<3> x21 = x2 - x1;
            Real l = x21.norm();
    
            if (l < 1e-8) {
                stress.setZero();
               
            }
    
            Real strain = (l - data.rest_length) / data.rest_length;
            Real stress_magnitude = data.k * strain;
            stress = stress_magnitude * (x21 / l);
        }
    };
    

    struct StVKModel
    {
        static inline Real getEnergy(
            const VecS<9>& positions,
            const FEMTriData& data)
        {
            // Compute Deformation Gradient
            MatS<3,2> F;
            MatS<3,2> B;
            computeDeformationGradient(positions, data.invRestMat, B, F);

            // Compute Green Strain
            const MatS<2,2> greenStrain = (F.transpose() * F - MatS<2,
            2>::Identity()) * 0.5;

            // Compute Energy
            const Real trace = greenStrain.trace();
            const Real energy = data.restVolume * (data.mu * greenStrain.squaredNorm() +
            data.lambda * 0.5 * trace * trace);

            return energy;
        }

        static inline void getStress(
                const MatS<3, 2>& F, const FEMTriData& data, MatS<3, 2>& sigma)
        {
            const MatS<2,2> greenStrain = (F.transpose() * F - MatS<2, 2>::Identity()) * 0.5;
            sigma = F * (2 * data.mu * greenStrain + data.lambda * greenStrain.trace() * MatS<2,2>::Identity());
        }
        static inline void getForce(
            const VecS<9>& x,
            const FEMTriData& data,
                VecS<9>& force)
        {

            MatS<3,2> B; B.setZero();
            MatS<3,2> F; F.setZero();

            computeDeformationGradient(x, data.invRestMat, B, F);

            MatS<3, 2> sigma;
            getStress(F, data, sigma);

            auto force_map = Eigen::Map<MatS<3, 3>>(reinterpret_cast<Real*>(force.data()));
            force_map = sigma * B.transpose() * data.restVolume;
            force = -force;
        }

        static inline void getForceJacobian(
            const VecS<9>& positions,
            const FEMTriData& data,
            MatS<9, 9>& forcesPositionDeriv)
        {
            MatS<3, 2> B;
            MatS<3, 2> F;
            computeDeformationGradient(positions, data.invRestMat, B, F);
            const MatS<2, 3> BT = B.transpose();
            const MatS<2, 2> greenStrain = (F.transpose() * F - MatS<2, 2>::Identity()) * 0.5;

            const Real greenStrainTrace = greenStrain.trace();
            const Real lambdaGreenTrace = data.lambda * greenStrainTrace;

            for (int s = 0; s < 2; s++) {
                const Eigen::Block<MatS<3, 2>, 3, 1, true> delF_rsdelx_j = B.col(s);

                for (int l = 0; l < 2; l++) {
                    Real delP_ildelF_rs_p1_cached = 2 * data.mu * greenStrain(s, l);
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
                                delP_ildelF_rs_p2 += data.mu * F.row(i).dot(delE_kldelF_rs_p1);
                            }

                            delP_ildelF_rs_p2 += F(i, s) * data.mu * F(r, l);
                            delP_ildelF_rs_p2 += F(i, l) * data.lambda * deltrEdelF_rs;

                            Real delP_ildelF_rs = delP_ildelF_rs_p1 + delP_ildelF_rs_p2;

                            for (int x = 0; x < 3; x++) {
                                for (int y = 0; y < 3; y++) {
                                    forcesPositionDeriv(i + 3 * x, r + 3 * y) +=
                                        -data.restVolume * delP_ildelF_rs * delF_rsdelx_j(y) *
                                        B(x, l);
                                }
                            }
                        }
                    }
                }
            }
        }
};

    enum MaterialType
    {
        SPRING,
        STVK,
        NEOHOOKEAN
    };

    using MaterialVariant = std::variant<SpringData, FEMTriData>;

    struct Material
    {
        MaterialType type;
        MaterialVariant data;
    };

}

#endif //MATERIALS_H






