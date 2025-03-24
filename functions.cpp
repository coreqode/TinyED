#include "functions.h"
#include <iostream>
static Real delta(Real i, Real j)
{
	return i == j ? (Real)1 : (Real)0;
}

void FEM::CalculateDeformationGradientTriangles(const Vectorr<9>& positions, const Matrixr<2,2>& invRestMat, Matrixr<3,2>& deformedCoordsToDeformationGradient, Matrixr<3, 2>& F)
{
	Matrixr<3,2> B;

  B << -1, -1,
      1, 0,
      0, 1;

	deformedCoordsToDeformationGradient = B * invRestMat;
	Matrixr<3,3> positionsMat;

	for(unsigned int p = 0; p < 3; p++)
	{
		positionsMat.col(p) = positions.segment(p*3, 3);
	}
	F = positionsMat * deformedCoordsToDeformationGradient;
}

Real FEM::stvkEnergy(
		const Vectorr<9>& positions,
        const Matrixr<2,2>& invRestMat,
        const Real restVolume,
        const Real mu, const Real lambda)
{ std::cout<<"stvkEnergy\n";
    Matrixr<3,2> sigma;

    Matrixr<3,2> deformedCoordsToDeformationGradient;
    Matrixr<3,2> F;

    CalculateDeformationGradientTriangles(positions, invRestMat, deformedCoordsToDeformationGradient, F);
    Matrixr<2,2> greenStrain = (F.transpose() * F - Matrixr<2,2>::Identity()) * 0.5;
    Real trace = greenStrain.trace();
    return restVolume * (mu * greenStrain.squaredNorm() + lambda * 0.5 * trace * trace);
}

void FEM::stvkPk1Stress(
		const Real mu, const Real lambda, const Matrixr<3,2>& F, Matrixr<3,2>& sigma)
{
	Matrixr<2,2> greenStrain = (F.transpose() * F - Matrixr<2,2>::Identity()) * 0.5;
	// St. Venant-Kirchhoff model
	sigma = F * (2 * mu * greenStrain + lambda * greenStrain.trace() * Matrixr<2,2>::Identity());
}

void FEM::femStretchForceTriangles(
		const Vectorr<9>& positions,
		const Matrixr<2,2>& invRestMat,
		const Real restVolume,
		const Real mu, const Real lambda,
    Vectorr<9>& force)
{
	Matrixr<3, 2> sigma;

	Matrixr<3, 2> deformedCoordsToDeformationGradient;
	Matrixr<3,2> F;
	CalculateDeformationGradientTriangles(positions, invRestMat, deformedCoordsToDeformationGradient, F);
  FEM::stvkPk1Stress(mu, lambda, F, sigma);

	auto force_map = Eigen::Map<Matrixr<3, 3>>(reinterpret_cast<Real*>(force.data()));
	force_map = sigma * deformedCoordsToDeformationGradient.transpose() * restVolume;
	force = -force;
}


void FEM::stvkForceJacobians(const Vectorr<9>& positions, const Matrixr<2,2>& invRestMat, const Real volume, const Real mu, const Real lambda, Vectorr<9>& force, Matrixr<9,9>& forcesPositionDeriv)
{

	femStretchForceTriangles(positions, invRestMat, volume, mu, lambda, force);

  const int N = 3; 

	Matrixr<N,N-1> deformedCoordsToDeformationGradient;
	Matrixr<3,N-1> F;
	CalculateDeformationGradientTriangles(positions, invRestMat, deformedCoordsToDeformationGradient, F);
	const Matrixr<N-1, N>& BT = deformedCoordsToDeformationGradient.transpose();
	const Matrixr<N-1, N-1> greenStrain = (F.transpose() * F - Matrixr<N - 1, N - 1>::Identity()) * 0.5;

	const Real greenStrainTrace = greenStrain.trace();
	const Real lambdaGreenTrace = lambda * greenStrainTrace;

	for (int s = 0; s < N - 1; s++)
	{
		const Eigen::Block<Matrixr<N, N-1>, N, 1, true> delF_rsdelx_j = deformedCoordsToDeformationGradient.col(s);
		for (int l = 0; l < N - 1; l++)
		{
			Real delP_ildelF_rs_p1_cached = 2 * mu * greenStrain(s, l);
			if (s == l)
			{
				delP_ildelF_rs_p1_cached += lambdaGreenTrace;
			}
			for (int r = 0; r < 3; r++)
			{
				Real deltrEdelF_rs = F(r, s);
				Matrixr<N - 1, 1> delE_kldelF_rs_p1 = Matrixr<N - 1, 1>::Zero();
				if (l == s)
				{
					for (int k = 0; k < N - 1; k++)
					{
						delE_kldelF_rs_p1(k) = F(r, k);
					}
				}
				for (int i = 0; i < 3; i++)
				{
					Real delP_ildelF_rs_p1 = i == r ? delP_ildelF_rs_p1_cached : 0;
					Real delP_ildelF_rs_p2 = 0;
					if (l == s)
					{
						delP_ildelF_rs_p2 += mu * F.row(i).dot(delE_kldelF_rs_p1);
					}
					delP_ildelF_rs_p2 += F(i, s) * (mu * F(r, l));
					delP_ildelF_rs_p2 += F(i, l) * lambda * deltrEdelF_rs;
					Real delP_ildelF_rs = delP_ildelF_rs_p1 + delP_ildelF_rs_p2;

					for (int x = 0; x < N; x++)
					{
						for (int y = 0; y < N; y++)
						{
							forcesPositionDeriv(i + 3 * x, r + 3 * y) += -volume * delP_ildelF_rs * delF_rsdelx_j(y) * deformedCoordsToDeformationGradient(x, l);
						}
					}
				}
			}
		}
	}
}

// Real HYLCMaterial::calculateEnergy(
// 	const Vectorr<18> &positions,
// 	const bool nn0_exists, const bool nn1_exists, const bool nn2_exists,
// 	const Matrixr<2,2>& invRestMat,
// 	const Real restVolume,
// 	const hylc::Config *hylc_.getStVKEnergyDerivativesconfig
// ){

// 	namespace mm = hylc::mathematica;
// 	hylc::Vec18 hylc_xlocal;
// 	for (int i = 0; i < 18; i++) hylc_xlocal(i) = positions(i);

// 	Mat2x2 hylc_invDm;
// 	for (int i = 0; i < 2; i++)
// 		for (int j = 0; j < 2; j++)
// 			hylc_invDm(i, j) = invRestMat(i, j);

//   	hylc::Vec6 strain;
//   	strain = mm::strain(hylc_xlocal, hylc_invDm, Vec3(nn0_exists, nn1_exists, nn2_exists));
// 	strain[0] -= 1;
// 	strain[2] -= 1;
//   	for (int i = 3; i < 6; i++) strain(i) *= hylc_config->bend_scale;

//   	Real E = hylc_config->material->psi(strain);
//   	E *= restVolume;
//   	return E * hylc_config->stiffness_mult;
// }


// void HYLCMaterial::calculateForceJacobian(
// 	const Vectorr<18> &positions,
// 	const bool nn0_exists, const bool nn1_exists, const bool nn2_exists,
// 	const Matrix2r& invRestMat,
// 	const Real restVolume,
// 	hylc::Config *hylc_config,
// 	Vector18r& phi,
// 	Matrixr<18,18>& jacobianFull)
// {
// 	namespace mm = hylc::mathematica;

// 	hylc::Vec18 hylc_xlocal;
// 	for (int i = 0; i < 18; i++) hylc_xlocal(i) = positions(i);

// 	Mat2x2 hylc_invDm;
// 	for (int i = 0; i < 2; i++)
// 		for (int j = 0; j < 2; j++)
// 			hylc_invDm(i, j) = invRestMat(i, j);

// 	std::tuple<std::vector<hylc::Mat18x18>, hylc::Mat6x18, hylc::Vec6> strain_hgv;
//   	strain_hgv = mm::strain_valdrv(hylc_xlocal, hylc_invDm,
//                                  Vec3(nn0_exists, nn1_exists, nn2_exists));

// 	hylc::Vec6 &strain = std::get<2>(strain_hgv);
// 	strain[0] -= 1;
// 	strain[2] -= 1;
// 	hylc::Mat6x18 &strain_grad               = std::get<1>(strain_hgv);
// 	std::vector<hylc::Mat18x18> &strain_hess = std::get<0>(strain_hgv);
// 	for (int i = 3; i < 6; i++) strain(i) *= hylc_config->bend_scale;

//   	std::pair<hylc::Mat6x6, hylc::Vec6> psidrv = hylc_config->material->psi_drv(strain);
//   	hylc::Vec18 g = transpose(strain_grad) * psidrv.second;
//   	hylc::Mat18x18 H = transpose(strain_grad) * psidrv.first * strain_grad;
//   	for (int i = 0; i < 6; ++i) H += psidrv.second[i] * strain_hess[i];

// 	g = g * restVolume;
// 	H = H * restVolume;
// 	g  = g * hylc_config->stiffness_mult;

// 	for (int i = 0; i < 18; i++) phi(i) = -g(i);

// 	 // SVD clamp eigenvalues of local matrix
//   if (hylc_config->hess_proj)
//   {
//     typedef Eigen::Matrix<double, 18, 18> emat18;
//     typedef Eigen::Matrix<double, 18, 1> evec18;
//     emat18 eH;
//     for (int i = 0; i < 18; ++i) {
//       for (int j = 0; j < 18; ++j) {
//       eH(i, j) = H(i, j);
//       }
//     }

//       Eigen::SelfAdjointEigenSolver<emat18> slv;
//       slv.compute(eH);

//     double mineig = 1e-8;
//     if (slv.eigenvalues()(0) < mineig) {
//       evec18 D = slv.eigenvalues().cwiseMax(mineig);

//       eH = slv.eigenvectors() * D.asDiagonal() * slv.eigenvectors().transpose();

//       for (int i = 0; i < 18; ++i) {
//       for (int j = 0; j < 18; ++j) {
//         H(i, j) = eH(i, j);
//       }
//       }
//     }
//   }

// 	H  = H * hylc_config->stiffness_mult;
// 	for (int i = 0; i < 18; i++)
// 		for (int j = 0; j < 18; j++)
// 			jacobianFull(i, j) = -H(i, j);
// }

void FEM::CalculateDeformationGradientTets(const Vectorr<12>& positions, const Matrixr<3,3>& invRestMat, Matrixr<4,3>& deformedCoordsToDeformationGradient, Matrixr<3, 3>& F)
{
	Matrixr<4, 3> B;

	 B << 1, 0, 0,
				0, 1, 0,
				0, 0, 1,
				-1, -1, -1;

	deformedCoordsToDeformationGradient = B * invRestMat;
	Matrixr<3,4> positionsMat;

	for(unsigned int p = 0; p < 4; p++)
	{
		positionsMat.col(p) = positions.segment(p*3, 3);
	}
	F = positionsMat * deformedCoordsToDeformationGradient;
}


Real FEM::stableNeohookeanEnergy( const Vectorr<12>& positions, const Matrixr<3, 3>& invRestMat, const Real restVolume, const Real mu, const Real lambda)
{
    Matrixr<3, 3> sigma;
    Matrixr<4, 3> deformedCoordsToDeformationGradient;
    Matrixr<3,3> F;

    CalculateDeformationGradientTets(positions, invRestMat, deformedCoordsToDeformationGradient, F);

    Matrix3r FTF = F.transpose() * F;
    Real I2 = FTF.trace();
    Real J = F.determinant();

    return restVolume * (mu * 0.5 * (I2 - 3) - mu * (J - 1) + lambda * 0.5 * (J - 1) * (J - 1));
}

void FEM::stableNeohookeanPk1Stress( const Real mu, const Real lambda, const Matrix3r& F, Matrix3r& sigma)
{
	// Stable Neo Hookean
	// PK = mu * F + (lambda + mu )* (J - 1 - mu/lambda) * (dJ/df)
	Real J = F.determinant();
	Matrix3r dJ_dF;
	dJ_dF.col(0) = F.col(1).cross(F.col(2));
	dJ_dF.col(1) = F.col(2).cross(F.col(0));
	dJ_dF.col(2) = F.col(0).cross(F.col(1));
	sigma = mu * F - mu * dJ_dF + lambda * (J - 1) * dJ_dF;
}

void FEM::femStretchForceTets(
		const Vectorr<12>& positions,
		const Matrixr<3,3>& invRestMat,
		const Real restVolume,
		const Real mu, const Real lambda,
    Vectorr<12>& force)
{
	Matrixr<3, 3> sigma;

	Matrixr<4, 3> deformedCoordsToDeformationGradient;
	Matrixr<3,3> F;
	CalculateDeformationGradientTets(positions, invRestMat, deformedCoordsToDeformationGradient, F);
  FEM::stableNeohookeanPk1Stress(mu, lambda, F, sigma);

	auto force_map = Eigen::Map<Matrixr<3, 4>>(reinterpret_cast<Real*>(force.data()));
	force_map = sigma * deformedCoordsToDeformationGradient.transpose() * restVolume;
	force = -force;
}


void FEM::stableNeohookeanForceJacobians(const Vectorr<12>& positions, const Matrixr<3, 3>& invRestMat, const Real volume, const Real mu, const Real lambda, const bool get_force_jacobian, Vector12r &force, Matrix1212r& forcesPositionDeriv)
{
	femStretchForceTets(positions, invRestMat, volume, mu, lambda, force);
  if (get_force_jacobian){
    Matrixr<4,3> B;
    Matrixr<3,3> F;
    CalculateDeformationGradientTets(positions, invRestMat, B, F);

        Matrix3r dJ_dF;
      dJ_dF.col(0) = F.col(1).cross(F.col(2));
      dJ_dF.col(1) = F.col(2).cross(F.col(0));
      dJ_dF.col(2) = F.col(0).cross(F.col(1));

      const Real J = F.determinant();
      const Real scale = lambda * (J - 1);

      //lambda function
      auto cross_product = [](const Matrix3r F, const int i)
      {
          Matrix3r fhat;
          fhat <<  0, -F(2,i),  F(1,i),
                    F(2,i),       0, -F(0,i),
                    -F(1,i),  F(0,i),  0;
          return fhat;
      };

      Matrix3r f0hat = cross_product(F, 0);
      Matrix3r f1hat = cross_product(F, 1);
      Matrix3r f2hat = cross_product(F, 2);

    // build the fractal cross-product
      Matrix9r hessJ;
      hessJ.setZero();
      for (int j = 0; j < 3; j++)
          for (int i = 0; i < 3; i++)
          {
            hessJ(i, j + 3) = -f2hat(i,j);
            hessJ(i + 3, j) =  f2hat(i,j);

            hessJ(i, j + 6) =  f1hat(i,j);
            hessJ(i + 6, j) = -f1hat(i,j);

            hessJ(i + 3, j + 6) = -f0hat(i,j);
            hessJ(i + 6, j + 3) =  f0hat(i,j);
          }

    for (int i = 0; i < 3; i++){
          for (int j = 0; j < 4; j++){
              for (int a = 0; a < 3; a++){
                  for (int b = 0; b < 4; b++){
                      Real term = 0;
                      for (int k = 0; k < 3; k++){
                          Real dPik_dXab = 0;
                          for (int m = 0; m < 3; m++){
                              for (int n = 0; n < 3; n++){
                                  Real T1 = mu * (delta(i,m) * delta(k,n) - hessJ(i+3*k, m+3*n));
                                  Real T2 = lambda * dJ_dF(i,k) * dJ_dF(m,n);
                                  Real T3 = lambda * (J - 1) * hessJ(i+3*k, m+3*n);

                                  Real dPik_dFmn = T1 + T2 + T3;
                                  Real dFmn_dXab = 0;
                                  for (int t = 0; t < 4; t++){
                                      dFmn_dXab += B(t,n) * delta(m,a) * delta(t,b);
                                  }
                                  dPik_dXab += dPik_dFmn * dFmn_dXab;
                              }
                          }
                          term += -volume * dPik_dXab * B(j,k);
                      }
                      forcesPositionDeriv(i + 3*j, a + 3*b) = term;
                  }
              }
          }
      }
  }
}

Real FEM::neohookeanEnergy(
		const Vectorr<12>& positions,
        const Matrixr<3, 3>& invRestMat,
        const Real restVolume,
        const Real mu, const Real lambda)
{
    Matrixr<3, 3> sigma;

    Matrixr<4, 3> deformedCoordsToDeformationGradient;
    Matrixr<3,3> F;

    CalculateDeformationGradientTets(positions, invRestMat, deformedCoordsToDeformationGradient, F);
	Real j = F.determinant();
	if (j < 1e-10)
	{
		return std::numeric_limits<Real>::infinity();
	}
	Real logJ = log(j);
    Real I1 = (F.transpose() * F).trace();

    return restVolume * (mu * 0.5 * (I1 - 3) - mu * logJ + lambda * 0.5 * logJ * logJ);
}

void FEM::neohookeanPk1Stress(
		const Real mu, const Real lambda, const Matrix3r& F, Matrix3r& sigma)
{
	Matrix3r FInvT = F.inverse().transpose();
	Real j = F.determinant();
	// Neohookean Model
	sigma = mu * (F - FInvT) + lambda * log(j) * FInvT;
}

void FEM::neohookeanForceJacobians(const Vectorr<12>& positions, const Matrixr<3, 3>& invRestMat, const Real volume, const Real mu, const Real lambda, const bool get_force_jacobian, Vector12r& force, Matrix1212r& forcesPositionDeriv)
{
	femStretchForceTets(positions, invRestMat, volume, mu, lambda, force);

  if (get_force_jacobian)
  {
    Matrixr<4,3> deformedCoordsToDeformationGradient;
    Matrixr<3,3> F;
    CalculateDeformationGradientTets(positions, invRestMat, deformedCoordsToDeformationGradient, F);

    Matrix3r FInv = F.inverse();
    Real J = F.determinant();
    Real term = (mu - lambda * log(J));

    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        Matrix3r delPij_delFkl_coalesced = term * FInv.col(i) * FInv.row(j) + lambda * FInv(j, i) * FInv;
        delPij_delFkl_coalesced(j, i) += mu;
        for (int l = 0; l < 3; l++)
        {
          for (int k = 0; k < 3; k++)
          {
            forcesPositionDeriv(Eigen::seqN(i, 4, 3), Eigen::seqN(k, 4, 3)) += -volume * delPij_delFkl_coalesced(l, k) * deformedCoordsToDeformationGradient.col(j) * deformedCoordsToDeformationGradient.col(l).transpose();
          }
        }
      }
    }
  }
}

