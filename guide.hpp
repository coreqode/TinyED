
class EnergyModel
    // Single Elements
    @abstract
    getEnergy(x) 

    @abstract
    getForceandJacobian(x, want_jacobian = True/False) 

    // Will loop over all elements
    @abstract
    getTotalEnergy(x)

    @abstract
    getTotalForceAndJacobian(x, want_jacobian = True/False)


class SpringModel: EnergyModel
    private: 
      k: springCoeff

    getEnergy(x)

    getForceandJacobian(x)

clas FemModel: EnergyModel
      mu
      lambda


class StVK: FemModel

    getEnergy(x)
      invRestMat
      deformationGradient

    getForceandJacobian
      invRestMat
      deformationGradient



main()
{
  Material mat(StVK); 
  energy = mat.getTotalEnergy(x)
}



