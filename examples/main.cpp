#include <iostream>
#include "materials.h"

using namespace TinyED;
int main()
{
    VecS<6> x;
    x << 1, 0, 0,
         0, 0, 0;

    Real l0 = 1.0;  // rest length
    Real k = 100.0; // stiffness constant

    const Real energy = SpringModel::getEnergy(x, l0, k);
    std::cout << "Computed energy: " << energy << std::endl;

}