#include <iostream>
#include "materials.h"
#include "examples_util.h"

using namespace TinyED;
int main()
{
    constexpr int TOTAL_FRAMES = 500;
    constexpr Real timeStep = 0.01;
    constexpr Real density = 1.0;
    constexpr int tesselation = 8;
    const Vec3 gravity = {0, -9.81, 0};

    auto [triangles, positions, uvs] = generateTesselatedSquare({0, 0, 0}, tesselation, 1.0);
    const std::vector<Mat2> restMatrices = computeRestShape(positions, uvs, triangles);
    const std::vector<Real> restVolumes = computeRestVolume(uvs, triangles);

    std::vector<Material> matData(triangles.size());
    std::vector<Element> elements(triangles.size());

    for (size_t i = 0; i < triangles.size(); ++i)
    {
        matData[i].type = MaterialType::STVK;
        matData[i].data = FEMTriData{};
        std::get<FEMTriData>(matData[i].data).mu = 100.0;
        std::get<FEMTriData>(matData[i].data).lambda = 150.0;
        std::get<FEMTriData>(matData[i].data).restVolume = restVolumes[i];
        std::get<FEMTriData>(matData[i].data).invRestMat = restMatrices[i];

        elements[i].elementData = triangles[i];
    }

    std::vector<unsigned int> constraintIndices = {0, tesselation};

    const VecN masses = computeVertexMasses(positions, triangles, density);
    VecN invMasses = 1.0 / masses.array();

    const int N = positions.cols();

    Mat3N velocities(3, N);
    velocities.setZero();


    VecN invMlumped(3*N);
    for (size_t i = 0; i < N; i++)
    {
        if (std::find(constraintIndices.begin(), constraintIndices.end(), i) != constraintIndices.end())
            invMasses(i) = 0;
        invMlumped.segment<3>(i*3).setConstant(invMasses(i));
    }

    for (int frame = 0; frame < TOTAL_FRAMES; frame++)
    {
        advanceStep(positions, velocities, timeStep, elements, matData, invMlumped, masses, gravity);
        const std::string filename = "frame_" + std::to_string(frame);
        saveSheet(filename, positions, triangles);
    }
}