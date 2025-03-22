#include "energy.h"

namespace tinyed
{
    struct SpringData {
        Real rest_length;
        Real stiffness;
    };

    // TODO: Replace Dynamic matrices with static matrices with templates
    class SpringModel : public EnergyModel {
        private:
            std::vector<SpringData> data;
            std::vector<std::pair<int, int>> constraints;
            int dim;

        public:
            // TODO: replace std::pair with VectorXi
            SpringModel(const std::vector<VectorXd>& rest_vertices, const std::vector<std::pair<int, int>>& constraints, Real k)
                : constraints(constraints), dim(rest_vertices[0].size()) {
                for (const auto& edge : constraints) {
                    Real rest_length = (rest_vertices[edge.first] - rest_vertices[edge.second]).norm();
                    data.push_back({rest_length, k});
                }
            }

            Real getEnergy(int idx, const std::vector<VectorXd>& x) const override;

            // TODO: pass by reference for force and jacobians
            // TODO: replace std::vector  with flattened VectorXd
            // void getForce(int idx, const std::vector<VectorXd>& x, VectorXd& force) const override;
            VectorXd getForce(int idx, const std::vector<VectorXd>& x) const override;
            MatrixXd getForceJacobian(int idx, const std::vector<VectorXd>& x) const override;
        };
}