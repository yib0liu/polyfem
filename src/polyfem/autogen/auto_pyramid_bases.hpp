#pragma once

#include <Eigen/Dense>
#include <cassert>

namespace polyfem {
namespace autogen {
void pyramid_nodes_3d(const int o, Eigen::MatrixXd &val);

void pyramid_basis_value_3d(const int o, const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &val);

void pyramid_grad_basis_value_3d(const int o, const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &val);



static const int MAX_Pyramid_BASES = 2;

}}
