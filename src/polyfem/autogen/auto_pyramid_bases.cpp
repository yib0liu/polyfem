#include "auto_pyramid_bases.hpp"
#include "auto_b_bases.hpp"
#include "p_n_bases.hpp"


namespace polyfem {
namespace autogen {
namespace {
void pyramid_1_basis_value_3d(const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &result_0){

auto x=uv.col(0).array();
auto y=uv.col(1).array();
auto z=uv.col(2).array();

switch(local_index){
	case 0: {const auto helper_0 = z - 1;
result_0 = -(helper_0*(helper_0 + x + y) + x*y)/helper_0;} break;
	case 1: {const auto helper_0 = z - 1;
result_0 = x*(helper_0 + y)/helper_0;} break;
	case 2: {result_0 = -x*y/(z - 1);} break;
	case 3: {const auto helper_0 = z - 1;
result_0 = y*(helper_0 + x)/helper_0;} break;
	case 4: {result_0 = z;} break;
	default: assert(false);
}}
void pyramid_1_basis_grad_value_3d(const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &val){

auto x=uv.col(0).array();
auto y=uv.col(1).array();
auto z=uv.col(2).array();

val.resize(uv.rows(), uv.cols());
 Eigen::ArrayXd result_0(uv.rows());
switch(local_index){
	case 0: {{const auto helper_0 = z - 1;
result_0 = -(helper_0 + y)/helper_0;val.col(0) = result_0; }{const auto helper_0 = z - 1;
result_0 = -(helper_0 + x)/helper_0;val.col(1) = result_0; }{const auto helper_0 = 2*z;
const auto helper_1 = pow(z, 2);
result_0 = (helper_0 - helper_1 + x*y - 1)/(-helper_0 + helper_1 + 1);val.col(2) = result_0; }} break;
	case 1: {{const auto helper_0 = z - 1;
result_0 = (helper_0 + y)/helper_0;val.col(0) = result_0; }{result_0 = x/(z - 1);val.col(1) = result_0; }{result_0 = -x*y/pow(z - 1, 2);val.col(2) = result_0; }} break;
	case 2: {{result_0 = -y/(z - 1);val.col(0) = result_0; }{result_0 = -x/(z - 1);val.col(1) = result_0; }{result_0 = x*y/pow(z - 1, 2);val.col(2) = result_0; }} break;
	case 3: {{result_0 = y/(z - 1);val.col(0) = result_0; }{const auto helper_0 = z - 1;
result_0 = (helper_0 + x)/helper_0;val.col(1) = result_0; }{result_0 = -x*y/pow(z - 1, 2);val.col(2) = result_0; }} break;
	case 4: {{result_0.setZero();val.col(0) = result_0; }{result_0.setZero();val.col(1) = result_0; }{result_0.setOnes();val.col(2) = result_0; }} break;
	default: assert(false);
}}


void pyramid_1_nodes_3d(Eigen::MatrixXd &res) {
 res.resize(5, 3); res << 
0, 0, 0,
1, 0, 0,
1, 1, 0,
0, 1, 0,
0, 0, 1;
}


void pyramid_2_basis_value_3d(const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &result_0){

auto x=uv.col(0).array();
auto y=uv.col(1).array();
auto z=uv.col(2).array();

switch(local_index){
	case 0: {const auto helper_0 = z - 1;
const auto helper_1 = pow(z, 2);
const auto helper_2 = helper_1 - 2*z + 1;
const auto helper_3 = pow(x, 2);
const auto helper_4 = pow(y, 2);
const auto helper_5 = 4*z;
result_0 = (helper_0*helper_2*(2*helper_1 + 2*helper_3 + 2*helper_4 + helper_5*x + helper_5*y - 3*x - 3*y - 3*z + 1) + 4*helper_0*helper_3*helper_4 + helper_2*x*y*(6*x + 6*y + 10*z - 9))/(helper_0*helper_2);} break;
	case 1: {const auto helper_0 = z - 1;
const auto helper_1 = 2*z;
const auto helper_2 = -helper_1 + pow(z, 2) + 1;
result_0 = x*(helper_0*helper_2*(2*x - 1) + 4*helper_0*x*pow(y, 2) + helper_2*y*(helper_1 + 6*x + 2*y - 3))/(helper_0*helper_2);} break;
	case 2: {const auto helper_0 = z - 1;
const auto helper_1 = 2*z;
const auto helper_2 = -helper_1 + pow(z, 2) + 1;
const auto helper_3 = x*y;
result_0 = helper_3*(4*helper_0*helper_3 + helper_2*(helper_1 + 2*x + 2*y - 1))/(helper_0*helper_2);} break;
	case 3: {const auto helper_0 = z - 1;
const auto helper_1 = 2*z;
const auto helper_2 = -helper_1 + pow(z, 2) + 1;
result_0 = y*(helper_0*helper_2*(2*y - 1) + 4*helper_0*pow(x, 2)*y + helper_2*x*(helper_1 + 2*x + 6*y - 3))/(helper_0*helper_2);} break;
	case 4: {result_0 = z*(2*z - 1);} break;
	case 5: {const auto helper_0 = z - 1;
const auto helper_1 = pow(z, 2) - 2*z + 1;
result_0 = -4*x*(helper_0*helper_1*(helper_0 + x) + 2*helper_0*x*pow(y, 2) + helper_1*y*(3*x + 2*y + 3*z - 3))/(helper_0*helper_1);} break;
	case 6: {const auto helper_0 = 2*z;
const auto helper_1 = -helper_0 + pow(z, 2) + 1;
const auto helper_2 = 2*x;
result_0 = -4*x*y*(helper_0*x + helper_1 + helper_2*y - helper_2 + y*z - y)/helper_1;} break;
	case 7: {const auto helper_0 = 2*z;
const auto helper_1 = -helper_0 + pow(z, 2) + 1;
const auto helper_2 = 2*y;
result_0 = -4*x*y*(helper_0*y + helper_1 + helper_2*x - helper_2 + x*z - x)/helper_1;} break;
	case 8: {const auto helper_0 = z - 1;
const auto helper_1 = pow(z, 2) - 2*z + 1;
result_0 = -4*y*(helper_0*helper_1*(helper_0 + y) + 2*helper_0*pow(x, 2)*y + helper_1*x*(2*x + 3*y + 3*z - 3))/(helper_0*helper_1);} break;
	case 9: {const auto helper_0 = z - 1;
result_0 = -4*z*(helper_0*(helper_0 + x + y) + x*y)/helper_0;} break;
	case 10: {const auto helper_0 = z - 1;
result_0 = 4*x*z*(helper_0 + y)/helper_0;} break;
	case 11: {result_0 = -4*x*y*z/(z - 1);} break;
	case 12: {const auto helper_0 = z - 1;
result_0 = 4*y*z*(helper_0 + x)/helper_0;} break;
	case 13: {const auto helper_0 = pow(z, 2) - 2*z + 1;
const auto helper_1 = x*y;
result_0 = 16*helper_1*(helper_0 + helper_1 + x*z - x + y*z - y)/helper_0;} break;
	default: assert(false);
}}
void pyramid_2_basis_grad_value_3d(const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &val){

auto x=uv.col(0).array();
auto y=uv.col(1).array();
auto z=uv.col(2).array();

val.resize(uv.rows(), uv.cols());
 Eigen::ArrayXd result_0(uv.rows());
switch(local_index){
	case 0: {{const auto helper_0 = z - 1;
const auto helper_1 = pow(z, 2) - 2*z + 1;
const auto helper_2 = 6*x;
const auto helper_3 = helper_1*y;
result_0 = (helper_0*helper_1*(4*x + 4*z - 3) + 8*helper_0*x*pow(y, 2) + helper_2*helper_3 + helper_3*(helper_2 + 6*y + 10*z - 9))/(helper_0*helper_1);val.col(0) = result_0; }{const auto helper_0 = z - 1;
const auto helper_1 = pow(z, 2) - 2*z + 1;
const auto helper_2 = 6*x;
result_0 = (helper_0*helper_1*(4*y + 4*z - 3) + 8*helper_0*pow(x, 2)*y + helper_1*helper_2*y + helper_1*x*(helper_2 + 6*y + 10*z - 9))/(helper_0*helper_1);val.col(1) = result_0; }{const auto helper_0 = pow(z, 3);
const auto helper_1 = pow(z, 2);
const auto helper_2 = 4*x;
const auto helper_3 = 4*y;
const auto helper_4 = x*y;
const auto helper_5 = 12*z;
const auto helper_6 = pow(y, 2);
const auto helper_7 = 6*helper_6*x;
const auto helper_8 = 12*helper_1;
const auto helper_9 = pow(x, 2);
const auto helper_10 = 6*helper_9*y;
result_0 = -(-helper_0*helper_2 - helper_0*helper_3 + 15*helper_0 - 21*helper_1 + helper_10*z - helper_10 + helper_2 + helper_3 + helper_4*z - helper_4 - helper_5*x - helper_5*y + 8*helper_6*helper_9 + helper_7*z - helper_7 + helper_8*x + helper_8*y - 4*pow(z, 4) + 13*z - 3)/(helper_0 - 3*helper_1 + 3*z - 1);val.col(2) = result_0; }} break;
	case 1: {{const auto helper_0 = 2*z;
const auto helper_1 = pow(z, 2);
const auto helper_2 = 4*x;
const auto helper_3 = 12*x*y;
const auto helper_4 = 8*x;
const auto helper_5 = pow(y, 2);
result_0 = (helper_0*helper_5 + helper_0 + helper_1*helper_2 + 2*helper_1*y - helper_1 + helper_2 + helper_3*z - helper_3 + helper_4*helper_5 - helper_4*z - 2*helper_5 - 5*y*z + 3*y - 1)/(-helper_0 + helper_1 + 1);val.col(0) = result_0; }{const auto helper_0 = pow(z, 2);
const auto helper_1 = 6*x;
const auto helper_2 = 4*y;
result_0 = x*(2*helper_0 + helper_1*z - helper_1 + helper_2*z - helper_2 + 8*x*y - 5*z + 3)/(helper_0 - 2*z + 1);val.col(1) = result_0; }{const auto helper_0 = 6*x;
const auto helper_1 = 2*y;
const auto helper_2 = x*y;
result_0 = -helper_2*(helper_0*z - helper_0 + helper_1*z - helper_1 + 8*helper_2 - z + 1)/(pow(z, 3) - 3*pow(z, 2) + 3*z - 1);val.col(2) = result_0; }} break;
	case 2: {{const auto helper_0 = 2*z;
const auto helper_1 = pow(z, 2);
const auto helper_2 = 4*x;
result_0 = y*(helper_0*y + 2*helper_1 + helper_2*z - helper_2 + 8*x*y - 2*y - 3*z + 1)/(-helper_0 + helper_1 + 1);val.col(0) = result_0; }{const auto helper_0 = 2*z;
const auto helper_1 = pow(z, 2);
const auto helper_2 = 4*y;
result_0 = x*(helper_0*x + 2*helper_1 + helper_2*z - helper_2 + 8*x*y - 2*x - 3*z + 1)/(-helper_0 + helper_1 + 1);val.col(1) = result_0; }{const auto helper_0 = 2*x;
const auto helper_1 = 2*y;
const auto helper_2 = x*y;
result_0 = -helper_2*(helper_0*z - helper_0 + helper_1*z - helper_1 + 8*helper_2 + z - 1)/(pow(z, 3) - 3*pow(z, 2) + 3*z - 1);val.col(2) = result_0; }} break;
	case 3: {{const auto helper_0 = pow(z, 2);
const auto helper_1 = 4*x;
const auto helper_2 = 6*y;
result_0 = y*(2*helper_0 + helper_1*z - helper_1 + helper_2*z - helper_2 + 8*x*y - 5*z + 3)/(helper_0 - 2*z + 1);val.col(0) = result_0; }{const auto helper_0 = 2*z;
const auto helper_1 = pow(z, 2);
const auto helper_2 = 4*y;
const auto helper_3 = 12*x*y;
const auto helper_4 = 8*y;
const auto helper_5 = pow(x, 2);
result_0 = (helper_0*helper_5 + helper_0 + helper_1*helper_2 + 2*helper_1*x - helper_1 + helper_2 + helper_3*z - helper_3 + helper_4*helper_5 - helper_4*z - 2*helper_5 - 5*x*z + 3*x - 1)/(-helper_0 + helper_1 + 1);val.col(1) = result_0; }{const auto helper_0 = 2*x;
const auto helper_1 = 6*y;
const auto helper_2 = x*y;
result_0 = -helper_2*(helper_0*z - helper_0 + helper_1*z - helper_1 + 8*helper_2 - z + 1)/(pow(z, 3) - 3*pow(z, 2) + 3*z - 1);val.col(2) = result_0; }} break;
	case 4: {{result_0.setZero();val.col(0) = result_0; }{result_0.setZero();val.col(1) = result_0; }{result_0 = 4*z - 1;val.col(2) = result_0; }} break;
	case 5: {{const auto helper_0 = 2*z;
const auto helper_1 = pow(z, 2);
const auto helper_2 = 2*x;
const auto helper_3 = 6*y;
const auto helper_4 = helper_3*x;
const auto helper_5 = 4*x;
const auto helper_6 = pow(y, 2);
const auto helper_7 = 3*helper_1;
result_0 = -4*(helper_0*helper_6 + helper_1*helper_2 + helper_2 - helper_3*z + helper_4*z - helper_4 + helper_5*helper_6 - helper_5*z - 2*helper_6 + helper_7*y - helper_7 + 3*y + pow(z, 3) + 3*z - 1)/(-helper_0 + helper_1 + 1);val.col(0) = result_0; }{const auto helper_0 = pow(z, 2);
const auto helper_1 = 3*x;
const auto helper_2 = 4*y;
result_0 = -4*x*(3*helper_0 + helper_1*z - helper_1 + helper_2*x + helper_2*z - helper_2 - 6*z + 3)/(helper_0 - 2*z + 1);val.col(1) = result_0; }{const auto helper_0 = 3*z;
const auto helper_1 = pow(z, 3);
const auto helper_2 = 3*pow(z, 2);
const auto helper_3 = x*y;
const auto helper_4 = pow(y, 2);
const auto helper_5 = 2*helper_4;
const auto helper_6 = 4*x;
result_0 = helper_6*(helper_0*helper_3 - helper_0 - helper_1 + helper_2 - 3*helper_3 + helper_4*helper_6 + helper_5*z - helper_5 + 1)/(helper_0 + helper_1 - helper_2 - 1);val.col(2) = result_0; }} break;
	case 6: {{const auto helper_0 = pow(z, 2) - 2*z + 1;
const auto helper_1 = 4*x;
result_0 = -4*y*(helper_0 + helper_1*y + helper_1*z - helper_1 + y*z - y)/helper_0;val.col(0) = result_0; }{const auto helper_0 = 2*z;
const auto helper_1 = -helper_0 + pow(z, 2) + 1;
const auto helper_2 = 4*x;
result_0 = -helper_2*(helper_0*x + helper_0*y + helper_1 + helper_2*y - 2*x - 2*y)/helper_1;val.col(1) = result_0; }{const auto helper_0 = 2*x;
const auto helper_1 = 4*x*y;
result_0 = helper_1*(helper_0*z - helper_0 + helper_1 + y*z - y)/(pow(z, 3) - 3*pow(z, 2) + 3*z - 1);val.col(2) = result_0; }} break;
	case 7: {{const auto helper_0 = 2*z;
const auto helper_1 = -helper_0 + pow(z, 2) + 1;
const auto helper_2 = 4*y;
result_0 = -helper_2*(helper_0*x + helper_0*y + helper_1 + helper_2*x - 2*x - 2*y)/helper_1;val.col(0) = result_0; }{const auto helper_0 = pow(z, 2) - 2*z + 1;
const auto helper_1 = 4*y;
result_0 = -4*x*(helper_0 + helper_1*x + helper_1*z - helper_1 + x*z - x)/helper_0;val.col(1) = result_0; }{const auto helper_0 = 2*y;
const auto helper_1 = 4*x*y;
result_0 = helper_1*(helper_0*z - helper_0 + helper_1 + x*z - x)/(pow(z, 3) - 3*pow(z, 2) + 3*z - 1);val.col(2) = result_0; }} break;
	case 8: {{const auto helper_0 = pow(z, 2);
const auto helper_1 = 4*x;
const auto helper_2 = 3*y;
result_0 = -4*y*(3*helper_0 + helper_1*y + helper_1*z - helper_1 + helper_2*z - helper_2 - 6*z + 3)/(helper_0 - 2*z + 1);val.col(0) = result_0; }{const auto helper_0 = 2*z;
const auto helper_1 = pow(z, 2);
const auto helper_2 = 2*y;
const auto helper_3 = 6*x;
const auto helper_4 = helper_3*y;
const auto helper_5 = 4*y;
const auto helper_6 = pow(x, 2);
const auto helper_7 = 3*helper_1;
result_0 = -4*(helper_0*helper_6 + helper_1*helper_2 + helper_2 - helper_3*z + helper_4*z - helper_4 + helper_5*helper_6 - helper_5*z - 2*helper_6 + helper_7*x - helper_7 + 3*x + pow(z, 3) + 3*z - 1)/(-helper_0 + helper_1 + 1);val.col(1) = result_0; }{const auto helper_0 = 3*z;
const auto helper_1 = pow(z, 3);
const auto helper_2 = 3*pow(z, 2);
const auto helper_3 = x*y;
const auto helper_4 = pow(x, 2);
const auto helper_5 = 2*helper_4;
const auto helper_6 = 4*y;
result_0 = helper_6*(helper_0*helper_3 - helper_0 - helper_1 + helper_2 - 3*helper_3 + helper_4*helper_6 + helper_5*z - helper_5 + 1)/(helper_0 + helper_1 - helper_2 - 1);val.col(2) = result_0; }} break;
	case 9: {{const auto helper_0 = z - 1;
result_0 = -4*z*(helper_0 + y)/helper_0;val.col(0) = result_0; }{const auto helper_0 = z - 1;
result_0 = -4*z*(helper_0 + x)/helper_0;val.col(1) = result_0; }{const auto helper_0 = z - 1;
const auto helper_1 = x + y;
const auto helper_2 = helper_0*(helper_0 + helper_1) + x*y;
result_0 = 4*(-helper_0*(helper_2 + z*(helper_1 + 2*z - 2)) + helper_2*z)/pow(helper_0, 2);val.col(2) = result_0; }} break;
	case 10: {{const auto helper_0 = z - 1;
result_0 = 4*z*(helper_0 + y)/helper_0;val.col(0) = result_0; }{result_0 = 4*x*z/(z - 1);val.col(1) = result_0; }{const auto helper_0 = 2*z;
const auto helper_1 = pow(z, 2);
result_0 = -4*x*(helper_0 - helper_1 + y - 1)/(-helper_0 + helper_1 + 1);val.col(2) = result_0; }} break;
	case 11: {{result_0 = -4*y*z/(z - 1);val.col(0) = result_0; }{result_0 = -4*x*z/(z - 1);val.col(1) = result_0; }{result_0 = 4*x*y/pow(z - 1, 2);val.col(2) = result_0; }} break;
	case 12: {{result_0 = 4*y*z/(z - 1);val.col(0) = result_0; }{const auto helper_0 = z - 1;
result_0 = 4*z*(helper_0 + x)/helper_0;val.col(1) = result_0; }{const auto helper_0 = 2*z;
const auto helper_1 = pow(z, 2);
result_0 = -4*y*(helper_0 - helper_1 + x - 1)/(-helper_0 + helper_1 + 1);val.col(2) = result_0; }} break;
	case 13: {{const auto helper_0 = 2*z;
const auto helper_1 = -helper_0 + pow(z, 2) + 1;
const auto helper_2 = 2*x;
result_0 = 16*y*(helper_0*x + helper_1 + helper_2*y - helper_2 + y*z - y)/helper_1;val.col(0) = result_0; }{const auto helper_0 = 2*z;
const auto helper_1 = -helper_0 + pow(z, 2) + 1;
const auto helper_2 = 2*y;
result_0 = 16*x*(helper_0*y + helper_1 + helper_2*x - helper_2 + x*z - x)/helper_1;val.col(1) = result_0; }{const auto helper_0 = x*y;
result_0 = -16*helper_0*(2*helper_0 + x*z - x + y*z - y)/(pow(z, 3) - 3*pow(z, 2) + 3*z - 1);val.col(2) = result_0; }} break;
	default: assert(false);
}}


void pyramid_2_nodes_3d(Eigen::MatrixXd &res) {
 res.resize(14, 3); res << 
0.0, 0.0, 0.0,
1.0, 0.0, 0.0,
1.0, 1.0, 0.0,
0.0, 1.0, 0.0,
0.0, 0.0, 1.0,
0.5, 0.0, 0.0,
1.0, 0.5, 0.0,
0.5, 1.0, 0.0,
0.0, 0.5, 0.0,
0.0, 0.0, 0.5,
0.5, 0.0, 0.5,
0.5, 0.5, 0.5,
0.0, 0.5, 0.5,
0.5, 0.5, 0.0;
}


}

void pyramid_nodes_3d(const int o, Eigen::MatrixXd &val){
switch(o){
	case 1: pyramid_1_nodes_3d(val); break;
	case 2: pyramid_2_nodes_3d(val); break;
	default: assert(false);
}}
void pyramid_basis_value_3d(const int o, const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &val){

switch(o){
	case 1: pyramid_1_basis_value_3d(local_index, uv, val); break;
	case 2: pyramid_2_basis_value_3d(local_index, uv, val); break;
	default: assert(false);
}}

void pyramid_grad_basis_value_3d(const int o, const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &val){

switch(o){
	case 1: pyramid_1_basis_grad_value_3d(local_index, uv, val); break;
	case 2: pyramid_2_basis_grad_value_3d(local_index, uv, val); break;
	default: assert(false);
}}

}}
