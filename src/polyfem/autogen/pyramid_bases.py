# https://raw.githubusercontent.com/sympy/sympy/master/examples/advanced/fem.py
from sympy import *
import os
import numpy as np
import argparse
from sympy.printing import ccode
import pretty_print
import sympy as sp

x, y, z = symbols('x,y,z')

def shuffle(a,b):
    return [a[i] for i in b]

if __name__ == "__main__":
    orders = [1,2]
    swapmap1 = [0,1,3,2,4]
    swapmap2 = [0,1,3,2,4,5,8,10,6,7,9,12,11,13]
    real_nodes=[
        np.array([
            [0,0,0],
            [1,0,0],
            [0,1,0],
            [1,1,0],
            [0,0,1]
        ])[swapmap1],
        np.array([
            [0,0,0],
            [1,0,0],
            [0,1,0],
            [1,1,0],
            [0,0,1],
            [1/2.0, 0, 0],
            [0, 1/2.0, 0],
            [0, 0, 1/2.0],
            [1, 1/2.0, 0],
            [1/2.0, 0, 1/2.0],
            [1/2.0, 1, 0],
            [0, 1/2.0, 1/2.0],
            [1/2.0, 1/2.0, 1/2.0],
            [1/2.0, 1/2.0, 0],
        ])[swapmap2]
    ]

    phi0 = (-x*y + (z - 1)*(-x - y - z + 1)) / (z - 1)
    phi1 = (x*(y + z - 1)) / (z - 1)
    phi2 = (y*(x + z - 1)) / (z - 1)
    phi3 = -(x*y) / (z - 1)
    phi4 = z

    basis_1 = [sp.simplify(p) for p in shuffle([phi0, phi1, phi2, phi3, phi4],swapmap1)]

    phi0  = (4*x**2*y**2*(z-1) + x*y*(z**2-2*z+1)*(6*x + 6*y + 10*z - 9) +
            (z-1)*(z**2-2*z+1)*(2*x**2 + 4*x*z - 3*x + 2*y**2 + 4*y*z - 3*y + 2*z**2 - 3*z + 1)) / ((z-1)*(z**2-2*z+1))

    phi1  = x*(4*x*y**2*(z-1) + y*(z**2-2*z+1)*(6*x + 2*y + 2*z - 3) + (2*x-1)*(z-1)*(z**2-2*z+1)) / ((z-1)*(z**2-2*z+1))

    phi2  = y*(4*x**2*y*(z-1) + x*(z**2-2*z+1)*(2*x + 6*y + 2*z - 3) + (2*y-1)*(z-1)*(z**2-2*z+1)) / ((z-1)*(z**2-2*z+1))

    phi3  = x*y*(4*x*y*(z-1) + (z**2-2*z+1)*(2*x + 2*y + 2*z - 1)) / ((z-1)*(z**2-2*z+1))

    phi4  = z*(2*z-1)

    phi5  = 4*x*(2*x*y**2*(1-z) + y*(z**2-2*z+1)*(-3*x -2*y -3*z +3)
                + (z-1)*(-x - z +1)*(z**2-2*z+1)) / ((z-1)*(z**2-2*z+1))

    phi6  = 4*y*(2*x**2*y*(1-z) + x*(z**2-2*z+1)*(-2*x -3*y -3*z +3)
                + (z-1)*(-y - z +1)*(z**2-2*z+1)) / ((z-1)*(z**2-2*z+1))

    phi7  = 4*z*(-x*y + (z-1)*(-x - y - z + 1)) / (z-1)

    phi8  = 4*x*y*(-2*x*y - 2*x*z + 2*x - y*z + y - z**2 + 2*z - 1) / (z**2 - 2*z +1)

    phi9  = 4*x*z*(y + z - 1) / (z-1)

    phi10 = 4*x*y*(-2*x*y - x*z + x - 2*y*z + 2*y - z**2 + 2*z - 1) / (z**2 - 2*z +1)

    phi11 = 4*y*z*(x + z - 1) / (z-1)

    phi12 = -4*x*y*z / (z-1)

    phi13 = 16*x*y*(x*y + x*z - x + y*z - y + z**2 - 2*z + 1) / (z**2 - 2*z +1)

    basis_2 = [sp.simplify(p) for p in shuffle([
        phi0, phi1, phi2, phi3, phi4,
        phi5, phi6, phi7, phi8, phi9,
        phi10, phi11, phi12, phi13
    ],swapmap2)]

    bases=[basis_1, basis_2]

    bases_grad = []

    #######
    bletter="pyramid"
    dim=3
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=str, help="path to the output folder")
    args=parser.parse_args()

    cpp = f"#include \"auto_pyramid_bases.hpp\""
    cpp = cpp + "\n#include \"auto_b_bases.hpp\""  # ?
    cpp = cpp + "\n#include \"p_n_bases.hpp\""  # ?
    cpp = cpp + "\n\n\n" \
        "namespace polyfem {\nnamespace autogen " + "{\nnamespace " + "{\n"

    hpp = "#pragma once\n\n#include <Eigen/Dense>\n#include <cassert>\n"
    hpp = hpp + "\nnamespace polyfem {\nnamespace autogen " + "{\n"

    suffix = "3d"
    unique_nodes = f"void {bletter}_nodes_{suffix}" + \
        f"(const int o, Eigen::MatrixXd &val)"

    unique_fun = f"void {bletter}_basis_value_{suffix}" + \
        f"(const int o, const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &val)"
    dunique_fun = f"void {bletter}_grad_basis_value_{suffix}" + \
        f"(const int o, const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &val)"

    hpp = hpp + unique_nodes + ";\n\n"

    hpp = hpp + unique_fun + ";\n\n"
    hpp = hpp + dunique_fun + ";\n\n"

    unique_nodes = unique_nodes + f"{{\nswitch(o)" + "{\n"

    unique_fun = unique_fun + "{\n"
    dunique_fun = dunique_fun + "{\n"

    unique_fun = unique_fun + f"\nswitch(o)" + "{\n"
    dunique_fun = dunique_fun + f"\nswitch(o)" + "{\n"

    for j, order in enumerate(orders):
        print(j,order)
        ni = real_nodes[j]
        bi = bases[j]
        print(ni)

        # nodes code gen
        nodes = f"void {bletter}_{order}_nodes_{suffix}(Eigen::MatrixXd &res) {{\n res.resize(" + str(
            len(bi)) + ", " + str(dim) + "); res << \n"
        unique_nodes = unique_nodes + f"\tcase {order}: " + \
            f"{bletter}_{order}_nodes_{suffix}(val); break;\n"

        
        for ii in range(len(bi)):
            print(ii)
            nodes = nodes + ccode(ni[ii, 0]) + ", " + ccode(ni[ii, 1]) + (
                (", " + ccode(ni[ii, 2])) if dim == 3 else "") + ",\n"
        nodes = nodes[:-2]
        nodes = nodes + ";\n}"

        # bases code gen
        func = f"void {bletter}_{order}_basis_value_{suffix}" + \
            "(const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &result_0)"
        dfunc = f"void {bletter}_{order}_basis_grad_value_{suffix}" + \
            "(const int local_index, const Eigen::MatrixXd &uv, Eigen::MatrixXd &val)"

        unique_fun = unique_fun + \
            f"\tcase {order}: {bletter}_{order}_basis_value_{suffix}(local_index, uv, val); break;\n"
        dunique_fun = dunique_fun + \
            f"\tcase {order}: {bletter}_{order}_basis_grad_value_{suffix}(local_index, uv, val); break;\n"

        # hpp = hpp + func + ";\n"
        # hpp = hpp + dfunc + ";\n"

        base = "auto x=uv.col(0).array();\nauto y=uv.col(1).array();\nauto z=uv.col(2).array();"
        base = base + "\n\n"
        dbase = base

        base = base + "switch(local_index){\n"
        dbase = dbase + \
            "val.resize(uv.rows(), uv.cols());\n Eigen::ArrayXd result_0(uv.rows());\n" + \
            "switch(local_index){\n"

        for i in range(len(bi)):
            base = base + "\tcase " + str(i) + ": {" + pretty_print.C99_print(
                simplify(bi[i])).replace(" = 1;", ".setOnes();") + "} break;\n"
            dbase = dbase + "\tcase " + str(i) + ": {" + \
                "{" + pretty_print.C99_print(simplify(diff(bi[i], x))).replace(" = 0;", ".setZero();").replace(" = 1;", ".setOnes();").replace(" = -1;", ".setConstant(-1);") + "val.col(0) = result_0; }" \
                "{" + pretty_print.C99_print(simplify(diff(bi[i], y))).replace(" = 0;", ".setZero();").replace(
                    " = 1;", ".setOnes();").replace(" = -1;", ".setConstant(-1);") + "val.col(1) = result_0; }"
            dbase = dbase + "{" + pretty_print.C99_print(simplify(diff(bi[i], z))).replace(" = 0;", ".setZero();").replace(
                " = 1;", ".setOnes();").replace(" = -1;", ".setConstant(-1);") + "val.col(2) = result_0; }"

            dbase = dbase + "} break;\n"

        base = base + "\tdefault: assert(false);\n}"
        dbase = dbase + "\tdefault: assert(false);\n}"

        cpp = cpp + func + "{\n\n"
        cpp = cpp + base + "}\n"

        cpp = cpp + dfunc + "{\n\n"
        cpp = cpp + dbase + "}\n\n\n"

        cpp = cpp + nodes + "\n\n\n"


    unique_nodes = unique_nodes + "\tdefault: assert(false);\n}}"

    unique_fun = unique_fun + "\tdefault: assert(false);\n}}"
    dunique_fun = dunique_fun + "\tdefault: assert(false);\n}}"

    cpp = cpp + "}\n\n" + unique_nodes + "\n" + unique_fun + \
        "\n\n" + dunique_fun + "\n" 
    hpp = hpp + "\n"

    hpp = hpp + \
        f"\nstatic const int MAX_{bletter.capitalize()}_BASES = {max(orders)};\n"

    cpp = cpp + "\n}}\n"
    hpp = hpp + "\n}}\n"

    path = os.path.abspath(args.output)

    print("saving...")
    with open(os.path.join(path, f"auto_{bletter}_bases.cpp"), "w") as file:
        file.write(cpp)

    with open(os.path.join(path, f"auto_{bletter}_bases.hpp"), "w") as file:
        file.write(hpp)

    print("done!")

