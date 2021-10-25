// Standard
#include <iostream> // cout, endl
#include <cmath> // exp, abs
#include <vector>
// Thirdparties
#include <Eigen/Dense> // Eigen
// Lib
#include "s0s/solve.h" // s0s
// Simple
#include "func.h" // Func

constexpr unsigned int DIM = 1000;
using TypeScalar = double;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;

using TypeFunction = Func;

template<typename ...Args>
using TypeContainer = std::vector<Args...>;

int main () { 
    // Init
    TypeVector x0 = TypeVector::Constant(1.0);
    TypeScalar t0 = 0.0;
    TypeScalar dt = 1e-4;
    TypeScalar tMax = 1e0;
    TypeFunction f;
    // computation
    TypeContainer<TypeVector> x;
    std::vector<TypeScalar> t;
    std::tie(x, t) = s0s::solve(f, x0, t0, tMax, dt);
    // out
    std::cout << "\n";
    std::cout << "Solver solved exp(" << t0 << " -> " << tMax << ") = " << std::endl;
    std::cout << "\n";
    for(const TypeVector& v : x) {
        std::cout << v(0) << " ";
    }
    std::cout << "\n";
    std::cout << std::endl;
}
