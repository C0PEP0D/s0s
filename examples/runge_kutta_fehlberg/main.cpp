// Standard
#include <iostream> // cout, endl
#include <cmath> // exp, abs
#include <vector>
// Thirdparties
#include <Eigen/Dense> // Eigen
// Lib
#include "s0s/euler.h" // s0s
// Simple
#include "func.h" // Func

constexpr unsigned int DIM = 10;
using TypeScalar = double;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;
template<typename ...Args>
using TypeView = Eigen::Map<Args...>;

using TypeFunction = Func<TypeVector>;

int main () { 
    // Parameters
    TypeScalar dt = 1e-4;
    TypeScalar tMax = 1e0;
    std::size_t nt = tMax / dt;
    TypeFunction f;
    s0s::SolverRungeKuttaFehlberg<TypeVector, TypeView> solver;
    // Init
    TypeVector x = TypeVector::Constant(1.0);
    TypeScalar t = 0.0;
    // computation
    for(unsigned int i = 0; i < nt; i++) {
        solver(f, x.data(), x.size(), t, dt);
        t += dt;
    }
    // out
    std::cout << "\n";
    std::cout << "Solver solved exp(" << 0.0 << " -> " << tMax << ") = " << std::endl;
    std::cout << "\n";
    std::cout << x.transpose() << "\n";
    std::cout << "\n";
    std::cout << std::endl;
}
