#ifndef S0S_PARALLEL_H
#define S0S_PARALLEL_H
#pragma once

// std
#include <algorithm> // for_each
#include <execution> // execution (parallel algorithms)

namespace s0s {

template<typename TypeSolver>
class SolverParallel {
    public:
        SolverParallel() {
        }

        template<typename TypeFunction, typename TypeVector>
        TypeVector operator()(const TypeFunction& f, const TypeVector& x, const double& t, const double& dt) {
            TypeVector y = x;
            std::for_each(std::execution::par_unseq, y.begin(), y.end(), [this, f, t, dt](double& value) { 
                value = solver(f, value, t, dt);
            });
            return y;
        }
    public:
        TypeSolver solver;
};

}

#endif
