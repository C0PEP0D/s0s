#ifndef S0S_EULER_H
#define S0S_EULER_H
#pragma once

#include <cassert> // assert

namespace s0s {

class SolverEuler {
    public:
        template<typename TypeFunction, typename TypeVector>
        TypeVector operator()(const TypeFunction& f, TypeVector& x, const double& t,  const double& dt) {
            assert(("dt should be greater than 0.0. Did you initialize dt ?", dt > 0.0));
            return x + f(x, t) * dt;
        }
};

}

#endif
