#ifndef S0S_EULER_H
#define S0S_EULER_H
#pragma once

#include <cassert> // assert

namespace s0s {

template<typename TypeVector, template<typename...> typename TypeView>
class SolverEuler {
    public:
        template<typename TypeFunction>
        static void step(const TypeFunction& f, double* pX, std::size_t xSize, const double& t,  const double& dt) {
            assert(("dt should be greater than 0.0. Did you initialize dt ?", dt > 0.0));
            TypeView<TypeVector> x(pX, xSize);
            x += f(pX, t) * dt;
        }
};

}

#endif
