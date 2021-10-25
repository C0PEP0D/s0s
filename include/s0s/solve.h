#ifndef S0S_SOLVE_H
#define S0S_SOLVE_H
#pragma once

#include <vector>
#include <cassert> // assert
#include "s0s/runge_kutta_fehlberg.h"
#include "s0s/parallel.h"

namespace s0s {

template<typename TypeFunction, typename TypeVector, template<typename...> class TypeContainer = std::vector>
std::tuple<TypeContainer<TypeVector>, TypeContainer<double>> solve(const TypeFunction& f, const TypeVector& x0, const double& t0, const double& tMax, const double& dt) {
    assert(("tMax should be bigger than t0", t0 < tMax));
    assert(("dt should be greater than 0.0", dt > 0.0));
    // Init solve variable
    TypeVector xSolve = x0;
    double tSolve = t0;
    // Init result array
    TypeContainer<TypeVector> x;
    x.emplace_back(x0);
    TypeContainer<double> t;
    t.emplace_back(t0);
    // Init solver
    SolverParallel<SolverRungeKuttaFehlberg> solver;
    // Computation
    while(tSolve < tMax) {
        // solve
        xSolve = solver(f, xSolve, tSolve, dt);
        tSolve += dt;
        // append results
        x.emplace_back(xSolve);
        t.emplace_back(tSolve);
    }
    return std::make_tuple(x, t);
}

}

#endif
