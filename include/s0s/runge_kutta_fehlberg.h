#ifndef S0S_RUNGE_KUTTA_FEHLBERG_H
#define S0S_RUNGE_KUTTA_FEHLBERG_H
#pragma once

#include <array>
#include <vector>
#include <limits> // numeric limits
#include <cassert> // assert

namespace s0s {

class SolverRungeKuttaFehlberg {
    // Source : https://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method
    public:
        constexpr static const std::array<double, 6> c = {
            0.0, // c(1)
            1.0/4.0, // c(2)
            3.0/8.0, // c(3)
            12.0/13.0, // etc ...
            1.0,
            1.0/2.0
        };
        constexpr static const std::array<std::array<double, 5>, 5> a = {{
            {1.0/4.0, std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN()}, // a(2,1)
            {3.0/32.0     , 9.0/32.0, std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN()}, // a(3, 1), a(3, 2)
            {1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN()}, // a(4, 1), a(4, 2), a(4, 3)
            {439.0/216.0  , -8.0          , 3680.0/513.0  , -845.0/4104.0, std::numeric_limits<double>::signaling_NaN()}, // etc ...
            {-8.0/27.0    , 2.0           , -3544.0/2565.0, 1859.0/4104.0 , -11.0/40.0}
        }};
        constexpr static const std::array<double, 6> b5 = {
            16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0
        };
        constexpr static const std::array<double, 5> b4 = {
            25.0/216.0, 0.0, 1408.0/2565.0 , 2197.0/4104.0  , -1.0/5.0
        };
    public:
        template<typename TypeFunction, typename TypeVector>
        TypeVector operator()(const TypeFunction& f, const TypeVector& x, const double& t,  const double& dt) {
            assert(("dt should be greater than 0.0. Did you initialize dt ?", dt > 0.0));
            std::array<TypeVector, 5> k;
            k[0] = f(x, t);
            TypeVector s = k[0] * b4[0];
            for(size_t i = 1; i < k.size(); i++) {
                TypeVector v = k[0] * a[i - 1][0];
                for(size_t j = 1; j < i; j++) {
                    v += k[j] * a[i - 1][j];
                }
                v = x + v * dt;
                k[i] = f(v, t + c[i] * dt);
                s += k[i] * b4[i];
            }
            return x + s * dt;
        }
};

}

#endif
