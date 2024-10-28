#ifndef S0S_NEWTON_RAPHSON_H
#define S0S_NEWTON_RAPHSON_H
#pragma once

#include <cmath> // abs

namespace s0s {

template<typename TypeVector, template<typename...> typename TypeView>
class SolverNewtonRaphson {
	// Source : https://en.wikipedia.org/wiki/Newton%27s_method
	public:
		template<typename TypeFunction>
		static void step(const TypeFunction& f, const TypeFunction& df_dx, double* pX, std::size_t xSize) {
			TypeView<TypeVector> x(pX, xSize);
			x -= f(pX, xSize) / df_dx(pX, xSize);
		}
		template<typename TypeFunction>
		static void solve(const TypeFunction& f, const TypeFunction& df_dx, double* pX, std::size_t xSize, const double tolerance) {
			TypeView<TypeVector> x(pX, xSize);
			while(std::abs(f(pX, xSize)) > tolerance) { // TODO: max iterations ?
				x -= f(pX, xSize) / df_dx(pX, xSize);
			}
		}
};

}

#endif
