#ifndef S0S_GRADIENT_DESCENT_H
#define S0S_GRADIENT_DESCENT_H
#pragma once

#include <cmath>

namespace s0s {

template<typename TypeVector, template<typename...> typename TypeView>
class SolverGradientDescent {
	// Source : https://en.wikipedia.org/wiki/Gradient_descent
	public:
		template<typename TypeFunction>
		static void solve(const TypeFunction& grad_x, double* pX, std::size_t xSize, const double tolerance) {
			// view
			TypeView<TypeVector> x1(pX, xSize);
			// initial step
			double step = 1.0;
			TypeVector x0 = x1;
			TypeVector grad_x0 = grad_x(x0.data(), x0.size());
			x1 -= step * grad_x0;
			TypeVector grad_x1 = grad_x(x1.data(), x1.size());
			// iterations
			while(grad_x1.norm() > tolerance) { // TODO: max iterations ?
				step = std::abs((grad_x1 - grad_x0).dot(x1 - x0)) / (grad_x1 - grad_x0).squaredNorm();
				x0 = x1;
				grad_x0 = grad_x1;
				x1 -= step * grad_x0;
				grad_x1 = grad_x(x1.data(), x1.size());
			}
		}
};

}

#endif
