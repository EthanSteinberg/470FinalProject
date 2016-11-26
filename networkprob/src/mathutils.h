#pragma once

#include <cmath>

/**
 * Compute the factorial.
 */
inline int factorial(int n) {
	int result = 1;

	for (int i = n; i > 1; i --) {
		result *= i;
	}

	return result;
}

/**
 * Compute -1^n.
 */
inline int pow1(int num) {
	if (num % 2 == 0) {
		return 1;
	} else {
		return -1;
	}
}

static double puvArray[8][8][8] = {};

/**
 * Compute all the puv values.
 */
inline int initPuvArray() {
	for (int u = 1; u < 8; u++) {
		for (int v = u; v > 0; v--) {
			for (int k = v; k <= u; k++) {
				double term1 = (2*k -1) * std::pow(-1, k-v)/((double)factorial(v) * factorial(k-v) * (v + k -1));


				double term2 = 1;

				for (int y =0; y <= k-1; y++) {
					term2 *= (v + y)*(u -y)/( (double) (u + y));
				}

				puvArray[u][v][k] = term1 * term2;
			}
		}
	}
	return 1;
}

static int dummy = initPuvArray();

/**
 * Compute the puv function.
 */
inline double puv(int u, int v, double T) {
	if (std::isinf(T)) {
		if (v == 1) {
			return 1.0;
		} else {
			return 0.0;
		}
	}

	if (v == 0 && u == 0) {
		return 1.0;
	}

	double sum = 0;

	for (int k = v; k <=u; k++) {
		sum += std::exp(-k*(k-1) * T/2.0) * puvArray[u][v][k];
	}

	return sum;
}

/**
 * Compute the derivative of the puv function.
 */
inline double derivatePuv(int u, int v, double T) {
	if (std::isinf(T) || (v == 0 && u == 0)) {
		return 0.0;
	}

	double sum = 0;

	for (int k = v; k <=u; k++) {
		sum += -k*(k-1) * 1.0/2.0 * std::exp(-k*(k-1) * T/2.0) * puvArray[u][v][k];
	}

	return sum;
}

static double numberOfOptionsArray[8][8] = {};

/**
 * Precompute the options array.
 */
inline int initNumberOfOptionsArray() {
	for (int starting = 0; starting < 8; starting++) {
		for (int ending = 0; ending < 8; ending++) {
			double product = 1.0;

			for (int i=0; i < (starting - ending); i++) {
				product *= factorial(starting -i)/(2.0 * factorial(starting -i - 2));
			}

			numberOfOptionsArray[starting][ending] = product;
		}
	}

	return 0;
}

static int dummy2 = initNumberOfOptionsArray();

/**
 * Get the number of options for a path with starting and ending genelogies.
 */
inline double getNumberOfOptions(int starting, int ending) {
	return numberOfOptionsArray[starting][ending];
}