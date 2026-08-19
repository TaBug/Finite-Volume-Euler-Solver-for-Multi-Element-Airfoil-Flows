#pragma once
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
using namespace std;

/*Find the maximum value of a double vector
INPUT: vector<double>
OUTPUT: double 
*/
inline double max(vector<double>& v) {
	double output = v[0];
	for (int i = 1; i < v.size(); i++) {
		if (v[i] > output) {
			output = v[i];
		}
	}
	return output;
}

inline vector<double> absolute(const vector<double>& v) {
	vector<double> output(v.size());
	for (int i = 0; i < v.size(); i++) {
		if (v[i] < 0) {
			output[i] = -v[i];
		}
		else {
			output[i] = v[i];
		}
	}
	return output;
}

/*Element-wise sum of two equal-length vectors
INPUT: vector<double> a, b (must be the same size)
OUTPUT: vector<double> a + b
*/
inline std::vector<double> addVectors(const std::vector<double>& a, const std::vector<double>& b) {
	// Callers index the returned vector without checking its length, so a size
	// mismatch must abort here rather than return an empty vector.
	assert(a.size() == b.size() && "addVectors: vectors have different sizes");

	std::vector<double> out(a.size());
	for (std::size_t i = 0; i < a.size(); i++) {
		out[i] = a[i] + b[i];
	}
	return out;
}

/*Element-wise difference of two equal-length vectors
INPUT: vector<double> a, b (must be the same size)
OUTPUT: vector<double> a - b
*/
inline std::vector<double> subtractVectors(const std::vector<double>& a, const std::vector<double>& b) {
	// Callers index the returned vector without checking its length, so a size
	// mismatch must abort here rather than return an empty vector.
	assert(a.size() == b.size() && "subtractVectors: vectors have different sizes");

	std::vector<double> out(a.size());
	for (std::size_t i = 0; i < a.size(); i++) {
		out[i] = a[i] - b[i];
	}
	return out;
}