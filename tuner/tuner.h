/*
Copyright (C) 2014 Daniel Dyer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef TUNER_H
#define TUNER_H

#include <stdint.h>
#include <cmath>

#include "ceres/ceres.h"
#include "types.h"
#include "state.h"

#define TUNER_HORIZON_LENGTH 100

class DynamicsCostFunctor {
	const StateVector* reference_state;
	const ControlVector* reference_control;

public:
	DynamicsCostFunctor(const StateVector* s, const ControlVector* c):
		reference_state(s),
		reference_control(c) {}

	bool operator()(const double* const params, double* residuals) const;
};

#ifdef __cplusplus
extern "C" {
#endif

void run_tuning(
	const real_t (&s)[TUNER_HORIZON_LENGTH*(UKF_STATE_DIM+1)],
	const real_t (&c)[TUNER_HORIZON_LENGTH*UKF_CONTROL_DIM]);

uint32_t get_horizon_length() { return TUNER_HORIZON_LENGTH; }

#ifdef __cplusplus
}
#endif

#endif
