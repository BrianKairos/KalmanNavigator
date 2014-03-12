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

#include <cmath>

#include "ceres/ceres.h"
#include "types.h"
#include "state.h"
#include "integrator.h"
#include "dynamics.h"
#include "debug.h"
#include "tuner.h"

using namespace ceres;

Problem problem;

double params[16];

bool DynamicsCostFunctor::operator()(
const double* const params, double* residuals) const {
	uint32_t i, j;

	real_t params_real[16];

	for(i = 0; i < 16; i++) {
		params_real[i] = params[i];
	}

	X8DynamicsModel dynamics_model = X8DynamicsModel(params_real);
	State current;
	IntegratorRK4 integrator;

	/* Set initial state. */
	current << reference_state[0];

	/* Calculate cost function residuals. */
	for(i = 0; i < TUNER_HORIZON_LENGTH; i++) {
		Eigen::Map<Eigen::Matrix<double, 13, 1> > residual_map(
			&residuals[i*13]);

		/* Calculate residual. */
		residual_map.segment<6>(0) = (current.segment<6>(0)
			- reference_state[i].segment<6>(0)).cast<double>();
		residual_map.segment<7>(6) = (current.segment<7>(9)
			- reference_state[i].segment<7>(9)).cast<double>();

		/* Evaluate dynamics model and integrate. */
		AccelerationVector temp = dynamics_model.evaluate(
			current, reference_control[i]);
		current.segment<3>(6) = temp.segment<3>(0);
		current.segment<3>(16) = temp.segment<3>(3);
		current = integrator.integrate(State(current), (1.0/50.0));
	}

	return true;
}

/* Initialise the tuner with a reference trajectory. */
void run_tuning(
const real_t (&s)[TUNER_HORIZON_LENGTH*(UKF_STATE_DIM+1)],
const real_t (&c)[TUNER_HORIZON_LENGTH*UKF_CONTROL_DIM]) {
	uint32_t i;
	StateVector reference_state[TUNER_HORIZON_LENGTH];
	ControlVector reference_control[TUNER_HORIZON_LENGTH];

	for(i = 0; i < TUNER_HORIZON_LENGTH; i++) {
		Eigen::Map<const StateVector> reference_state_map(
			&s[i*(UKF_STATE_DIM+1)]);
		Eigen::Map<const Vector3r> reference_control_map(
			&c[i*UKF_CONTROL_DIM]);

		reference_state[i] = reference_state_map;
		reference_control[i] = reference_control_map;
	}

	/*
	Dimensions of residuals includes just position, velocity, attitude,
	angular velocity.
	*/
	CostFunction* cost_function =
		new NumericDiffCostFunction<
		DynamicsCostFunctor,
		CENTRAL,
		TUNER_HORIZON_LENGTH*13,
		16>(new DynamicsCostFunctor(reference_state, reference_control));

	/* Set default dynamics model parameters. */
    params[0] = 0.8;
    params[1] = 0.15;
    params[2] = 0.05;
    params[3] = 0.7;
    params[4] = 0.3;
    params[5] = 0.0;
    params[6] = 0.015;
    params[7]= -0.003;
    params[8] = -0.005;
    params[9] = 0.02;
    params[10] = 0.08;
    params[11] = 0.01;
    params[12] = -0.02;
    params[13] = -0.05;
    params[14] = -0.001;
    params[15] = 0.0025;

	problem.AddResidualBlock(cost_function, NULL, params);

	Solver::Options options;
	options.max_num_iterations = 1000;
	options.linear_solver_type = DENSE_QR;
	options.use_nonmonotonic_steps = true;
	options.max_consecutive_nonmonotonic_steps = 20;
	options.function_tolerance=1e-12;
	options.parameter_tolerance=1e-12;
	options.numeric_derivative_relative_step_size = 1e0;
	options.min_relative_decrease = 1e-1;
	options.minimizer_progress_to_stdout = true;

	Solver::Summary summary;
	Solve(options, &problem, &summary);
	std::cout << summary.BriefReport() << "\n";
	std::cout << "Final params:\n";

	for(int i = 0; i < 16; i++) {
		std::cout << params[i] << "\n";
	}
}
