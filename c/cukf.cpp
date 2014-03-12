/*
Copyright (C) 2013 Daniel Dyer

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

#include "types.h"
#include "state.h"
#include "integrator.h"
#include "dynamics.h"
#include "sensors.h"
#include "ukf.h"

#include "cukf.h"

static real_t dynamics_params[16] = {
    0.8, 0.15,
    0.05, 0.7,
    0.3,
    0.0, 0.015, -0.003, -0.005,
    0.02, 0.08, 0.01,
    -0.02, -0.05, -0.001,
    0.0025
};

static IOBoardModel model = IOBoardModel(
    Quaternionr(1, 0, 0, 0), /* NOTE: W, x, y, z */
    Vector3r(0, 0, 0),
    Quaternionr(1, 0, 0, 0), /* NOTE: W, x, y, z */
    Quaternionr(1, 0, 0, 0), /* NOTE: W, x, y, z */
    Vector3r(1, 0, 0));
static UnscentedKalmanFilter ukf = UnscentedKalmanFilter(model);
static CentripetalModel centripetal_model = CentripetalModel();
static CustomDynamicsModel custom_model = CustomDynamicsModel();
static X8DynamicsModel x8_model = X8DynamicsModel(dynamics_params);

void ukf_init(void) {

}

void ukf_set_position(real_t lat, real_t lon, real_t alt) {
    State temp = ukf.get_state();
    temp.position() << lat, lon, alt;
    ukf.set_state(temp);
}

void ukf_set_velocity(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.velocity() << x, y, z;
    ukf.set_state(temp);
}

void ukf_set_acceleration(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.acceleration() << x, y, z;
    ukf.set_state(temp);
}

void ukf_set_attitude(real_t w, real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.attitude() << x, y, z, w;
    ukf.set_state(temp);
}

void ukf_set_angular_velocity(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.angular_velocity() << x, y, z;
    ukf.set_state(temp);
}

void ukf_set_angular_acceleration(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.angular_acceleration() << x, y, z;
    ukf.set_state(temp);
}

void ukf_set_wind_velocity(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.wind_velocity() << x, y, z;
    ukf.set_state(temp);
}

void ukf_set_gyro_bias(real_t x, real_t y, real_t z) {
    State temp = ukf.get_state();
    temp.gyro_bias() << x, y, z;
    ukf.set_state(temp);
}

void ukf_get_state(struct ukf_state_t *in) {
    State current = ukf.get_state();
    in->position[0] = current.position()[0];
    in->position[1] = current.position()[1];
    in->position[2] = current.position()[2];
    in->velocity[0] = current.velocity()[0];
    in->velocity[1] = current.velocity()[1];
    in->velocity[2] = current.velocity()[2];
    in->acceleration[0] = current.acceleration()[0];
    in->acceleration[1] = current.acceleration()[1];
    in->acceleration[2] = current.acceleration()[2];
    /*
    Eigen stores quaternions in {x, y, z, W}, matching the order of the
    output array
    */
    in->attitude[0] = current.attitude()[0];
    in->attitude[1] = current.attitude()[1];
    in->attitude[2] = current.attitude()[2];
    in->attitude[3] = current.attitude()[3];
    in->angular_velocity[0] = current.angular_velocity()[0];
    in->angular_velocity[1] = current.angular_velocity()[1];
    in->angular_velocity[2] = current.angular_velocity()[2];
    in->angular_acceleration[0] = current.angular_acceleration()[0];
    in->angular_acceleration[1] = current.angular_acceleration()[1];
    in->angular_acceleration[2] = current.angular_acceleration()[2];
    in->wind_velocity[0] = current.wind_velocity()[0];
    in->wind_velocity[1] = current.wind_velocity()[1];
    in->wind_velocity[2] = current.wind_velocity()[2];
    in->gyro_bias[0] = current.gyro_bias()[0];
    in->gyro_bias[1] = current.gyro_bias()[1];
    in->gyro_bias[2] = current.gyro_bias()[2];
}

void ukf_set_state(struct ukf_state_t *in) {
    State temp = ukf.get_state();
    temp <<
        in->position[0],
        in->position[1],
        in->position[2],
        in->velocity[0],
        in->velocity[1],
        in->velocity[2],
        in->acceleration[0],
        in->acceleration[1],
        in->acceleration[2],
        /*
        Eigen stores quaternions in {x, y, z, W}, matching the order of the
        input array
        */
        in->attitude[0],
        in->attitude[1],
        in->attitude[2],
        in->attitude[3],
        in->angular_velocity[0],
        in->angular_velocity[1],
        in->angular_velocity[2],
        in->angular_acceleration[0],
        in->angular_acceleration[1],
        in->angular_acceleration[2],
        in->wind_velocity[0],
        in->wind_velocity[1],
        in->wind_velocity[2],
        in->gyro_bias[0],
        in->gyro_bias[1],
        in->gyro_bias[2];
    ukf.set_state(temp);
}

void ukf_get_state_covariance(
real_t state_covariance[UKF_STATE_DIM*UKF_STATE_DIM]) {
    Eigen::Map<StateCovariance> covariance_map(state_covariance);
    covariance_map = ukf.get_state_covariance();
}

void ukf_get_state_covariance_diagonal(
real_t state_covariance_diagonal[UKF_STATE_DIM]) {
    Eigen::Map< Eigen::Matrix<real_t, UKF_STATE_DIM, 1> > covariance_map =
        Eigen::Map< Eigen::Matrix<real_t, UKF_STATE_DIM, 1>
            >(state_covariance_diagonal);
    covariance_map = ukf.get_state_covariance().diagonal();
}

void ukf_get_state_error(real_t state_error[UKF_STATE_DIM]) {
    Eigen::Map< Eigen::Matrix<real_t, UKF_STATE_DIM, 1> > error_map =
        Eigen::Map< Eigen::Matrix<real_t, UKF_STATE_DIM, 1>
            >(state_error);
    error_map =
        ukf.get_state_covariance().cwiseAbs().rowwise().sum().cwiseSqrt();
}

void ukf_sensor_clear() {
    model.clear();
}

void ukf_sensor_set_accelerometer(real_t x, real_t y, real_t z) {
    model.set_accelerometer(Vector3r(x, y, z));
}

void ukf_sensor_set_gyroscope(real_t x, real_t y, real_t z) {
    model.set_gyroscope(Vector3r(x, y, z));
}

void ukf_sensor_set_magnetometer(real_t x, real_t y, real_t z) {
    model.set_magnetometer(Vector3r(x, y, z));
}

void ukf_sensor_set_gps_position(real_t lat, real_t lon, real_t alt) {
    model.set_gps_position(Vector3r(lat, lon, alt));
}

void ukf_sensor_set_gps_velocity(real_t x, real_t y, real_t z) {
    model.set_gps_velocity(Vector3r(x, y, z));
}

void ukf_sensor_set_pitot_tas(real_t tas) {
    model.set_pitot_tas(tas);
}

void ukf_sensor_set_barometer_amsl(real_t amsl) {
    model.set_barometer_amsl(amsl);
}

void ukf_set_params(struct ukf_ioboard_params_t *in) {
    model = IOBoardModel(
        /*
        accel_orientation, gyro_orientation and mag_orientation are
        {x, y, z, W}; Eigen's Quaternionr constructor is {W, x, y, z}.
        */
        Quaternionr(
            in->accel_orientation[3],
            in->accel_orientation[0],
            in->accel_orientation[1],
            in->accel_orientation[2]),
        Vector3r(
            in->accel_offset[0],
            in->accel_offset[1],
            in->accel_offset[2]),
        Quaternionr(
            in->gyro_orientation[3],
            in->gyro_orientation[0],
            in->gyro_orientation[1],
            in->gyro_orientation[2]),
        Quaternionr(
            in->mag_orientation[3],
            in->mag_orientation[0],
            in->mag_orientation[1],
            in->mag_orientation[2]),
        Vector3r(
            in->mag_field[0],
            in->mag_field[1],
            in->mag_field[2]));

    MeasurementVector covariance(17);
    covariance <<
        in->accel_covariance[0], in->accel_covariance[1],
            in->accel_covariance[2],
        in->gyro_covariance[0], in->gyro_covariance[1],
            in->gyro_covariance[2],
        in->mag_covariance[0], in->mag_covariance[1], in->mag_covariance[2],
        in->gps_position_covariance[0], in->gps_position_covariance[1],
            in->gps_position_covariance[2],
        in->gps_velocity_covariance[0], in->gps_velocity_covariance[1],
            in->gps_velocity_covariance[2],
        in->pitot_covariance,
        in->barometer_amsl_covariance;
    model.set_covariance(covariance);
}

void ukf_choose_dynamics(enum ukf_model_t t) {
    switch(t) {
        case UKF_MODEL_NONE:
            ukf.set_dynamics_model((DynamicsModel *)NULL);
            break;
        case UKF_MODEL_CENTRIPETAL:
            ukf.set_dynamics_model(&centripetal_model);
            break;
        case UKF_MODEL_CUSTOM:
            ukf.set_dynamics_model(&custom_model);
            break;
        case UKF_MODEL_X8:
            ukf.set_dynamics_model(&x8_model);
            break;
        default:
            assert(false);
            break;
    }
}

void ukf_set_custom_dynamics_model(ukf_model_function_t func) {
    assert(func);
    custom_model.set_function(func);
    ukf.set_dynamics_model(&custom_model);
}

void ukf_iterate(float dt, real_t control_vector[UKF_CONTROL_DIM]) {
    ukf.iterate(dt,
        Eigen::Matrix<real_t, UKF_CONTROL_DIM, 1>(control_vector));
}

void ukf_set_process_noise(real_t process_noise_covariance[UKF_STATE_DIM]) {
    Eigen::Map< Eigen::Matrix<real_t, UKF_STATE_DIM, 1> > covariance_map =
        Eigen::Map< Eigen::Matrix<real_t, UKF_STATE_DIM, 1>
            >(process_noise_covariance);
    ProcessCovariance covariance = covariance_map;
    ukf.set_process_noise(covariance);
}

uint32_t ukf_config_get_state_dim() {
    return UKF_STATE_DIM;
}

uint32_t ukf_config_get_measurement_dim() {
    return UKF_MEASUREMENT_DIM;
}

uint32_t ukf_config_get_control_dim() {
    return UKF_CONTROL_DIM;
}

enum ukf_precision_t ukf_config_get_precision() {
#ifdef UKF_SINGLE_PRECISION
    return UKF_PRECISION_FLOAT;
#else
    return UKF_PRECISION_DOUBLE;
#endif
}
