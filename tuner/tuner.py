import os
import sys
import copy
import math
from ctypes import *

import ukf
ukf.init()

_tuner = None
_HORIZON_LENGTH = None

def init():
    global _tuner, _HORIZON_LENGTH

    # Load the requested library and determine configuration parameters
    lib = os.path.join(
        os.path.dirname(__file__), "libtuner.dylib")

    _tuner = cdll.LoadLibrary(lib)

    _tuner.get_horizon_length.argtype = []
    _tuner.get_horizon_length.restype = c_long

    _HORIZON_LENGTH = int(_tuner.get_horizon_length())

    _tuner.run_tuning.argtypes = [
        POINTER(ukf._REAL_T
            * ((ukf._STATE_DIM + 1) * _HORIZON_LENGTH)),
        POINTER(ukf._REAL_T
            * (ukf._CONTROL_DIM * _HORIZON_LENGTH))]
    _tuner.run_tuning.restype = None

#WGS84 reference ellipsoid constants
A = 6378137.0
B = 6356752.314245
E2 = 0.0066943799901975848
A2 = A**2 #to speed things up a bit
B2 = B**2

def lla2ecef(lat, lon, alt):
    # Convert geodetic lat/long/alt to ECEF
    sin_lat = math.sin(lat)
    sin_lon = math.sin(lon)
    cos_lat = math.cos(lat)
    cos_lon = math.cos(lon)

    ntheta = A / math.sqrt(1.0 - E2 * sin_lat * sin_lat)
    x = (ntheta + alt) * cos_lat * cos_lon
    y = (ntheta + alt) * cos_lat * sin_lon
    z = (ntheta * (1.0 - E2) + alt) * sin_lat

    return x, y, z

def ecef2ned(x, y, z, lat, lon):
    sin_lat = math.sin(lat)
    sin_lon = math.sin(lon)
    cos_lat = math.cos(lat)
    cos_lon = math.cos(lon)

    # Convert ECEF to NED
    n = (-sin_lat * cos_lon * x) + (-sin_lat * sin_lon * y) + (cos_lat * z)
    e = (-sin_lon * x) + (cos_lon * y)
    d = (-cos_lat * cos_lon * x) + (-cos_lat * sin_lon * y) + (-sin_lat * z)

    return n, e, d

def optimise(horizon):
    state_horizon = []
    control_horizon = []

    for s in horizon:
        state = (
            float(s["lat"]),
            float(s["lon"]),
            float(s["alt"]),
            float(s["v_n"]),
            float(s["v_e"]),
            float(s["v_d"]),
            0, 0, 0,
            float(s["att_x"]),
            float(s["att_y"]),
            float(s["att_z"]),
            float(s["att_w"]),
            float(s["angular_v_x"]),
            float(s["angular_v_y"]),
            float(s["angular_v_z"]),
            0, 0, 0,
            float(s["wind_n"]),
            float(s["wind_e"]),
            float(s["wind_d"]),
            0, 0, 0
        )
        state_horizon.extend(state)

        control = (
            float(s["pwm_t"]),
            float(s["pwm_l"]),
            float(s["pwm_r"])
        )
        control_horizon.extend(control)

    _tuner.run_tuning(
        (ukf._REAL_T
            * ((ukf._STATE_DIM + 1) * _HORIZON_LENGTH))
            (*state_horizon),
        (ukf._REAL_T
            * (ukf._CONTROL_DIM * _HORIZON_LENGTH))
            (*control_horizon))

if __name__ == "__main__":
    init()

    fields = [
        "time",
        "video_roll", "video_pitch",
        "roll", "pitch",
        "lat", "lon", "alt",
        "v_n", "v_e", "v_d",
        "att_x", "att_y", "att_z", "att_w",
        "angular_v_x", "angular_v_y", "angular_v_z",
        "wind_n", "wind_e", "wind_d",
        "pwm_t", "pwm_l", "pwm_r"]
    horizons = []
    current_horizon = []
    last_time = 0
    for line in sys.stdin:
        try:
            readings = dict(zip(fields, map(float, line.strip("\n").split("\t"))))
        except ValueError:
            current_horizon = []
            continue

        # Get every 20th millisecond.
        if int(readings["time"])/20 > last_time/20:
            # Ensure video roll, pitch exist and are within 5 degrees of
            # UKF values.
            if abs(float(readings["video_roll"]) - float(readings["roll"])) > 15 or \
            abs(float(readings["video_pitch"]) - float(readings["pitch"])) > 15:
                print "Ignoring sample at %d ms" % readings["time"]
                current_horizon = []
                last_time = int(readings["time"])
                continue

            # Add to the current horizon, or start a new one if it's full.
            current_horizon.append(readings)
            if len(current_horizon) == _HORIZON_LENGTH:
                horizons.append(copy.deepcopy(current_horizon))
                current_horizon = []
                print "Finished horizon ending at %d ms" % readings["time"]

        last_time = int(readings["time"])

    # TODO: Pick a horizon to optimise.
    optimise(horizons[len(horizons)/2])