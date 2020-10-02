from flexWecDesignOpt.analysis import boundary_condition_frequency_solver
import math
import numpy as np


def linear_equation(x, a):
    return 2 * (x - a)


def quadratic_equation(x):
    return (x - 0.5) * (x - 1.25)


def cubic_equation(x):
    return (x - 0.5) * (x - 1.25) * (x - 9.4)


def sine_wave_equation(x):
    return math.sin(x)


def free_free_barge_equation(x):
    return math.cos(x) * math.cosh(x) - 1


def flexible_tube_mode_type_1(w, L, Di, Ts, rho):
    lk = math.sqrt(((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) - 1))
    uk = math.sqrt(((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) + 1))
    return (lk * L / 2) * math.tanh(uk * L / 2) - (uk * L / 2) * math.tan(lk * L / 2)
    # TODO: look into if fabs from Jennifer's code is ABSOLUTEly necessary


def flexible_tube_mode_type_2(w, L, Di, Ts, rho, r_s, M, K_a):
    S_s = math.pi * (r_s ** 2)
    lk = math.sqrt(((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) - 1))
    uk = math.sqrt(((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) + 1))
    return (uk * L / 2) * math.tanh(uk * L / 2) + (lk * L / 2) * math.tan(lk * L / 2) \
           - (((w ** 2) * rho * S_s * L) / (-M * (w ** 2) + 2 * K_a)) * (uk / lk + lk / uk) \
           * math.tanh(uk * L / 2) * math.tan(lk * L / 2)


def test_0_roots_linear():
    exp = []
    obs = boundary_condition_frequency_solver(linear_equation, 0)
    assert np.allclose(exp, obs)


def test_1_roots_linear__1():
    exp = [1]
    obs = boundary_condition_frequency_solver(linear_equation, 1, extra_args=1)
    assert np.allclose(exp, obs)


def test_1_roots_linear__2():
    exp = [math.pi]
    obs = boundary_condition_frequency_solver(linear_equation, 1, extra_args=math.pi)
    assert np.allclose(exp, obs)


def test_2_roots_quadratic():
    exp = [0.5, 1.25]
    obs = boundary_condition_frequency_solver(quadratic_equation, 2, h=0.3)
    assert np.allclose(exp, obs)


def test_3_roots_quadratic():
    exp = [0.5, 1.25]
    obs = boundary_condition_frequency_solver(quadratic_equation, 3, h=0.3, freq_limit=15)
    assert np.allclose(exp, obs)


def test_3_roots_cubic():
    exp = [0.5, 1.25, 9.4]
    obs = boundary_condition_frequency_solver(cubic_equation, 3, h=0.3)
    assert np.allclose(exp, obs)


def test_10_roots_sine_wave():
    exp = [math.pi, 2 * math.pi]
    obs = boundary_condition_frequency_solver(sine_wave_equation, 2)
    assert np.allclose(exp, obs)


def test_4_roots_flexible_barge():
    exp = [4.7300, 7.8532, 10.9956, 14.1372]
    obs = boundary_condition_frequency_solver(free_free_barge_equation, 4)
    assert np.allclose(exp, obs)


def test_5_roots_flexible_tube__mode_type_1():
    exp = [2.1114360764249511, 4.5933899504176336, 7.7089354146509326, 11.6086922395909120, 16.3723249238653743]
    obs = boundary_condition_frequency_solver(flexible_tube_mode_type_1, 5, extra_args=(10, 1.12e-4, 1.8770e4, 1000))
    assert np.allclose(exp, obs)


def test_5_roots_flexible_tube__mode_type_2():
    exp = [0.5489434835683205, 1.4244105181770936, 1.9121722822149714, 3.9035022567988560, 6.6024612335389854]
    obs = boundary_condition_frequency_solver(flexible_tube_mode_type_2, 5, h=0.5, freq_limit=8,
                                              extra_args=(10, 1.12e-4, 1.8770e4, 1000, 0.274, 311.7, 510))
    assert np.allclose(exp, obs, rtol=1e-2)
