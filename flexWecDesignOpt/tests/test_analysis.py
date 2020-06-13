from flexWecDesignOpt.analysis import eigenvalue_solver
import math
import numpy as np


def linear_equation(x):
    return 2 * (x - 1)


def quadratic_equation(x):
    return (x - 0.5) * (x - 1.25)


def cubic_equation(x):
    return (x - 0.5) * (x - 1.25) * (x - 9.4)


def sine_wave_equation(x):
    return math.sin(x)


def free_free_barge_equation(x):
    return math.cos(x) * math.cosh(x) - 1


def test_0_roots_linear():
    exp = []
    obs = eigenvalue_solver(linear_equation, 0)
    assert np.allclose(exp, obs)


def test_1_roots_linear():
    exp = [1]
    obs = eigenvalue_solver(linear_equation, 1)
    assert np.allclose(exp, obs)


def test_2_roots_quadratic():
    exp = [0.5, 1.25]
    obs = eigenvalue_solver(quadratic_equation, 2, h=0.3)
    assert np.allclose(exp, obs)


def test_3_roots_quadratic():
    exp = [0.5, 1.25]
    obs = eigenvalue_solver(quadratic_equation, 3, h=0.3, freq_limit=15)
    assert np.allclose(exp, obs)


def test_3_roots_cubic():
    exp = [0.5, 1.25, 9.4]
    obs = eigenvalue_solver(cubic_equation, 3, h=0.3)
    assert np.allclose(exp, obs)


def test_10_roots_sine_wave():
    exp = [math.pi, 2 * math.pi]
    obs = eigenvalue_solver(sine_wave_equation, 2)
    assert np.allclose(exp, obs)


def test_4_roots_flexible_barge():
    exp = [4.7300, 7.8532, 10.9956, 14.1372]
    obs = eigenvalue_solver(free_free_barge_equation, 4)
    assert np.allclose(exp, obs)
