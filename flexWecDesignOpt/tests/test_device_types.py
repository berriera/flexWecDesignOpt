import numpy as np


def test_device():
    from flexWecDesignOpt.device_types import Device
    exp = {'var1name': 0, 'var2name': 0, 'var3name': 0}
    device = Device([0, 0, 0])
    obs = device.substitutions()
    device.geometry()
    assert exp == obs


def test_cylinder():
    from flexWecDesignOpt.device_types import Cylinder
    cylinder = Cylinder([1, 1])
    exp = np.array([1, 1, 0, (1 / 3) ** (1 / 2), (1 / 3) ** (1 / 2), (1 / 2 ** (1 / 2))])
    obs = cylinder.substitutions()
    cylinder.geometry()
    assert np.allclose(exp, obs)


def test_barge_higher_order():
    from flexWecDesignOpt.device_types import FlexibleBarge
    barge = FlexibleBarge([80, 10, 10])
    barge.substitutions()  # TODO: compare dictionary values
    barge.geometry()


def test_flexible_barge():
    from flexWecDesignOpt.device_types import FlexibleBarge
    barge = FlexibleBarge([80, 10, 10])
    barge.substitutions()  # TODO: compare dictionary values
    barge.geometry()


def test_reference_model_3():
    from flexWecDesignOpt.device_types import RM3
    barge = RM3([20, 20, 30, 30])
    barge.substitutions()  # TODO: compare dictionary values
    barge.geometry()
