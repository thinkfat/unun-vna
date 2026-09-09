"""Synthetic RF-model sanity tests.

These tests intentionally do not access a VNA.  They verify that the equations
used by analysis.py reproduce two exact reference cases for which the expected
answer is known analytically.
"""

import unittest

import numpy as np

from unun_vna.analysis import (
    ideal_transformer_reference,
    series_fixture_impedance,
    unun_efficiency,
    vswr_from_s11,
)


class MathTests(unittest.TestCase):
    def test_series_fixture_2400(self):
        """A perfect 2400-ohm series element must be recovered exactly."""
        z0 = 50.0
        z = 2400.0 + 0j

        # Closed-form S-parameters of one series impedance between equal Z0
        # ports.  These are synthetic, noise-free reference data.
        s11 = np.array([z / (2 * z0 + z)])
        s21 = np.array([2 * z0 / (2 * z0 + z)])

        result = series_fixture_impedance(s11, s21, z0)

        self.assertAlmostEqual(
            result["z_series"][0].real, 2400.0, places=9
        )
        self.assertAlmostEqual(
            result["z_series"][0].imag, 0.0, places=9
        )
        self.assertAlmostEqual(
            abs(result["identity_error"][0]), 0.0, places=12
        )

    def test_ideal_1_to_49_efficiency(self):
        """A lossless 1:49 transformer must calculate as 100% efficient."""
        # With 2400 ohms in series with the 50-ohm receiving port, the complete
        # high-side load is 2450 ohms.  A lossless 1:49 transformer is perfectly
        # matched from 50 ohms to that load, therefore S11=0.
        s11 = np.array([0 + 0j])

        # The 7:1 voltage ratio means only 1/49 of the high-side load power is
        # dissipated in the 50-ohm Port-2 termination, hence |S21|=1/7.
        s21 = np.array([1 / 7 + 0j])
        z = np.array([2400 + 0j])

        result = unun_efficiency(s11, s21, z, 50.0)

        self.assertAlmostEqual(result["efficiency"][0], 1.0, places=12)
        self.assertAlmostEqual(result["loss_db"][0], 0.0, places=12)

    def test_ideal_1_to_9_efficiency(self):
        """An ideal 1:9 transformer also evaluates to 100%."""
        s11 = np.array([0 + 0j])
        s21 = np.array([1 / 3 + 0j])
        z = np.array([400 + 0j])

        result = unun_efficiency(s11, s21, z, 50.0)

        self.assertAlmostEqual(result["efficiency"][0], 1.0, places=12)
        self.assertAlmostEqual(result["loss_db"][0], 0.0, places=12)

    def test_ideal_reference_1_to_9(self):
        """The ideal-reference model reproduces a matched 1:9 case."""
        ref = ideal_transformer_reference(
            np.array([400 + 0j]),
            impedance_ratio=9.0,
            z0=50.0,
        )

        self.assertAlmostEqual(abs(ref["gamma_in"][0]), 0.0, places=12)
        self.assertAlmostEqual(ref["s21_magnitude"][0], 1 / 3, places=12)
        self.assertAlmostEqual(ref["target_series_resistance"], 400.0, places=12)
        self.assertAlmostEqual(ref["turns_ratio"], 3.0, places=12)

    def test_ideal_reference_1_to_64(self):
        """The ideal-reference model is generic beyond 1:49."""
        ref = ideal_transformer_reference(
            np.array([3150 + 0j]),
            impedance_ratio=64.0,
            z0=50.0,
        )

        self.assertAlmostEqual(abs(ref["gamma_in"][0]), 0.0, places=12)
        self.assertAlmostEqual(ref["s21_magnitude"][0], 1 / 8, places=12)
        self.assertAlmostEqual(ref["target_series_resistance"], 3150.0, places=12)
        self.assertAlmostEqual(ref["turns_ratio"], 8.0, places=12)

    def test_non_nominal_fixture_is_used_as_measured(self):
        """Efficiency accepts an arbitrary complex measured fixture value."""
        s11 = np.array([0.12 + 0.03j])
        s21 = np.array([0.18 - 0.02j])
        z = np.array([731.0 + 13.0j])

        result = unun_efficiency(s11, s21, z, 50.0)
        expected = (
            abs(s21[0]) ** 2
            * (1.0 + z[0].real / 50.0)
            / (1.0 - abs(s11[0]) ** 2)
        )

        self.assertAlmostEqual(result["efficiency"][0], expected, places=15)

    def test_vswr_perfect_match(self):
        """S11=0 must give VSWR=1."""
        self.assertAlmostEqual(vswr_from_s11(np.array([0j]))[0], 1.0, places=12)

    def test_vswr_known_reflection(self):
        """|Gamma|=1/3 gives VSWR=2."""
        self.assertAlmostEqual(
            vswr_from_s11(np.array([1 / 3 + 0j]))[0], 2.0, places=12
        )

    def test_vswr_total_reflection(self):
        """|Gamma|=1 gives infinite VSWR."""
        self.assertTrue(np.isinf(vswr_from_s11(np.array([1 + 0j]))[0]))


if __name__ == "__main__":
    unittest.main()
