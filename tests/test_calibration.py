"""Synthetic tests for the host-side enhanced-response calibration."""

import cmath
import sys
import types
import unittest

import numpy as np

# cli.py imports pynanovna for hardware access; calibration math itself does not
# need hardware.  Provide a stub so these tests remain host-independent.
sys.modules.setdefault("pynanovna", types.ModuleType("pynanovna"))

from unun_vna.cli import apply_cal, compute_cal


class EnhancedResponseCalibrationTests(unittest.TestCase):
    @staticmethod
    def raw_reflection(gamma, e00, e11, etrack):
        return e00 + etrack * gamma / (1.0 - e11 * gamma)

    def test_high_reflection_dut_is_recovered(self):
        """Enhanced response must recover S21 for a strongly reflecting DUT."""
        # Synthetic but realistic complex error terms.
        e00 = 0.015 + 0.008j
        e11 = 0.085 - 0.030j
        etrack = 0.93 * cmath.exp(0.12j)
        leak0 = 0.0012 + 0.0005j
        leak_r = 0.0040 - 0.0015j
        forward_k = 0.31 * cmath.exp(-0.20j)  # e.g. path including a 10 dB pad

        def m11(gamma):
            return self.raw_reflection(gamma, e00, e11, etrack)

        def leakage(raw11):
            return leak0 + leak_r * raw11

        def gain(gamma):
            return 1.0 / (1.0 - e11 * gamma)

        # SOL standards.
        short11 = m11(-1.0 + 0j)
        open11 = m11(1.0 + 0j)
        load11 = m11(0.0 + 0j)

        # SHORT and OPEN have zero true S21; measured S21 is leakage only.
        short21 = leakage(short11)
        open21 = leakage(open11)

        # Non-perfect THRU reflection deliberately exercises THRU-S11 handling.
        true_thru_s11 = 0.025 - 0.018j
        thru11 = m11(true_thru_s11)
        thru21 = forward_k * gain(true_thru_s11) + leakage(thru11)

        cal = compute_cal(
            np.array([short11]),
            np.array([short21]),
            np.array([open11]),
            np.array([open21]),
            np.array([load11]),
            np.array([thru11]),
            np.array([thru21]),
        )

        # Fixture-like DUT: nearly total reflection, small forward transmission.
        true_s11 = 0.958 + 0.012j
        true_s21 = 0.0405 * cmath.exp(-0.37j)
        raw11 = m11(true_s11)
        raw21 = forward_k * true_s21 * gain(true_s11) + leakage(raw11)

        got_s11, got_s21 = apply_cal(
            np.array([raw11]), np.array([raw21]), *cal
        )

        self.assertAlmostEqual(got_s11[0].real, true_s11.real, places=12)
        self.assertAlmostEqual(got_s11[0].imag, true_s11.imag, places=12)
        self.assertAlmostEqual(got_s21[0].real, true_s21.real, places=12)
        self.assertAlmostEqual(got_s21[0].imag, true_s21.imag, places=12)

    def test_calibration_thru_returns_unity(self):
        """Applying a calibration to its own THRU must give S21=1."""
        e00 = 0.01 + 0.004j
        e11 = 0.06 + 0.02j
        etrack = 0.95 * cmath.exp(-0.08j)
        leak0 = 0.0008 - 0.0003j
        leak_r = 0.003 + 0.001j
        k = 0.30 * cmath.exp(0.1j)

        def m11(g):
            return self.raw_reflection(g, e00, e11, etrack)
        def leak(m):
            return leak0 + leak_r*m
        def gain(g):
            return 1/(1-e11*g)

        short11, open11, load11 = m11(-1), m11(1), m11(0)
        short21, open21 = leak(short11), leak(open11)
        thru_gamma = 0.015 + 0.006j
        thru11 = m11(thru_gamma)
        thru21 = k*gain(thru_gamma) + leak(thru11)

        cal = compute_cal(
            np.array([short11]), np.array([short21]),
            np.array([open11]), np.array([open21]),
            np.array([load11]), np.array([thru11]), np.array([thru21]),
        )
        s11, s21 = apply_cal(np.array([thru11]), np.array([thru21]), *cal)
        self.assertAlmostEqual(s11[0].real, thru_gamma.real, places=12)
        self.assertAlmostEqual(s11[0].imag, thru_gamma.imag, places=12)
        self.assertAlmostEqual(s21[0].real, 1.0, places=12)
        self.assertAlmostEqual(s21[0].imag, 0.0, places=12)


if __name__ == "__main__":
    unittest.main()
