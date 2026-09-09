"""Numerical models used by :mod:`unun_vna.cli`.

This module intentionally contains no VNA I/O.  It translates calibrated
complex S-parameters into the two quantities that are specific to the
one-transformer UnUn measurement method:

* the complex impedance of the series-load fixture, and
* the power efficiency of the transformer itself.

The efficiency calculation intentionally does **not** require a transformer
ratio or a nominal load resistor. It uses the load fixture that was actually
measured. A transformer impedance ratio is used only by
``ideal_transformer_reference`` for matching/reference diagnostics.

The separation is deliberate: the measurement code can be changed or replaced
without changing the RF model, and the model can be unit-tested with synthetic
S-parameters.

Important limitation
--------------------
The S-A-A-2 is a T/R VNA, not a full bidirectional 2-port VNA.  We therefore do
*not* attempt a general 2-port/ABCD de-embedding of the resistor fixture.  Such
an operation would require the complete set S11, S21, S12 and S22.  Instead we
use the physical fact that the fixture was deliberately constructed as a
single two-terminal *series impedance* between the two VNA centre conductors.
The code verifies how well the measurement obeys that model before using it.
"""

from __future__ import annotations

import numpy as np


def vswr_from_s11(s11):
    """Convert complex S11 values to voltage standing-wave ratio (VSWR).

    For reflection-coefficient magnitude ``rho = |S11|``,

        VSWR = (1 + rho) / (1 - rho).

    A perfect match has VSWR = 1.0. Total reflection (rho = 1) gives infinite
    VSWR. Values with rho > 1 can occur because of residual calibration/noise
    error; these are also returned as infinity instead of a non-physical
    negative VSWR.
    """
    rho = np.abs(np.asarray(s11, dtype=complex))
    out = np.full_like(rho, np.inf, dtype=float)
    valid = rho < 1.0
    out[valid] = (1.0 + rho[valid]) / (1.0 - rho[valid])
    return out



def series_fixture_impedance(s11, s21, z0: float = 50.0):
    """Extract the complex series impedance of the resistor fixture.

    Parameters
    ----------
    s11, s21:
        Calibrated complex S-parameter arrays of the fixture.  The arrays must
        refer to the same frequency grid.
    z0:
        Real reference impedance of both VNA ports.  The S-A-A-2 normally uses
        50 ohms.

    RF model
    --------
    The fixture is assumed to be a *single series impedance* ``Zs`` inserted
    between two equal ``Z0`` ports::

        Port 1 ---- Zs ---- Port 2
          Z0                  Z0

    For this special network, elementary two-port analysis gives

        S11 = Zs / (2*Z0 + Zs)
        S21 = 2*Z0 / (2*Z0 + Zs)

    Therefore ``Zs`` can be recovered in two independent ways:

        Zs = 2*Z0 * S11/S21

    and

        Zs = 2*Z0 * (1/S21 - 1)

    With an ideal series-only fixture and an ideal calibration, both estimates
    are identical.  Their disagreement is useful because it exposes effects
    that cannot be represented by a single series impedance, for example a
    relevant shunt capacitance to a connector shell or residual VNA
    load-match/calibration error.

    A second identity follows immediately from the two equations above:

        S11 + S21 = 1

    ``identity_error`` stores the complex departure from this identity at every
    frequency.  It is a convenient, dimensionless fixture-quality diagnostic.

    Why average the two impedance estimates?
    ----------------------------------------
    Neither estimator is universally superior in the presence of measurement
    noise.  ``2*Z0*S11/S21`` uses both measured channels, while the expression
    based on S21 alone is especially sensitive when S21 is very small.  Taking
    their mean is a simple symmetric choice; importantly, the code also stores
    both estimates separately so their disagreement remains visible.

    Returns
    -------
    dict of numpy arrays
        ``z_series``
            Mean complex series impedance used by the subsequent efficiency
            calculation.
        ``z_from_ratio`` / ``z_from_s21``
            The two independent impedance estimates described above.
        ``identity_error``
            ``S11 + S21 - 1``.
        ``model_error``
            RMS complex S-parameter error after reconstructing S11 and S21 from
            the extracted series impedance.
        ``pred_s11`` / ``pred_s21``
            S-parameters predicted by the extracted pure-series model.
    """
    # Convert input sequences to vectorised complex numpy arrays.  Keeping the
    # full complex values is essential: phase contains the fixture reactance.
    s11 = np.asarray(s11, dtype=complex)
    s21 = np.asarray(s21, dtype=complex)

    # Both extraction formulae divide by S21. A vanishing transmission would
    # make the inverse problem numerically singular. No nominal fixture
    # resistance is assumed here.
    if np.any(np.abs(s21) < 1e-15):
        raise ValueError("S21 too close to zero for series-impedance extraction")

    # Estimator 1: divide the reflected wave by the transmitted wave.  For a
    # pure series element the common denominator (2*Z0 + Zs) cancels exactly.
    z_ratio = 2.0 * z0 * s11 / s21

    # Estimator 2: invert the closed-form S21 expression directly.
    z_s21 = 2.0 * z0 * (1.0 / s21 - 1.0)

    # Use the average as the nominal fixture impedance, but keep both originals
    # in the return value so that their disagreement can be inspected.
    z_series = 0.5 * (z_ratio + z_s21)

    # For a pure series impedance between equal reference impedances:
    #     S11 + S21 = (Zs + 2*Z0)/(2*Z0 + Zs) = 1.
    # A non-zero result points to parasitics or residual measurement errors.
    identity_error = s11 + s21 - 1.0

    # Reconstruct the S-parameters from the extracted series model.  Comparing
    # these to the measurements gives a second, directly S-parameter-domain
    # measure of how well the one-element model describes the fixture.
    pred_s21 = 2.0 * z0 / (2.0 * z0 + z_series)
    pred_s11 = z_series / (2.0 * z0 + z_series)
    model_error = np.sqrt(
        (np.abs(s11 - pred_s11) ** 2 + np.abs(s21 - pred_s21) ** 2) / 2.0
    )

    return {
        "z_series": z_series,
        "z_from_ratio": z_ratio,
        "z_from_s21": z_s21,
        "identity_error": identity_error,
        "model_error": model_error,
        "pred_s11": pred_s11,
        "pred_s21": pred_s21,
    }


def unun_efficiency(s11, s21, z_series, z0: float = 50.0):
    """Calculate UnUn efficiency for the one-transformer measurement method.

    Measurement topology
    --------------------
    The calibrated VNA sees the complete chain::

        Port 1 ---- UnUn ---- Z_series ---- Port 2
          Z0                                 Z0

    ``Z_series`` is the independently measured series-load fixture. Port 2
    itself contributes the reference resistance ``Z0``. The total complex load
    presented to the transformer's high-impedance side is therefore

        Z_high_load = Z_series + Z0.

    No nominal resistor value and no transformer ratio are assumed by the
    efficiency calculation.

    What S21 measures
    -----------------
    By definition ``|S21|^2`` is the fraction of *incident Port-1 power* that
    reaches the matched 50-ohm receiver load at Port 2.  It does **not** include
    the real power dissipated in the series resistor.

    The resistor and the 50-ohm Port-2 load are in series, hence the same RF
    current flows through both.  Real power is proportional to resistance, so
    the total real power delivered by the transformer to the combined
    high-side load is

        P_high / P_incident
            = |S21|^2 * (Z0 + Re{Z_series}) / Z0
            = |S21|^2 * (1 + Re{Z_series}/Z0).

    Only the *real* part of ``Z_series`` appears here.  Its imaginary part can
    alter S11/S21 and therefore the operating point, but an ideal reactance does
    not dissipate average power.

    What S11 corrects
    -----------------
    Some incident power can be reflected at Port 1.  Transformer efficiency
    should be referred to power actually accepted by the DUT, not to the
    generator's incident power.  For a calibrated 50-ohm port,

        P_accepted / P_incident = 1 - |S11|^2.

    The resulting efficiency estimate is therefore

        eta = |S21|^2 * (1 + Re{Z_series}/Z0)
              -----------------------------------
                         1 - |S11|^2

    The multiplicative load factor is determined entirely by the *measured*
    real part of the series fixture. For example, if Re(Z_series)=2400 ohms
    and Z0=50 ohms, the factor happens to be 49. Other fixture values work
    without changing the equation.

    Assumptions and limitations
    ---------------------------
    * The load fixture must be adequately described by a *series* impedance.
      ``series_fixture_impedance`` provides diagnostics for this.
    * S11 and S21 must be calibrated on the same reference planes used for the
      fixture measurement.  In practice the 10 dB Port-2 pad, if used to tame
      receiver load-match error, must remain present for calibration, fixture
      measurement and DUT measurement.
    * This is a small-signal VNA efficiency measurement.  It does not measure
      power-dependent core heating/non-linearity that may appear at 50-100 W.
    * Residual T/R-VNA load-match errors can still bias S21.  The external pad
      is intended to reduce that error; efficiencies above 100% are therefore
      treated by the CLI as a useful warning sign rather than physical results.

    Returns
    -------
    dict of numpy arrays
        ``accepted_fraction``
            Fraction of incident Port-1 power accepted by the DUT.
        ``delivered_fraction``
            Fraction of incident Port-1 power reconstructed at the complete
            high-side resistive load (series resistor + Port 2).
        ``load_factor``
            ``1 + Re(Z_series)/Z0``.
        ``efficiency``
            ``delivered_fraction / accepted_fraction``.
        ``loss_db``
            Transformer insertion loss referred to accepted power,
            ``-10*log10(efficiency)``.
    """
    s11 = np.asarray(s11, dtype=complex)
    s21 = np.asarray(s21, dtype=complex)
    z_series = np.asarray(z_series, dtype=complex)

    # Fraction of the incident forward power that actually enters the UnUn.
    # |S11|^2 is reflected *power*, hence no 20*log convention here.
    accepted = 1.0 - np.abs(s11) ** 2

    # Port 2 measures only the power dissipated in its own Z0 termination.  The
    # same current also dissipates power in Re(Z_series); scale by the resistance
    # ratio to reconstruct the total real load power on the transformer's high
    # side.
    load_factor = 1.0 + np.real(z_series) / z0
    delivered = np.abs(s21) ** 2 * load_factor

    # Avoid dividing by zero or a negative accepted-power value.  The latter
    # would be non-physical and normally indicates a corrupt calibration/data
    # set.  NaN makes such points visible in subsequent analysis rather than
    # silently producing nonsense.
    eta = np.full_like(accepted, np.nan, dtype=float)
    good = accepted > 0.0
    eta[good] = delivered[good] / accepted[good]

    # Express the same result as dissipative loss.  eta>1 yields a negative loss
    # value; the CLI separately warns when eta is implausibly above unity.
    loss_db = np.full_like(eta, np.nan, dtype=float)
    good_eta = eta > 0.0
    loss_db[good_eta] = -10.0 * np.log10(eta[good_eta])

    return {
        "accepted_fraction": accepted,
        "delivered_fraction": delivered,
        "load_factor": load_factor,
        "efficiency": eta,
        "loss_db": loss_db,
    }


def ideal_transformer_reference(
    z_series,
    impedance_ratio: float = 49.0,
    z0: float = 50.0,
):
    """Calculate lossless-transformer reference values for the measured fixture.

    This function is for diagnostics only. It does not participate in the
    measured efficiency calculation.

    ``impedance_ratio`` is the high-side/low-side impedance ratio

        N = Z_high / Z_low.

    A transformer commonly described as 1:49 therefore uses ``N=49`` and has a
    7:1 turns ratio.

    The measured high-side load is

        Z_load = Z_series + Z0.

    A lossless ideal transformer reflects that load to Port 1 as

        Z_in = Z_load / N.

    The input reflection caused only by the selected ratio and the actual
    measured fixture is

        Gamma_in = (Z_in - Z0) / (Z_in + Z0).

    Since the transformer is lossless, all accepted real power is dissipated in
    Re(Z_series) and the matched Port-2 resistance Z0. The fraction reaching
    Port 2 is therefore

        Z0 / (Z0 + Re(Z_series)).

    Thus

        |S21_ideal|^2 =
            (1 - |Gamma_in|^2) * Z0/(Z0 + Re(Z_series)).

    The nominal purely resistive series fixture required for a perfect match is

        R_series,target = (N - 1) * Z0.

    That target is reported only as a reference. It is never substituted for
    the measured fixture impedance.
    """
    z_series = np.asarray(z_series, dtype=complex)

    if not np.isfinite(impedance_ratio) or impedance_ratio <= 0.0:
        raise ValueError("impedance_ratio must be finite and > 0")
    if not np.isfinite(z0) or z0 <= 0.0:
        raise ValueError("z0 must be finite and > 0")

    n = float(impedance_ratio)

    z_load = z_series + z0
    z_in = z_load / n
    gamma_in = (z_in - z0) / (z_in + z0)

    accepted = 1.0 - np.abs(gamma_in) ** 2
    real_total_load = z0 + np.real(z_series)

    s21_power = np.full_like(accepted, np.nan, dtype=float)
    valid = (accepted >= 0.0) & (real_total_load > 0.0)
    s21_power[valid] = accepted[valid] * z0 / real_total_load[valid]

    s21_magnitude = np.full_like(s21_power, np.nan, dtype=float)
    valid_power = s21_power >= 0.0
    s21_magnitude[valid_power] = np.sqrt(s21_power[valid_power])

    return {
        "z_load": z_load,
        "z_in": z_in,
        "gamma_in": gamma_in,
        "accepted_fraction": accepted,
        "s21_power": s21_power,
        "s21_magnitude": s21_magnitude,
        "target_high_side_resistance": n * z0,
        "target_series_resistance": (n - 1.0) * z0,
        "turns_ratio": np.sqrt(n),
    }
