# unun-vna 0.5.1

`unun-vna` is a command-line tool for measuring the **small-signal efficiency**
of HF UnUn transformers with an S-A-A-2 / NanoVNA V2 using the
**one-transformer method**.

The tool is generalized:

- no nominal load/series resistor is assumed;
- the actual series-load fixture is measured as a complex impedance versus
  frequency;
- transformer efficiency is calculated from measured S11/S21 and the measured
  fixture;
- transformer impedance ratio is configurable;
- 1:49 remains the default only for convenience.

The lower-impedance side of the transformer is connected to Port 1. The
higher-impedance side is connected to the series fixture and then Port 2.

---

## 1. Transformer ratio

`--ratio` means the **impedance** ratio

\[
N = \frac{Z_\text{high}}{Z_\text{low}}.
\]

These are equivalent:

```bash
--ratio 49
--ratio 1:49
--ratio 1/49
```

All mean `N = 49`.

Examples:

```text
1:4   -> N=4,  turns ratio 1:2
1:9   -> N=9,  turns ratio 1:3
1:16  -> N=16, turns ratio 1:4
1:49  -> N=49, turns ratio 1:7
1:64  -> N=64, turns ratio 1:8
```

The default is `N=49`.

Ratios below 1 are rejected because this measurement topology expects the
low-impedance side at Port 1. Reverse the transformer connection instead.

---

## 2. Generic measurement topology

```text
Port 1 ---- UnUn ---- measured series fixture ---- Port 2
  Z0        1:N                                  Z0
```

`Z0` is normally 50 Ω.

The fixture may have any resistance value. It is not inferred from the
transformer ratio and it is not assumed to equal its nominal/DC value.

It is measured as

\[
Z_s(f)=R_s(f)+jX_s(f).
\]

---

## 3. Ideal matching load for a selected ratio

For an ideal transformer with ratio `N`, the target total high-side resistance
for a perfect resistive match is

\[
R_{\text{high,target}} = N Z_0.
\]

Port 2 itself contributes `Z0`, so the ideal purely resistive series fixture
would be

\[
R_{\text{series,target}} = (N-1)Z_0.
\]

At `Z0=50 Ω`:

| UnUn | N | Total high-side target | Ideal series fixture |
|---|---:|---:|---:|
| 1:4 | 4 | 200 Ω | 150 Ω |
| 1:9 | 9 | 450 Ω | 400 Ω |
| 1:16 | 16 | 800 Ω | 750 Ω |
| 1:25 | 25 | 1250 Ω | 1200 Ω |
| 1:49 | 49 | 2450 Ω | 2400 Ω |
| 1:64 | 64 | 3200 Ω | 3150 Ω |

These numbers are **diagnostics only**. They are never substituted for the
measured fixture impedance.

---

## 4. Fixture characterization

For a pure series impedance between two equal `Z0` ports:

\[
S_{11}=\frac{Z_s}{2Z_0+Z_s}
\]

\[
S_{21}=\frac{2Z_0}{2Z_0+Z_s}.
\]

Therefore:

\[
Z_{s,1}=2Z_0\frac{S_{11}}{S_{21}}
\]

and

\[
Z_{s,2}=2Z_0\left(\frac{1}{S_{21}}-1\right).
\]

A pure series element also satisfies

\[
S_{11}+S_{21}=1.
\]

The `fixture` command reports both impedance estimates, their disagreement, and
`|S11+S21-1|`.

This checks whether the physical fixture can safely be represented by one
complex series impedance.

---

## 5. Efficiency does not depend on the configured ratio

By S-parameter definition,

\[
|S_{21}|^2
\]

is the incident Port-1 power fraction delivered to the matched Port-2
termination `Z0`.

The same current flows through the Port-2 resistance and the series fixture.
Thus total real high-side load power is

\[
\frac{P_\text{high}}{P_\text{incident}}
=
|S_{21}|^2
\left(1+\frac{\Re\{Z_s\}}{Z_0}\right).
\]

Power accepted at Port 1 is

\[
\frac{P_\text{accepted}}{P_\text{incident}}
=
1-|S_{11}|^2.
\]

Therefore

\[
\boxed{
\eta=
\frac{
|S_{21}|^2
\left(1+\frac{\Re\{Z_s\}}{Z_0}\right)
}{
1-|S_{11}|^2
}
}
\]

and

\[
L_\text{dB}=-10\log_{10}(\eta).
\]

### VSWR

The analysis also reports Port-1 VSWR derived directly from calibrated S11.
With

\[
\rho = |S_{11}|,
\]

\[
\boxed{
\mathrm{VSWR}=\frac{1+\rho}{1-\rho}
}
\]

for \(\rho<1\).

A perfect match gives VSWR = 1.0. Total reflection gives infinite VSWR.
Values with \(|S_{11}|\ge 1\), which can arise from residual calibration or
measurement error, are reported as infinite rather than as a non-physical
negative VSWR.

The tool also reports the **ideal lossless VSWR** for the configured
transformer ratio with the actually measured fixture. This helps separate
load mismatch from transformer loss.

There is no `N` in the measured efficiency equation. That is intentional.

The measured power flow already contains the transformer's actual behaviour;
the ratio is not needed to reconstruct efficiency.

---

## 6. What `--ratio` is used for

The configured ratio is used to calculate a **lossless ideal reference** for
the actual measured fixture.

The real measured high-side load is

\[
Z_\text{load}=Z_s+Z_0.
\]

An ideal lossless transformer with impedance ratio `N` reflects that to Port 1
as

\[
Z_\text{in}=\frac{Z_\text{load}}{N}.
\]

The input reflection expected from load mismatch alone is

\[
\Gamma_\text{ideal}
=
\frac{Z_\text{in}-Z_0}
     {Z_\text{in}+Z_0}.
\]

Therefore

\[
|S_{11,\text{ideal}}|=|\Gamma_\text{ideal}|.
\]

The accepted incident power of that ideal network is

\[
1-|\Gamma_\text{ideal}|^2.
\]

Since an ideal transformer is lossless, the fraction of that real power
dissipated in Port 2 is

\[
\frac{Z_0}{Z_0+\Re\{Z_s\}}.
\]

Hence

\[
|S_{21,\text{ideal}}|^2
=
\left(1-|\Gamma_\text{ideal}|^2\right)
\frac{Z_0}{Z_0+\Re\{Z_s\}}.
\]

The analysis output shows measured S11/S21 beside these ideal values.

This is useful when the fixture is not exactly the perfect matching resistance
for the selected ratio.

---

## 7. S-A-A-2 Port-2 load-match limitation

The S-A-A-2 is a T/R VNA, not a full bidirectional two-port VNA.

It measures S11 and forward S21, but cannot perform a complete independent
Port-2 load-match correction.

High-ratio transformer/fixture combinations can be strongly reflecting when
viewed from Port 2 even when Port 1 is well matched.

A THRU verification can therefore look excellent while a highly reflective DUT
still has biased S21.

### Recommended Port-2 pad

Use a fixed 50 Ω attenuator, typically 10 dB, immediately before Port 2:

```text
Port 2 ---- 10 dB pad ---- calibration plane ---- DUT
```

The same pad must remain in place during:

1. calibration;
2. fixture characterization;
3. transformer measurement.

The THRU calibration removes its forward loss and phase. Reflections from Port
2 pass through it twice, so a 10 dB pad reduces the load-match influence by
roughly 20 dB.

---

## 8. Calibration

```bash
unun-vna calibrate \
    --start 1M \
    --stop 35M \
    --points 501 \
    -o hf35-pad
```

The program requests:

1. SHORT on Port 1;
2. OPEN on Port 1;
3. 50 Ω LOAD on Port 1;
4. THRU between the calibration planes.

Files:

```text
hf35-pad.npz
hf35-pad.npz.json
hf35-pad.npz.through.csv
```

---

## 9. Characterize a fixture

The fixture extraction itself is independent of ratio.

For a 1:9 reference:

```bash
unun-vna fixture \
    --cal hf35-pad.npz \
    --ratio 1:9 \
    -o fixture-1to9
```

For 1:64:

```bash
unun-vna fixture \
    --cal hf35-pad.npz \
    --ratio 64 \
    -o fixture-1to64
```

`--ratio` changes only the target-load diagnostics printed by this command.

---

## 10. Measure a transformer

Measurement acquisition itself requires no ratio:

```bash
unun-vna measure \
    --cal hf35-pad.npz \
    -o transformer
```

The measured chain is:

```text
Port 1 ---- transformer ---- characterized fixture ---- Port 2 pad ---- Port 2
```

---

## 11. Analyze

Default 1:49:

```bash
unun-vna analyze \
    --dut transformer.npz \
    --fixture fixture.npz \
    -o efficiency.csv
```

Explicit 1:49:

```bash
unun-vna analyze \
    --dut transformer.npz \
    --fixture fixture.npz \
    --ratio 1:49 \
    -o efficiency.csv
```

1:9:

```bash
unun-vna analyze \
    --dut transformer.npz \
    --fixture fixture-1to9.npz \
    --ratio 1:9 \
    -o efficiency-1to9.csv
```

1:64:

```bash
unun-vna analyze \
    --dut transformer.npz \
    --fixture fixture-1to64.npz \
    --ratio 64 \
    -o efficiency-1to64.csv
```

---

## 12. Analysis CSV

The analysis CSV contains:

```text
frequency_hz
s11_db
vswr
s21_db
fixture_r_ohm
fixture_x_ohm
impedance_ratio_high_over_low
target_high_side_resistance_ohm
target_series_resistance_ohm
ideal_lossless_s11_db_for_measured_fixture
ideal_lossless_vswr_for_measured_fixture
ideal_lossless_s21_db_for_measured_fixture
accepted_fraction
delivered_fraction
efficiency
efficiency_percent
transformer_loss_db
```

The `target_*` and `ideal_lossless_*` columns depend on `--ratio`.

The actual efficiency columns do not.

---

## 13. Ideal sanity checks

### 1:9

At `Z0=50 Ω`, use an ideal 400 Ω series fixture:

\[
400+50=450\ \Omega
\]

and

\[
450/9=50\ \Omega.
\]

Thus an ideal transformer has

\[
S_{11}=0
\]

and

\[
|S_{21}|=1/3
\]

or approximately

\[
S_{21}=-9.54\ \text{dB}.
\]

Efficiency evaluates to 100%.

### 1:49

With 2400 Ω series:

\[
|S_{21}|=1/7
\]

or

\[
S_{21}\approx-16.90\ \text{dB}.
\]

### 1:64

With 3150 Ω series:

\[
|S_{21}|=1/8
\]

or

\[
S_{21}\approx-18.06\ \text{dB}.
\]

---

## 14. Remaining assumptions

The tool does **not** assume:

- a 1:49 transformer;
- a 2400 Ω fixture;
- any nominal fixture resistance;
- that the fixture is purely resistive.

It does assume:

1. the lower-impedance transformer side is connected to Port 1;
2. the fixture is adequately represented by one complex series impedance;
3. Port 2 is a matched real `Z0` termination at the calibrated reference plane;
4. calibration, fixture measurement and DUT measurement use the same RF setup;
5. a fixed Port-2 pad, if used, remains in place;
6. the measurement is small-signal.

A fully arbitrary two-port fixture would require S12 and S22 and therefore a
full bidirectional VNA.

---

## 15. Full workflow example

```bash
pipx install --editable . --force

unun-vna info

unun-vna calibrate \
    --start 1M \
    --stop 35M \
    --points 501 \
    -o hf35-pad

unun-vna fixture \
    --cal hf35-pad.npz \
    --ratio 1:49 \
    -o fixture

unun-vna measure \
    --cal hf35-pad.npz \
    -o transformer

unun-vna analyze \
    --dut transformer.npz \
    --fixture fixture.npz \
    --ratio 1:49 \
    -o efficiency.csv
```

---

## 16. Source structure

```text
unun-vna/
├── pyproject.toml
├── README.md
├── src/
│   └── unun_vna/
│       ├── __init__.py
│       ├── cli.py
│       └── analysis.py
└── tests/
    └── test_math.py
```

`analysis.py` deliberately separates:

- ratio-independent measured efficiency;
- ratio-dependent ideal-reference diagnostics.

This prevents a nominal transformer ratio from silently changing the measured
power calculation.
