"""Tests for the analysis-CSV visualization layer."""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from unun_vna.plotting import plot_analysis_files, read_analysis_csv


class PlottingTests(unittest.TestCase):
    def _write_csv(self, path: Path, offset: float = 0.0) -> None:
        with path.open("w", newline="", encoding="utf-8") as fp:
            w = csv.writer(fp)
            w.writerow(["frequency_hz", "vswr", "efficiency_percent"])
            w.writerow([7_000_000, 1.25 + offset, 87.0 - offset])
            w.writerow([14_000_000, 1.45 + offset, 84.0 - offset])
            w.writerow([21_000_000, 1.55 + offset, 80.0 - offset])
            w.writerow([28_500_000, 1.70 + offset, 72.0 - offset])

    def test_read_analysis_csv(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "unun1-efficiency.csv"
            self._write_csv(path)
            s = read_analysis_csv(path)
            self.assertEqual(s.label, "unun1")
            self.assertEqual(len(s.frequency_mhz), 4)
            self.assertAlmostEqual(s.frequency_mhz[0], 7.0)
            self.assertAlmostEqual(s.efficiency_percent[-1], 72.0)

    def test_single_png_plot(self):
        with tempfile.TemporaryDirectory() as td:
            csv_path = Path(td) / "unun1.csv"
            png_path = Path(td) / "plot.png"
            self._write_csv(csv_path)
            result = plot_analysis_files([csv_path], png_path)
            self.assertEqual(result, png_path)
            self.assertTrue(png_path.is_file())
            self.assertGreater(png_path.stat().st_size, 1000)

    def test_multi_svg_plot(self):
        with tempfile.TemporaryDirectory() as td:
            a = Path(td) / "a.csv"
            b = Path(td) / "b.csv"
            out = Path(td) / "compare.svg"
            self._write_csv(a)
            self._write_csv(b, offset=0.1)
            plot_analysis_files(
                [a, b],
                out,
                labels=["Transformer A", "Transformer B"],
                title="Comparison",
            )
            self.assertTrue(out.is_file())
            text = out.read_text(encoding="utf-8")
            self.assertIn("Comparison", text)
            self.assertIn("Transformer A", text)
            self.assertIn("Transformer B", text)

    def test_label_count_must_match(self):
        with tempfile.TemporaryDirectory() as td:
            a = Path(td) / "a.csv"
            b = Path(td) / "b.csv"
            self._write_csv(a)
            self._write_csv(b)
            with self.assertRaises(ValueError):
                plot_analysis_files([a, b], Path(td) / "x.png", labels=["A"])


if __name__ == "__main__":
    unittest.main()
