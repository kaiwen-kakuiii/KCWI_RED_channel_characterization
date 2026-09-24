"""File-contract tests for the throughput steps. They do not run PypeIt."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import common

PYPEIT = """# Auto-generated PypeIt input file using PypeIt version: 2.0.1
# UTC 2026-08-24T03:43:15.030+00:00

# User-defined execution parameters
[rdx]
    spectrograph = keck_kcrm

# Setup
setup read
Setup A:
  binning: 1,1
  cenwave: 6679.9814453
  decker: Small
  dispname: RH1
setup end

# Data block
data read
 path /tmp/raw
                   filename |            frametype |        target |       exptime | ra_off | dec_off |
  KR.fe.fits |                  arc |    DOME FLATS |           1.0 |    0.0 |     0.0 |
  KR.th.fits |                 tilt |    DOME FLATS |           1.0 |    0.0 |     0.0 |
  KR.tr.fits |      illumflat,trace |    DOME FLATS |          10.0 |    0.0 |     0.0 |
  KR.pf.fits | pixelflat,scattlight |    DOME FLATS |          10.0 |    0.0 |     0.0 |
  KR.al.fits |                align |    DOME FLATS |          10.0 |    0.0 |     0.0 |
  KR.s5.fits |              science |       G191B2B |           5.0 |    0.0 |     0.0 |
  KR.s1.fits |              science |       G191B2B |         200.0 |    0.0 |     0.0 |
  KR.s2.fits |              science |       feige34 |         180.0 |    3.3 |     0.0 |
data end
"""


def frame(name, star, exptime, mjd, airmass, ra=0.0, dec=0.0, rate=None):
    row = {"filename": name, "star": star, "exptime": exptime, "mjd": mjd,
           "airmass": airmass, "ra": ra, "dec": dec, "spec2d": name}
    if rate is not None:
        row["rate"] = rate
    return row


class ConfigTests(unittest.TestCase):
    def test_every_grating_file_resolves(self):
        for name in ("RL", "RH1", "RH2", "RH3", "RH4", "RM1", "RM2"):
            cfg = common.load_config(name)
            self.assertEqual(cfg["grating"], name)
            self.assertEqual(cfg["boxcar_arcsec"], 3.4)
            self.assertEqual(cfg["exptime_floor_s"], 10)
            self.assertEqual(cfg["min_frames"], 2)
            self.assertEqual(cfg["max_airmass_spread"], 0.10)
            for entry in cfg["templates"].values():
                common.reid_value(name, entry)
            common.lamps_value(cfg["lamps"])

    def test_published_orders_and_the_rl_window(self):
        self.assertEqual(common.load_config("RH3")["polyorder"], 5)
        self.assertEqual(common.load_config("RH3")["bridge"], 5)
        self.assertEqual(common.load_config("RH3")["lamps"], "ThArRH3")
        self.assertEqual(common.load_config("RL")["slice_span_A"], [3400, 4500])
        self.assertEqual(common.load_config("RL")["shift_tol_A"], 400)
        self.assertEqual(common.load_config("RM2")["nline_min"], 8)
        self.assertEqual(common.load_config("RM2")["count_rate_floor"], 0.20)

    def test_rh1_templates_and_the_dead_column(self):
        cfg = common.load_config("RH1")
        self.assertEqual(common.template_entry(cfg, 6520, "Large"),
                         "templates/keck_kcrm_RH1_6520.fits")
        self.assertEqual(common.template_entry(cfg, 6680, "Small"),
                         "templates/keck_kcrm_RH1_small.fits")
        self.assertEqual(cfg["exclude_regions"]["2023-11-08/A"], "1:652:655,")
        with self.assertRaises(common.StepError):
            common.template_entry(cfg, 7000, "Large")

    def test_shift_tol_does_not_cross_slicers(self):
        rm1 = common.load_config("RM1")
        self.assertEqual(common.shift_tol(rm1, "Large"), 100)
        self.assertEqual(common.shift_tol(rm1, "Small"), 120)
        with self.assertRaises(common.StepError):
            common.shift_tol(rm1, "Medium")
        self.assertEqual(common.reid_value("RM1", "shipped"), "keck_kcrm_RM1.fits")

    def test_drop_and_flag_are_the_published_cubes(self):
        drops, flags = common.configured_lists()
        self.assertEqual(set(drops), {
            "g191b2b_2023-11-07_B",
            "feige34_2024-03-13_B",
            "feige110_2023-09-23_B",
        })
        self.assertEqual(set(flags), {"feige34_2024-03-15_B"})

    def test_exclude_without_a_comma_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "RH1.yaml")
            with open(path, "w") as fh:
                fh.write(
                    "grating: RH1\npolyorder: 15\nbridge: 0\nlamps: null\n"
                    "boxcar_arcsec: 3.4\nexptime_floor_s: 10\nmin_frames: 2\n"
                    "max_airmass_spread: 0.10\ncount_rate_floor: null\n"
                    "slice_span_A: [400, 900]\nshift_tol_A: 50\nnline_min: 20\n"
                    "templates:\n  \"6200 Large\": shipped\n"
                    "exclude_regions:\n  \"2023-11-08/A\": \"1:652:655\"\n")
            old = common.CONFIG_DIR
            common.CONFIG_DIR = tmp
            try:
                with self.assertRaises(common.StepError):
                    common.load_config("RH1")
            finally:
                common.CONFIG_DIR = old


class HeaderTests(unittest.TestCase):
    def test_wavelengths_stay_out_of_skysub_and_the_data_block_is_kept(self):
        reid = "/tmp/keck_kcrm_RH1_small.fits"
        rendered = common.render_header(
            PYPEIT, reid=reid, lamps="/tmp/ThArRH3",
            exclude="1:652:655,", length_range=0.3, user_regions=":35,64:")
        head, body = common.split_header(rendered)
        self.assertIn("[[wavelengths]]", head)
        self.assertLess(head.index("lamps"), head.index("[[skysub]]"))
        self.assertLess(head.index("exclude_regions"), head.index("[[wavelengths]]"))
        self.assertIn("user_regions = :35,64:", head)
        self.assertIn("KR.s1.fits |              science |       G191B2B", body)
        self.assertIn("# UTC 2026-08-24", rendered)
        values = common.header_values(rendered)
        self.assertEqual(values["reid_arxiv"], reid)
        self.assertEqual(values["exclude_regions"], "1:652:655,")
        self.assertEqual(values["length_range"], "0.3")

    def test_exclude_without_a_comma_is_not_written(self):
        with self.assertRaises(common.StepError):
            common.render_header(PYPEIT, exclude="1:652:655")

    def test_thar_retype_is_idempotent(self):
        lamps = {"KR.fe.fits": "FeAr", "KR.th.fits": "ThAr"}
        once, note = common.retype_thar(PYPEIT, lamps.get)
        self.assertIn("arc,tilt", once)
        self.assertIn("# FeAr, unused", once)
        self.assertNotIn("\n  KR.fe.fits |                  arc |", once)
        self.assertIn("1 ThAr", note)
        twice, again = common.retype_thar(once, lamps.get)
        self.assertEqual(again, "already retyped")
        self.assertEqual(twice, once)

    def test_activate_comments_only_the_other_star(self):
        text = common.activate(PYPEIT, {"KR.s1.fits"})
        rows = {row["filename"]: row for row in common.parse_rows(text)}
        self.assertFalse(rows["KR.s1.fits"]["commented"])
        self.assertTrue(rows["KR.s2.fits"]["commented"])
        self.assertFalse(rows["KR.al.fits"]["commented"])
        rows = {row["filename"]: row for row in common.parse_rows(common.activate(text, {"KR.s1.fits", "KR.s2.fits"}))}
        self.assertFalse(rows["KR.s2.fits"]["commented"])

    def test_cenwave_rounds_to_10(self):
        self.assertEqual(common.round_cenwave("6679.9814453"), 6680)
        self.assertEqual(common.round_cenwave("7149.8920898"), 7150)
        self.assertEqual(common.round_cenwave("8849.6"), 8850)
        self.assertEqual(common.round_cenwave("8900.4"), 8900)

    def test_science_floor_and_missing_calibrations(self):
        meta = common.setup_meta(PYPEIT)
        self.assertEqual(meta["letter"], "A")
        self.assertEqual(meta["decker"], "Small")
        cfg = {"exptime_floor_s": 10, "min_frames": 2}
        kept = common.science_frames(PYPEIT, cfg)
        self.assertEqual([row["filename"] for row in kept], ["KR.s1.fits", "KR.s2.fits"])
        self.assertEqual(kept[0]["star"], "g191b2b")


class GroupTests(unittest.TestCase):
    def test_exposure_floor_yields_for_the_only_two_frames(self):
        frames = [frame("a", "feige34", 7, 1, 1.1), frame("b", "feige34", 14, 2, 1.1)]
        kept, dropped = common.apply_exptime_floor(frames, 10, 2)
        self.assertEqual(dropped, [])
        self.assertEqual(len(kept), 2)

    def test_exposure_floor_drops_an_acquisition_when_two_remain(self):
        frames = [frame("a", "feige110", 5, 1, 1.1),
                  frame("b", "feige110", 30, 2, 1.1),
                  frame("c", "feige110", 30, 3, 1.1)]
        kept, dropped = common.apply_exptime_floor(frames, 10, 2)
        self.assertEqual([item["filename"] for item in dropped], ["a"])
        self.assertEqual(len(kept), 2)

    def test_airmass_and_pointing(self):
        # 1.2 arcsec nod stays one visit. 3.3 arcsec repoint splits.
        nod = [frame("a", "feige110", 60, 1, 1.20, ra=0.0),
               frame("b", "feige110", 60, 2, 1.21, ra=-0.7, dec=1.0)]
        groups, _, _ = common.group_frames(nod, 0.10, None, 2, 10)
        self.assertEqual(len(groups), 1)
        self.assertFalse(groups[0]["refused"])
        repoint = [frame("a", "g191b2b", 60, 1, 1.24, ra=0.0),
                   frame("b", "g191b2b", 60, 2, 1.60, ra=3.3)]
        groups, _, _ = common.group_frames(repoint, 0.10, None, 2, 10)
        self.assertEqual(len(groups), 2)

    def test_two_stars_and_a_one_frame_cube(self):
        frames = [frame("a", "feige110", 50, 1, 1.20),
                  frame("b", "feige110", 240, 2, 1.21),
                  frame("c", "feige34", 180, 3, 1.15)]
        groups, _, _ = common.group_frames(frames, 0.10, None, 2, 10)
        by_star = {group["star"]: group for group in groups}
        self.assertEqual(len(by_star["feige110"]["frames"]), 2)
        self.assertTrue(by_star["feige34"]["refused"])
        self.assertEqual(common.cube_tag("feige110", "2024-12-28", "B", 1, 1),
                         "feige110_2024-12-28_B")
        self.assertEqual(common.cube_tag("g191b2b", "2023-10-15", "B", 2, 2),
                         "g191b2b_2023-10-15_B_v2")

    def test_count_rate_floor_and_its_yield(self):
        bright = frame("a", "feige34", 60, 1, 1.1, rate=4.7e6)
        faint = frame("b", "feige34", 66, 2, 1.1, rate=2.1e4)
        more = [frame("c", "feige34", 60, 3, 1.1, rate=4.6e6),
                frame("d", "feige34", 60, 4, 1.1, rate=4.8e6)]
        groups, _, dropped = common.group_frames([bright, faint] + more, 0.10, 0.20, 2, 10)
        self.assertEqual([item["filename"] for item in dropped], ["b"])
        self.assertEqual(len(groups[0]["frames"]), 3)
        groups, _, dropped = common.group_frames([bright, faint], 0.10, 0.20, 2, 10)
        self.assertEqual(dropped, [])
        self.assertEqual(len(groups[0]["frames"]), 2)

    def test_coadd_file_lists_only_the_visit(self):
        frames = [frame("spec2d_a.fits", "g191b2b", 120, 1, 1.60),
                  frame("spec2d_b.fits", "g191b2b", 120, 2, 1.62)]
        text = common.render_coadd3d(
            "g191b2b_2023-10-15_B_v1", "g191b2b", "2023-10-15", "B", 1, 2, frames)
        self.assertIn("combine = True", text)
        self.assertIn("Science/spec2d_a.fits", text)
        self.assertNotIn("spec2d_c", text)

    def test_sky_windows_cluster_only_when_they_agree(self):
        separate = common.cluster_regions([
            (("feige", 0.0, 0.0), ":38,65:", ["a"]),
            (("feige", -0.7, 1.0), ":31,59:", ["b"]),
        ])
        self.assertEqual(len(separate), 2)
        merged = common.cluster_regions([
            (("a", 0, 0), ":35,64:", ["a"]),
            (("b", 0, 0), ":36,65:", ["b"]),
        ])
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0]["frames"], ["a", "b"])

    def test_sens_file_records_order_and_refuses_a_zero_red_trim(self):
        text = common.render_sens(120, 126, 15)
        self.assertIn("algorithm = IR", text)
        self.assertIn("extr = BOX", text)
        self.assertIn("polyorder = 15", text)
        self.assertIn("trim_std_pixs = 120, 126", text)
        sens, fits = common.sens_names("feige110_2024-12-04_B", 5)
        self.assertEqual(fits, "feige110_2024-12-04_B_sens_IR_p5cov.fits")
        with self.assertRaises(common.StepError):
            common.render_sens(10, 0, 15)


if __name__ == "__main__":
    unittest.main()
