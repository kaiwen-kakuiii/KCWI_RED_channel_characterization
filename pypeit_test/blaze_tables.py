#!/usr/bin/env python
"""Render GRATING_SUMMARY.md: hand-written prose, machine-written tables.

The prose lives in `grating_summary.template.md` beside this file, with a
`{{NAME}}` placeholder wherever a table goes.  Every table comes from the CSVs
`plot_blaze_summary.py` writes, so a re-reduction that changes a number changes
the document and the figures together rather than leaving the document quietly
stale.  **Edit the template, never GRATING_SUMMARY.md** -- the latter is output
and a hand edit to it is lost on the next render.

Tables available: COMPOSITE (composite vs envelope), OVERALL (the seven-grating
overview), and per grating <G>_ROLL (blaze rollup) and <G>_CUBES (one row per
cube).

Usage:
    python blaze_tables.py --render          # rewrite ../GRATING_SUMMARY.md
    python blaze_tables.py --check           # verify it is in sync, exit 1 if not
    python blaze_tables.py [--section RM2_ROLL]
"""
import argparse
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ORDER = ["RL", "RM1", "RM2", "RH1", "RH2", "RH3", "RH4"]


def build():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "blaze_summary.csv"))))
    peaks = {(r["grating"], int(r["cenwave"]), r["slicer"]): r
             for r in csv.DictReader(open(os.path.join(ROOT,
                                                       "blaze_peaks.csv")))}
    comp = {r["grating"]: r for r in csv.DictReader(
        open(os.path.join(ROOT, "composite_summary.csv")))}
    trim = {r["grating"]: r for r in csv.DictReader(
        open(os.path.join(ROOT, "trim_summary.csv")))}
    out = {}

    # Per config: the three cuts, in Angstroms and as a fraction of what the
    # reduction delivered.  Reported per config because that is the unit the
    # cuts act on; the union table below is the same accounting per grating.
    L = ["| grating | delivered | after `measure_trim` | fit range | cut by "
         "`measure_trim` | cut by the fit | kept |",
         "|---|--:|--:|--:|--:|--:|--:|"]
    for g in ORDER:
        t = trim[g]
        d = float(t["delivered_A"])
        pc = lambda k: f"{float(t[k]):.0f} ({100 * float(t[k]) / d:.0f}%)"
        L.append(f"| **{g}** | {t['delivered_A']} A | {t['after_trim_std_A']} A "
                 f"| {t['fit_range_A']} A | −{pc('cut_trim_std_A')} | "
                 f"−{pc('cut_fit_A')} | **{t['kept_pct']}%** |")
    out["TRIM"] = "\n".join(L)

    # Per grating, stacking every central wavelength: what a config gives up at
    # its own ends, the next config's centre usually covers.
    L = ["| grating | delivered by all configs | span | usable after both cuts "
         "| span | kept | gaps left |",
         "|---|---|--:|---|--:|--:|---|"]
    for g in ORDER:
        t = trim[g]
        L.append(f"| **{g}** | {t['union_lo_A']}–{t['union_hi_A']} A | "
                 f"{t['union_span_A']} A | "
                 f"{t['union_fit_lo_A']}–{t['union_fit_hi_A']} A | "
                 f"{t['union_fit_span_A']} A | {t['union_kept_pct']}% | "
                 f"{t['interior_gaps_A'] or '—'} |")
    out["TRIM_UNION"] = "\n".join(L)

    edge = {r["grating"]: r for r in csv.DictReader(
        open(os.path.join(ROOT, "edge_check.csv")))}
    L = ["| grating | fit span | edge RMS | edge max | interior RMS | interior "
         "max | RMS ratio | max ratio |",
         "|---|--:|--:|--:|--:|--:|--:|--:|"]
    for g in ORDER:
        e = edge[g]
        if not e["interior_rms_mag"]:
            L.append(f"| **{g}** | {e['mean_fit_span_A']} A | "
                     f"{e['edge_rms_mag']} | {e['edge_max_mag']} | "
                     f"*no interior zone* | — | — | — |")
            continue
        L.append(f"| **{g}** | {e['mean_fit_span_A']} A | {e['edge_rms_mag']} | "
                 f"{e['edge_max_mag']} | {e['interior_rms_mag']} | "
                 f"{e['interior_max_mag']} | **{e['rms_ratio']}** | "
                 f"**{e['max_ratio']}** |")
    out["TRIM_EDGE"] = "\n".join(L)

    # Composite vs envelope, side by side.  The gap between the two is the whole
    # point of having both figures: it is how much the central wavelengths that
    # happened to be scheduled cost you against tuning the grating to the
    # wavelength you care about.
    L = ["| grating | composite max | best blaze | gap | blazes | n over "
         "range | Small-slicer drag (mean / worst) |",
         "|---|--:|--:|--:|--:|---|---|"]
    for g in ORDER:
        c = comp[g]
        best = max((float(p["combined_peak_pct"]) for p in peaks.values()
                    if p["grating"] == g))
        gap = float(c["composite_max_pct"]) - best
        drag = ("—" if not c["small_drag_mean_pt"]
                else f"{c['small_drag_mean_pt']} / {c['small_drag_worst_pt']} pt")
        L.append(f"| **{g}** | {c['composite_max_pct']}% @ {c['at_wave_A']} A | "
                 f"{best:.1f}% | {gap:+.1f} pt | {c['nblazes']} | "
                 f"{c['n_min']}–{c['n_max']} | {drag} |")
    out["COMPOSITE"] = "\n".join(L)

    L = ["| grating | nights | blazes | cubes | frames | central wavelengths "
         "(A) | slicers | standards | best blaze |",
         "|---|--:|--:|--:|--:|---|---|---|---|"]
    for g in ORDER:
        # Excluded cubes are counted nowhere: a night whose photons did not
        # arrive is not a night this grating was measured on.
        ok = [r for r in rows if r["grating"] == g and r["status"] != "excluded"]
        bl = sorted({(int(r["cenwave"]), r["slicer"]) for r in ok})
        cws = sorted({int(r["cenwave"]) for r in ok})
        best = max((peaks[(g, cw, s)] for cw, s in bl),
                   key=lambda p: float(p["combined_peak_pct"]))
        L.append(f"| **{g}** | {len({r['night'] for r in ok})} | {len(bl)} | "
                 f"{len(ok)} | {sum(int(r['nframes']) for r in ok)} | "
                 f"{', '.join(str(c) for c in cws)} | "
                 f"{', '.join(sorted({r['slicer'] for r in ok}))} | "
                 f"{', '.join(sorted({r['star'] for r in ok}))} | "
                 f"**{best['combined_peak_pct']}%** @ {best['peak_wave_A']} A "
                 f"({best['cenwave']} {best['slicer']}) |")
    out["OVERALL"] = "\n".join(L)

    for g in ORDER:
        rr = [r for r in rows if r["grating"] == g]
        L = ["| blaze (cenwave) | slicer | nights | cubes | frames | common "
             "fit range (A) | combined peak | at | top within 0.5 pt |",
             "|--:|---|--:|--:|--:|---|--:|--:|---|"]
        for cw, sl in sorted({(int(r["cenwave"]), r["slicer"]) for r in rr}):
            ok = [r for r in rr if int(r["cenwave"]) == cw
                  and r["slicer"] == sl and r["status"] != "excluded"]
            if not ok:
                L.append(f"| {cw} | {sl} | — | 0, all excluded | — | — | — | "
                         f"— | — |")
                continue
            p = peaks[(g, cw, sl)]
            # The common fit range, not the union: the combined mean is only
            # defined where every cube of the blaze contributes.
            lo = max(int(r["fit_lo_A"]) for r in ok)
            hi = min(int(r["fit_hi_A"]) for r in ok)
            mark = " *(lower limit)*" if p["status"] != "ok" else ""
            L.append(f"| {cw} | {sl} | {len({r['night'] for r in ok})} | "
                     f"{len(ok)} | {sum(int(r['nframes']) for r in ok)} | "
                     f"{lo}–{hi} | {p['combined_peak_pct']}%{mark} | "
                     f"{p['peak_wave_A']} A | "
                     f"{p['plateau_lo_A']}–{p['plateau_hi_A']} A |")
        out[g + "_ROLL"] = "\n".join(L)

        L = ["| night | setup | blaze | slicer | bin | standard | frames | "
             "exp (s) | airmass | fit range (A) | peak | at | interior mean | |",
             "|---|---|--:|---|---|---|--:|--:|--:|---|--:|--:|--:|---|"]
        for r in sorted(rr, key=lambda r: (int(r["cenwave"]), r["night"],
                                           float(r["airmass"]))):
            st = {"ok": "", "lower limit": "lower limit",
                  "excluded": "**excluded**"}[r["status"]]
            L.append(f"| {r['night']} | {r['setup']} | {r['cenwave']} | "
                     f"{r['slicer']} | {r['binning'].replace(',', '×')} | "
                     f"{r['star']} | {r['nframes']} | {r['exptime_s']} | "
                     f"{r['airmass']} | {r['fit_lo_A']}–{r['fit_hi_A']} | "
                     f"{r['peak_pct']}% | {r['peak_wave_A']} A | "
                     f"{r['interior_mean_pct']}% | {st} |")
        out[g + "_CUBES"] = "\n".join(L)
    return out


TEMPLATE = os.path.join(HERE, "grating_summary.template.md")
DOC = os.path.join(ROOT, "GRATING_SUMMARY.md")


def render(tables):
    body = open(TEMPLATE).read()
    for k, v in tables.items():
        body = body.replace("{{%s}}" % k, v)
    left = re.findall(r"\{\{(\w+)\}\}", body)
    if left:
        raise SystemExit(f"template asks for tables that do not exist: {left}")
    return body


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", default=None,
                    help="print one table only, e.g. OVERALL or RM2_ROLL")
    ap.add_argument("--render", action="store_true",
                    help="rewrite GRATING_SUMMARY.md from the template")
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if GRATING_SUMMARY.md is not what --render "
                         "would write")
    args = ap.parse_args()
    tables = build()
    if args.section:
        print(tables[args.section])
        return
    if args.render:
        open(DOC, "w").write(render(tables))
        print(f"wrote {DOC}")
        return
    if args.check:
        want = render(tables)
        have = open(DOC).read() if os.path.exists(DOC) else ""
        if want == have:
            print("GRATING_SUMMARY.md is in sync")
            return
        raise SystemExit("GRATING_SUMMARY.md is STALE -- run --render "
                         "(and put prose edits in the template, not the doc)")
    for k, v in tables.items():
        print(f"\n<!-- {k} -->\n{v}")


if __name__ == "__main__":
    main()
