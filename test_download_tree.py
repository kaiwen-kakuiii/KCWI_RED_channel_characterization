"""Test the SHIPPED step-7 tree code from download_std_red.py (exec'd verbatim,
not a copy) against the trap that produced RM2/RH4 folders:

  night A (2024-12-22-like): RL science + bias/dark calibs whose own RGRATNAM
     says RM2 (wheel parked) -> must land in RL/, no RM2/ folder
  night B (2023-12-07-like): RL + RH1 science sharing one bias -> bias in both
"""
import os, re, shutil, collections, types
import pandas as pd

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "download_std_red.py")
code = open(SRC).read()
start = code.index("    grat = dict(zip(df.koaid")
end = code.index("absent += 1", start) + len("absent += 1")
snippet = "\n".join(l[4:] for l in code[start:end].split("\n"))  # dedent

work = os.path.join(os.path.dirname(os.path.abspath(__file__)), "step7_test")
shutil.rmtree(work, ignore_errors=True)
lev0 = os.path.join(work, "lev0")
calib_lev0 = os.path.join(work, "calib", "lev0")
os.makedirs(lev0), os.makedirs(calib_lev0)

SCI = {"KR.A.SCI1.fits": "RL", "KR.B.SCI2.fits": "RL", "KR.B.SCI3.fits": "RH1"}
CAL = ["KR.A.BIAS_RM2.fits", "KR.A.DARK_RH4.fits", "KR.B.BIAS_SHARED.fits"]
for k in SCI:
    open(os.path.join(lev0, k), "w").write("x")
for k in CAL:
    open(os.path.join(calib_lev0, k), "w").write("x")

ns = dict(
    os=os, collections=collections,
    KOAID_RE=re.compile(r"^KR\..*\.fits$"),
    args=types.SimpleNamespace(outdir=work),
    lev0=lev0, calib_lev0=calib_lev0,
    df=pd.DataFrame({"koaid": list(SCI), "rgratnam": list(SCI.values())}),
    nights={"KR.A.SCI1.fits": "2024-12-22", "KR.B.SCI2.fits": "2023-12-07",
            "KR.B.SCI3.fits": "2023-12-07"},
    # calib rows carry wrong/parked rgratnam on purpose -- must be ignored
    pairs={"KR.A.SCI1.fits": [{"koaid": "KR.A.BIAS_RM2.fits", "rgratnam": "RM2"},
                              {"koaid": "KR.A.DARK_RH4.fits", "rgratnam": "RH4"}],
           "KR.B.SCI2.fits": [{"koaid": "KR.B.BIAS_SHARED.fits", "rgratnam": "RM2"}],
           "KR.B.SCI3.fits": [{"koaid": "KR.B.BIAS_SHARED.fits", "rgratnam": "RM2"}]},
    linked=0, absent=0, cleared=0,
)
exec(snippet, ns)

got = {os.path.relpath(os.path.join(dp, f), os.path.join(work, "by_night"))
       for dp, _, fs in os.walk(os.path.join(work, "by_night")) for f in fs}
want = {"RL/2024-12-22/KR.A.SCI1.fits",
        "RL/2024-12-22/KR.A.BIAS_RM2.fits",
        "RL/2024-12-22/KR.A.DARK_RH4.fits",
        "RL/2023-12-07/KR.B.SCI2.fits",
        "RL/2023-12-07/KR.B.BIAS_SHARED.fits",
        "RH1/2023-12-07/KR.B.SCI3.fits",
        "RH1/2023-12-07/KR.B.BIAS_SHARED.fits"}
gratdirs = set(os.listdir(os.path.join(work, "by_night")))

assert got == want, f"layout wrong:\n extra={got-want}\n missing={want-got}"
assert gratdirs == {"RL", "RH1"}, f"stray grating folders: {gratdirs - {'RL','RH1'}}"
assert ns["absent"] == 0 and ns["linked"] == 7
# shared bias is the same inode in both grating dirs
a = os.stat(os.path.join(work, "by_night", "RL", "2023-12-07", "KR.B.BIAS_SHARED.fits"))
b = os.stat(os.path.join(work, "by_night", "RH1", "2023-12-07", "KR.B.BIAS_SHARED.fits"))
assert a.st_ino == b.st_ino
print("PASS: no RM2/RH4 dirs created; calibs follow science grating; shared bias hardlinked into both")

# second run = rebuild idempotence (stale clear + relink)
ns.update(linked=0, absent=0, cleared=0)
exec(snippet, ns)
got2 = {os.path.relpath(os.path.join(dp, f), os.path.join(work, "by_night"))
        for dp, _, fs in os.walk(os.path.join(work, "by_night")) for f in fs}
assert got2 == want and ns["cleared"] == 7
print("PASS: rebuild idempotent (7 stale cleared, same layout)")
shutil.rmtree(work)
