"""BD-rates between runs of tools/eval/run.sh (or x264.sh).

usage: bd.py BASE NAME... [--ratio]

For each NAME, per clip: the BD-rate against BASE on luma PSNR, YUV PSNR
(Y:U:V weighted 6:1:1) and luma SSIM (in dB, -10 log10(1 - SSIM)), then
the mean over clips, encoder/decoder mismatches and the slowest encode
and decode per frame. With --ratio, BD-rates print as bitrate ratios
(1.24x), for comparisons against x264. Classic Bjontegaard: cubic fit of
log rate against quality, over the overlapping quality range.
"""
import sys, os, math
import numpy as np

EVAL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "build", "eval")

def load(name):
    pts = {}
    for line in open(os.path.join(EVAL, name + ".txt")):
        f = line.split()
        clip, rate, y, u, v, s = f[0], int(f[3]), float(f[7]), float(f[8]), float(f[9]), float(f[10])
        pts.setdefault(clip, []).append(dict(rate=rate, y=y, yuv=(6 * y + u + v) / 8,
                                             ssim=-10 * math.log10(1 - s), mm=int(f[11]),
                                             enc=float(f[12]), dec=float(f[13])))
    return pts

def bd(a, t, key):
    ra, qa = np.log([p["rate"] for p in a]), np.array([p[key] for p in a])
    rt, qt = np.log([p["rate"] for p in t]), np.array([p[key] for p in t])
    lo, hi = max(qa.min(), qt.min()), min(qa.max(), qt.max())
    if hi <= lo:
        return float("nan")
    fa, ft = np.polyint(np.polyfit(qa, ra, 3)), np.polyint(np.polyfit(qt, rt, 3))
    d = ((np.polyval(ft, hi) - np.polyval(ft, lo)) - (np.polyval(fa, hi) - np.polyval(fa, lo))) / (hi - lo)
    return math.exp(d)

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    ratio = "--ratio" in sys.argv
    base = load(args[0])
    fmt = (lambda r: f"{r:.3f}x") if ratio else (lambda r: f"{100 * (r - 1):+.2f}%")
    for name in args[1:]:
        t = load(name)
        clips = [c for c in base if c in t]
        print(f"{name} vs {args[0]}   (Y / YUV / SSIM)")
        sums = np.zeros(3)
        for c in clips:
            r = [bd(base[c], t[c], k) for k in ("y", "yuv", "ssim")]
            sums += r
            print(f"  {c:<8} " + "  ".join(fmt(x) for x in r))
        mm = sum(p["mm"] for c in t for p in t[c])
        enc = max(p["enc"] for c in t for p in t[c])
        dec = max(p["dec"] for c in t for p in t[c])
        print(f"  {'mean':<8} " + "  ".join(fmt(x) for x in sums / len(clips))
              + f"   mismatches {mm}, slowest encode {enc:.0f} ms, decode {dec:.1f} ms per frame")

main()
