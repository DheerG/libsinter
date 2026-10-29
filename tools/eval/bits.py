"""Where the bits go, from a run with BITS=1 (tools/eval/run.sh).

usage: bits.py build/eval/NAME.bits [GROUP...]

Per clip and QP: each frame type's share of the bytes, then each model
group's share of all bits (ideal code lengths, from snte_bits_dump), and
headers with escape bits. GROUP names select the columns (default: all
groups of SntCdfSet, as harness.c lists them).
"""
import sys, os, glob

def main():
    files = sorted(glob.glob(os.path.join(sys.argv[1], "*_q*.txt")),
                   key=lambda p: (os.path.basename(p).split("_q")[0], int(p.rsplit("_q", 1)[1][:-4])))
    want = sys.argv[2:]
    header = False
    for path in files:
        clip, qp = os.path.basename(path)[:-4].rsplit("_q", 1)
        groups, bits, other, ft = {}, {}, 0.0, [0, 0, 0]
        order = []
        for line in open(path):
            p = line.split()
            if p[0] == "group":
                groups[p[1]] = range(int(p[2]), int(p[2]) + int(p[3])); order.append(p[1])
            elif p[0] == "bits":
                bits[int(p[2])] = bits.get(int(p[2]), 0) + float(p[3])
            elif p[0] == "other":
                other += float(p[2])
            elif p[0] == "frametypes":
                ft = [int(x) for x in p[1:4]]
        cols = want or order
        if not header:
            print(f"{'clip':<8} {'QP':>3} {'key%':>5} {'P%':>5} {'B%':>5} |"
                  + "".join(f" {c:>6}" for c in cols) + "  hdr+esc  (% of all bits)")
            header = True
        tot = sum(bits.values()) + other
        T = sum(ft) or 1
        print(f"{clip:<8} {qp:>3} {100 * ft[0] / T:5.1f} {100 * ft[1] / T:5.1f} {100 * ft[2] / T:5.1f} |"
              + "".join(f" {100 * sum(bits.get(k, 0) for k in groups[c]) / tot:6.1f}" for c in cols)
              + f"  {100 * other / tot:7.1f}")

main()
