"""voltageCalib per P-File, header units, amplitude statistics and clipping plateaus of nm000436.

Usage:
    python 2026-10-22_units-na-voltagecalib.py <nm000436 root> <git clone> <P-File dir>
"""

import hashlib
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REV = "origin/main"
CHANNELS = ("Ant", "Cen", "Pos")
CLIP_FROM = 25000  # |value| above which clipping plateaus are searched
BIN = 100
MIN_MODE = 100
PEAK = 2
CELL = r"(base|stim)\((\d+),1:3\)(?: trial (\d) of \d)?"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def git(clone: Path, *args: str) -> str:
    out = subprocess.run(["git", "-C", clone, *args], capture_output=True, text=True)
    return out.stdout.strip() if out.returncode == 0 else "n/a"


def annex_ok(path: Path, clone: Path, rel: str) -> bool:
    key = git(clone, "cat-file", "-p", f"{REV}:{rel}").rsplit("--", 1)[-1]
    return key.startswith(sha256(path))


def cell(f: h5py.File, group: str, field: str, row: int, col: int) -> np.ndarray:
    return np.asarray(f[f[group][field][row, col]]).ravel()


def plateaus(x: np.ndarray) -> tuple[float, int, float, int]:
    """Modal BIN-wide bin (centre, count) beyond +CLIP_FROM and beyond -CLIP_FROM.

    count is 0 unless the bin holds >= MIN_MODE samples and >= PEAK x the first bin
    (a peak, not the decaying tail).
    """
    edges = np.arange(CLIP_FROM, 50000 + BIN, BIN)
    out = []
    for sign in (1, -1):
        counts, _ = np.histogram(sign * x, edges)
        i = counts.argmax()
        peak = counts[i] >= MIN_MODE and counts[i] >= PEAK * counts[0]
        out += [sign * (edges[i] + BIN / 2), int(counts[i]) if peak else 0]
    return tuple(out)


def main(root: Path, clone: Path, pfiles: Path) -> None:
    print(f"h5py {h5py.__version__}, numpy {np.__version__}, pandas {pd.__version__}")
    print(f"{REV} of the clone: {git(clone, 'rev-parse', REV)}")

    calib = []
    print(
        "\n## voltageCalib read from each P-File (fields /base/voltageCalib and /stim/voltageCalib)"
    )
    print(
        "| P-File | sha256 = annex key | dtype | non-empty cells (base + stim) | Ant | Cen | Pos |"
    )
    print("|---|---|---|---:|---|---|---|")
    for path in sorted(pfiles.glob("P-CH1052-*.mat")):
        rel = f"sourcedata/P-Files/{path.name}"
        cells = []
        with h5py.File(path, "r") as f:
            for group in ("base", "stim"):
                refs = f[group]["voltageCalib"]
                for row in range(refs.shape[0]):
                    for col in range(refs.shape[1]):
                        dataset = f[refs[row, col]]
                        v = np.asarray(dataset).ravel()
                        if v.size == 1:  # empty cells hold [0, 0]
                            cells.append(
                                (
                                    path.name,
                                    group,
                                    CHANNELS[row],
                                    float(v[0]),
                                    str(dataset.dtype),
                                )
                            )
        df = pd.DataFrame(
            cells, columns=["pfile", "group", "channel", "voltageCalib", "dtype"]
        )
        values = " | ".join(
            "; ".join(
                repr(float(v)) for v in df[df.channel == ch].voltageCalib.unique()
            )
            for ch in CHANNELS
        )
        n = df.group.value_counts()
        print(
            f"| {rel} | {annex_ok(path, clone, rel)} | {', '.join(df.dtype.unique())} "
            f"| {n.get('base', 0)} + {n.get('stim', 0)} | {values} |"
        )
        calib.append(df)
    calib = pd.concat(calib)
    per_ch = calib.groupby("channel").agg(
        distinct=("voltageCalib", "nunique"),
        cells=("voltageCalib", "size"),
        pfiles=("pfile", "nunique"),
    )
    print(f"per electrode over all P-Files: {per_ch.to_dict('index')}")
    factor = calib.groupby("channel").voltageCalib.first().to_dict()

    vhdrs = sorted(root.glob("sub-*/ieeg/*_ieeg.vhdr"))
    units = []
    for vhdr in vhdrs:
        lines = vhdr.read_text(encoding="utf-8").splitlines()
        units += [
            line.rsplit(",", 1)[-1] for line in lines if re.match(r"Ch\d+=", line)
        ]
    headers_ok = [annex_ok(v, clone, v.relative_to(root).as_posix()) for v in vhdrs]
    print(
        f"\nBrainVision headers: {len(vhdrs)}, sha256 = git-annex key: {sum(headers_ok)}; "
        f"channel records {len(units)}, units {pd.Series(units).value_counts().to_dict()}"
    )

    rows, files = [], []
    for scans in sorted(root.glob("sub-*/sub-*_scans.tsv")):
        sub = scans.parent.name[4:]
        pfile = pfiles / f"P-CH1052-{sub}.mat"
        with h5py.File(pfile, "r") as f:
            for run in pd.read_csv(scans, sep="\t").itertuples():
                match = re.fullmatch(CELL, run.source_cell)
                group, col, trial = match[1], int(match[2]), int(match[3] or 1)
                vhdr = scans.parent / run.filename
                rel = vhdr.with_suffix(".eeg").relative_to(root).as_posix()
                files.append(annex_ok(root / rel, clone, rel))
                # MULTIPLEXED float32, channels Ant, Cen, Pos
                stored = np.fromfile(root / rel, dtype="<f4").reshape(-1, 3).T
                for k, ch in enumerate(CHANNELS):
                    # HDF5 holds the MATLAB (samples x trials) array transposed
                    x = np.asarray(f[f[group]["raw"][k, col - 1]])[trial - 1]
                    pos, n_pos, neg, n_neg = plateaus(x)
                    rows.append(
                        {
                            "sub": sub,
                            "file": vhdr.name.removesuffix("_ieeg.vhdr"),
                            "cell": f"{group}({col},{k + 1}) trial {trial}",
                            "channel": ch,
                            "n": x.size,
                            "bids_eq_float32": np.array_equal(np.float32(x), stored[k]),
                            "max_rounding": np.abs(stored[k] - x).max(),
                            "max_rel_rounding": (np.abs(stored[k] - x) / np.abs(x))[
                                x != 0
                            ].max(),
                            "integer": np.mean(x == np.round(x)),
                            "calib_multiple": np.mean(
                                np.abs(x / factor[ch] - np.round(x / factor[ch])) < 1e-6
                            ),
                            "mean": x.mean(),
                            "std": x.std(),
                            "pos": pos,
                            "n_pos": n_pos,
                            "neg": neg,
                            "n_neg": n_neg,
                        }
                    )
    df = pd.DataFrame(rows)
    print(f"\nBrainVision files: {len(files)}, sha256 = git-annex key: {sum(files)}")
    print(
        f"traces: {len(df)}, all {df.n.unique()} samples; BIDS = float32(P-File raw): "
        f"{df.bids_eq_float32.sum()}; max |BIDS - raw| {df.max_rounding.max():.5f}, "
        f"max |BIDS - raw| / |raw| {df.max_rel_rounding.max():.1e}"
    )
    print(
        f"max fraction of samples per trace that are integers: {df.integer.max():.1e}, "
        f"integer multiples of voltageCalib: {df.calib_multiple.max():.1e}"
    )
    print(
        f"|mean| {df['mean'].abs().min():.1e} to {df['mean'].abs().max():.1e}; "
        f"std {df['std'].min():.0f} to {df['std'].max():.0f}; "
        f"max |mean|/std {(df['mean'].abs() / df['std']).max():.1e}"
    )

    clipped = df[(df.n_pos > 0) & (df.n_neg > 0)]
    print(
        f"\n## traces with a plateau beyond +{CLIP_FROM} and beyond -{CLIP_FROM}: "
        f"modal {BIN}-wide bin with >= {MIN_MODE} samples and >= {PEAK} x the first bin "
        f"({len(clipped)} of {len(df)})"
    )
    print("| file | cell | channel | + plateau (n) | - plateau (n) | span |")
    print("|---|---|---|---:|---:|---:|")
    for r in clipped.itertuples():
        print(
            f"| {r.file} | {r.cell} | {r.channel} | {r.pos:.0f} ({r.n_pos}) "
            f"| {r.neg:.0f} ({r.n_neg}) | {r.pos - r.neg:.0f} |"
        )


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
