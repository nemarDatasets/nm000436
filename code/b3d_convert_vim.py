"""d3-dryad-vim-taps (Dryad doi:10.5061/dryad.cjsxksngp, CC0-1.0) -> BIDS-iEEG DERIVATIVE dataset (Bruno 2026-10-08 08:38).

The release has no unprocessed recording (FHC .apm / 'R files' never released). The P-File field `raw` is the authors'
software-filtered (steep digital low-pass near 8-10 kHz), offset-removed (zero-mean), 60 s-trimmed 48 kHz bipolar
microelectrode signal (evidence: ledger/d3-dryad-vim-taps/result.json, Jobs ieeg-b3w10-vim-insp-1008005934 and
ieeg-b3w10-vim-quant-1008010513). That field is repackaged here, one BIDS run per 60 s trial:
  P-File cell <cond>(row, col) -> sub-<NN>_task-taps_acq-<off|on>_run-<k>_ieeg (3 channels Ant/Cen/Pos = MATLAB columns 1-3)
  cond base -> acq-off (TAPS OFF), stim -> acq-on (TAPS ON); run k = chronological order (trialTimes) within the subject.
Data format: BrainVision IEEE_FLOAT_32 (BrainVision has no float64); units n/a (the release documents no unit for `raw`).
The float64 P-File values are cast to float32 (relative rounding <= 6e-8); the exact float64 originals are in sourcedata/.
Round trip: every .eeg is read back and compared bit for bit with float32(P-File raw); max |float32 - float64| is reported.
Usage: python b3d_convert_vim.py <work_dir> <bids_root>
"""
import glob
import json
import os
import re
import shutil
import sys
import zipfile
from collections import defaultdict

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b3w3_common import scrub_mat_header, wjson, wtsv  # noqa: E402
from b3w10_common import flag_strings, mat_strings, write_run  # noqa: E402

W, OUT = sys.argv[1], sys.argv[2]
PDIR = os.path.join(W, "extract", "P-Files")
RFDIR = os.path.join(W, "extract", "PNSxET Root Folder ")
DL = os.path.join(W, "sourcedata_dl")
SKIP_SRC = os.environ.get("SKIP_SOURCEDATA") == "1"
CH = ["Ant", "Cen", "Pos"]
TRACK = {"Ant": "anterior", "Cen": "center", "Pos": "posterior"}
DOI = "10.5061/dryad.cjsxksngp"
rep = {"recordings": [], "notes": [], "empty_cells": [], "mat_privacy": {}, "sourcedata": []}
VCAL = defaultdict(set)  # channel -> P-File voltageCalib values over all subjects


def val(f, ref):
    ds = f[ref]
    if ds.attrs.get("MATLAB_empty", 0):
        return None
    cls = ds.attrs.get("MATLAB_class", b"").decode()
    o = ds[()]
    if cls == "char":
        return "".join(chr(c) for c in np.asarray(o).flatten(order="F"))
    return np.asarray(o)


def num(v, k=None):
    if v is None:
        return None
    a = np.asarray(v, dtype=float).ravel()
    if a.size == 0:
        return None
    return a if k is None else float(a[k])


def orient(a, n):
    """HDF5 array of a MATLAB matrix -> (ntrials, n) with the time axis last."""
    a = np.asarray(a)
    if a.ndim == 1:
        a = a[None, :]
    if a.shape[-1] == n:
        return a
    if a.shape[0] == n:
        return a.T
    raise ValueError(f"cannot orient {a.shape} to n={n}")


def trial_pairs(tt, ntr, n):
    """trialTimes (MATLAB 1x2 per trial, possibly ntr x 2 or 2 x ntr) -> list of (start, end) ticks, or Nones."""
    if tt is None:
        return [(None, None)] * ntr
    m = np.asarray(tt, dtype=float).T  # MATLAB orientation
    cands = []
    if m.shape == (ntr, 2):
        cands.append(m)
    if m.shape == (2, ntr):
        cands.append(m.T)
    if m.size == 2 * ntr:
        cands += [m.reshape(ntr, 2), m.reshape(2, ntr).T, m.ravel(order="F").reshape(ntr, 2)]
    for c in cands:
        d = c[:, 1] - c[:, 0]
        if np.all((d >= n - 1) & (d <= n + 1000)):
            return [(float(a), float(b)) for a, b in c]
    a = np.asarray(tt, dtype=float).ravel()
    if ntr > 1 and a.size == 2 and n - 1 <= a[1] - a[0] <= n + 1000:
        rep["notes"].append(f"trialTimes holds one pair {a.tolist()} for a {ntr}-trial cell: assigned to trial 1; later trials n/a")
        return [(float(a[0]), float(a[1]))] + [(None, None)] * (ntr - 1)
    rep["notes"].append(f"trialTimes {np.asarray(tt).tolist()} not matched to {ntr} trials of {n} samples")
    return [(None, None)] * ntr


def mask_segments(mk, fs):
    """boolean vector -> list of (onset_sample, n_samples) of True runs."""
    m = np.asarray(mk).astype(bool).astype(np.int8)
    if not m.any():
        return []
    d = np.diff(np.concatenate([[0], m, [0]]))
    on, off = np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]
    return list(zip(on.tolist(), (off - on).tolist()))


files = sorted(glob.glob(os.path.join(PDIR, "P-CH1052-*.mat")))
assert len(files) == 9, files
if os.environ.get("ONLY"):  # smoke test on a subset
    files = [f for f in files if re.search(r"-(%s)\.mat$" % os.environ["ONLY"], f)]
tet = {}
for row in open(os.path.join(RFDIR, "Data", "T-Files", "TETRAS Data Log.csv"), encoding="utf-8-sig").read().splitlines()[1:]:
    c = row.split(",")
    tet[c[0].split("-")[-1]] = c[1]
parts = []
for fn in files:
    nn = re.search(r"P-CH1052-(\d+)\.mat$", fn).group(1)
    sub = f"sub-{nn}"
    D = os.path.join(OUT, sub, "ieeg")
    recs = []  # each: dict with cond,row,trial,chans{name: vector},meta
    sides, freqs, amps, implant, therapy = set(), set(), set(), set(), set()
    with h5py.File(fn, "r") as f:
        for cond in ("base", "stim"):
            g = f[cond]
            R = g["raw"][()]
            ncol, nrow = R.shape
            assert ncol == 3, (fn, cond, R.shape)
            for i in range(nrow):
                cells = {}
                for j in range(3):
                    ref = R[j, i]
                    ds = f[ref]
                    if ds.attrs.get("MATLAB_empty", 0) or ds.size < 10:
                        rep["empty_cells"].append(f"{os.path.basename(fn)} {cond}({i + 1},{j + 1})")
                        continue
                    meta = {}
                    for k in ("fs", "voltageCalib", "timeCalib", "recDepth", "chan", "chanCoord", "implantCoord", "therapyCoord",
                              "implantSide", "stimFreq", "stimMiliAmp", "trialTimes", "spkSTD", "blankingMask"):
                        if k in g:
                            meta[k] = val(f, g[k][()][j, i])
                    fs = num(meta["fs"], 0)
                    raw = np.asarray(ds[()], dtype=np.float64)
                    n = max(raw.shape) if raw.ndim == 2 and min(raw.shape) <= 4 else raw.shape[-1]
                    raw = orient(raw, n)
                    cells[j] = (raw, meta)
                if not cells:
                    continue
                if len(cells) < 3:
                    rep["notes"].append(f"{sub} {cond}({i + 1},:) only columns {[c + 1 for c in cells]} present")
                ns = {c[0].shape for c in cells.values()}
                assert len(ns) == 1, (fn, cond, i, ns)
                ntr, n = next(iter(ns))
                j0 = min(cells)
                m0 = cells[j0][1]
                fs = num(m0["fs"], 0)
                assert fs == 48000.0, fs
                pairs = trial_pairs(m0.get("trialTimes"), ntr, n)
                for jj, (_, mm) in cells.items():
                    assert num(mm.get("chan"), 0) == jj + 1, (fn, cond, i, jj, mm.get("chan"))
                    assert mm.get("recDepth") == m0.get("recDepth")
                if cond == "stim":
                    sides.add(m0.get("implantSide"))
                    sf, sa = num(m0.get("stimFreq"), 0), num(m0.get("stimMiliAmp"), 0)
                    freqs.add(sf if sf is not None and not np.isnan(sf) else None)
                    amps.add(sa)
                    if m0.get("implantCoord") is not None:
                        implant.add(tuple(num(m0["implantCoord"]).tolist()))
                    if m0.get("therapyCoord") is not None:
                        therapy.add(tuple(num(m0["therapyCoord"]).tolist()))
                for t in range(ntr):
                    masks = {}
                    for jj, (_, mm) in cells.items():
                        bm = mm.get("blankingMask")
                        if bm is not None and np.asarray(bm).size == ntr * n:
                            masks[CH[jj]] = orient(bm, n)[t].astype(bool)
                        elif bm is not None:
                            rep["notes"].append(f"{sub} {cond}({i + 1},{jj + 1}) blankingMask size {np.asarray(bm).shape} != {ntr}x{n}; not used")
                    recs.append(dict(cond=cond, row=i, trial=t, ntr=ntr, n=n, fs=fs, pair=pairs[t], cellstart=pairs[0][0],
                                     chans={CH[jj]: c[0][t] for jj, c in sorted(cells.items())},
                                     coords={CH[jj]: num(c[1].get("chanCoord")).tolist() if c[1].get("chanCoord") is not None else None for jj, c in cells.items()},
                                     vcal={CH[jj]: num(c[1].get("voltageCalib"), 0) for jj, c in cells.items()},
                                     depth=m0.get("recDepth"), tcal=num(m0.get("timeCalib"), 0),
                                     sf=num(m0.get("stimFreq"), 0) if cond == "stim" else None,
                                     sa=num(m0.get("stimMiliAmp"), 0) if cond == "stim" else None, masks=masks))
    # chronological order (trialTimes start); fallback to release order
    recs.sort(key=lambda r: (r["cellstart"] if r["cellstart"] is not None else 1e18, r["trial"], r["cond"], r["row"]))
    scans, elec = [], {}
    side = sorted(s for s in sides if s)
    side = side[0] if len(side) == 1 else ("n/a" if not side else "/".join(side))
    for k, r in enumerate(recs, start=1):
        acq = "off" if r["cond"] == "base" else "on"
        stem = f"{sub}_task-taps_acq-{acq}_run-{k:02d}"
        names = list(r["chans"].keys())
        X = np.stack([r["chans"][c] for c in names], axis=1)  # (n, nch) float64
        fs, n = r["fs"], r["n"]
        dl = r["depth"] or "n/a"
        z = r["coords"][names[0]][2] if r["coords"][names[0]] else None
        ch_extra = []
        for c in names:
            ename = f"{c}_d{re.sub(r'[^0-9A-Za-z]', '', dl)}"
            elec[ename] = (r["coords"][c], c, dl)
            VCAL[c].add(r["vcal"][c])
            ch_extra.append({"type": "DBS", "units": "n/a", "reference": "bipolar: macro contact 1 cm dorsal to this microelectrode tip",
                             "group": f"BenGun_{side}", "description": f"{TRACK[c]} microelectrode (P-File column {CH.index(c) + 1}) at depth label '{dl}'; position = electrodes.tsv row {ename}; authors' processed 'raw' field; P-File voltageCalib {r['vcal'][c]:.8f} (unit not documented)"})
        cell_txt = f"{r['cond']}({r['row'] + 1},{'1:3' if len(names) == 3 else '[' + ' '.join(str(CH.index(c) + 1) for c in names) + ']'})" + (f" trial {r['trial'] + 1} of {r['ntr']}" if r["ntr"] > 1 else "")
        sidecar = {
            "TaskName": "taps",
            "TaskDescription": "Intraoperative microelectrode recording at a fixed depth in/near the ventral intermediate thalamus (VIM) during DBS surgery for essential tremor, "
                               "without (acq-off, P-File 'base') or during (acq-on, P-File 'stim') transcutaneous afferent patterned stimulation (TAPS) of the treated hand.",
            "Instructions": "n/a",
            "SamplingFrequency": fs, "PowerLineFrequency": 60,
            "SoftwareFilters": {"AuthorsProcessing": {"Description": "Authors' MER_preprocess output (P-File field 'raw'): digital low-pass (not documented; measured: power falls 13-20 orders of magnitude from 6-8 kHz to 20-24 kHz), offset removed (zero mean), trimmed to exactly 60 s per trial. Not the unprocessed recording."}},
            "HardwareFilters": {"BandLimit": {"Description": "Guideline 4000 LP+ (FHC) band limit 20 kHz, digitised at 48 kHz (article)"}},
            "Manufacturer": "FHC", "ManufacturersModelName": "Guideline 4000 LP+",
            "iEEGReference": "Bipolar: each of 3 microelectrodes (Ben-gun, 2 mm apart antero-posterior) against the low-impedance macro contact 1 cm dorsal to its tip (article)",
            "iEEGGround": "n/a", "iEEGPlacementScheme": "Three parallel Ben-gun microelectrode tracks (anterior/center/posterior) advanced towards the planned VIM target",
            "ElectrodeManufacturer": "FHC", "DBSChannelCount": len(names), "SEEGChannelCount": 0, "ECOGChannelCount": 0,
            "RecordingType": "continuous", "RecordingDuration": round(n / fs, 6),
            "ElectricalStimulation": acq == "on",
            "ElectricalStimulationParameters": ("TAPS (transcutaneous afferent patterned stimulation, wearable peripheral nerve stimulator) on the treated wrist: "
                                                f"frequency {r['sf'] if r['sf'] is not None and not np.isnan(r['sf']) else 'n/a (NaN in P-File)'} Hz, amplitude {r['sa'] if r['sa'] is not None else 'n/a'} mA (P-File stimFreq/stimMiliAmp)") if acq == "on" else "n/a (TAPS off)",
            "SourceCell": f"{os.path.basename(fn)} {cell_txt}",
        }
        res = write_run(D, stem, X, fs, names, ["n/a"] * len(names), ["DBS"] * len(names), sidecar, ch_extra=ch_extra,
                        comment=f"b3d_convert_vim.py: {os.path.basename(fn)} {cell_txt} field 'raw' (float64, unit not documented) as float32")
        back = np.fromfile(os.path.join(D, stem + "_ieeg.eeg"), dtype="<f4").reshape(-1, len(names))
        exact32 = bool(back.shape == X.shape and np.array_equal(back, X.astype(np.float32)))
        err64 = float(np.max(np.abs(back.astype(np.float64) - X)))
        relmax = float(np.max(np.abs(back.astype(np.float64) - X) / np.maximum(np.abs(X), 1e-30)))
        # events: whole trial + authors' motion-artefact blanking segments
        ev = [[0.0, round(n / fs, 6), f"taps_{acq}_trial", 0, "n/a", dl]]
        if r["masks"]:
            ms = list(r["masks"].values())
            same = all(np.array_equal(ms[0], m) for m in ms[1:])
            groups = [("all", ms[0])] if same else [(c, m) for c, m in r["masks"].items()]
            for cname, m in groups:
                for s0, ln in mask_segments(m, fs):
                    ev.append([round(s0 / fs, 6), round(ln / fs, 6), "artifact_blanked_by_authors", s0, cname, dl])
            afrac = float(np.mean(np.any(np.stack(ms), axis=0)))
        else:
            afrac = None
        ev.sort(key=lambda e: (e[0], e[2] != f"taps_{acq}_trial"))
        wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "trial_type", "sample", "channel", "depth_label"], ev)
        t0 = r["pair"][0]
        scans.append([f"ieeg/{stem}_ieeg.vhdr", "n/a", acq, dl, z, cell_txt, t0, r["pair"][1],
                      round(t0 * r["tcal"] * 1e-6, 4) if t0 is not None and r["tcal"] else None,
                      (r["sf"] if r["sf"] is not None and not np.isnan(r["sf"]) else None) if acq == "on" else None,
                      r["sa"] if acq == "on" else None, round(afrac, 6) if afrac is not None else None])
        rr = dict(sub=sub, stem=stem, cell=cell_txt, depth=dl, acq=acq, n_samples=int(n), n_ch=len(names), sfreq=fs,
                  exact_float32_bitwise=exact32, max_abs_err_vs_float64_uV=err64, max_rel_err=relmax, mne_ok=res["roundtrip_ok"],
                  n_artifact_events=len(ev) - 1, artifact_fraction=afrac)
        rr["roundtrip_ok"] = bool(exact32 and res["roundtrip_ok"] and relmax < 1e-7)
        rep["recordings"].append(rr)
        print(stem, cell_txt, dl, n, "exact32", exact32, "err64", f"{err64:.3g}", "mne", res["roundtrip_ok"], "ev", len(ev) - 1, flush=True)
    wtsv(os.path.join(OUT, sub, f"{sub}_scans.tsv"),
         ["filename", "acq_time", "taps", "depth_label", "depth_mm_above_target", "source_cell", "trial_start_tick", "trial_end_tick",
          "trial_start_s_acquisition_clock", "taps_frequency_hz", "taps_amplitude_ma", "artifact_blanked_fraction"], scans)
    # electrodes (space-Other: offsets from the planned target along the trajectory, frame not stated by the authors)
    rows = []
    for en, (xyz, c, dl) in sorted(elec.items(), key=lambda kv: (CH.index(kv[1][1]), -(kv[1][0][2] if kv[1][0] else 0), kv[0])):
        x, y, zz = (xyz if xyz else ["n/a"] * 3)
        rows.append([en, x, y, zz, "n/a", f"BenGun_{side}", "microelectrode", TRACK[c], dl])
    for lab, S in (("DBS_implant", implant), ("DBS_therapy", therapy)):
        for q, xyz in enumerate(sorted(S)):
            rows.append([lab if len(S) == 1 else f"{lab}{q + 1}", *xyz, "n/a", "DBS_lead", "DBS lead position" if lab == "DBS_implant" else "therapeutic DBS contact", "n/a", "n/a"])
        if len(S) > 1:
            rep["notes"].append(f"{sub}: {len(S)} distinct {lab} coordinates {sorted(S)}")
    wtsv(os.path.join(D, f"{sub}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "group", "type", "track", "depth_label"], rows)
    wjson(os.path.join(D, f"{sub}_space-Other_electrodes.json"), {
        "name": {"Description": "<track>_d<depth label>: microelectrode track (Ant/Cen/Pos) at one recording depth; DBS_implant / DBS_therapy: final DBS lead / therapeutic contact (P-File implantCoord / therapyCoord)"},
        "x": {"Description": "P-File coordinate x (mm); frame not stated by the authors (see coordsystem.json)", "Units": "mm"},
        "y": {"Description": "P-File coordinate y (mm): +2 anterior, 0 center, -2 posterior Ben-gun track", "Units": "mm"},
        "z": {"Description": "P-File coordinate z (mm) = recording depth label, mm above the planned target along the trajectory", "Units": "mm"},
        "type": {"Description": "microelectrode, DBS lead position or therapeutic DBS contact"},
        "track": {"Description": "Ben-gun track (anterior/center/posterior)"},
        "depth_label": {"Description": "P-File recDepth label (e.g. '10a' = 10 mm above target)"},
        "SpatialReference": "none: no image; coordinates are offsets (mm) from the planned surgical target as stored in the P-Files (see coordsystem.json)"})
    wjson(os.path.join(D, f"{sub}_space-Other_coordsystem.json"), {
        "iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "mm",
        "iEEGCoordinateSystemDescription": "P-File chanCoord/implantCoord/therapyCoord in mm. The frame is NOT stated by the authors: the values (x 0, y +2/0/-2 for the anterior/center/posterior Ben-gun microelectrodes, z = recording depth in mm above the planned target) look like offsets from the planned VIM target along the surgical trajectory. Not MNI, not Talairach, not a scanner/ACPC frame; do not use them for atlas lookup."})
    sf = sorted(x for x in freqs if x is not None)
    sa = sorted(x for x in amps if x is not None)
    parts.append([sub, "n/a", "n/a", "essential tremor", side, tet.get(nn, "n/a"), "/".join(f"{x:g}" for x in sf) or None, "/".join(f"{x:g}" for x in sa) or None])

wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "diagnosis", "implant_side", "taps_treated_hand", "taps_frequency_hz", "taps_amplitude_ma"], parts)
wjson(os.path.join(OUT, "participants.json"), {
    "age": {"Description": "Not released per subject (article: group 65 +/- 7 years)", "Units": "year"},
    "sex": {"Description": "Not released per subject (article: 2 of 9 female)"},
    "diagnosis": {"Description": "Medication-refractory essential tremor (article)"},
    "implant_side": {"Description": "Hemisphere of the DBS implant / recording (P-File implantSide)"},
    "taps_treated_hand": {"Description": "Hand that received TAPS (TETRAS Data Log.csv treated_hand)"},
    "taps_frequency_hz": {"Description": "TAPS frequency (P-File stimFreq); n/a where the P-File stores NaN", "Units": "Hz"},
    "taps_amplitude_ma": {"Description": "TAPS amplitude (P-File stimMiliAmp)", "Units": "mA"}})
wjson(os.path.join(OUT, "task-taps_events.json"), {
    "trial_type": {"Description": "Event type", "Levels": {"taps_off_trial": "whole 60 s trial without TAPS (P-File base)", "taps_on_trial": "whole 60 s trial during TAPS (P-File stim)",
                                                         "artifact_blanked_by_authors": "segment marked by the authors' motion-artefact blanking mask (P-File blankingMask); the samples in the data files are NOT blanked"}},
    "sample": {"Description": "Onset sample (0-based, 48 kHz)"},
    "channel": {"Description": "Channels the blanking segment applies to ('all' = identical mask on every channel)"},
    "depth_label": {"Description": "P-File recDepth label of the recording"}})
for s in parts:
    wjson(os.path.join(OUT, s[0], f"{s[0]}_scans.json"), {
        "taps": {"Description": "TAPS state: off (P-File base) or on (P-File stim)"},
        "depth_label": {"Description": "P-File recDepth label ('10a' = 10 mm above the planned target; suffixes such as '_3030' or '_filt' as released, meaning not documented)"},
        "depth_mm_above_target": {"Description": "P-File chanCoord z", "Units": "mm"},
        "source_cell": {"Description": "P-File struct cell <base|stim>(row, columns) (MATLAB 1-based) and trial index for multi-trial cells"},
        "trial_start_tick": {"Description": "P-File trialTimes start (acquisition clock ticks of timeCalib = 20.833 us)"},
        "trial_end_tick": {"Description": "P-File trialTimes end"},
        "trial_start_s_acquisition_clock": {"Description": "trial_start_tick x timeCalib: seconds on the FHC acquisition clock (relative, not a date)", "Units": "s"},
        "taps_frequency_hz": {"Description": "P-File stimFreq (n/a if NaN or TAPS off)", "Units": "Hz"},
        "taps_amplitude_ma": {"Description": "P-File stimMiliAmp (TAPS on only)", "Units": "mA"},
        "artifact_blanked_fraction": {"Description": "Fraction of samples in the authors' blanking mask (any channel)"}})

# ---- privacy: strings in the P-Files ----
for fn in files:
    s = mat_strings(fn)
    rep["mat_privacy"][os.path.basename(fn)] = {"strings": sorted(s)[:100], "flagged": flag_strings(s), "header_scrubbed_view": scrub_mat_header(open(fn, "rb").read(116))[0].decode("latin-1").strip()}

# ---- sourcedata: P-Files (MAT header creation date -> day 01), Dryad README, root folder zip (MAT headers + member dates -> day 01) ----
S = os.path.join(OUT, "sourcedata")
os.makedirs(os.path.join(S, "P-Files"), exist_ok=True)
if not SKIP_SRC:
    for fn in files:
        dst = os.path.join(S, "P-Files", os.path.basename(fn))
        shutil.copyfile(fn, dst)
        with open(dst, "r+b") as fh:
            head = fh.read(128)
            nh, ch = scrub_mat_header(head)
            fh.seek(0)
            fh.write(nh[:128])
        new = open(dst, "rb").read(116).decode("latin-1").strip()
        rep["sourcedata"].append({"file": f"sourcedata/P-Files/{os.path.basename(fn)}", "header_scrubbed": ch, "header_now": new, "bytes": os.path.getsize(dst)})
        print("sourcedata", dst, ch, new, flush=True)
shutil.copyfile(os.path.join(DL, "README.md"), os.path.join(S, "README.md"))
zin = os.path.join(DL, "PNSxET_Root_Folder_.zip")
zout = os.path.join(S, "PNSxET_Root_Folder_.zip")
changed = []
with zipfile.ZipFile(zin) as zi, zipfile.ZipFile(zout, "w", compression=zipfile.ZIP_DEFLATED) as zo:
    for it in zi.infolist():
        b = zi.read(it.filename)
        if it.filename.lower().endswith(".mat"):
            b2, ch = scrub_mat_header(b)
            if ch:
                changed.append(it.filename)
            b = b2
        y, mo = it.date_time[0], it.date_time[1]
        ni = zipfile.ZipInfo(it.filename, date_time=(y, mo, 1, 0, 0, 0))
        ni.compress_type = zipfile.ZIP_DEFLATED if not it.is_dir() else zipfile.ZIP_STORED
        ni.external_attr = it.external_attr
        zo.writestr(ni, b)
rep["sourcedata"].append({"file": "sourcedata/PNSxET_Root_Folder_.zip", "members": len(zipfile.ZipFile(zout).infolist()), "mat_headers_scrubbed": changed, "member_dates": "day set to 01"})
open(os.path.join(S, "README_sourcedata.txt"), "w").write(
    "Original Dryad release doi:10.5061/dryad.cjsxksngp (CC0-1.0), downloaded 2026-10 and SHA-256 verified against Dryad.\n"
    "- P-Files/P-CH1052-NN.mat: the authors' processed MATLAB v7.3 files (fields raw, lfp, spk, rawLPCheby, rawBPCheby, raw_blanked, spk_blanked,\n"
    "  blankingMask, pxx, f, spike waveforms, coordinates, TAPS settings, ...). Identical to Dryad except the 116-byte MAT text header,\n"
    "  whose file-creation date was set to day 01 (privacy rule); the HDF5 content is unchanged.\n"
    "- PNSxET_Root_Folder_.zip: the authors' analysis scripts, TETRAS scores and stimVbaseSPK.mat; repacked with member dates set to day 01\n"
    "  and MAT header dates set to day 01; file contents otherwise unchanged.\n"
    "- README.md: the Dryad README (unchanged).\n"
    "- E-Files.zip (P-Files + binned spectra) is not included; get it from Dryad.\n")

# ---- code, metadata ----
os.makedirs(os.path.join(OUT, "code"), exist_ok=True)
here = os.path.dirname(os.path.abspath(__file__))
for x in ("b3d_convert_vim.py", "b3w10_common.py", "b3w3_common.py"):
    shutil.copy(os.path.join(here, x), os.path.join(OUT, "code", x))
n_rec = len(rep["recordings"])
hrs = sum(r["n_samples"] / r["sfreq"] for r in rep["recordings"]) / 3600
rt_all = all(r["roundtrip_ok"] for r in rep["recordings"])
rep["summary"] = {"subjects": [p[0] for p in parts], "runs": n_rec, "hours": hrs, "rt_all": rt_all,
                  "max_abs_err_vs_float64_uV": max(r["max_abs_err_vs_float64_uV"] for r in rep["recordings"]),
                  "max_rel_err": max(r["max_rel_err"] for r in rep["recordings"])}
json.dump(rep, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1, default=str)
AUTH = ["Cuong P. Luu", "Jordan Ranum", "Youngwon Youn", "Jennifer L. Perrault", "Bryan M. Krause", "Matthew I. Banks",
        "Laura Buyan-Dent", "Kip A. Ludwig", "Wendell B. Lake", "Aaron J. Suminski"]
wjson(os.path.join(OUT, "dataset_description.json"), {
    "Name": "Intraoperative VIM thalamus microelectrode recordings during transcutaneous afferent patterned stimulation (TAPS) for essential tremor (Luu et al.; derivative: authors' processed 48 kHz signals)",
    "BIDSVersion": "1.10.0", "DatasetType": "derivative", "License": "CC0-1.0", "Authors": AUTH,
    "Acknowledgements": "Data and code released on Dryad by Cuong P. Luu et al. (University of Wisconsin-Madison). Repackaged to BIDS for NEMAR.",
    "HowToAcknowledge": "Cite Luu et al. (2025) Brain Stimulation, doi:10.1016/j.brs.2025.06.004, and the Dryad dataset doi:10.5061/dryad.cjsxksngp.",
    "Funding": ["University of Wisconsin-Madison School of Medicine and Public Health", "Partial funding of AJS and WBL from Cala Health, Inc. (preprint)"],
    "EthicsApprovals": ["University of Wisconsin-Madison Institutional Review Board, protocol 2018-1052"],
    "ReferencesAndLinks": ["https://doi.org/10.5061/dryad.cjsxksngp", "https://doi.org/10.1016/j.brs.2025.06.004", "https://doi.org/10.1101/2024.12.02.24317799"],
    "Keywords": ["essential tremor", "ventral intermediate nucleus", "thalamus", "microelectrode recording", "deep brain stimulation", "peripheral nerve stimulation", "TAPS", "intraoperative"],
    "SourceDatasets": [{"DOI": "doi:10.5061/dryad.cjsxksngp", "URL": "https://datadryad.org/dataset/doi:10.5061/dryad.cjsxksngp", "Version": "Dryad version 412969 (published 2025-12-10)"}],
    "GeneratedBy": [
        {"Name": "Authors' processing (MER_loadData -> MER_preprocess)", "Description": "FHC Guideline 4000 LP+ recordings (48 kHz, .apm) loaded into 'R files' (not released) and processed into the released P-Files by MER_preprocess.m (not released): field 'raw' is software low-pass filtered (steep digital roll-off near 8-10 kHz, measured), offset-removed (zero mean) and trimmed to exactly 60 s per trial; motion-artefact masks computed (blankingMask).",
         "CodeURL": "https://doi.org/10.5061/dryad.cjsxksngp"},
        {"Name": "b3d_convert_vim.py", "Description": "Repackaging of the P-File field 'raw' (float64, unit not documented) to BrainVision IEEE float32 (relative rounding <= 6e-8, bit-exact to float32 of the P-File values), one run per 60 s trial; TAPS state as acq-off/on; depth, source cell and TAPS settings in scans.tsv; authors' blanking masks as events; coordinates as released (space-Other). Code in code/.", "Version": "2026-10-08"}],
})
VC = {c: v for c, (v,) in VCAL.items()}  # one voltageCalib per electrode in all P-Files
open(os.path.join(OUT, "README"), "w").write(f"""# VIM thalamus microelectrode recordings during TAPS for essential tremor (Luu et al. 2025) - DERIVATIVE

Intraoperative microelectrode recordings (MER) from the ventral intermediate nucleus of the thalamus (VIM) region of
9 adults with medication-refractory essential tremor undergoing DBS surgery (University of Wisconsin-Madison, IRB 2018-1052).
At each recording depth, 60 s were recorded without and with transcutaneous afferent patterned stimulation (TAPS) of the
treated hand, a wearable peripheral nerve stimulation therapy. Remifentanil was stopped 15 min before mapping (article).

## THIS IS A DERIVATIVE DATASET: the original recordings were never released
- The original FHC recordings (.apm files) and the authors' intermediate 'R files' (output of MER_loadData.m) are NOT part
  of the Dryad release. The release contains only the authors' processed MATLAB files ('P-' = "Processed neural data after
  signal processing and artifact removal", Dryad README) and derived files.
- The signals here are the P-File field `raw` exactly as released (cast to float32): the authors' **software-filtered,
  offset-removed (zero-mean), 60 s-trimmed** signals. Evidence (measured on all files in Voyager Jobs
  ieeg-b3w10-vim-insp-1008005934 and ieeg-b3w10-vim-quant-1008010513):
  - every trace is exactly 2,880,000 samples = 60.0 s at 48 kHz per trial, while the acquisition timestamps (trialTimes)
    span 2,880,013-2,880,061 ticks: trimmed to exactly 60 s;
  - the mean of each 60 s trace is 4.9e-8..7.2e-4 against a standard deviation of 472..7285 (|mean| / std <= 6.5e-7;
    code/receipts/2026-10-22_units-na-voltagecalib.md): offset removed in software;
  - the power spectrum falls by 13-20 orders of magnitude between 6-8 kHz and 20-24 kHz, i.e. a steep double-precision
    digital low-pass near 8-10 kHz, whereas the article states a 20 kHz hardware band limit;
  - the values are not on an ADC quantisation grid: no trace is on the integer grid or on the voltageCalib grid
    (code/receipts/2026-10-22_units-na-voltagecalib.md), whereas the authors' FHC reader (APMReadData.m) reads .apm
    samples as integers and stores voltageCalib beside them without applying it; so they are not a scaled copy of the
    acquired samples.
- Further derived fields of the P-Files (lfp, spk, rawLPCheby, rawBPCheby, raw_blanked, spk_blanked, spike waveforms, PSDs)
  are not converted; they remain available in sourcedata/P-Files.

## Contents
- {len(parts)} participants (`sub-<NN>` = study subject number 02, 03, 06-12 of IRB study CH1052), {n_rec} runs, {hrs:.2f} h at 48 kHz.
- `task-taps_acq-off` = P-File `base` (TAPS off); `task-taps_acq-on` = P-File `stim` (TAPS on). Runs are numbered in
  chronological order of acquisition within each subject (P-File trialTimes). A few P-File cells hold two consecutive 60 s
  trials; each trial is a run, placed right after the cell's first trial (trialTimes gives only the first trial's times,
  so later trials have n/a trial ticks in scans.tsv).
- `sub-<NN>_scans.tsv`: TAPS state, depth label and depth (mm above target), source P-File cell, acquisition-clock trial
  times, TAPS frequency/amplitude, fraction of samples in the authors' artefact mask.
- Channels `Ant`, `Cen`, `Pos`: the anterior, center and posterior Ben-gun microelectrodes (2 mm apart), each recorded
  bipolar against the macro contact 1 cm dorsal to its tip (article); type DBS, units n/a: the release gives no unit for
  `raw` (Dryad README: '-'); the authors' figure code (Figures_Paper.m) plots it unscaled with a uV axis label. Each P-File
  cell also gives voltageCalib: {', '.join(f'{c} {VC[c]:.8f}' for c in CH)}, the same in all {len(parts)} P-Files
  (Dryad README: 'Voltage calibration factor', unit '-'; not applied here).
- `events.tsv`: the whole trial (`taps_off_trial` / `taps_on_trial`) and the segments of the authors' motion-artefact
  blanking mask (`artifact_blanked_by_authors`). The data samples themselves are not blanked.
- Data format: BrainVision IEEE float32. BrainVision has no float64; the float64 P-File values are stored as float32
  (relative rounding <= 6e-8, absolute <= {rep['summary']['max_abs_err_vs_float64_uV']:.1g}).
  Every file was read back and is bit-identical to float32 of the P-File values. The float64 originals are in sourcedata/.
- Electrodes: `space-Other` per subject, with the P-File coordinates in mm (chanCoord for each track and depth, implantCoord
  and therapyCoord for the final DBS lead and therapeutic contact). CAVEAT: the authors do not state the frame. The values
  (x 0, y +2/0/-2 for the anterior/center/posterior tracks, z = depth in mm above the planned target) look like offsets
  from the planned VIM target along the trajectory. They are NOT MNI/Talairach coordinates; do not use them for atlas lookups.
- participants.tsv: age and sex are n/a (the article reports only group values: 65 +/- 7 years, 2 of 9 female);
  implant side (P-File), TAPS-treated hand (TETRAS log), TAPS settings (P-File).
- Empty P-File cells (no recording) are skipped; they are listed in code/conversion_report.json.
- sourcedata/: P-Files (MAT header creation date set to day 01, HDF5 content unchanged), the Dryad README, and the root
  folder zip (analysis scripts, TETRAS scores, stimVbaseSPK.mat; member and MAT header dates set to day 01).

## Source
- Dryad: {', '.join(AUTH)}. "Real-time VIM thalamus recordings during peripheral nerve stimulation treatment for essential
  tremor: DBS intraoperative dataset". doi:10.5061/dryad.cjsxksngp (CC0-1.0).
- Article: Luu et al. (2025) "Wearable peripheral nerve stimulator reduces essential tremor symptoms through targeted brain
  modulation", Brain Stimulation, doi:10.1016/j.brs.2025.06.004. Preprint doi:10.1101/2024.12.02.24317799.
- Funding: University of Wisconsin-Madison School of Medicine and Public Health; AJS and WBL were partly funded by Cala Health, Inc. (preprint).

## Privacy
The release contains study codes only (CH1052-NN). Strings in the P-Files are depth labels and implant side. MAT header
file-creation dates (June 2024, processing dates) were set to day 01. No names or recording dates were found.
""")
open(os.path.join(OUT, "CHANGES"), "w").write("1.0.0 2026-10-08\n  - Derivative repackaging of the P-File 'raw' signals of Dryad doi:10.5061/dryad.cjsxksngp for NEMAR.\n")
open(os.path.join(OUT, ".bidsignore"), "w").write("code/\n")
print(json.dumps(rep["summary"]))
print("empty_cells", len(rep["empty_cells"]), "notes", rep["notes"][:20])
print("mat_flags", {k: v["flagged"] for k, v in rep["mat_privacy"].items() if v["flagged"]})
print("CONVERT_DONE")
