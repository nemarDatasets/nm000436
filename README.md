[![DOI](https://img.shields.io/badge/DOI-10.82901%2Fnemar.nm000436-blue)](https://doi.org/10.82901/nemar.nm000436)

# VIM thalamus microelectrode recordings during TAPS for essential tremor (Luu et al. 2025) - DERIVATIVE

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
- 9 participants (`sub-<NN>` = study subject number 02, 03, 06-12 of IRB study CH1052), 84 runs, 1.40 h at 48 kHz.
- `task-taps_acq-off` = P-File `base` (TAPS off); `task-taps_acq-on` = P-File `stim` (TAPS on). Runs are numbered in
  chronological order of acquisition within each subject (P-File trialTimes). A few P-File cells hold two consecutive 60 s
  trials; each trial is a run, placed right after the cell's first trial (trialTimes gives only the first trial's times,
  so later trials have n/a trial ticks in scans.tsv).
- `sub-<NN>_scans.tsv`: TAPS state, depth label and depth (mm above target), source P-File cell, acquisition-clock trial
  times, TAPS frequency/amplitude, fraction of samples in the authors' artefact mask.
- Channels `Ant`, `Cen`, `Pos`: the anterior, center and posterior Ben-gun microelectrodes (2 mm apart), each recorded
  bipolar against the macro contact 1 cm dorsal to its tip (article); type DBS, units n/a: the release gives no unit for
  `raw` (Dryad README: '-'); the authors' figure code (Figures_Paper.m) plots it unscaled with a uV axis label. Each P-File
  cell also gives voltageCalib: Ant 0.02195717, Cen 0.02259168, Pos 0.02262028, the same in all 9 P-Files
  (Dryad README: 'Voltage calibration factor', unit '-'; not applied here). The µV in the current BrainVision headers
  is a repackaging label, not a documented unit; the converter in code/ writes n/a there.
- `events.tsv`: the whole trial (`taps_off_trial` / `taps_on_trial`) and the segments of the authors' motion-artefact
  blanking mask (`artifact_blanked_by_authors`). The data samples themselves are not blanked.
- Data format: BrainVision IEEE float32. BrainVision has no float64; the float64 P-File values are stored as float32
  (relative rounding <= 6e-8, absolute <= 0.002).
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
- Dryad: Cuong P. Luu, Jordan Ranum, Youngwon Youn, Jennifer L. Perrault, Bryan M. Krause, Matthew I. Banks, Laura Buyan-Dent, Kip A. Ludwig, Wendell B. Lake, Aaron J. Suminski. "Real-time VIM thalamus recordings during peripheral nerve stimulation treatment for essential
  tremor: DBS intraoperative dataset". doi:10.5061/dryad.cjsxksngp (CC0-1.0).
- Article: Luu et al. (2025) "Wearable peripheral nerve stimulator reduces essential tremor symptoms through targeted brain
  modulation", Brain Stimulation, doi:10.1016/j.brs.2025.06.004. Preprint doi:10.1101/2024.12.02.24317799.
- Funding: University of Wisconsin-Madison School of Medicine and Public Health; AJS and WBL were partly funded by Cala Health, Inc. (preprint).

## Privacy
The release contains study codes only (CH1052-NN). Strings in the P-Files are depth labels and implant side. MAT header
file-creation dates (June 2024, processing dates) were set to day 01. No names or recording dates were found.
