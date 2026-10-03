# Radar signal processing chain (MATLAB)

```
video ─► matched filter ─► 3-pulse canceler ─► Doppler FFT ─► non-coherent integration ─► SO-CFAR
         (decoder)          (canceler)                          (integral)                   (cfar)
```

Every block is a separate function. Every value is a parameter in `rsp_default_params.m`.
The output of every block is returned, so it can be compared with the matching lane of
the MCPS log (`decoder`, `canceler`, `integral`, `cfar`). No toolbox is needed. The code runs in
MATLAB R2016b+ and in GNU Octave 7+.

## Quick start

```matlab
% 1. load the log (your script)
main_mcps                 % -> video, decoder, canceler, integral, cfar

% 2. set your values in section 2 of main_rsp.m (fs, p1, p2, CFAR ...), then
main_rsp                  % runs the chain, compares every block, plots

% without the log: set  source = 'simulate'  in main_rsp.m
rsp_selftest              % checks every block against a direct computation
```

Programmatic use:

```matlab
P   = rsp_default_params();
P.fs = 4;  P.pulse(1).widthUs = 2;  P.pulse(2).widthUs = 20;   % fs [MHz], p1, p2 [us]
out = rsp_chain(video, P);         % out.decoder, out.canceler, out.doppler, out.integral, out.cfar
rsp_compare(decoder, out.decoder, 'decoder', 'plot', true);
rsp_plot(out, video);
```

## Files

| File | Block |
|---|---|
| `rsp_default_params.m` | all parameters, with comments |
| `rsp_chain.m` | runs the whole chain and returns every block output, its row indices and timing |
| `rsp_matched_filter.m` | pulse compression of the two pulses, range alignment, stitching |
| `rsp_canceler.m`, `rsp_canceler_taps.m` | N-pulse MTI canceler (default 3-pulse `[1 -2 1]`) |
| `rsp_doppler_fft.m` | slow-time FFT, sliding (`hop = 1`) or block CPI |
| `rsp_nci.m` | non-coherent integration of each bin over `n` consecutive FFT outputs |
| `rsp_cfar.m` | SO / CA / GO range CFAR, all cells at once (cumulative sums) |
| `rsp_cfar_factor.m` | threshold factor for a target Pfa (exact Gamma model) |
| `rsp_cfar_looks.m` | effective looks / reference cells of the real noise, for the Pfa design |
| `rsp_waveform.m`, `rsp_window.m` | pulse replicas (LFM / unmodulated / phase code / user samples), windows |
| `rsp_compare.m` | finds shift and gain, then gives correlation and NMSE (or hits / misses for detections) |
| `rsp_plot.m`, `rsp_db_limits.m` | figures |
| `rsp_simulate.m` | synthetic video: targets, clutter, noise, eclipsing |
| `rsp_selftest.m` | automatic tests |
| `main_rsp.m` | main script |

## Blocks

**Matched filter.** Each pulse has its own width (`p1`, `p2`) and LFM bandwidth (0 = unmodulated).
It also has a centre frequency, chirp slope, window and transmit delay. A user replica (`samples`)
or a phase code (`code`) can replace the generated waveform. Each pulse is correlated in the
frequency domain: one forward FFT of the data, then one inverse FFT per pulse, with no circular
wrap-around. All pulses are aligned to the same range grid: output cell `r` holds a scatterer
whose echo of pulse k starts at sample `r + delay_k`. With `P.mf.combine = 'stitch'`, cells
`1..switchCell` come from the short pulse and the rest from the long pulse. By default
`switchCell` is the blind zone of the long pulse, `max(delay + length) - delay_long`.
`P.mf.norm = 'noise'` keeps the noise floor equal for both pulses, so the stitch is seamless.

**Canceler.** `y(m) = x(m) - 2x(m-1) + x(m-2)` along slow time. You can change the order or set
your own taps. `output = 'same'` keeps every pulse (the first 2 rows are 0); `'valid'` drops them.

**Doppler FFT.** `nPulses` pulses per FFT, `nfft` bins and a window. Frames are `hop` pulses apart:
`hop = 1` gives one FFT per pulse, `hop = nPulses` gives non-overlapping CPIs. The output is
`frames x range x nfft`, the same layout as the log. The work is split into blocks of frames
(`chunk`) so memory stays bounded.

**Non-coherent integration.** For each range cell and bin, the sum or mean over a buffer of `n`
consecutive FFT outputs. The detector law is `square` (|x|²), `linear` (|x|) or `log` (dB).

**SO-CFAR.** It uses `nRef` reference cells and `nGuard` guard cells on each side.
The noise estimate is `min(mean(lead), mean(lag))`; `CA` and `GO` can also be selected.
The threshold is either

* `thresholdMode = 'factor'`: the factor is `factorDb`; or
* `thresholdMode = 'pfa'`: the factor is computed for `pfa`.

In `'pfa'` mode the design accounts for two effects:
1. Matched-filter noise is correlated in range, so 16 oversampled cells hold only about 7
   independent ones.
2. Overlapping FFT windows and the canceler correlate the frames.

`rsp_cfar_looks` computes both effects from the real filters, using the eigenvalues of the
noise covariance. On simulated noise the measured Pfa matches the design (≈1.1e-6 for 1e-6).
Ignoring effect 1 would give 1e-4.
`edge` sets what happens at the ends of the range: `'oneSided'` uses the complete window,
`'partial'` uses the available cells, `'none'` skips the cell.
The output fields are `det` (logical), `threshold`, `map` (the value on detections, 0 elsewhere,
like the log lane) and `list` (frame, range, bin, value, SNR).

## Comparing with the log

`main_rsp.m` does two comparisons:

* **a) full chain from video.** Every output of our chain is compared with the log lane.
* **b) block by block.** Each block gets the **previous lane of the log** as input, so a
  mismatch points at exactly one block.

`rsp_compare` first finds the pulse and range shift, then the least-squares gain, then reports
the correlation (1 = identical up to a gain) and the NMSE. For CFAR lanes it reports matched,
missed and extra detections, with a tolerance in range cells. `'plot', true` shows the reference,
the test and their difference.

Typical causes of a mismatch:

* a constant range shift → `P.mf.rangeOffset` or the pulse `delayUs`
* a constant gain → `P.mf.norm`, `P.fft.norm`, `P.nci.average`
* the canceler matches but the integral does not → `nfft`, `hop`, window, `fft.shift`, `nci.law`

## Notes

* Rows: `out.idx.<block>` gives, for every output row, the index of the newest input pulse.
  Use it to align with `videoInfo.pulseSeq` or to pick log rows.
* In the last `length_long + delay_long` cells the long pulse is only partly received, so the
  noise falls there. SO-CFAR then raises false alarms at that edge; limit the tested range with
  `P.cfar.rangeCells`.
* The Pfa design assumes thermal noise. Strong clutter residue near zero Doppler stays correlated
  from frame to frame, so expect extra alarms there. You can exclude bins (`P.cfar.bins`) or set
  `P.cfar.nIntEffective = 1`.
