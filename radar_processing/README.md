# Radar signal processing chain (MATLAB)

```
video ─► decoder ─► 3-pulse canceler ─► Doppler FFT ─► non-coherent ─► SO-CFAR ─► max over ─► plots / PPI
        (code/LFM,                                      integration               Doppler bins
         matched or your taps)
```

Each block is a separate function, and every value is a parameter. Each block's output is
returned so it can be compared with the matching MCPS log lane (`decoder`, `canceler`,
`integral`, `cfar`). No toolbox is needed. The code runs in MATLAB R2016b+ and GNU Octave 7+.

## Quick start

1. **Put your radar's values in `radar_params.m`:** fs, the codes and decoder taps of the
   two pulses, PRF, rpm and antenna. Every field is documented in `rsp_default_params.m`.
2. **Your data:** run `main_mcps`, then `main_rsp`. This runs the chain and compares every
   block with the log. It then prints the plots (m / deg / m/s) and draws the block and a PPI sector.
3. **Simulated rotating radar:** run `main_scan`. It first prints a link budget (the expected
   margin of every target over the CFAR threshold). It then simulates two revolutions while
   the PPI sweep turns at the antenna speed. Finally it prints every target echo against its
   plot and replays the revolutions in real time.
4. **Self-test:** run `rsp_selftest` (45 checks).

### Test scenario (`main_scan`)

| Target | Range | Azimuth | Velocity | Power |
|---|---|---|---|---|
| 1 | 10 km | 40° | 20 m/s | 10 dB |
| 2 | 50 km | 20° | 60 m/s | 5 dB |
| 3 | 80 km | 80° | 100 m/s | 0 dB |

* **Levels:** noise 0 dB, clutter −50 dB up to 25 km. Every level is a per-sample power at
  the receiver on one dB scale (`S.noiseDb`, `targets.powerDb`, `clutter.powerDb`).
* **Waveform:** pulse 1 is Barker 13 with a 52-tap mismatched decoder; pulse 2 is an LFM of
  100 µs and 1 MHz with Taylor weighting.
* **Results:** all three targets pass the CFAR by 27–32 dB, as `rsp_budget` predicts. Errors
  are about 0.01° in azimuth, about 1 m in range and below 0.5 m/s in velocity. No other plots
  were found.
* **Target 1** is the hardest echo. At 10 km it lies in the long pulse's blind zone
  (0–16.9 km), so only the 13 µs Barker pulse sees it. At 20 m/s it is close to the
  3-pulse MTI notch, but the Doppler filter after the canceler recovers most of that loss.

```matlab
P    = radar_params();
out  = rsp_chain(video, P);              % every block: out.decoder ... out.cfar, out.max
scan = rsp_scan(P, S, 'degrees', 360);   % simulated scan + live PPI
rsp_ppi_replay(scan, 'speed', 1, 'bins', [6 7 8]);   % replay at 6 rpm, chosen bins only
```

Two revolutions of the test scenario: the sweep turns, each target is painted as it passes
and gets its plot marker.

![Rotating PPI](docs/ppi_rotating.gif)

Earlier example output (Octave rendering; MATLAB draws the PPI on black):

| PPI after one simulated revolution | Blocks of one 350-pulse block |
|---|---|
| ![PPI](docs/ppi_full_scan.png) | ![Blocks](docs/chain_blocks.png) |

## Pulses and decoders

Each pulse is either `type = 'code'` or `type = 'lfm'`. The bandwidth `bwMHz = 1` gives
1 µs chips, or a 1 MHz LFM sweep.

| Field | Meaning |
|---|---|
| `code` | your chip values, any length, binary or polyphase (`rsp_code` has Barker, MLS, Frank, P3, P4) |
| `decoder` | your decoder taps; `[]` = matched filter |
| `decoderForm` | `'fir'`: FIR taps applied by convolution (matched = `conj(fliplr(code))`); `'reference'`: correlation reference (matched = `code`) |
| `decoderRate` | `'chip'` (one tap per chip) or `'sample'` (taps already at fs) |
| `decoderLag` | `[]` aligns on the main peak automatically, so mismatched filters longer than the code also land on the right range cell |
| `fcMHz`, `delayUs` | frequency offset and transmit delay of the pulse |

For every pulse the program prints the PSL, the ISL and the loss against a matched filter.
If the loss is above 3 dB it warns you, because then `decoderForm` or `decoderRate` is
probably wrong.
`rsp_code_mmf(code, len)` designs a least-squares mismatched filter. For Barker 13 with
39 taps it gives −38.5 dB PSL for a 0.2 dB loss.

The short and long outputs are stitched: the short pulse covers the blind zone of the long
pulse. The CFAR runs only on cells whose echo is received in full (`P.cfar.validOnly`).

## Radar, antenna and units

* **Radar:** `P.radar.prfHz = 1000`, `rpm = 6` (36°/s, 0.036° per pulse, 10 s per revolution),
  `fcMHz = 600` (λ = 0.5 m).
  With these values a range cell is c/2fs = 25 m at 6 MHz, and the unambiguous velocity is
  ±λ·PRF/4 = ±125 m/s, so 20, 60 and 100 m/s are all measured directly.
  At L band (1300 MHz) the limit would be ±58 m/s, and 100 m/s would read as −15 m/s.
* **Antenna:** `P.antenna`: 4° 3-dB beam, −25 dB sidelobes, −40 dB back lobe. The pattern is
  one-way; echoes use it two-way.
* **Targets (simulation):** given as `rangeM`, `azDeg`, `velocityMps` (positive = approaching)
  and `powerDb`. The echo is the transmitted signal (both pulses) at a fractional delay
  2R(t)/c with carrier phase −2π f₀ τ, so Doppler, its sign and range walk come from the
  geometry. Its amplitude follows the two-way antenna pattern as the antenna turns.
  Clutter is ground clutter at all azimuths, with a Gaussian spectrum and log-normal texture.
  Its mean power equals `powerDb`.
* **Simulator check:** `rsp_selftest` compares the simulator with an independent computation.
  The echo matches the delayed transmission × carrier phase × amplitude exactly. The
  pulse-to-pulse phase equals 2v/λ/PRF. Noise and clutter powers come out at the requested dB.
* **Plots (`rsp_extract_plots`):** range [m], azimuth [deg], radial velocity [m/s] and
  CFAR margin. Azimuth is corrected for the chain delay: an output row comes from the
  pulses before it. Velocity uses the Doppler bin refined between bins by a 3-point fit, so
  it is not limited to the 15.6 m/s bin spacing.

## Output after the CFAR

`rsp_cfar_max` takes the maximum over the Doppler bins. This is the only max in the chain.
`P.output.bins` selects which bins enter it:

* `[]`: all bins
* `5`: only bin 5
* `2:16`: every bin except zero Doppler

The output is in `out.max`:

* `.value`, `.bin`: the max value and its Doppler bin
* `.binFrac`: the refined Doppler bin
* `.marginDb`: how far the detection is above the CFAR threshold

`P.cfar.mapValue` sets what the CFAR output carries:

* `'value'`: the integrated value, like the log lane
* `'snr'`: the value divided by its noise estimate. Bins with different noise floors
  (MTI-shaped) then compare fairly in the max.

## PPI

`rsp_ppi_init` and `rsp_ppi_update` draw the PPI: north up, range rings in km, a rotating
sweep, afterglow (`P.ppi.fadeDb` per revolution), and plot markers labelled with velocity.
`P.ppi.source` selects what is painted:

* `'max'`: the output after CFAR (default)
* `'video'`, `'decoder'` or `'canceler'`: raw views

The sweep turns in `P.ppi.stepDeg` steps (0.5°) at `P.ppi.speed` times the real antenna
speed (6 rpm), painting each sector as it passes. Detections are shown in dB above the CFAR
threshold, plots are marked as they are extracted, and plots from earlier revolutions stay as a
trail (`P.ppi.trailScans`).
`rsp_scan` paints live while it processes; if processing is slower than real time, the sweep
waits for the data. `rsp_ppi_replay(scan, 'speed', 1, 'loops', n)` replays at exactly the
antenna speed.

## Processing blocks (summary)

* **Decoder:** frequency-domain correlation (one FFT of the data, one inverse FFT per pulse),
  with range alignment and normalisation (`'noise'` keeps the same noise floor for both pulses).
* **Canceler:** `x(m) − 2x(m−1) + x(m−2)`. The order or the taps can be changed.
* **Doppler FFT:** `nPulses` pulses per FFT, `nfft` bins, a window, and `hop` pulses between FFTs.
* **Integration:** sum or mean of `n` consecutive FFT outputs per bin, with a square, linear or log law.
* **CFAR:** SO (or CA/GO). The threshold is either `factorDb` or computed from `pfa`. The Pfa
  design uses the exact noise correlation of the chain (`rsp_cfar_looks`). Measured on
  simulated noise: 2–3e-6 for a design of 1e-6. With `nGuard = []` / `nRef = []` the window
  follows the waveform (`rsp_cfar_window`):
  * guard = widest compressed main lobe to its first nulls + 2 cells (14 here)
  * reference = 16 resolution cells (96 here)
* **Integration length:** n = 32 consecutive FFT outputs. Together with the FFT window this
  spans 47 pulses (1.7°), which fits inside the 111-pulse beam dwell. The extra looks lower
  the SO-CFAR threshold from 14 to 8.5 dB compared with n = 4.
* **Link budget (`rsp_budget`):** the expected SNR after every block (decoder gain, the
  canceler + Doppler-bin filter at the target's Doppler) and the margin over the actual CFAR
  threshold. It agrees with the simulated margins within about 1 dB.

## Comparing with your log (`main_rsp`)

* **a)** The full chain from video is compared with every log lane.
* **b)** Each block is fed the **previous log lane**, so a mismatch points at exactly one block.

`rsp_compare` finds the pulse and range shift, then the gain, correlation and NMSE. For CFAR
lanes it reports matched, missed and extra detections.

## Files

| File | Purpose |
|---|---|
| `radar_params.m` | **your radar's values** |
| `rsp_default_params.m` | every parameter with documentation |
| `main_rsp.m`, `main_scan.m` | log data and comparison / simulated rotating radar with PPI |
| `rsp_chain.m` | the whole chain for one block |
| `rsp_waveform.m`, `rsp_decoder.m`, `rsp_code.m`, `rsp_code_mmf.m`, `rsp_window.m` | pulses, codes, decoders |
| `rsp_matched_filter.m` | pulse compression, alignment, stitching |
| `rsp_canceler.m`, `rsp_canceler_taps.m` | MTI canceler |
| `rsp_doppler_fft.m`, `rsp_nci.m` | Doppler FFT, non-coherent integration |
| `rsp_cfar.m`, `rsp_cfar_cells.m`, `rsp_cfar_factor.m`, `rsp_cfar_looks.m` | CFAR and its design |
| `rsp_cfar_max.m` | max over the selected bins after the CFAR (refined Doppler, margin) |
| `rsp_cfar_window.m`, `rsp_budget.m` | CFAR window from the waveform, link budget |
| `rsp_extract_plots.m` | plot extraction |
| `rsp_simulate.m`, `rsp_truth.m`, `rsp_levels.m`, `rsp_antenna_pattern.m`, `rsp_geometry.m` | simulation and units |
| `rsp_scan.m` | block-by-block scan processing (output independent of block size) |
| `rsp_ppi_init.m`, `rsp_ppi_update.m`, `rsp_ppi_replay.m` | PPI display |
| `rsp_plot.m`, `rsp_compare.m`, `rsp_db_limits.m` | figures and comparison |
| `rsp_selftest.m` | tests |

## Notes

* **Frequency diversity and channel filters:** in the simulation the short and long pulses
  sit at −1.5 / +1.5 MHz. A target in the blind zone, like target 1, has a strong long echo
  partly received after the receiver turns on. Without channel separation that echo leaks into
  the short decoder and makes a ghost plot about 3 km closer to the radar. `P.mf.channelFilter`
  puts a band-pass FIR in front of each decoder, as a dual-frequency receiver does.
  It lowers the leakage from −23 to −39 dB, and the ghost disappears. Set `channelFilter = false`
  to match a plain decoder (e.g. your log), and `fcMHz = 0` if your radar uses one frequency.
* **Velocity folding:** beyond ±λ·PRF/4 the velocity folds. Plots report the measured value;
  `rsp_truth` gives the folded truth for comparison.
* **Zero-Doppler clutter residue:** it can still pass the CFAR (the Pfa design assumes thermal
  noise). Use `P.output.bins = 2:16` to leave bin 1 out of the output.
