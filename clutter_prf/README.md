# Clutter-map MTI with a scan-to-scan PRF set (MATLAB)

```
baseband after pulse compression ─► platform Doppler removed ─► per-cell MTI filter ─► Doppler FFT ...
  (online or simulated, same struct)    (moving radar, TACCAR)      (from the clutter map)
                                              │
                                              └─► clutter map (learned over scans, in m/s) ─► PRF set (P PRFs)
```

* `clutter_prf/`: the processor. It does not know whether the data is online or simulated.
* `clutter_sim/`: an independent simulator that writes test data in the same format.
* `main_clutter_prf.m` (repo root): runs both together.
* `cpf_selftest`: 13 checks on synthetic data.

No toolbox is needed. The code runs in MATLAB R2016b+ and GNU Octave 7+.

## Input: one block of pulses

The data is complex baseband (I + jQ) after the decoder / pulse compression. The carrier
frequency is used only to convert Hz to m/s. A block may hold any number of pulses (1, one
CPI, or a whole scan), but all pulses of a block must use the same PRF.

| Field | Meaning |
|---|---|
| `blk.x` | pulses × range, complex |
| `blk.prfHz` | PRF of these pulses |
| `blk.azDeg` | antenna azimuth of every pulse, 0 = north, clockwise |
| `blk.rangeM` | range of every range cell [m] |
| `blk.velEN` | platform velocity [east north] in m/s, 1×2 or pulses×2 (optional, default stationary) |
| `blk.posEN` | platform position [m] (optional; integrated from `velEN` when missing) |

```matlab
cfg = cpf_config();                 % edit: fc, beamwidth, PRF range, P, max length, velocity band
st  = cpf_init(cfg);                % empty map + first PRF set
for each scan
    [prf, st, sel] = cpf_next_prf(st, speed);   % PRF of this scan (set re-selected now and then)
    for each block of the scan
        [y, st] = cpf_process(st, blk);         % y -> your Doppler FFT / CFAR
    end
end
```

From the MCPS log, for example:
`blk.x = decoder; blk.azDeg = decoderInfo.azimuthDeg; blk.prfHz = 1e6/median(diff(decoderInfo.timestampUs)); blk.rangeM = (decoderInfo.rangeCell - 1)*cellM;`

## Settings you give (`cpf_config`)

| Field | Default | Meaning |
|---|---|---|
| `radar.fcHz` | 600e6 | carrier frequency (λ) |
| `radar.beamwidthDeg`, `radar.rpm` | 4, 6 | for the scan and platform spread of the clutter |
| `prf.minHz`, `prf.maxHz`, `prf.stepHz` | 800, 1500, 5 | candidate PRFs |
| `prf.count` | 4 | P |
| `mti.maxLength` | 8 | maximum filter length |
| `target.vMinMps`, `target.vMaxMps` | 20, 350 | radial speed band over ground, both signs |

## How the PRF set is chosen (`cpf_select_prfs`)

For every candidate PRF the filter is designed for the design clutter. The set of P PRFs is
then chosen by these criteria:

| Criterion | How it is measured | Default use |
|---|---|---|
| Worst-case clutter attenuation | CA of the worst PRF of the set; each dB below the needed CA costs (`prf.wCa`). The needed CA is the design CNR: clutter brought down to the noise. More CA is not rewarded. | cost |
| Ripple of the mean response | max/min over ±[vMin vMax] of the response averaged over the P PRFs (`prf.wRipple`) | cost |
| Blind PRFs at once | at any velocity at most `prf.maxBlind` of the P PRFs may be blind (default 1). Blind means the folded Doppler lands in the clutter region: clutter edge + half a Doppler-FFT bin (`prf.blindMode='clutter'`). Alternatively, the filter gain is below `prf.blindDb` (`'notch'`). | limit; if it cannot be met, the smallest over-limit part of the band wins |
| Cumulative Pd over the P scans | Swerling 1 scan to scan, SNR `prf.snrDb` × filter gain | reported (`prf.wPd > 0` adds it to the cost) |

The search tries every subset when there are at most `prf.maxExhaustive` subsets. Otherwise it
runs a local search from several starts, exchanging one PRF at a time. `cpf_report(sel)` prints
the set, and `cpf_prf_eval(cfg, design, prfs, lambda)` scores a set of your own.

## The filter (`cpf_design_filter`)

* **Clutter:** the taps stay in the span of the clutter eigenvectors whose power, CNR × λₖ,
  is below `mti.residueDb` (an eigencanceler). Any filter in that span leaves at most that
  much clutter.
* **Passband:** the passband is every target velocity whose folded Doppler is outside the
  clutter edge. Inside the span, the taps maximise the mean gain in dB over the passband. This
  keeps the notch narrow and the passband free of holes.
* **Why not the plain mean gain:** maximising the plain (linear) mean gain gives the widest
  notch, because for unit-norm taps the mean over all Doppler is always 1.
* **Normalisation:** every filter has unit noise gain.
* **Length:** `mti.lengthMode = 'auto'` uses the shortest length that meets the residue.

Each PRF has a bank of filters over clutter width / PRF × CNR. A filter is designed the first
time a cell needs it. Every sample takes the filter of its map cell, rounded up to the next grid
value. The filter is shifted to the cell's mean
velocity. A cell with CNR below `mti.bypassCnrDb` passes unfiltered.

## The clutter map

* **Where it lives:** the map is fixed to the ground, as a polar grid around a map origin
  (`map.azCellDeg` = half the beamwidth, `map.rangeCellM` = a few range cells).
* **What it holds:** each cell holds the clutter power, the intrinsic spread σv and the mean
  velocity. σv and the mean velocity are in m/s, so a cell learned at one PRF serves every PRF.
  Changing the PRF set needs no new learning.
* **Spread estimate:** σ comes from a line fit of log|R(k)| against k², over lags
  1…`map.maxLag`. Narrow clutter hardly decorrelates in one pulse, so lag 1 alone is far too
  noisy. The known antenna-scan spread (0.265·ω/θ) and platform spread are subtracted.
* **Updates:** a cell is folded in once the beam has left it. The update is recursive with
  weight `map.alpha`. Power is limited to ±`map.clipDb` per update, so a passing target does
  not wipe out what was learned.
* **Moving radar:** when the radar moves more than `map.recenterM` from the origin, the map
  is re-centred and keeps its values.

## Moving radar (ground vehicle)

* **Clutter Doppler:** fixed ground appears at 2·v·u/λ, where u is the beam direction. This
  changes with azimuth as the antenna turns. `cpf_process` removes it pulse by pulse
  (TACCAR), so the clutter is back at 0 Hz.
* **Same filters:** the weights stay the same family. The map, which was learned over ground,
  stays valid.
* **Wider spectrum:** the clutter is wider by 2/λ · 0.300·θ · v⊥, largest broadside. Every
  cell's filter includes this, and the PRF design uses it with `update.motionCase`.
* **Output frame:** by default the output keeps the platform Doppler removed, so the FFT
  shows velocities over ground. Set `mti.output='radar'` to put it back.

## Periodic update (`cpf_next_prf`)

* **Map:** updated every scan by `cpf_process`.
* **PRF set:** selected again every `update.everyScans` scans, or when the platform speed
  changed by `update.speedDeltaMps`. The design clutter comes from the map: the 95th
  percentile of CNR and the 90th percentile of σv, plus the scan and platform spread.
* **Hysteresis:** the new set replaces the old one only when its cost is lower by
  `update.hysteresis`, or when the old set no longer meets the limits.

## Demo result (`main_clutter_prf`, 600 MHz, 4°, 6 rpm, PRF 800–1500 Hz, P = 4, N ≤ 8)

* **First PRF set:** 815 / 930 / 1080 / 1490 Hz.
  * Worst CA is 75.6 dB against the 66 dB needed.
  * At most one PRF is blind at any velocity.
  * The ripple of the mean response is 11.2 dB, and the lowest mean gain is at ±20 m/s.
* **Stationary radar:** clutter goes from 48 dB over noise to 0.0 dB after the filter.
  * The map is learned after 2 scans; median intrinsic σv is 0.26 m/s (simulated: 0.3).
  * Every target's power after the filter matches the gain expected for its velocity
    (about ±1 dB away from the notch).
  * A target near the blind speed of one PRF (−220 m/s at 815 / 930 Hz) is seen in the
    other scans.
* **Moving at 15 m/s (from scan 9):**
  * The PRF set is selected again (805 / 930 / 1130 / 1390 Hz) because the speed changed.
  * The residue stays at 0.0 dB with motion compensation. Without it the residue is 30.6 dB,
    and then 19, 14 and 5 dB as the map learns the shifted clutter.
* **Limit of these numbers:** with N = 8 the notch is about 0.1·PRF wide, so 20 m/s sits on
  its edge. vMin or a longer filter moves this; with N = 12 the −10 dB notch is also met at
  1 blind PRF.
