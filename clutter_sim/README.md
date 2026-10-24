# Test-data simulator for the clutter-map MTI

The simulator is independent of `clutter_prf`: neither one calls the other. It writes
complex baseband after pulse compression, block by block, in the same struct the processor
reads from online data (`x`, `prfHz`, `azDeg`, `rangeM`, `velEN`, `timeS`).

```matlab
S   = csim_config();                 % radar, clutter, platform, targets
sim = csim_init(S);
[blk, sim, truth] = csim_block(sim, prfHz);   % next block; truth.scanDone at north
sim.vel = [10 5];                     % platform starts to move (m/s east, north)
```

* **Echo model:** every echo is amplitude × two-way Gaussian beam × exp(−j4πR(t)/λ). The
  Doppler of the clutter and the targets therefore comes from the geometry, both their own
  motion and the platform's.
* **Ground clutter:** one scatterer per range cell every `clutter.azStepDeg`, fixed on the
  ground.
  * **Power:** mean CNR `clutter.cnrDb` at `refRangeM`, falling as R^−`rangeExp`, with a
    log-normal texture and stronger patches (`clutter.patches`).
  * **Motion:** each scatterer has an intrinsic motion, a sum of `components` lines with a
    Gaussian Doppler spread of σv = `clutter.sigmaVMps`.
* **Targets:** each target has a start range and azimuth, a radial and a tangential speed
  over ground, and an SNR. Set `S.swerling = true` for Swerling 1 (amplitude drawn again
  every revolution), or `false` for constant amplitude.
* **Scans:** a scan always ends at north, so the PRF can change between scans.
