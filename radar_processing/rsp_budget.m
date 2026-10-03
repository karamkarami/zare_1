function B = rsp_budget(P, S, quiet)
%RSP_BUDGET Expected SNR of every target through the chain and its CFAR margin.
%
%   B = rsp_budget(P, S)
%   B = rsp_budget(P, S, true)       no printing
%
%   For each target of the scenario S (rsp_simulate) at the beam peak:
%     .snrInDb     SNR per sample at the receiver (powerDb - noiseDb)
%     .pulse       pulse that covers its range (1 = short, 2 = long ...)
%     .decoderDb   decoder (pulse compression) gain, mismatch loss included
%     .mtiDb       canceler gain for its Doppler against the canceler's total
%                  noise gain (negative near zero Doppler: blind speed)
%     .fftDb       what the Doppler bin adds: the bin filter (canceler * window
%                  * DFT) gain minus mtiDb. The canceler also lowers the noise
%                  inside the bin, so mtiDb + fftDb is the real gain of the bin
%     .snrOutDb    SNR after the integration
%     .bin         Doppler bin where it lands
%     .thresholdDb CFAR threshold factor of that bin
%     .marginDb    snrOutDb against the threshold (> 0: passes the CFAR at the
%                  beam peak; noise makes single looks scatter by a few dB)
%     .widthDeg    azimuth extent where the margin stays > 0 (beam pattern)
%     .foldedMps   velocity the radar will measure (+-lambda*PRF/4)
%
%   The prediction uses the same decoders, canceler taps, FFT window and
%   CFAR design as rsp_chain, so it explains (before any simulation) why a
%   target passes the CFAR or not.

if nargin < 3
    quiet = false;
end
S  = rsp_levels(S);
R  = S.nRange;
g  = rsp_geometry(P, R);
tr = rsp_truth(P, S);

[~, mf] = rsp_matched_filter(zeros(1, R, 'single'), P);
[P, ~]  = rsp_cfar_window(P, mf);
[nInt, nRefEff] = rsp_cfar_looks(P, mf.replica{mf.longPulse});

c    = rsp_canceler_taps(P.canceler);
if ~P.canceler.enable, c = 1; end
w    = rsp_window(P.fft.window, P.fft.nPulses, P.fft.windowParam);
nfft = max(P.fft.nfft, P.fft.nPulses);
m    = (0:P.fft.nPulses-1)';
kb   = (0:nfft-1) - P.fft.shift * floor(nfft/2);

B = struct('target', {}, 'rangeM', {}, 'velocityMps', {}, 'snrInDb', {}, 'pulse', {}, ...
           'decoderDb', {}, 'mtiDb', {}, 'fftDb', {}, 'snrOutDb', {}, 'bin', {}, ...
           'thresholdDb', {}, 'marginDb', {}, 'widthDeg', {}, 'foldedMps', {});
for k = 1:numel(tr)
    t = tr(k);
    if t.cell <= mf.switchCell
        ip = mf.shortPulse;
    else
        ip = mf.longPulse;
    end
    dec  = mf.dec{ip};
    gDec = max(abs(dec.response))^2 / sum(abs(dec.h).^2);       % decoder SNR gain
    fd   = 2 * t.velocityMps / g.lambdaM / P.radar.prfHz;        % cycles per pulse
    % canceler + windowed DFT of bin b form one slow-time filter g_b; the
    % noise in that bin is the white noise through g_b (the canceler also
    % suppresses the noise there), so the SNR gain of the bin is
    % |G_b(fd)|^2 / ||g_b||^2
    gB = zeros(1, nfft);
    for b = 1:nfft
        gb    = conv(w .* exp(-1j*2*pi*kb(b)*m/nfft), c(:));
        gB(b) = abs(sum(gb .* exp(1j*2*pi*fd*(0:numel(gb)-1)')))^2 / sum(abs(gb).^2);
    end
    [gBin, bin] = max(gB);
    gMti = abs(sum(c(:).' .* exp(-1j*2*pi*fd*(0:numel(c)-1))))^2 / sum(c.^2);
    gFft = gBin / gMti;                                          % share of the bin filter
    snrOut = t.snrDb + 10*log10(gDec * gBin);
    alpha  = rsp_cfar_factor(P.cfar.type, nRefEff(bin), nInt(bin), P.cfar.pfa);
    if strcmpi(P.cfar.thresholdMode, 'factor')
        alpha = 10^(P.cfar.factorDb/10);
    end
    margin = 10*log10(1 + 10^(snrOut/10)) - 10*log10(alpha);

    % azimuth extent with a positive margin (two-way pattern)
    az = 0:0.01:3*P.antenna.beamwidthDeg;
    Gw = 20*log10(rsp_antenna_pattern(az, P.antenna));
    ok = 10*log10(1 + 10.^((snrOut + Gw)/10)) - 10*log10(alpha) > 0;
    last = find(~ok, 1) - 1;
    if isempty(last), last = numel(az); end

    B(k).target      = k;
    B(k).rangeM      = t.rangeM;
    B(k).velocityMps = t.velocityMps;
    B(k).snrInDb     = t.snrDb;
    B(k).pulse       = ip;
    B(k).decoderDb   = 10*log10(gDec);
    B(k).mtiDb       = 10*log10(gMti);
    B(k).fftDb       = 10*log10(gFft);
    B(k).snrOutDb    = snrOut;
    B(k).bin         = bin;
    B(k).thresholdDb = 10*log10(alpha);
    B(k).marginDb    = margin;
    B(k).widthDeg    = 2 * az(max(last, 1)) * (last > 0);
    B(k).foldedMps   = t.foldedMps;
end

if ~quiet
    fprintf(['\nLink budget at the beam peak (CFAR %s, Pfa %.0e, %d guard + %d reference ' ...
             'cells, unambiguous +-%.0f m/s)\n'], P.cfar.type, P.cfar.pfa, P.cfar.nGuard, ...
             P.cfar.nRef, g.vUnambMps/2);
    fprintf(['  target  range[km]  vel[m/s]  SNRin  pulse  decoder   MTI    FFT  | SNRout  bin  ' ...
             'thresh  margin  width[deg]  measured vel\n']);
    for k = 1:numel(B)
        b = B(k);
        fprintf('  %6d  %9.1f  %8.1f  %5.1f  %5d  %7.1f %6.1f %6.1f  | %6.1f  %3d  %6.1f  %6.1f  %10.1f  %8.1f\n', ...
            b.target, b.rangeM/1e3, b.velocityMps, b.snrInDb, b.pulse, b.decoderDb, b.mtiDb, ...
            b.fftDb, b.snrOutDb, b.bin, b.thresholdDb, b.marginDb, b.widthDeg, b.foldedMps);
    end
    fprintf('  (dB; margin > 0 = above the CFAR threshold; width = azimuth extent above it)\n\n');
end
end
