function out = rsp_chain(video, P, cache)
%RSP_CHAIN Full processing chain, returns the output of every block.
%
%   out = rsp_chain(video, P)
%   out = rsp_chain(video, P, out.cache)     reuse the CFAR design (block processing)
%
%   video : pulses x range samples, complex baseband (e.g. 350 x 5469 single)
%   P     : parameters, see rsp_default_params
%
%   out.decoder   pulses x range               matched filter / decoder
%   out.canceler  pulses x range               3-pulse canceler
%   out.doppler   frames x range x nfft        Doppler FFT (complex)
%   out.integral  frames x range x nfft        non-coherent integration
%   out.cfar      struct (det, threshold, map, factor, list), see rsp_cfar
%                 (only on out.mf.validCells when P.cfar.validOnly)
%   out.max       output stage after the CFAR (rsp_cfar_max over P.output.bins):
%                 .value frames x range, .bin frames x range (uint8)
%   out.mf        decoder info (replicas, alignment, PSL, loss, switch cell ...)
%   out.idx       pulse index of each output row, per block, to align the
%                 rows with the log lanes (e.g. videoInfo.pulseSeq(idx))
%   out.nInt      effective looks per integrated cell, per bin (CFAR design)
%   out.nRef      effective independent cells per reference window, per bin
%   out.cache     CFAR design, pass it back for the next block
%   out.time      processing time of each block [s]
%   out.P         parameters used

if ~isfloat(video)
    video = single(video);
end
M = size(video, 1);

t = tic;
[out.decoder, out.mf] = rsp_matched_filter(video, P);
out.time.decoder = toc(t);
out.idx.decoder  = (1:M)';

t = tic;
[out.canceler, out.idx.canceler] = rsp_canceler(out.decoder, P.canceler, out.idx.decoder);
out.time.canceler = toc(t);

t = tic;
[out.doppler, out.idx.doppler] = rsp_doppler_fft(out.canceler, P.fft, out.idx.canceler);
out.time.doppler = toc(t);

t = tic;
[out.integral, out.idx.integral] = rsp_nci(out.doppler, P.nci, out.idx.doppler);
out.time.integral = toc(t);
if ~P.output.keepDoppler
    out.doppler = [];
end

t = tic;
if nargin < 3 || isempty(cache)
    [cache.nInt, cache.nRef] = rsp_cfar_looks(P, out.mf.replica{out.mf.longPulse});
end
out.nInt  = cache.nInt;
out.nRef  = cache.nRef;
out.cache = cache;
if P.cfar.validOnly
    out.cfar = rsp_cfar_cells(out.integral, P, out.mf.validCells, out.nInt, out.nRef);
else
    out.cfar = rsp_cfar(out.integral, P.cfar, P.nci.law, out.nInt, out.nRef);
end
out.time.cfar = toc(t);
out.idx.cfar  = out.idx.integral;

t = tic;
[out.max.value, out.max.bin] = rsp_cfar_max(out.cfar.map, P.output.bins);
out.time.max = toc(t);
out.idx.max  = out.idx.cfar;

out.P = P;
end
