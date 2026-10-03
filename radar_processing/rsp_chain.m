function out = rsp_chain(video, P)
%RSP_CHAIN Full processing chain, returns the output of every block.
%
%   out = rsp_chain(video, P)
%
%   video : pulses x range samples, complex baseband (e.g. 350 x 5469 single)
%   P     : parameters, see rsp_default_params
%
%   out.decoder   pulses x range               matched filter (pulse compression)
%   out.canceler  pulses x range               3-pulse canceler
%   out.doppler   frames x range x nfft        Doppler FFT (complex)
%   out.integral  frames x range x nfft        non-coherent integration
%   out.cfar      struct (det, threshold, map, factor, list), see rsp_cfar
%   out.mf        matched-filter info (replicas, delays, switch cell, each pulse)
%   out.idx       pulse index of each output row, per block, to align the
%                 rows with the log lanes (e.g. videoInfo.pulseSeq(idx))
%   out.nInt      effective looks per integrated cell, per bin (CFAR design)
%   out.nRef      effective independent cells per reference window
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

t = tic;
[out.nInt, out.nRef] = rsp_cfar_looks(P, out.mf.replica{out.mf.longPulse});
out.cfar = rsp_cfar(out.integral, P.cfar, P.nci.law, out.nInt, out.nRef);
out.time.cfar = toc(t);
out.idx.cfar  = out.idx.integral;

out.P = P;
end
