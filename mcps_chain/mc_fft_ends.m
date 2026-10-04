function endPulses = mc_fft_ends(rows, refInfo, nRef, P)
%MC_FFT_ENDS Pulses on which the Doppler FFTs end, taken from the log "integral" rows.
%
%   endPulses = mc_fft_ends(rows, refInfo, nRef, P)
%
%   rows    : row bookkeeping of the FFT input (canceler output)
%   refInfo : info of the log integral lane, nRef its number of rows
%   endPulses : rows of the FFT input whose pulse number equals the pulse number
%               of a log integral row (+ P.doppler.logPulseShift), so both chains
%               compute the same FFTs. [] (= every P.doppler.hop pulses) when
%               P.doppler.alignToLog is false or there is no pulseSeq.

endPulses = [];
refSeq = mc_info(refInfo, 'pulseSeq');
if ~P.doppler.alignToLog || ~rows.hasSeq || isempty(refSeq) || numel(refSeq) ~= nRef
    return
end
[found, loc] = ismember(double(refSeq) + P.doppler.logPulseShift, rows.seq);
endPulses = loc(found);
end
