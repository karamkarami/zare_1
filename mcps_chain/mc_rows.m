function rows = mc_rows(info, n)
%MC_ROWS Row bookkeeping of a lane, carried through the blocks of the chain.
%
%   rows = mc_rows(info, n)      n = number of pulses (rows) of the lane
%
%   rows.pulse   pulse (row of the input lane) that each output row ends on
%   rows.seq     pulse number of that pulse (info.pulseSeq), or the row index
%                when the info has no pulseSeq
%   rows.hasSeq  true when rows.seq comes from the info
%   rows.valid   false while a row still needs pulses from before the record
%                (first rows of the canceler, FFT and integration)

seq = mc_info(info, 'pulseSeq');
rows.pulse  = (1:n)';
rows.hasSeq = numel(seq) == n;
if rows.hasSeq
    rows.seq = double(seq);
else
    rows.seq = (1:n)';
end
rows.valid = true(n, 1);
end
