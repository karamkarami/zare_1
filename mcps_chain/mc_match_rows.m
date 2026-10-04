function [im, ir, how] = mc_match_rows(rows, refInfo, nRef)
%MC_MATCH_ROWS Pair my rows with the rows of a log lane.
%
%   [im, ir, how] = mc_match_rows(rows, refInfo, nRef)
%
%   rows    : my row bookkeeping (mc_rows / output of a block)
%   refInfo : info of the log lane, nRef its number of rows
%   im, ir  : my rows im(i) and log rows ir(i) belong to the same pulse
%   how     : 'pulseSeq' (same pulse number) or 'order' (no pulseSeq: row by row)
%
%   My rows with an incomplete history (rows.valid = false) are left out.

refSeq = mc_info(refInfo, 'pulseSeq');
if rows.hasSeq && numel(refSeq) == nRef
    [found, loc] = ismember(rows.seq, double(refSeq));
    im  = find(found & rows.valid);
    ir  = loc(im);
    how = 'pulseSeq';
else
    n   = min(numel(rows.seq), nRef);
    im  = find(rows.valid(1:n));
    ir  = im;
    how = 'order';
end
end
