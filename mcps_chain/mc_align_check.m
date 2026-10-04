function A = mc_align_check(each, rows, ref, refInfo, P, filt)
%MC_ALIGN_CHECK Range delay and gain of each pulse that best fit the log decoder.
%
%   A = mc_align_check(each, rows, ref, refInfo, P, filt)
%
%   each : {short output, long output} over the whole range (mc_decoder)
%   rows : row bookkeeping of the decoder output
%   ref  : log decoder lane, refInfo its info
%   filt : filter info of each pulse (third output of mc_decoder)
%
%   For each pulse, only the cells that the pulse feeds are used (short: before
%   P.decoder.switchCell, long: from it on). Prints the txDelayUs and gain that
%   line my output up with the log; put them in mc_params.m and run again.

[im, ir] = mc_match_rows(rows, refInfo, size(ref, 1));
if numel(im) > P.compare.maxRows
    pick = round(linspace(1, numel(im), P.compare.maxRows));
    im = im(pick);
    ir = ir(pick);
end
nc = min(size(each{1}, 2), size(ref, 2));
pB = mean(abs(double(ref(ir, 1:nc))).^2, 1);
sw = min(P.decoder.switchCell, nc + 1);
regions = {1:sw-1, sw:nc};

fprintf('  Range alignment of each pulse against the log decoder:\n');
A = struct('shift', {}, 'rho', {}, 'gainDb', {}, 'txDelayUs', {}, 'gain', {});
for k = 1:numel(each)
    cells = regions{min(k, 2)};
    if isempty(each{k}) || numel(cells) < 20
        continue
    end
    X  = abs(double(each{k}(im, 1:nc)));
    [s, rho] = mc_range_shift(mean(X.^2, 1), pB, cells, P.compare.maxShift);
    c  = cells(cells + s >= 1 & cells + s <= nc);
    a  = X(:, c + s);
    b  = abs(double(ref(ir, c)));
    g  = sum(a(:) .* b(:)) / sum(a(:).^2);
    A(k).shift     = s;
    A(k).rho       = rho;
    A(k).gainDb    = 20*log10(g);
    A(k).txDelayUs = P.pulse(k).txDelayUs + s / P.fsMHz;
    A(k).gain      = P.pulse(k).gain * g;
    fprintf(['    pulse %d (%s), cells %d..%d: shift %+d cells (profile corr %.3f), ' ...
             'gain %+.2f dB  ->  P.pulse(%d).txDelayUs = %g;  P.pulse(%d).gain = %.6g;\n'], ...
            k, P.pulse(k).name, cells(1), cells(end), s, rho, A(k).gainDb, ...
            k, A(k).txDelayUs, k, A(k).gain);
    if s ~= 0 && ~strcmp(filt(k).lagFrom, 'code')
        fprintf('      (the peak lag %d was %s; a code in P.pulse(%d).code measures it)\n', ...
                filt(k).lag, filt(k).lagFrom, k);
    end
end
end
