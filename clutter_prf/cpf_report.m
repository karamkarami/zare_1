function cpf_report(sel)
%CPF_REPORT Print a PRF set and its criteria.

fprintf('PRF set (design CNR %.1f dB, clutter sigma %.2f Hz, platform %.1f m/s)\n', ...
        sel.design.cnrDb, sel.design.sigmaHz, getf(sel.design, 'speedMps'));
fprintf('  PRF [Hz]  blind speed lambda*PRF/2 [m/s]  CA [dB]\n');
for i = 1:numel(sel.prfHz)
    fprintf('  %7.1f   %8.1f                       %6.1f\n', sel.prfHz(i), ...
            sel.prfHz(i) * 299792458 / (2 * sel.design.fcHz), sel.caDb(i));
end
ok = {'NOT met', 'met'};
fprintf('  worst-case CA       %6.1f dB (needed %.1f dB)\n', sel.caWorstDb, sel.caReqDb);
fprintf('  ripple of mean resp.%6.2f dB (lowest mean gain at %g m/s)\n', sel.rippleDb, sel.vWorstMean);
fprintf('  blind PRFs at once  %d (limit %d, %s; over the limit on %.1f%% of the band)\n', ...
        sel.blindMax, sel.maxBlind, ok{(sel.blindMax <= sel.maxBlind) + 1}, sel.overPct);
fprintf('  cumulative Pd over the set: min %.3f (at %g m/s), mean %.3f\n', sel.pcMin, sel.vWorstPd, sel.pcMean);
end

function v = getf(s, f)
if isfield(s, f), v = s.(f); else, v = 0; end
end
