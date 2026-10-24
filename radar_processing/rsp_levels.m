function S = rsp_levels(S)
%RSP_LEVELS Scenario levels in dB (accepts the older linear / SNR field names).
%
%   S = rsp_levels(S)
%
%   Fills S.noiseDb, S.targets(k).powerDb and S.clutter.powerDb from
%   S.noisePower (linear), S.targets(k).snrDb and S.clutter.cnrDb when the
%   dB fields are missing.
if ~isfield(S, 'noiseDb') || isempty(S.noiseDb)
    if isfield(S, 'noisePower') && ~isempty(S.noisePower)
        S.noiseDb = 10*log10(S.noisePower);
    else
        S.noiseDb = 0;
    end
end
for k = 1:numel(S.targets)
    if ~isfield(S.targets, 'powerDb') || isempty(S.targets(k).powerDb)
        S.targets(k).powerDb = S.targets(k).snrDb + S.noiseDb;
    end
end
if isfield(S, 'clutter') && isfield(S.clutter, 'cnrDb')
    for k = 1:numel(S.clutter)
        if (~isfield(S.clutter, 'powerDb') || isempty(S.clutter(k).powerDb)) && ...
                ~isempty(S.clutter(k).cnrDb)
            S.clutter(k).powerDb = S.clutter(k).cnrDb + S.noiseDb;
        end
    end
end
end
