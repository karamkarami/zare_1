function v = mc_info(info, name)
%MC_INFO One field of a log info (struct, struct array or table) as a column, [] if missing.
%
%   seq = mc_info(videoInfo, 'pulseSeq')

v = [];
if isempty(info)
    return
end
try
    if isstruct(info) && numel(info) > 1
        v = [info.(name)];
    else
        v = info.(name);
    end
    v = v(:);
catch
    v = [];
end
end
