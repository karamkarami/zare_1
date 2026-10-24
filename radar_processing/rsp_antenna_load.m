function A = rsp_antenna_load(A, file)
%RSP_ANTENNA_LOAD Use a measured / designed azimuth pattern (table or file).
%
%   A = rsp_antenna_load(A, file)    read the pattern from a file
%   A = rsp_antenna_load(A)          use A.file, or the table already in
%                                    A.patternDeg / A.patternDb
%
%   A    : P.antenna
%   file : .csv / .txt  two columns: angle from the beam axis [deg], one-way
%                       power gain [dB]; a header line is allowed
%          .mat         variables angleDeg and gainDb (or the first two
%                       numeric vectors in the file)
%          Angles may be -180..180, 0..360 (converted) or 0..180 (symmetric).
%          A voltage pattern can be given in dB as well: 20*log10(|E|) is the
%          same power gain.
%
%   A.type = 'table', A.patternDeg / A.patternDb sorted and normalised to a
%   0 dB peak, and A.beamwidthDeg is measured from the table (3 dB width),
%   so everything that uses the beamwidth (plot extraction, clutter, budget)
%   follows the loaded pattern. Without a file or table, A is unchanged.

if nargin >= 2 && ~isempty(file)
    A.file = file;
end
if isfield(A, 'file') && ~isempty(A.file)
    [a, d] = readPattern(A.file);
elseif isfield(A, 'patternDeg') && ~isempty(A.patternDeg)
    a = A.patternDeg;
    d = A.patternDb;
else
    return
end
a = double(a(:));
d = double(d(:));
if numel(a) ~= numel(d) || numel(a) < 3
    error('rsp_antenna_load:table', 'The pattern needs at least 3 (angle, gain) pairs.');
end
if max(a) > 180                                        % 0..360 -> -180..180
    a = mod(a + 180, 360) - 180;
end
[a, i] = unique(a);
d = d(i) - max(d);                                     % 0 dB peak
A.type       = 'table';
A.patternDeg = a;
A.patternDb  = d;
A.beamwidthDeg = beamwidth(a, d);
end

% ======================================================================================
function [a, d] = readPattern(file)
if ~exist(file, 'file')
    error('rsp_antenna_load:file', 'Pattern file "%s" not found.', file);
end
[~, ~, ext] = fileparts(file);
if strcmpi(ext, '.mat')
    s = load(file);
    if isfield(s, 'angleDeg') && isfield(s, 'gainDb')
        a = s.angleDeg;  d = s.gainDb;
    else
        f = fieldnames(s);
        v = {};
        for i = 1:numel(f)
            if isnumeric(s.(f{i})) && isvector(s.(f{i}))
                v{end+1} = s.(f{i});                   %#ok<AGROW>
            end
        end
        if numel(v) < 2
            error('rsp_antenna_load:mat', '%s: need angleDeg and gainDb.', file);
        end
        a = v{1};  d = v{2};
    end
else
    fid = fopen(file, 'r');
    txt = fread(fid, Inf, '*char')';
    fclose(fid);
    lines = regexp(txt, '[\r\n]+', 'split');
    m = zeros(0, 2);
    for i = 1:numel(lines)
        x = sscanf(strrep(lines{i}, ',', ' '), '%f')';
        if numel(x) >= 2
            m(end+1, :) = x(1:2);                      %#ok<AGROW>
        end
    end
    a = m(:, 1);  d = m(:, 2);
end
end

function bw = beamwidth(a, d)
% 3 dB width around the peak, interpolated between table points.
[~, k] = max(d);
if min(a) >= 0                                         % one-sided table
    j = k + find(d(k:end) < -3, 1) - 1;
    bw = 2 * cross(a, d, j, -1);
else
    j1 = k + find(d(k:end) < -3, 1) - 1;
    j0 = k - find(d(k:-1:1) < -3, 1) + 1;
    bw = cross(a, d, j1, -1) - cross(a, d, j0, 1);
end
end

function x = cross(a, d, j, s)
% angle where the gain crosses -3 dB between table point j (below -3 dB) and
% its neighbour j + s towards the peak
if isempty(j)
    error('rsp_antenna_load:bw', 'The pattern never falls 3 dB below its peak.');
end
i = j + s;
x = a(i) + (-3 - d(i)) * (a(j) - a(i)) / (d(j) - d(i));
end
