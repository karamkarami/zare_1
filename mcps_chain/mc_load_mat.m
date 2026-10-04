function [video, videoInfo, ref] = mc_load_mat(file)
%MC_LOAD_MAT Lanes of one scan from the .mat file written by log.export.
%
%   [video, videoInfo, ref] = mc_load_mat(file)      file '' = choose in a dialog
%
%   ref.decoder, ref.canceler, ref.integral, ref.cfar and their ...Info
%   ([] when the file does not hold the lane).

if nargin < 1 || isempty(file)
    [f, p] = uigetfile('*.mat', 'File written by log.export');
    if isequal(f, 0)
        error('mc_load_mat:file', 'No file chosen.');
    end
    file = fullfile(p, f);
end
S = load(file);
if ~isfield(S, 'video')
    error('mc_load_mat:video', ['"%s" has no variable "video". Variables: %s\n' ...
          'Adapt mc_load_mat.m to these names.'], file, strjoin(fieldnames(S)', ', '));
end
video     = S.video;
videoInfo = field(S, 'videoInfo');
lanes = {'decoder', 'canceler', 'integral', 'cfar'};
for i = 1:numel(lanes)
    ref.(lanes{i})            = field(S, lanes{i});
    ref.([lanes{i} 'Info'])   = field(S, [lanes{i} 'Info']);
end
end

function v = field(S, name)
v = [];
if isfield(S, name)
    v = S.(name);
end
end
