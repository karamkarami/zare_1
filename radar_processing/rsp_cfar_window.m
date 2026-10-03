function [P, w] = rsp_cfar_window(P, mf)
%RSP_CFAR_WINDOW Choose the CFAR guard and reference cells from the decoders.
%
%   [P, w] = rsp_cfar_window(P, mf)
%
%   P  : parameters; P.cfar.nGuard / P.cfar.nRef left empty are filled in
%   mf : decoder info (rsp_matched_filter second output)
%   w  : .nGuard .nRef .mainlobe (samples) .resolution (samples per
%        resolution cell) .auto (which values were chosen here)
%
%   Guard cells: the widest compressed main lobe of the pulses (out to its
%   first nulls, so no main-lobe energy leaks into the noise estimate) plus
%   2 cells for range straddle and range walk during the integration.
%   Reference cells: 16 resolution cells (1/bandwidth each) per side. The
%   decoder output is oversampled (fs / bandwidth samples per resolution
%   cell), so fewer cells would hold too few independent noise samples and
%   the SO-CFAR loss would grow quickly; many more would span clutter
%   edges and neighbouring targets.

lobe = 0;
res  = 1;
for k = 1:numel(mf.dec)
    lobe = max(lobe, max(abs(mf.dec{k}.mainlobe)));
    res  = max(res, mf.dec{k}.spc);
end
w.mainlobe   = lobe;
w.resolution = res;
w.auto       = {};
if isempty(P.cfar.nGuard)
    P.cfar.nGuard = ceil(lobe) + 2;
    w.auto{end+1} = 'nGuard';
end
if isempty(P.cfar.nRef)
    P.cfar.nRef = max(16, ceil(16 * res));
    w.auto{end+1} = 'nRef';
end
w.nGuard = P.cfar.nGuard;
w.nRef   = P.cfar.nRef;
end
