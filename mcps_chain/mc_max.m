function [M, bin] = mc_max(D, P)
%MC_MAX Maximum over the Doppler bins after the CFAR.
%
%   [M, bin] = mc_max(D, P)
%
%   D   : rows x range cells x bins (CFAR output)
%   M   : rows x range cells, max over the bins P.output.bins ([] = all)
%   bin : Doppler bin of the max (0 where nothing was detected)

bins = P.output.bins;
if isempty(bins)
    bins = 1:size(D, 3);
end
[M, i] = max(D(:, :, bins), [], 3);
bin = reshape(bins(i), size(i));
bin(M == 0) = 0;
end
