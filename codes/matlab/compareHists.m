%
% s = compareHists(h1,h2)
%       returns a histogram similarity in the range 0..1
%
% Compares 2 normalised histograms using the Bhattacharyya coefficient.
% Assumes that sum(h1) == sum(h2) == 1
%
function s = compareHists(h1,h2)

s = sum(sum(sqrt(h1).*sqrt(h2)));

