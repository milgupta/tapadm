function [BestStart, quality_strength, quality_index] = FindBestStart(Amps, Derivatives, channel);
% FINDBESTSTART analyzes AG500 amplitude data to find a auspicious sample to start with.
%           When calculating positions, results are sometimes depending on
%           the first sample, since the first result yields the start value
%           for the 2nd sample and so on. Former strategies were to start
%           with the first or last sample of an trial and traverse the data
%           to the other end. This function tries to estimate the 'quality'
%           of all samples with respect to the first order derivatives of the amplitudes.
%           It returns the number of the sample, which is expected to be most 
%           easy to compute. 
%
%           BestStart = FindBestStart(Amps, Derivatives);
%
%           The additional return values quality_strength and quality_index
%           are providing the complete rating of all samples.
%
%           [BestStart, quality_strength, quality_index] = FindBestStart(...
%
%           (needs the Statistics Toolboox)
%
%   SEE ALSO DERIVATIVE

%---------------------------------------------------------------------
% Copyright © 2005-2007 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
% 
% This program is free software; you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your option) 
% any later version.
%
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with 
% this program; if not, write to the Free Software Foundation, Inc., 
% 51 Franklin Street, Fifth Floor, Boston, MA 
%---------------------------------------------------------------------

PERCENTILE_THRESHOLD = 60;					% How many percent of the samples must be smaller, to be rated as 'relatively big'?
NumPoints = size(Amps, 1);

% Find points where both amplitudes and derivatives of most transmitters are relatively big.
% We expect those positions to be uncomplicated points for position-calculation to begin with.

best_amplitude_bound = prctile(abs(Amps(:, :, channel)), PERCENTILE_THRESHOLD);	% For the 6 transmitters, get the amplitude values, where x% of the samples are smaller 
best_derivative_bound = prctile(abs(Derivatives(:, :, channel)), PERCENTILE_THRESHOLD); % same for the 1st order derivatives of the amplitudes

% Build a (NumPoints x 6) quality-matrix with ones where both amplitude and derivative are above the %-bound.
% Multiply with it's line-sum as weight-parameter, favor points where many transmitters exceeds the bounds
quality_matrix = (abs(Derivatives(:, :, channel)) > repmat(best_derivative_bound, NumPoints, 1)) & (abs(Amps(:, :, channel)) > repmat(best_amplitude_bound, NumPoints, 1));
quality_matrix = quality_matrix .* repmat(sum(quality_matrix, 2), 1, 6);

% Now add up- and down-shifted matrices to favor points where predecessor/successor also exceeds the 75%-bound. 
% Build the line-sum of the resulting matrix to get a quality-vector which honors 
% steep strong amplitudes for many transmitters which are clustered in time - that is what we expect to be 
% a 'good' starting point for the position calculation of the trial.
quality_matrix_down_shift = [quality_matrix(2:end, :); zeros(1, 6)]; 
quality_matrix_up_shift = [zeros(1, 6); quality_matrix(1:(end-1), :)]; 
quality_vector = sum(quality_matrix_down_shift+quality_matrix+quality_matrix_up_shift, 2);

% additionally favor (a little) points in the center, the factor 1.3 is empirical 
med = median(1:NumPoints); quality_vector = quality_vector + 1.3/med*(med-abs((1:NumPoints)-med))';

% finally sort the quality-vector in ascending order, so last one is best one
[quality_strength, quality_index] = sort(quality_vector); 
BestStart = quality_index(end);
