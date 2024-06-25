function ReCalibrationFactors = discaleval(DisCalData, ReCalibrationFactors); 

% DISCALEVAL evaluates DisCal calibration data.
%
%           S = discaleval(Filename);
%           S = discaleval((DisCalData, ReCalibrationFactors);
%           
%           see discal.m

%---------------------------------------------------------------------
% Copyright © 2007 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
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
error(nargchk(1, 2, nargin));

if (nargin < 2)
	ReCalibrationFactors = ones(6, 12);
end

if ischar(DisCalData)						% if DisCalData contains a file name, load that file
	load(DisCalData); DisCalData = result; clear result;
end

ChanIdx = DisCalData.ChanIdx; NChannels = length(ChanIdx);
NCombinations = factorial(NChannels); SNC = sqrt(NCombinations);
N = DisCalData.cnt;
% 7 x NumChan x N 
AllReCalibrationFactors = double(DisCalData.ReCalibrationFactors);
OriginalErr = AllReCalibrationFactors(7, :, 1); % first page contains values for original calibration
disp('Estimated spatial accuracy based on calibration:');
for (i_chan = 1:length(ChanIdx))
	channel = ChanIdx(i_chan);
	disp(['Channel ' num2str(channel, '%02d') ' ' num2str(OriginalErr(1, i_chan)/SNC, '%4.2f') ' mm']);
end

disp(' '); disp(['Result of the Monte-Carlo Simulation (' num2str(N) ' trials):']);

for (i_chan=1:length(ChanIdx))
	channel = ChanIdx(i_chan);
	
	RCF = squeeze(AllReCalibrationFactors(:, i_chan, 1:N))';	
	[RCF, idx] = sortrows(RCF, 7); 	
	better = find(RCF(:, 7) <= OriginalErr(i_chan));
	if (length(better) > 1)
		RFC = RCF(better, :);
		RCF(:, 7) = RCF(:, 7) / SNC;	% Estimated position error in mm (rms over all combinations)	
		disp(['Possible calibration factors for channel ' num2str(channel)])
		disp(num2str([RCF idx]))
	else
		disp(['No better solution found for channel ' num2str(channel)])
	end
	ReCalibrationFactors(:, channel) = RCF(1, 1:6)';
end

RCF = squeeze(AllReCalibrationFactors(1:6, 1, 1:N))';
CN = repmat('00', NChannels, 1);
for (i_chan = 1 : NChannels)
	CN(i_chan, :) = num2str(ChanIdx(i_chan), '%02d');
end

hist(RCF)
xlabel('re-calibration factors'); ylabel('N'); legend(CN);
title(['Distribution of RCFs, N=' num2str(N)])