function rescaledata(BasePath, InputPath, OutputPath, TrialIdx, ChanIdx, ScaleFactors);
% RESCALEDATA rescales AG500 amplitude data, i.e. multiplies it with the provided
%           factors and saves the altered data in a MAT file.
%           The function writes a appropriate comment to the new file to
%           document the change.
%           
%           rescaledata(BasePath, InputPath, OutputPath, TrialIdx, ChanIdx, ScaleFactors);
%
%           ScaleFactors can be 6x1 vector (which will be expanded) or a
%           6x12 matrix with individual scaling factors for each channel.
%
%   SEE ALSO DISCAL, PREPDATA

%---------------------------------------------------------------------
% Copyright © 2005 - 2007 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
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
error(nargchk(6, 6, nargin));

if (size(ScaleFactors, 1) ~= 6) | ((size(ScaleFactors, 2) ~= 1) & (size(ScaleFactors, 2) ~= 12))
	error('size(ScaleFactors) has to be 6x1 or 6x12')
end
if (size(ScaleFactors, 2) == 1)
	ScaleFactors = repmat(ScaleFactors, 1, 12);
end
ScaleFactors = shiftdim(ScaleFactors, -1);

for i_trial = 1:length(TrialIdx)
   trial = TrialIdx(i_trial);
   
   InFileName = fullfile(BasePath, InputPath, trialfile(trial)); % load original data
	[data, comment, descriptor, dimension, private, samplerate, unit]...
					= loadsdata(InFileName); NumPoints = size(data, 1);
	if (isempty(data))
		disp(['Skipping trial ' int2str(trial) '; no such file: ' InFileName]); continue;
	end
	
	SF = repmat(ScaleFactors, NumPoints, 1);	data = data .* SF; % scale the data
	
	ScaleFactorsmsg = sprintf('%6.5f ', ScaleFactors);		% write a comment
	ScaleFactorsmsg = frametext('rescaledata', ScaleFactorsmsg, {'written', datestr(now,0)});
	comment = sprintf('%s\n%s\n', comment, ScaleFactorsmsg);
	
   OutFileName = fullfile(BasePath, OutputPath, trialfile(trial)); % save the rescaled data
   if (exist(OutFileName, 'file'))		% backup a existing output file
	   OutFileBakName = fullfile(BasePath, OutputPath, trialfile(trial, '', '_bak'));
		copyfile(OutFileName, OutFileBakName);
	end

	v_tmp = version; Matlab_Is_Newer_Than_6 = str2num(v_tmp(1)) > 6; clear v_tmp;
	if Matlab_Is_Newer_Than_6	
		save(OutFileName, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit', '-v6');
	else
		save(OutFileName, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit');
	end
end	% of  for i_trial = 1:length(TrialIdx)


%------------- subfunctions -----------
function savedata(filename, data, samplerate, sensor_names, comment, DParams)
	descriptor = ['Trans1'; 'Trans2'; 'Trans3'; 'Trans4'; 'Trans5'; 'Trans6'];	unit = repmat('NormalizedAmp',[6 1]);
	dimension.descriptor = char(ones(3, 11)); 
	dimension.descriptor(1, :) = 'Time       ';
	dimension.descriptor(2, :) = 'Transmitter';
	dimension.descriptor(3, :) = 'Sensor     ';
	dimension.unit = char(repmat(32, 3, 1));
	dimension.axis = cell(3, 1); 
	dimension.axis(2) = {descriptor}; 
	dimension.axis(3) = {sensor_names}; 
	private.DerivativParameters = DParams;

function refreshDisplay(trial)
	clc; disp(['ReScaleData Version ' tapadversion]); disp(' ');
	disp(['trial: ', num2str(trial)]);

