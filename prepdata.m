function prepdata(BasePath, AmpPath, ResPath, TrialIdx, FilterType, DownFac, SensorNames, Comment);
% PREPDATA stores AG500 amplitude data as MAT files.
%           The function also performs data pre-processing, if desired.
%           BasePath is somethig like 'C:\MyStuf\EMA', i.e. the base 
%           for all data. 
%           AmpPath and ResPath defines the relative path for input and
%           output data, e.g. AmpPath = 'Data\Amp'; ResPath = 'Data\Mat'; 
%            
%           TrialIdx defines the trials (files) to be processed. Input data is 
%           expected from one directory, with a 4-digits file pattern, like '0001.amp'.
%            
%           FilterType and DownFac can be used to apply a low-pass filter or to perform
%           a initial downsampling of the amplitude data. See DECFIX.m for  
%           further information on downsampling. 
%
%           SensorNames has to be a (12x1) cell array, containing strings
%           with the sensor names. If omitted, default is 'Sensor1' to 'Sensor12'.
%           
%           Comment is a optional string with additional remarks, which are to be 
%           stored with the data.
%           
%           prepdata(BasePath, AmpPath, ResPath, TrialIdx, FilterType, DownFac, SensorNames, Comment);
%
%
%   SEE ALSO TAPADM, DECIFIX

%---------------------------------------------------------------------
% Copyright © 2005 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
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

AMP_NAMEBASE = ''; 			AMP_FILEEXT = 'amp'; 
RESULT_NAMEBASE = ''; 		RESULT_FILEEXT = 'mat';

VERSIONINFO = ['PREPDATA ' tapadversion ' Result File. MATLAB Version ' version];
clc; 								disp(['Prepdata ']);
errormsg = nargchk(4, 8, nargin);
if (~isempty(errormsg))
	error(errormsg);
end

if (nargin < 8) | isempty(Comment)
	Comment = '';
end

if (nargin < 7) | isempty(SensorNames)
	SensorNames = {'Sensor1', 'Sensor2', 'Sensor3', 'Sensor4', 'Sensor5', 'Sensor6', 'Sensor7', 'Sensor8', 'Sensor9', 'Sensor10', 'Sensor11', 'Sensor12'}';
end

if (nargin < 6) | isempty(DownFac)
	samplerate = 200; DownFac = 1;
else
	samplerate = 200/DownFac;
end

if (nargin < 5) | isempty(FilterType)
	FilterType = 'NONE';
end

for i_trial = 1:length(TrialIdx)
   trial = TrialIdx(i_trial);
	refreshDisplay(trial);
   AmpFileName = fullfile(BasePath, AmpPath, trialfile(trial, AMP_NAMEBASE, '', AMP_FILEEXT));
   ResultFileName = fullfile(BasePath, ResPath, trialfile(trial, RESULT_NAMEBASE));
	if (exist(ResultFileName, 'file'))
		disp(['Skipping trial ' int2str(i_trial) '; already existing MAT file: ' ResultFileName]);
		continue;
	end
	[MAmps, NumPoints] = loaddata(AmpFileName); % load the data
	if (isempty(MAmps))
		disp(['Skipping trial ' int2str(i_trial) '; no such file: ' AmpFileName]);
		continue;
	end
	
	newcomment = sprintf('%s\nAMP file: %s\n', VERSIONINFO, AmpFileName);
	fileinfo = dir(AmpFileName);
	if ~isempty(fileinfo)
		newcomment = sprintf('%sTimestamp of AMP file: %s\n', newcomment, fileinfo(1).date);
	end
	if (length(Comment) > 0)
		newcomment = sprintf('%sComment: %s\n', newcomment, Comment);
	end
	if (~strcmp(FilterType, 'NONE'))
		[MAmps, FilterDesc] = decifix(MAmps, FilterType);
		newcomment = sprintf('%sFilter: %s\n', newcomment, FilterDesc);
	end
	if (DownFac ~= 1)
		newcomment = sprintf('%sDownsampling: %d\n', newcomment, round(DownFac));
		Tmp = downsample(MAmps(:,:,1), DownFac);
		TAmps = zeros(size(Tmp, 1), 6, 12);	% rem(NumPoints,DownFac) has to be zero!
		TAmps(:,:,1) = Tmp;
		for i=2:12
			TAmps(:,:,i) = downsample(MAmps(:,:,i), DownFac);
		end
		MAmps = TAmps;
	end
	DParams = derivative(MAmps);
	newcomment = sprintf('%sSampling rate: %s Hz', newcomment, num2str(samplerate));
	fcomment = [frametext('prepdata',  newcomment, {'written', datestr(now,0)})];
	fcomment = sprintf('%s\n', fcomment);

	v_tmp = version; Matlab_Is_Newer_Than_6 = str2num(v_tmp(1)) > 6; clear v_tmp;
	savedata(ResultFileName, single(MAmps), samplerate,  SensorNames, fcomment, DParams);
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
	v_tmp = version; Matlab_Is_Newer_Than_6 = str2num(v_tmp(1)) > 6; clear v_tmp;
	if Matlab_Is_Newer_Than_6	
		save(filename, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit', '-v6');
	else
		save(filename, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit');
	end
	
function refreshDisplay(trial)
	clc; disp(['PrepData Version ' tapadversion]); disp(' ');
	disp(['trial: ', num2str(trial)]);
