function discal(BasePath, TrialNo, ChanIdx, SampleSize, DiscalOptions);
% DISCAL    Calculates calibration factors from arbitrary quantities of data. 
%           The function repeatedly calculates positions for the given samples and varies
%           re-calibration factors to harmonize the distances between
%           channels. The function can operate on any data set but may be
%           very slow when operating on large number of points. 
%
%           DisCal runs in an infinite loop. To stop it, you must delete the
%           semaphore file 'discal.lck'. (This way, the function be can easily
%           controlled remotely via a ssh or ftp connection or via a Windows share)
%
%           To save disk space and bandwith if loaded via a network connection, the
%           re-calibration data is written in a single precision binary format.
%           Use the function LoadDCData to read the binary file back into MATLAB.
%           
%           discal(BasePath, TrialNo, ChanIdx, SampleSize, DiscalOptions);
%
%           The function also creates a (text-) logfile, containing the calibration data
%           in a readable but truncated format.
%
%
%           see  CAL_PREPDATA, DISCALEVAL

%---------------------------------------------------------------------
% Copyright © 2006-2007 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
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
LOGFILENAME = fullfile(BasePath, 'discal_log.txt');
DATFILENAME = fullfile(BasePath, 'discal.mat');
LOCKFILENAME = fullfile(BasePath, 'discal.lck');
rand('state', sum(100*clock));			% random seed

if (nargin < 5)
	DiscalOptions = '';
end

if (nargin < 4)
	SampleSize = [];
end

if (nargin < 3)
	ChanIdx = 1:12;
end

debug_mode = ~isempty(findstr('-d', DiscalOptions));
quiet_mode = ~isempty(findstr('-q', DiscalOptions));
save_mode  = ~isempty(findstr('-s', DiscalOptions));

data.lsq_options = optimset('lsqnonlin');	% TAPAD-Options
data.lsq_options = optimset(data.lsq_options, 'LargeScale', 'off', 'LevenbergMarquardt', 'on');
if (debug_mode)
	data.lsq_options = optimset(data.lsq_options, 'TolX', 1E-6, 'TolFun', 1E-6, 'MaxIter', 3, 'MaxFunEvals', 200);
else
	data.lsq_options = optimset(data.lsq_options, 'TolX', 1E-8, 'TolFun', 1E-8, 'MaxIter', 30, 'MaxFunEvals', 2000);
end
data.lsq_options = optimset(data.lsq_options, 'Display', 'off', 'Diagnostics', 'off');

% load amplitude data
data.AmpFileName = fullfile(BasePath, trialfile(TrialNo)); 
[data.MAmps, data.amp_comment, data.amp_descriptor, data.amp_dimension,...
	data.amp_private, data.samplerate, data.amp_unit] = loadsdata(data.AmpFileName); 
if (isempty(data.MAmps))
	error(['no such file: ' data.AmpFileName]);
end


data.StartFileName = fullfile(BasePath, 'start_values.mat'); 
data.StartValues = loadsdata(data.StartFileName);

% TAPAD variables
data.NumPoints = size(data.MAmps, 1);

data.ReCalibrationFactors = zeros(6, 12);	% set unused channels to zero 
data.ReCalibrationFactors(:, ChanIdx) = ones(6, length(ChanIdx));

data.Result = NaN * ones(data.NumPoints, 7, 12);
data.Residuals = NaN * ones(data.NumPoints, 6, 12);	
data.TAPADIterations = zeros(data.NumPoints, 1, 12);

data.DerivativParameters = zeros(5*3+1, 6, 12);
if (isempty(SampleSize))
	data.PointIdx = 1:data.NumPoints;
elseif (SampleSize > data.NumPoints)
	error('SampleSize too large');
else
	idx = randperm(data.NumPoints); data.PointIdx = sort(idx(1:SampleSize));
end

if (length(ChanIdx) < 2)
	error('At least two channels required');
else
	data.ChanIdx = ChanIdx;
end
if (exist(DATFILENAME) & (~debug_mode))
	error(['first remove/rename ' DATFILENAME]);
end

% use this semaphore file (lock-file) as a simple method to remotely
% stop execution
[lockfile, msg] = fopen(LOCKFILENAME, 'w');
if (lockfile == -1)
	error(['can not create ' LOCKFILENAME ' ' msg]);
else
	fprintf(lockfile, '%s', 'Remove this file to stop discal');
	fclose(lockfile); clear lockfile;
end
	
if (~quiet_mode)
	if (debug_mode)
		disp(['Distance related re-calibration (DisCal) Version ' tapadversion ' started at ' datestr(now) ' DEBUG-MODE!']);	
	else
		disp(['Distance related re-calibration (DisCal) Version ' tapadversion ' started at ' datestr(now)]);	
	end
	disp(['Amplitude file: ' data.AmpFileName '  DisCal data file: ' DATFILENAME]);	
	disp(''); disp('Remove the file ''discal.lck'' to stop the calculation!');
end
% initialize result data structure
result.Version = ['Distance related re-calibration (DisCal) Version ' tapadversion];
result.data.AmpFileName = data.AmpFileName;
result.NumDataPoints = data.NumPoints;
result.ChanIdx = ChanIdx;
result.PointIdx = uint16(data.PointIdx);
BLOCKSIZE = 100;
result.ReCalibrationFactors = repmat(single(0), [6+1 length(ChanIdx) BLOCKSIZE]);
result.TAPAD = repmat(single(0), [7 length(ChanIdx) BLOCKSIZE]);

tic; data= do_tapad(data);  t = toc; result.cnt = 1;
[r, S] = SensorSpacing('', ChanIdx, data.Result(data.PointIdx, :, :), data.amp_comment, data.samplerate);
result.ReCalibrationFactors(1:6, :, result.cnt) = data.ReCalibrationFactors(:, ChanIdx);
result.ReCalibrationFactors(7, :, result.cnt) = calc_quality(ChanIdx, S);

if (save_mode)
	variant = 0;
	mkdir(BasePath, 'discal_data'); VariantPath = fullfile(BasePath, 'discal_data'); 	
	if (~quiet_mode)
		disp(['saving variants in: ' VariantPath])
	end
	mkdir(VariantPath, 'amps');	mkdir(VariantPath, 'pos');	
	save_intermediate_result(data, VariantPath, variant);
end

if (~quiet_mode)
	disp(['First calculation took ' num2str(round(t/60)) ' min and ' num2str(rem(t, 60)) ' s for ' num2str(length(data.PointIdx)) ' points']);
	disp(['Initial Quality: ' num2str(double(result.ReCalibrationFactors(7, :, result.cnt)) , '%6.3f ')]); 
	diary; diary; % flush diary
end

OZ = zeros(6, 12);	% set unused channels to zero
ON = ones(6, length(ChanIdx)); % and used channels to one
v_tmp = version; Matlab_Is_Newer_Than_6 = str2num(v_tmp(1)) > 6; clear v_tmp;

while exist(LOCKFILENAME)					% semi-infinite monte-carlo loop 
	if (Matlab_Is_Newer_Than_6)
		save(DATFILENAME, 'result', '-v6'); % save last result
	else
		save(DATFILENAME, 'result'); % save last result
	end
	% calculate randomly distributed re-calibration factors 
	% The variation about the neutral value of 1 is normally distributed with mean 0 an standard deviation 1 
	data.ReCalibrationFactors = OZ;	% set unused channels to zero 
	data.ReCalibrationFactors(:, ChanIdx) = ON + 0.01 * randn(6, length(ChanIdx)); % add random values (approx. 0 - 25 Digit, 5 Digit is System Noise)
	try
		data= do_tapad(data);	% calculate positions based on the altered ReCalibrationFactors
		result.cnt = result.cnt + 1;
		if (mod(result.cnt, BLOCKSIZE) == 0) % if buffer full, append memory
			result.ReCalibrationFactors = cat(3, result.ReCalibrationFactors, repmat(single(0), [6+1 length(ChanIdx) BLOCKSIZE]));
		end
		% now calculate the variance of the inter-channel distances
		[r, S] = SensorSpacing('', data.ChanIdx, data.Result(data.PointIdx, :, :), data.amp_comment, data.samplerate);
		result.ReCalibrationFactors(1:6, :, result.cnt) = data.ReCalibrationFactors(:, ChanIdx);
		result.ReCalibrationFactors(7, :, result.cnt) = calc_quality(ChanIdx, S);
		if (save_mode)
			variant = variant + 1 ;
			save_intermediate_result(data, VariantPath, variant);
		end
	catch
		disp(['TAPAD failed: ' lasterr]); 
	end
	if (debug_mode)
 		disp(num2str(result.cnt)); 
	end
end

result.ReCalibrationFactors = result.ReCalibrationFactors(:, :, 1:result.cnt); % discard unused entries
if (Matlab_Is_Newer_Than_6)
	save(DATFILENAME, 'result', '-v6'); 
else
	save(DATFILENAME, 'result'); 
end
return 

%------------- subfunctions -----------
% objective function to rate the 'quality' of the given calibration factors 
function Q = calc_quality(ChanIdx, S)
	Q = zeros(size(ChanIdx));
	for i_chan = 1:length(ChanIdx)		% iterate channels
		channel = ChanIdx(i_chan);
		idx = find(sum(S.Channels == channel)); % find all rows with this sensor in the spacing table
		Q(i_chan) = norm(S.SpatialDist.IQR(idx));
	end
return
		
% perform position calculation without History, e.g. calculate each sample
% independently (like 'tapad -h'). This is necessary because discal evaluates
% the data statistically and the order of the samples should not alter the result. 
function data = do_tapad(data)
	startpoint = [0 0 0 0 0];
	Derivatives = ones(size(data.MAmps));
	for i_chan = 1:length(data.ChanIdx)		% iterate channels
		channel = data.ChanIdx(i_chan);
		for i_point = 1:length(data.PointIdx)					% iterate samples
			point = data.PointIdx(i_point);
			[data.Result(point, 1:7, channel), data.Residuals(point, :, channel), data.TAPADIterations(point, 1, channel)] =...
					calcpos(data.MAmps(point, :, channel) .* data.ReCalibrationFactors(:, channel)', Derivatives(point, :, channel), data.StartValues(point, 1:5, channel), data.lsq_options); 
			data.Result(point, 4:5, channel) = normalizeangles(data.Result(point, 4:5, channel));
		end
	end
return

%
function save_intermediate_result(data, VariantPath, variant)
	AmpFileName = fullfile(VariantPath, 'amps', trialfile(variant));
	CF =	repmat(shiftdim(data.ReCalibrationFactors, -1), [data.NumPoints 1 1]);
	CAmps = data.MAmps .* CF; DParams = derivative(CAmps);
	comment = frametext('discal',  ['discal ' tapadversion ' intermediate amplitude file. MATLAB Version ' version], {'written', datestr(now,0)});
	SensorNames = {'Sensor1', 'Sensor2', 'Sensor3', 'Sensor4', 'Sensor5', 'Sensor6', 'Sensor7', 'Sensor8', 'Sensor9', 'Sensor10', 'Sensor11', 'Sensor12'}';
	save_amp_data(AmpFileName, single(CAmps), data.samplerate, SensorNames, comment, DParams)

	ResultFileName = fullfile(VariantPath, 'pos', trialfile(variant));
	dimension = struct('descriptor', '', 'unit', '', 'axis', ''); private = struct('Iterations', ''); unit = ''; descriptor = ''; 
	
	Residuals = NaN * ones(data.NumPoints, 6, 12);	res_private.Iterations = zeros(data.NumPoints, 1, 12);
	private.DerivativParameters = DParams; private.PointIdx = data.PointIdx; comment = data.amp_comment;
	comment = frametext('discal',  ['discal ' tapadversion ' intermediate result file. MATLAB Version ' version], {'written', datestr(now,0)});
	save_pos_data(ResultFileName, single(data.Result), comment, descriptor, dimension, private, data.samplerate, unit, 'tapad_result.mat');
	
	ResidualFileName = fullfile(VariantPath, 'pos', trialfile(variant, 'resi'));
	res_private.Iterations = data.TAPADIterations;	res_private.DerivativParameters = DParams;
	save_pos_data(ResidualFileName, single(data.Residuals), comment, descriptor, dimension, res_private, data.samplerate, unit, 'tapad_residual.mat');
return


function save_amp_data(filename, data, samplerate, sensor_names, comment, DParams)
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


function save_pos_data(filename, data, comment, descriptor, dimension, private, samplerate, unit, filename2)
v_tmp = version; Matlab_Is_Newer_Than_6 = str2num(v_tmp(1)) > 6; clear v_tmp;

try
	if Matlab_Is_Newer_Than_6	
		save(filename, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit', '-v6');
	else	
		save(filename, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit');
	end
catch
	if Matlab_Is_Newer_Than_6	
		save(filename2, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit', '-v6');
	else	
		save(filename2, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit');
	end
    error([lasterr ' saved as ' filename2]);
end
    
