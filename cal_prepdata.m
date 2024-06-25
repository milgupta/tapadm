function cal_prepdata(DataPath, ResultPath, UserComment);
% CAL_PREPDATA stores AG500 calibration data as MAT files.
%           The function reads the data from the files '\channel..\circal.txt' and also some information
%           from '...\ag500.ini'. It stores the data in the file 'calibration.mat' under the DataPath or
%           any given ResultPath. If a UserComment is given, it will be saved with the data.
%           The calibration data is organized in a structure array, one enty per channel.
%           
%           cal_prepdata(DataPath, ResultPath, UserComment);
%
%           alpha, the rotation angle of the calibration device is converted into radians, assuming
%           800 encoder increments per revolution. alpha rotates counter-clockwise in the X/Z-Plane   
%           of the AG500, relative to the first point, i.e. alpha0. 
%	         If the performed circal movement covers more than one revolution, data with angles alpha > 360°
%           results. Such data is unwrapped by this function, so that the resulting data is always in the 
%           range 0 <= alpha < 2pi. 

%---------------------------------------------------------------------
% Copyright © 2006 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
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
NUMCHAN = 12;
VERSIONINFO = ['CAL_PREPDATA ' tapadversion ' Result File. MATLAB Version ' version];

if ((nargin == 0) | (length(DataPath)==0))	
	DataPath = pwd;
end

if (nargin < 2)
	ResultFile = fullfile(DataPath, 'calibration.mat');
else
	ResultFile = fullfile(ResultPath, 'calibration.mat');
end

Comment = VERSIONINFO;
if (nargin > 2)
	Comment = sprintf('%s\n%s\n', Comment, UserComment);
end

for (chan=1:NUMCHAN)
	[a, Ap] = loaddata(fullfile(DataPath, ['channel' num2str(chan)], 'circal.txt')); n = size(Ap,1);
	[CFac, TransPol] = ReadProfile(fullfile(DataPath, ['channel' num2str(chan)], 'ag500.ini'));	
	CalData(chan) = struct('alphas', a, 'Amps', Ap, 'NumPoints', n, 'OriginalCalibrationFactor', CFac,...
		'TransmitterPolarity', TransPol, 'Channel', num2str(chan));
end

if (exist(ResultFile, 'file'))
	ResultFileBak = ResultFile(1:findstr(ResultFile, '.mat')-1);
	ResultFileBak = [ResultFileBak '_bak.mat'];
	disp('Renaming existing result file.')
	copyfile(ResultFile, ResultFileBak);	% backup existing file
end

v_tmp = version; Matlab_Is_Newer_Than_6 = str2num(v_tmp(1)) > 6; clear v_tmp;
if (Matlab_Is_Newer_Than_6)
	save(ResultFile, 'VERSIONINFO', 'DataPath', 'CalData', '-v6');
else
	save(ResultFile, 'VERSIONINFO', 'DataPath', 'CalData');
end
%------------- subfunctions -----------
function [alphas, Amps] = loaddata(FileName) 
	cdata = load(FileName);						% Load text-file with AG500 calibration data
	Amps = cdata(:, [2 4 6 8 10 12]);		% get signal amplitudes
	Amps = Amps .* repmat([-1 1 -1 1 -1 1], size(Amps, 1), 1);		% get signal amplitudes
	% get alpha, the rotation angle of the calibration device and convert it to radians, 
	% relative to the first point. That will be the abscissa of the calibration data.
	ENCODER_INCREMENTS_PER_REVOLUTION = 8000;
	alphas = (cdata(:,1) - cdata(1,1)) .* 2*pi/ENCODER_INCREMENTS_PER_REVOLUTION; 
	% map data to one revolution: 0 <= alpha < 2pi 
	idx = find(alphas < 0); alphas(idx) = alphas(idx) + 2*pi;
	[alphas, idx] = sort(alphas); Amps = Amps(idx, :); clear idx; clear cdata;
return

function [CFac, TransPol] = ReadProfile (FileName);
% READPROFILE Read various parameters from ini-file 

[fid, xmsg] = fopen(FileName, 'r');
if fid == -1
	error([xmsg ' Filename: ' FileName]);
end

TPSection = 0;

while 1
	line = fgetl(fid);						% read line
	if ~isstr(line), break, end
	cpos = findstr(line, ';');				% remove remarks
	if (isempty(cpos) == 0)
		cidx = cpos(1) : length(line);
		line(cidx) = [];
	end
	line = lower(line);

	if TPSection 								% process only [tapad]-section of file
		if (isempty(findstr(line, '[')) == 0)
			TPSection = 0;
		end

		if (isempty(findstr(line, 'calibration factors')) == 0)
			cpos = findstr(line, '=');
			cidx = 1 : cpos(1);
			line(cidx) = [];
			line = [ '[' line ']' ];
			CFac = str2num(line);
		end

		if (isempty(findstr(line, 'transmitter polarity')) == 0)
			cpos = findstr(line, '=');
			cidx = 1 : cpos(1);
			line(cidx) = [];
			deblank(line); j=1;
			for i=1:length(line)
				if line(i) == '+'
					TransPol(j) = 1; j = j+1;
				elseif line(i) == '-' 
					TransPol(j) = -1; j = j+1;
				end
			end
			if j~= 7
				warning(['Transmitter polarity tag is corrupt: ' line]);
			end
		end

	else
		TPSection = (isempty(findstr(line, '[tapad]')) == 0);
	end

end
fclose(fid); 

