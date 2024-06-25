function statistics = TrialStat(trial, NumOfTrials, ChanIdx, Result, Residuals, Iterations)
% TRIALSTAT calculates per trial statistics. It can be used in two ways.
%
%   TAPAD callback function 
%           StatMsg = TrialStat(trial, NumOfTrials, ChanIdx, Result,...
%                        Residuals, Iterations);
%
%           The Function returns a text message, containing the statistics for a
%           randomly choosen channel. Printing all channels would overflow
%           the command window, so output is confined when trialstat.m is 
%           used as TAPAD callback function.
%           
%   Stand-alone function 
%           statistics = TrialStat(BasePath, ResPath, ChanIdx, TrialIdx);
%
%				If the output argument (statistics) is omitted, the function displays
%           the statistics for all channels and trials. That can be a lot of output,
%           so use this option with care.
%
%           Statistics are organized as array (NumTrials x 15 x 12-Channels)
%
%     statistcs column format
%         1 trial number
%     2 - 4 median of sensor position (x, y, z) [mm] (ignoring NaNs)
%     5 - 7 estimated robust standard deviation of sensor position (x, y, z) [mm]
%     8 - 9 median of sensor orientation (phi, theta) [°] (ignoring NaNs)
%    10 -11 estimated robust standard deviation of sensor orientation (phi, theta) [°]
%    12 -13 median and deviation of the Carstens "RMS-Value"
%    14 -15 median of iteration depth and number of points per trial
%
%           TrialStat needs the Statistics Toolbox.


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

if (nargin == 6)								% are we called from tapad?
	statistics = TAPADTrialStat(trial, NumOfTrials, ChanIdx, Result, Residuals, Iterations);
	channel = ChanIdx(1 + floor(rand*(length(ChanIdx)-1)));
	statistics = sprintf(['channel %d, trial %d: centeroid at %3.1f %3.1f %3.1f, spatial deviation: %3.1f %3.1f %3.1f [mm]\n   '...
				'average orientation: %3.1f %3.1f (+/- %3.1f %3.1f) [°],  average residual rms %3.3f +/- %3.3f\n   '...
				'%d iterations on average, %d points of %d solved\n'], channel, trial, statistics(2, :, channel));
else
	BasePath = trial;								% Re-map argumets for stand alone version
	ResPath = NumOfTrials; TrialIdx = Result;

	% constant definitions, build file names according to this schema:
	% [BasePath FILESEP Path FILESEP NAMEBASE num2str(trial, SUFIX_FORMAT) EXT_DEL FILEEXT]
	RESULT_NAMEBASE = ''; 				RESULT_FILEEXT = 'mat';
	RESIDUAL_NAMEBASE = 'resi'; 		RESIDUAL_FILEEXT = 'mat';
	SUFIX_FORMAT = '%04d';				EXT_DEL = '.';
	statistics = NaN * ones(length(TrialIdx), 15, 12);

	for i_trial = 1:length(TrialIdx)			% iterate trials (and data files)
		trial = TrialIdx(i_trial);
		ResultFileName = fullfile(BasePath, ResPath, trialfile(trial, RESULT_NAMEBASE));
		ResidualFileName = fullfile(BasePath, ResPath, trialfile(trial, RESIDUAL_NAMEBASE));

		[Residuals, res_comment, res_descriptor, res_dimension, res_private] = loadsdata(ResidualFileName); 
		Iterations = res_private.Iterations;
		[Result, comment, descriptor, dimension, private, samplerate, unit] = loadsdata(ResultFileName); 
		statistics(i_trial, 1, :) = trial;
		statistics(i_trial, 2:15, :) = TAPADTrialStat(trial, length(TrialIdx), ChanIdx, Result, Residuals, Iterations);
		statistics(i_trial, 16, :) = size(Result,1);	
		if (nargout < 1)						% display results, if no output argument given
			for channel = ChanIdx(1): ChanIdx(end)
				disp(sprintf(['channel %d, trial %d: centeroid at %3.1f %3.1f %3.1f, spatial deviation: %3.1f %3.1f %3.1f [mm]\n   '...
						'average orientation: %3.1f %3.1f (+/- %3.1f %3.1f) [°],  average residual rms %3.3f +/- %3.3f\n   '...
						'%d iterations on average, %d points of %d solved.\n'], channel, trial, statistics(i_trial, 2:16, channel)))
			end
		end
	end
end

%------------- subfunctions -----------
function stat = TAPADTrialStat(trial, NumOfTrials, ChanIdx, Result, Residuals, Iterations)
	stat = NaN * ones(1, 14, 12);
	for i_chan = 1:length(ChanIdx)		% iterate channels
		channel = ChanIdx(i_chan); Result(:, 4:5, channel) = NormalizeAngles(Result(:, 4:5, channel));
		stat(1, [1:3 7:8 11], channel) = nanmedian(Result(:, 1:6, channel));
		if ( sum(sum(isnan(Result(:, 1:6, channel)))) < 6)
			stat(1, [4:6 9:10 12], channel) = 0.7413 * iqr(Result(:, 1:6, channel)); % iqr(X) computes the difference between the 75th and the 25th percentiles  
		end
			
		if (size(Iterations, 3) > 1)
			stat(1, 13, channel) = nanmedian(Iterations(:, 1, channel)); 
		else
			stat(1, 13, channel) = NaN;
		end
		stat(1, 14, channel) = nansum(Result(:, 7, channel)); 

	end




