function tapad(BasePath, AmpPath, ResPath, TrialIdx, ChanIdx, TapadOptions, StartValues, UserFunction);
% TAPAD Three-dimensional Artikulographic Position and Align Determination.
%           This is the core function of the TAPADM Toolbox, which performs 
%           the position calculation for a set of EMA measured data files.
%
%           tapad(BasePath, AmpPath, ResPath, TrialIdx, ChanIdx, Options,...
%                 StartValues, UserFunction);
%
%   Mandatory Arguments 
%           BasePath is somethig like 'C:\MyStuf\EMA', i.e. the base 
%           for all data. Hold your BasePath in an global Variable,
%           or use the pwd command.
%            
%           AmpPath and ResPath are defining the relative path for input and
%           output data, e.g. AmpPath = 'Data\Amp'; ResPath = 'Data\Pos'; 
%            
%           TrialIdx defines the trials (files) to be processed. Input data is 
%           expected from one directory, with a 4-digits file pattern, like '0001.mat'.
%           Perform the data pre-processing prior to tapad.m to generate
%           these input data from the AG500 '.amp' files
%            
%           ChanIdx names the channels to be computed, e.g. 1:12. If there is 
%           already an existing TAPADM result file and TAPADM is used to calculate
%           position data just for some channels, the function will keep existing 
%           position data which is not affected by this run. Thus, position calculating 
%           can be performed in a sequential order: channel by channel.
%            
%   Options
%           -a automatically starts the position calculation at the best suited 
%              sample in the trial. Continue in both directions from that
%              point on. (Disables -r and -f)
%           -d use amplitude derivatives to weight errors
%           -f flip time, i. e. process data onwards from the last point
%           -h don't use history (i.e. last result) as start point (significantly 
%              increases computation time!)
%           -l use Levenberg-Marquardt instead of Newton method
%           -r recursive mode, where individual start values for every sample
%              are taken from former calculated result files. The path to
%              this files is expected as 7th function-argument
%           -s user supplied initial start value for the first sample in
%              the trial. Start values are expected as a 7th
%              function-argument. (see below)
%           -t turbo mode for test purposes.
%           -q quiet mode
%
%   Optional Arguments
%           StartValues is a manifold argument: 
%           In recursive mode (-r), StartValues defines the path where
%           files with individual start positions can be found. (similar to
%           AmpPath)
%           StartValues in conjunction with -s defines the start position to 
%           be used for first sample in trial (arranged 5 x channels). 
%           Accessed by channel number, not position in ChanIdx list!!!
%
%           UserFunction: If supplied, tapad.m will call this function for
%           every trial with the folliwing arguments:
%           msg = UserFunction(trial, NumOfTrials, ChanIdx, Result, Residuals, Iterations)
%
%           If the function returns a character string, tapad will display
%           it during the calculation of the next trial.

% TAPAD's error policy is: If there ist something wrong with a mandatory
% argument, generate an error. If it's a optional argument/feature, display
% a warning and ignore it. On principle, try to avoid fatal error
% conditions and keep on going.

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

% constant definitions, build file names according to this schema:
AMP_NAMEBASE = ''; 					AMP_FILEEXT = 'mat'; 
RESULT_NAMEBASE = ''; 				RESULT_FILEEXT = 'mat';
RESIDUAL_NAMEBASE = 'resi'; 		RESIDUAL_FILEEXT = 'mat';
DERIVATIVES_NAMEBASE = 'drvt';   DERIVATIVES_FILEEXT = 'mat';
SUFIX_FORMAT = '%04d';				EXT_DEL = '.';
VERSIONINFO = ['TAPADM ' tapadversion ' Result File. MATLAB Version ' version];

% check function arguments and evaluate TapadOptions
errormsg = nargchk(5, 8, nargin);
if (~isempty(errormsg))
	error(errormsg);
end
if (nargin < 6)								
	TapadOptions = '';
end
use_best_start_point =		~isempty(findstr('-a', TapadOptions));
use_derivatives = 			~isempty(findstr('-d', TapadOptions));
use_flip_data = 				~isempty(findstr('-f', TapadOptions));
use_no_history = 				~isempty(findstr('-h', TapadOptions));
use_initial_start_point =	~isempty(findstr('-s', TapadOptions));
use_turbo =						~isempty(findstr('-t', TapadOptions));
use_levenberg =				~isempty(findstr('-l', TapadOptions));
use_idividual_start_points=~isempty(findstr('-r', TapadOptions));
use_display =		        	isempty(findstr('-q', TapadOptions));
if (use_display)
    clc; disp(['Three-dimensional Artikulographic Position and Align Determination Version ' tapadversion]); % say hello
end

% honour exclusive options
if (use_best_start_point)
	if (use_idividual_start_points)
		use_idividual_start_points = logical(0);	
		warning('Option -r was ignored due to -a mode.'); pause(3);
	end
	if (use_flip_data)
		use_flip_data = logical(0);	
		warning('Option -f was ignored due to -a mode.'); pause(3);
	end
end

if (use_initial_start_point)
	if (use_idividual_start_points)
		use_idividual_start_points = logical(0);	
		use_initial_start_point = logical(0);	
		warning('Conflicting options -s and -r were ignored.'); pause(3);
	elseif (nargin < 7)
		use_initial_start_point = logical(0);	
		warning('Option -s was ignored due to missing start-value.'); pause(3);
	elseif (size(StartValues, 1) ~= 5) | (size(StartValues, 2) ~= 12)
		use_initial_start_point = logical(0);	
		warning('StartValue must be 5x12 matrix! Option -s was ignored.'); pause(3);
	end
end

if (use_idividual_start_points)
	if (nargin < 7)
		use_idividual_start_points = logical(0);	
		warning('Option -r was ignored due to missing start-value.'); pause(3);
	elseif ~ischar(StartValues)
		use_idividual_start_points = logical(0);	
		warning('User start value must be path to input positions in recursive mode! Option -s was ignored.'); pause(3);
	end
end
	  
if ((nargin == 7) & ~use_idividual_start_points & ~use_initial_start_point)
	warning('Extra argument ''StartValues'' ignored!'); pause(3);
end

% miscellaneous arguments
if (isempty(BasePath))
	BasePath = '';
end
use_user_function = (nargin > 7); 	user_msg = '';

%---------------------------------------------------------------------
% m a i n   l o o p
%---------------------------------------------------------------------
for i_trial = 1:length(TrialIdx)			% iterate trials (and data files)
   trial = TrialIdx(i_trial);
   AmpFileName = 			fullfile(BasePath, AmpPath, trialfile(trial, AMP_NAMEBASE, 			'',	AMP_FILEEXT, 			EXT_DEL));
   ResultFileName = 		fullfile(BasePath, ResPath, trialfile(trial, RESULT_NAMEBASE,		'',	RESULT_FILEEXT,		EXT_DEL));
   ResidualFileName = 	fullfile(BasePath, ResPath, trialfile(trial, RESIDUAL_NAMEBASE, 	'',	RESIDUAL_FILEEXT,		EXT_DEL));
   DerivativesFileName= fullfile(BasePath, ResPath, trialfile(trial, DERIVATIVES_NAMEBASE, '',	DERIVATIVES_FILEEXT, EXT_DEL));
   ResultFileBak = 		fullfile(BasePath, ResPath, trialfile(trial, RESULT_NAMEBASE,	 '_bak',	RESULT_FILEEXT,		EXT_DEL));
   ResidualFileBak = 	fullfile(BasePath, ResPath, trialfile(trial, RESIDUAL_NAMEBASE, '_bak',	RESIDUAL_FILEEXT,		EXT_DEL));
   DerivativesFileBak = fullfile(BasePath, ResPath, trialfile(trial, DERIVATIVES_NAMEBASE,'_bak',DERIVATIVES_FILEEXT,EXT_DEL));
	if (use_idividual_start_points)
		IndividualStartPosFileName = fullfile(BasePath, StartValues, trialfile(trial, RESULT_NAMEBASE, '', RESULT_FILEEXT, EXT_DEL));
	else
		IndividualStartPosFileName	= '';
	end
   ParameterInfo = ['Calculated Positions from ' AmpFileName ' Channel Index: ' num2str(ChanIdx) ' TAPAD Options: ' TapadOptions];

	[MAmps, amp_comment, amp_descriptor, amp_dimension, amp_private, samplerate, amp_unit] = loadsdata(AmpFileName); NumPoints = size(MAmps, 1); % load the data
	% We do not want tapad to quit the whole job, if some files are missing.
	% It could be very annoying if e.g. trial 5 of 500 fails and one notices after hours, that
	% the whole job has aborted...   			
	if (isempty(MAmps))
        if (use_display)
		    disp(['Skipping trial ' int2str(trial) '; no such file: ' AmpFileName]);
        end
		continue;
	end
		
	% There are two output files: one for the standard result and the other for more sophisticated 
	% output like residuals, iterations and other private data. For both, the function tries to preserve 
	% any existing data, to enable serial computing of the channels.
	% One output file is called 'residual-file' and is considered to depend on the standard result-file.
	% (The Carstens 'CalcPos'-Program calls the residuals 'PosAmps')
	if (exist(ResidualFileName, 'file'))
		[Residuals, res_comment, res_descriptor, res_dimension, res_private] = loadsdata(ResidualFileName);
		copyfile(ResidualFileName, ResidualFileBak);	% backup existing file
	else
		[Residuals, res_comment, res_descriptor, res_dimension, res_private] = loadsdata; % just create the variables		
	end
	[N, M, L] = size(Residuals);
	if ((N~=NumPoints) | (M~=6) | (L~=12))	% do we already have valid residual data? (from the file)
		if (~isempty(Residuals))	
			warning('Existing tapad residual data discarded due to incompatible format!');
		end	
		Residuals = NaN * ones(NumPoints, 6, 12);  % If not, initialize the result variables
		res_private.Iterations = zeros(NumPoints, 1, 12);
		res_private.DerivativParameters = zeros(5*3+1, 6, 12);
	else											% if residual array is ok, check for the other results
		[N, M, L] = size(res_private.Iterations);
		if ((N~=NumPoints) | (M~=6) | (L~=12))
			res_private.Iterations = zeros(NumPoints, 1, 12);
			res_private.DerivativParameters = zeros(5*3+1, 6, 12);
		end
	end	

	% Now check for a result-file which contains pre-computed positions (e.g. from different channels). 
	% If there is no usable result-file, the residual-file data will be nulled whether there is 
	% old residual data, or not (since keeping it without results would make no sense). 
	if (exist(ResultFileName, 'file'))
 	  [Result, comment, descriptor, dimension, private, res_samplerate, unit] = loadsdata(ResultFileName); 
   else
 	  [Result, comment, descriptor, dimension, private, res_samplerate, unit] = loadsdata; 
   end	  
	[N, M, L] = size(Result);
   if (~isempty(Result))	
		copyfile(ResultFileName, ResultFileBak);
	end
	if ((N~=NumPoints) | (M~=7) | (L~=12))
		if (~isempty(Result))	
			warning('Existing tapad data discarded due to incompatible format!');
		end	
		Result = NaN * ones(NumPoints, 7, 12); data = single([]); comment = ''; descriptor = ''; 
		dimension = struct('descriptor', '', 'unit', '', 'axis', ''); private = struct('Iterations', ''); unit = '';
		Residuals = NaN * ones(NumPoints, 6, 12);	res_private.Iterations = zeros(NumPoints, 1, 12);
		res_private.DerivativParameters = zeros(5*3+1, 6, 12); comment = amp_comment;
	end	
	Iterations = res_private.Iterations;

	% load idividual start points if requested
	if (use_idividual_start_points)
		StartValues = loadsdata(IndividualStartPosFileName); [N, M, L] = size(Result);
		if ((N~=NumPoints) | (M>=5) | (L~=12))
			use_idividual_start_points = logical(0);	
			warning('Existing start value data file discarded due to incompatible format!');
		end
	end	
	
	% Derivative-Parameters are already stored in the AmpFile, if it has been
	% generated by 'prepdata'. This parameters allow to roughly estimate the
	% run of the amplitude and the gradient. If TAPAD is started with the
	% '-a'-Option to find the best start-point, Derivative-Parameters are
	% required. In that case, the script computes them belatedly.
	% (Derivative-Parameters are stored in the result file, wether or not, 
	% they already exist in the 'amp'-file.)
	DParams = [];
	if (exist('amp_private') == 1)
		if	(isfield(amp_private,'DerivativParameters')==1)
	 		[N, M, L] = size(amp_private.DerivativParameters);
			if ((N==16) | (M==6) | (L==12))
				DParams = amp_private.DerivativParameters;
			else
				warning('Existing derivative parameters discarded due to incompatible format!');
				DParams = derivative(MAmps);	
			end
		end
	end
	if (use_best_start_point & isempty(DParams))
		warning('No derivative parameters in amplitude file found, calculating them belatedly');
		DParams = derivative(MAmps);	
	end
	
	% Based on the Derivative-Parameters, TAPAD can generate a approximation of
	% the first order derivative, to estimate the relation between residuals
	% and position errors. (Option '-d')
   Derivatives = ones(size(MAmps));	
   if (use_derivatives)	
       if (use_display)
		    disp('computing derivatives...');
        end
		[L, N, M] = size(MAmps);
		t = (1/samplerate) * (1 : L)';					% time axis 
		for i=1:N
			for j=1:M
				[approx_amp, approx_vel] = trig_approx_func(DParams(:,i,j)', t);
				Derivatives(:, i, j) = approx_vel';
			end
		end
	end	

	% Tha data can be 'flipped' to perform the position calculation in reverse
	% order (last to first point). 
	if (use_flip_data)
		Result = flipdim(Result,1);
		Residuals = flipdim(Residuals,1);
		Iterations = flipdim(Iterations,1);
		MAmps = flipdim(MAmps,1);
		Derivatives = flipdim(Derivatives,1);
	end	

	% In 'turbo' mode, the iteration depth is clearly reduced, to speed up the position calculation
	if (use_turbo)
		MaxIter = 10;
	else
		MaxIter = 1000;
	end		
	% Set up the parameters for the position calculation (all points)
	lsq_options = optimset('lsqnonlin'); 
	if (use_levenberg)	% 	select method		
		lsq_options = optimset(lsq_options, 'LargeScale', 'off', 'LevenbergMarquardt', 'on', 'TolX', 1E-8, 'TolFun', 1E-8, 'MaxIter', MaxIter);
	else
		lsq_options = optimset(lsq_options, 'LargeScale', 'on', 'LevenbergMarquardt', 'off', 'TolX', 1E-5, 'TolFun', 1E-7, 'MaxIter', MaxIter);
	end
	lsq_options = optimset(lsq_options, 'Display', 'off', 'Diagnostics', 'off');

	% Set up individual parameters for the position calculation of the
	% startpoint of the trial
	lsq_options_start = optimset('lsqnonlin');
	lsq_options_start = optimset(lsq_options_start,...
							'LargeScale', 'off',		'LevenbergMarquardt', 'on',...
							'TolX', 1E-8,				'TolFun', 1E-8,...
							'MaxIter', MaxIter,		'MaxFunEvals', 2000);
	if (use_display)
		lsq_options_start = optimset(lsq_options_start, 'Display', 'iter', 'Diagnostics', 'on');	
 	else
		lsq_options_start = optimset(lsq_options_start, 'Display', 'off', 'Diagnostics', 'off');	
	end
	
	%================== Position Calculation ==============================
	for i_chan = 1:length(ChanIdx)		% inner loop, iterate channels
		channel = ChanIdx(i_chan);
		msg = sprintf('%s\nprocessed trial: %d channel: %d', user_msg, trial, channel);
   	
 		if (use_initial_start_point)		% initialize the startpoint with... 
			startpoint = StartValues(channel,:);	% an individual value,
		else
			startpoint = [0 0 0 0 0];		% or zeros.
		end       
         
		if (use_best_start_point)
			StartIdx = FindBestStart(MAmps, Derivatives, channel); % look for a promising point to start with, than calculate position for that point
			[Result, Residuals, Iterations] = calcpos4trial(Result, Residuals, Iterations, MAmps, StartIdx, StartIdx, channel, startpoint, lsq_options_start, [], []);
			if (Iterations(StartIdx, 1, channel)  > -1)			% Does the calculation yields a result? (that's what we expect)
				startpoint = Result(StartIdx, 1:5, channel);		% if so, take the result as new startpoint!
				best_start_point_exists = true;
			else
				msg = sprintf('%s\nNo best start point found, switched to normal mode!', msg);
				best_start_point_exists = false;	clear StartIdx;			
			end
		else
			best_start_point_exists = false;				
		end			
			
		if (best_start_point_exists)		% There is an optimal start point inside the trial, so position calculation
			% Starting from the optimal StartIdx, calculate positions backwards to the beginning of the trial
	    	if (use_display)
				[Result, Residuals, Iterations] = calcpos4trial(Result, Residuals, Iterations, MAmps, StartIdx, 1, channel, startpoint, lsq_options, @refreshDisplay, msg);
			else
				[Result, Residuals, Iterations] = calcpos4trial(Result, Residuals, Iterations, MAmps, StartIdx, 1, channel, startpoint, lsq_options);
			end
			% After that, calculate positions from that StartIdx on to the end of the trial
	    	if (use_display)
				[Result, Residuals, Iterations] = calcpos4trial(Result, Residuals, Iterations, MAmps, StartIdx, NumPoints, channel, startpoint, lsq_options, @refreshDisplay, msg);
			else
				[Result, Residuals, Iterations] = calcpos4trial(Result, Residuals, Iterations, MAmps, StartIdx, NumPoints, channel, startpoint, lsq_options);
			end
		else										% don't use FindBestStart()
			for i = 1:NumPoints					% iterate points from first to last
				if (use_idividual_start_points)
					startpoint = StartValues(i, 1:5, channel); 	% overwrite startpoint
				end			
 				[Result(i, 1:7, channel), Residuals(i, :, channel), Iterations(i, 1, channel)] =...
					calcpos(MAmps(i, :, channel), Derivatives(i, :, channel), startpoint, lsq_options); % perform position calculation
				Result(i, 4:5, channel) = NormalizeAngles(Result(i, 4:5, channel)); % make sure, that -180 < phi <= +180 and -90 < theta <= +90
			
				if use_no_history				% Are we going to make further use of the just calculated result?  
					if use_initial_start_point
						startpoint = StartValues(:, channel); 	% If not, reset to use_initial_startvalue if available
					else
						startpoint = [];								% if no initial value exists, void start value
					end
				else								% With history...
					startpoint = Result(i, 1:5, channel);	% use result as startpoint for next sample
				end
            
				if (use_display)
    				if (rem(i,100) == 0)     
	 	    			refreshDisplay(i, msg);
					end
				end
				
			end			% of  for i = 1:NumPoints
		end
		
	end				% of inner loop:  for i_chan = 1:length(ChanIdx)
	%============== Position Calculation done =============================
	
	if (use_flip_data)	% re-flip data			
		Result = flipdim(Result,1);
		Residuals = flipdim(Residuals,1);
		Iterations = flipdim(Iterations,1);
		MAmps = flipdim(MAmps,1);
		Derivatives = flipdim(Derivatives,1);
	end	
   
	if (use_user_function)
		try
			user_msg = feval(UserFunction, trial, length(TrialIdx), ChanIdx, Result, Residuals, Iterations);
		catch
			[errmsg, msgid] = lasterr;
			warning(['User function failed: ' errmsg]);
		end
	end
			
   % Now save the results. 
	newcomment = sprintf('%s\n%s\n', VERSIONINFO, ParameterInfo);
	comment = [comment frametext('tapad',  newcomment, {'written', datestr(now,0)})];
	if use_initial_start_point
 	   private.startposition = StartValues;
	end;
	res_private.Iterations = Iterations;	res_private.DerivativParameters = DParams;
	savedata(ResultFileName, single(Result), comment, descriptor, dimension, private, samplerate, unit, 'tapad_result.mat');
	savedata(ResidualFileName, single(Residuals), comment, descriptor, dimension, res_private, samplerate, unit, 'tapad_residual.mat');
   if (use_derivatives)	
		savedata(DerivativesFileName, single(Derivatives), comment, descriptor, dimension, res_private, samplerate, unit);
	end;	
   pause(1); 									% a short break to catch a possible <Ctrl-C>
end	% of  for i_trial = 1:length(TrialIdx)
return											% end of TAPAD

%------------- subfunctions -----------
function savedata(filename, data, comment, descriptor, dimension, private, samplerate, unit, filename2)
v_tmp = version; Matlab_Is_Newer_Than_6 = str2num(v_tmp(1)) > 6; clear v_tmp;

try				% try to save the data into the given file
	if Matlab_Is_Newer_Than_6	
		save(filename, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit', '-v6');
	else	
		save(filename, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit');
	end
catch				% if that fails (e.g. file is locked), write to filename2
	if Matlab_Is_Newer_Than_6	
		save(filename2, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit', '-v6');
	else	
		save(filename2, 'data', 'comment', 'descriptor', 'dimension', 'private', 'samplerate', 'unit');
	end
    error([lasterr ' saved as ' filename2]);
end
    
function refreshDisplay(i, msg)
	clc; disp(['Three-dimensional Artikulographic Position and Align Determination Version ' tapadversion]); disp(' ');
	disp(msg);
	disp([num2str(i), ' points done']);
return