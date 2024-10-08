function T = loadtrial(ResultFileName, ResidualFileName, AmpFileName);
% LOADTRIAL Load all trial-related data from the three 'mat'-files:  
%           ResultFile   contains the result from a TAPAD run
%           ResidualFile contains additional TAPAD data 
%           AmpFile      contains the measured signal amplitudes
%           The residual and amp files are optional.
%     
%           T = loadtrial(ResultFileName, ResidualFileName, AmpFileName);
%
%           T = loadtrial(fullfile(BasePath, 'pos', trialfile(201)),... 
%                         fullfile(BasePath, 'pos', trialfile(201, 'resi')),... 
%                         fullfile(BasePath, 'amps', trialfile(201)));

%
%           The function returns a struct with all collected data.

%---------------------------------------------------------------------
% Copyright � 2005-2007 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
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
errormsg = nargchk(1, 3, nargin);
if (~isempty(errormsg))
	error(errormsg);
end

if (nargin<3)
	AmpFileName = [];
end
if (nargin<2)
	ResidualFileName = [];
end

% Load amplitude data -----------------------------------------------
if (exist(AmpFileName, 'file'))
	[MAmps, amp_comment, amp_descriptor, amp_dimension, amp_private, samplerate, amp_unit] = loadsdata(AmpFileName);
else
	[MAmps, amp_comment, amp_descriptor, amp_dimension, amp_private, samplerate, amp_unit] = loadsdata;
	if (~isempty(AmpFileName))
		warning('No amplitude data loaded!');
	end
end	
T.Amps = MAmps; T.Samplerate = samplerate; T.NumPoints = size(MAmps, 1);	
	
% Derivative-Parameters are already stored in the AmpFile, if it has been generated by 'prepdata'. 
% This parameters allow to roughly estimate the run of the amplitude and the gradient. 
T.DerivativParameters = [];
if (exist('amp_private') == 1)
	if	(isfield(amp_private,'DerivativParameters')==1)
 		[N, M, L] = size(amp_private.DerivativParameters);
		if ((N==16) | (M==6) | (L==12))
			T.DerivativParameters = amp_private.DerivativParameters;
		else
			warning('Existing derivative parameters discarded due to incompatible format!');
		end
	end
end


% There are two output files: one for the standard result and the other for more sophisticated 
% output like residuals, iterations and other private data. 
% One output file is called 'residual-file' and is considered to depend on the standard result-file.
% (The Carstens 'CalcPos'-Program calls the residuals 'PosAmps')
if (exist(ResidualFileName, 'file'))
	[Residuals, res_comment, res_descriptor, res_dimension, res_private] = loadsdata(ResidualFileName);
	[N, M, L] = size(Residuals);
	if (T.NumPoints == 0)					% if NumPoints is not jet initialized, adopt it from residual data
		T.NumPoints = N;
	else
		if (N~=T.NumPoints)	% if sample-count of the amp and residual file do not match, abort
			error(['Amplitude file contains ' num2str(T.NumPoints) ' samples, while residual file contains ' num2str(N)])
		end
	end
	if ((M~=6) | (L~=12))					% is the size of the residual data as expected?
		Residuals = [];						% if not, discard the data
		warning('Existing tapad residual data discarded due to incompatible format!');
	end
	[N, M, L] = size(res_private.Iterations);
	if ((N~=T.NumPoints) | (M~=1) | (L~=12)) % is the size of the Iterations array as expected?
		Residuals = [];						% if not, discard the data
		warning('Existing tapad residual data discarded due to incompatible format! (Iterations)');
	end	
	clear N M L;
else												% the residual file does not exist
	[Residuals, res_comment, res_descriptor, res_dimension, res_private] = loadsdata; % just create the variables		
end

% Now check for a result-file which contains pre-computed positions (e.g. from different channels). 
% If there is no usable result-file, the residual-file data will be nulled whether there is 
% old residual data, or not (since keeping it without results would make no sense). 
if (exist(ResultFileName, 'file'))
  [Result, comment, descriptor, dimension, private, res_samplerate, unit] = loadsdata(ResultFileName); 
else
  error(['Can''t load result file: ''' ResultFileName '''']);
end	  
[N, M, L] = size(Result);
if (T.NumPoints == 0)						% if NumPoints is not jet initialized, adopt it from result data
	T.NumPoints = N;
else
	if (N~=T.NumPoints)							% if sample-count of the result file does not match, abort
		error(['Result file contains ' num2str(N) ' samples, while other data contains ' num2str(T.NumPoints)])
	end
end

if (isempty(Result) | (M~=7) | (L~=12))
	if (~isempty(Result))	
		warning('Existing tapad data discarded due to incompatible format!');
	else
		warning('Empty result!');
	end
	[Result, comment, descriptor, dimension, private, res_samplerate, unit] = loadsdata; 
	Result = NaN * ones(T.NumPoints, 7, 12); 	comment = amp_comment;
	Residuals = [];							% missing result data invalidates any residual data
end	

if (isempty(Residuals))						% if there is no valid residual data, initialize the result variables
	Residuals = NaN * ones(T.NumPoints, 6, 12);   
	res_private.Iterations = zeros(T.NumPoints, 1, 12);
	res_private.DerivativParameters = zeros(5*3+1, 6, 12);
end	

if (isempty(T.DerivativParameters))
	T.DerivativParameters = res_private.DerivativParameters;
end

T.ExitFlag = Result(:, 7, :);
T.Result = Result;
T.Residuals = Residuals;
T.Iterations = res_private.Iterations;
T.Comment = comment;
