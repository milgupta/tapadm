function [result, FilterDesc] = decifix(data, FilterType)
% DECIFIX is a convenience function stub to easily filter amplitude data.
%           All filters are of FIR lowpass type with 75 taps and 
%           60 dB dampling at start of the stopband.
%   FilterType (character string)
%           'FIR2030' is designed for articulators. It has passband edge 20 Hz, 
%              stopband edge 30 Hz and allows domsampling by a factor 3
%              with respect to the sampling theorem
%
%           'FIR0515' is designed for for reference sensors, with: 
%              passband edge 5 Hz, stopband edge 15 Hz, allows domsampling factor 6
%
%           'FIR0512' is designed for downsampling by factor 8
%
%           [result, FilterDesc] = decifix(data, FilterType);
%
%	See Also FIRFILTER DECIMATE 

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
if (nargin < 2)		
	FilterType = 'FIR2030';
end

persistent FIR2030;		
if (isempty(FIR2030)) 
	FIR2030 = load(fullfile('private', 'fir_20_30.mat')); 	
end

persistent FIR0515;		
if (isempty(FIR0515)) 
	FIR0515 = load(fullfile('private', 'fir_05_15.mat')); 	
end

persistent FIR0512;		
if (isempty(FIR0512)) 
	FIR0515 = load(fullfile('private', 'fir_05_12.mat')); 	
end

if (strcmpi(FilterType, 'FIR0515'))
	result = firfilter(data, FIR0515.data);
	FilterDesc = sprintf('%s\n%s\n%s', '75 taps FIR lowpass filter for ref. sensors',...
		'passband edge 5 Hz, stopband edge 15 Hz, 60 dB damping at start of stopband',...
		'allows downsampling by a factor 6 with respect to the sampling theorem');
elseif (strcmpi(FilterType, 'FIR0512'))
	result = firfilter(data, FIR0512.data);
	FilterDesc = sprintf('%s\n%s\n%s', '75 taps FIR lowpass filter for downsampling',...
		'passband edge 5 Hz, stopband edge 12.5 Hz, 60 dB damping at start of stopband',...
		'allows downsampling by a factor 8 with respect to the sampling theorem');
else
	result = firfilter(data, FIR2030.data);
	FilterDesc = sprintf('%s\n%s\n%s', '75 taps FIR lowpass filter for articulators',...
		'passband edge 20 Hz, stopband edge 30 Hz, 60 dB damping at start of stopband',...
		'allows downampling by a factor 3 with respect to the sampling theorem');
end
	