function result = firfilter(data, FilterCoeff)
% FIRFILTER performs smart, transient reduced FIR filtering of the given data.
%           Filtering is performed along the 1st dimension (column) of data, 
%           which can be multi-dimensional.
%
%           result = firfilter(data, FilterCoeff);
%
%	See Also FILTER FILTFILT

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

DataSize = size(data); NumSamples = DataSize(1); DataSize(1) = 1;
Dimensions = length(DataSize); FilterLength = length(FilterCoeff);
subs = '';
for (i=2:Dimensions)
	subs = [subs ', :'];	
end

% If th filter order ist even, make it odd to simplify the handling of the 
% mirrored sequences and allow for group delay.
if rem(FilterLength, 2) == 0		
	FilterCoeff= [0 FilterCoeff];
	FilterLength = FilterLength + 1;   
end;

MirrorLength = round((FilterLength-1) / 2);

if NumSamples <= MirrorLength
   error('Input sequence is too short');
end;

% To eliminate transients, the data is extended at both ends.
% A reflection about a point is used to perform this data extension
% in a 'natural' way, i.e. a short sequence is reversed not only in time,
% but also it's inverted with respect to the amplitude of the sequence start. 

preidx = (MirrorLength+1) : -1 : 2;		% index the first MirrorLength-samples in reverse order 
eval(['predata = 2*data(1) - data(preidx' subs ');']);		% and invert amplitudes

postidx = (NumSamples-1) : -1 : (NumSamples-MirrorLength); % index the last MirrorLength-samples in reverse order 
eval(['postdata = 2*data(end) - data(postidx' subs ');']);	% and invert amplitudes

exdata = cat(1, predata, data, postdata);
result = filter(FilterCoeff, 1, exdata, [], 1);
eval(['result = result((MirrorLength+1 : MirrorLength+NumSamples)' subs ');']);