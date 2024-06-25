function showresult(DataFileName, ChanIdx, PointIdx, Styles);
% SHOWRESULT plots single- or multi-channel position data.
%
%            ShowResult(DataFileName, ChanIdx, PointIdx, Styles);

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

if (nargin < 4) 
   	Styles = ['r-'; 'g-'; 'b-'; 'm-'; 'y-'; 'c-'; 'r:'; 'g:'; 'b:'; 'm:'; 'y:'; 'c:'];
end   
if ((nargin < 2) | isempty(ChanIdx))
   	ChanIdx = 1:12;
end   

P = loaddata(DataFileName);
if ((nargin < 3) | isempty(PointIdx)) 
   	PointIdx = 1:size(P, 1);
end   

MaxP = max(max(max(abs(P(PointIdx, 1:3, ChanIdx)))));
MaxP = ceil(MaxP / 10)*10;

drapefig(MaxP);
fh = gcf;

t = title(DataFileName); set(t, 'Interpreter', 'none');

chanel = 1; style = 'b-'; 
for i_chan = 1:length(ChanIdx)
	chanel = ChanIdx(i_chan);
	plotsc(P, chanel, PointIdx, fh, Styles(ChanIdx(i_chan), 1:2));
end	

rotate3D on;									% enable interactive rotation of 3-d plots

