function y = decround(x, p);
% DECROUND round towards the next decimalplace. 
%           The function rounds abs(x) up to the next decimal place, e.g.
%           decround(134) yields 200. If x is a vector, the function takes
%           the maximum (absolute) value. 
%           The optional 2nd parameter denotes the requested number of
%           decimal places, e.g. decround(134, 10) yields 140,
%           decround(134, 1) yields 134 and decround(134, 1000) 1000.
%
%           y = decround(x, p);
%
%           The function is usefull for scaling purposes, i.e. in plots.
%           It allways returns positive values. 

%---------------------------------------------------------------------
% Copyright © 2007 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
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

errormsg = nargchk(1, 2, nargin);
ax = max(abs(x)); 
rp = floor(log10(ax));	% required decimal places
if (nargin < 2)
	e = 10^(rp);
else
	if (p < 1)
		e = 10^(rp);
	else	
		dp = floor(log10(p));	% desired decimal places
		e = 10^(dp);
	end
end

y = ceil(ax/e)*e;