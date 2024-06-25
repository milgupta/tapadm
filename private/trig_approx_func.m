function [ydata, velocity, acceleration] = trig_approx_func(x, xdata);
% TRIG_APPROX_FUNC calculates ydata from xdata (ydata=f(xdata)), where
%           f is a taylor series (f = SUM cos + sin)  with parameters x.
%           The function can be used as objective e.g. with lsqcurvefit
%           to approximate a periodic function. It is particularly usefull
%           to approximate the derivative of a noisy data set.
%
%           [ydata, velocity, acceleration] = trig_approx_func(x, xdata);
%
%           To estimate a start vector, try:
%           x0 = ones(1, 16) .* mean(MyYData); x0(2:3:14) = 1:5; x0(1)=x0(1)*2;
%				
%           see lsqcurvefit, lsqnonlin (Optimization Toolbox) 


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


ORDER = 5;

ydata = x(1); 

for i=1:3:3*ORDER
	ydata = ydata + x(i+2) .* cos(x(i+1) .* xdata) + x(i+3) .* sin(x(i+1) .* xdata);
end

if (nargout > 1)
	velocity = 0; 
	for i=1:3:3*ORDER
		velocity = velocity - x(i+2) .* sin(x(i+1) .* xdata) .* x(i+1)  +  x(i+3) .* cos( x(i+1) .* xdata) .* x(i+1);
	end
end

if (nargout > 2)
	acceleration = 0; 
	for i=1:3:3*ORDER
		acceleration = acceleration - x(i+2) .* cos(x(i+1) .* xdata) .* x(i+1).^2  +  x(i+3) .* sin( x(i+1) .* xdata) .* x(i+1).^2;
	end
end


