function TVec = rotat(AVec);
% RotAT     rotates the Nx3 matrix AVec with N vectors
%           from the Autokal into the TAPAD co-ordinate system. 
%           The 'Autokal' coordinate system is a rough approximation for the 
%           skull-fixed coordinate system. 
%           Notice that the function does not rescale the vectors,
%           so to transform position-vectors use:
%
%           TAPAD_Positions = 1/1000 * RotAT(Autokal_Positions);     
%
%           see RotTA.m
 
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

persistent RotAT;		% Matrix for coordinate transformation Autokal -> TAPAD
if (isempty(RotAT)) 
   AT = [ -1/sqrt(3)  -1/sqrt(2)  -1/sqrt(6);...
           1/sqrt(3)  -1/sqrt(2)   1/sqrt(6);...
          -1/sqrt(3)     0         2/sqrt(6)]; 
end
TVec = (AT * AVec')';								
