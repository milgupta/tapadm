function xml_text = frametext(tag, content, attributes);
% FRAMETEXT frames text in an xml conform manner.
%           The function parenthesizes a given text xml-like with 
%           a Start and End-Tag. A optional attribute list can be used, 
%           which has to be a N x 2 cell-array, containing attribute 
%           Name/Value pairs. The content-text can be omitted to form
%           an Empty-Element Tag.     
%     
%           xml_text = frametext(tag, content, attributes);
%
%           This function is usefull, when dealing with log-data, or when
%           adding lots of comments from different sources to some data.
%                
%   Examples
%           >> frametext('tapad', 'This is my text', {'date', 'now'; 'place', 'here'})
%                
%           ans =
%                
%           <tapad date="now" place="here">
%           This is my text
%           </tapad>
%
%           >> frametext('tapad', [], {'date', 'now'; 'place', 'here'})
%
%           ans =
%
%           <tapad date="now" place="here" />
%
%   Link
%           Please notice, that the XML-Specification requires the keeping 
%           of some constraints to form valid xml, which are not covered by 
%           this function. For more information, see
%
%           http://www.w3.org/TR/2004/REC-xml-20040204/#sec-starttags

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

errormsg = nargchk(1, 3, nargin);
if (~isempty(errormsg))
	error(errormsg);
end

has_content = (nargin >= 2) & (~isempty(content)) & (length(content) > 0);
has_attrbutes = (nargin >= 3);

if (has_attrbutes & ((iscellstr(attributes) ~= 1) | (size(attributes, 2) ~= 2)))
   warning('Attributes ignored (must be a Nx2 cell array).');
	has_attrbutes = logical(0);
end   

xml_text = ['<' tag];

if (has_attrbutes)
	for (i = 1 : size(attributes, 1))
		xml_text = [xml_text ' ' char(attributes(i, 1)) '="' char(attributes(i, 2)) '"'];	
	end
end

if (has_content)
	xml_text = sprintf('%s>\n%s\n</%s>', xml_text, content, tag);	
else
	xml_text = [xml_text ' />'];	
end
