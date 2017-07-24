function [o]=sumsq(x,DIM)
% SUMSQ calculates the sum of squares.
% 
% [y] = sumsq(x [,  DIM])
% 
% DIM	dimension
% 	N STD of  N-th dimension 
%	default or []: first DIMENSION, with more than 1 element
%
% y	estimated standard deviation
%
% features:
% - can deal with NaN's (missing values)
% - dimension argument also in Octave
% - compatible to Matlab and Octave
%
% see also: RMS, SUMSKIPNAN, MEAN, VAR, MEANSQ,
%
%
% References(s):


%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; If not, see <http://www.gnu.org/licenses/>.

%	$Id$
%	Copyright (C) 2009,2010 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

if nargin<2,
	DIM = []; 
end;
if isempty(DIM), 
        DIM = find(size(x)>1,1);
        if isempty(DIM), DIM=1; end;
end;

[s,n,o] = sumskipnan(x,DIM);

