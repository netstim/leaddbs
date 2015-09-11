function [varargout] = read_vista_mesh(varargin)

% READ_VISTA_MESH is implemented as mex file
%
% Use as
%   [nodes,elements,labels] = read_vista_mesh(filename);
% where
%   filename = the name of the Vista mesh (with extension .v)
%
% $Id: read_vista_mesh.m 8776 2013-11-14 09:04:48Z roboos $

error('The mex file %s is missing', [mfilename '.' mexext]);
