% fid2niftihdr.m
%
%        $Id$ 
%      usage: hdr = fid2niftihdr(fidname,<verbose>,<movepro=0>)
%         by: justin gardner
%       date: 12/09/09
%    purpose: strips out of fid2nifti the functionality to create nifti headers from procpar. Fidname
%             is either the name of a fid directory *or* can be a procpar structure returned by readprocpar.
%             
%
function [hdr info] = fid2niftihdr(fidname,verbose,varargin)

hdr = [];info = [];
if ieNotDefined('verbose'),verbose=1;end

% parse arguments
movepro=0;loadArgs = [];movepss=0;
validArgs = {'movepro','movepss'};
getArgs(varargin,{validArgs{:},'loadArgs'});
% now evaluate the load args (these are passed from mlrImageLoad)
getArgs(loadArgs,validArgs,'suppressUnknownArgMessage=1');

% check arguments
if nargin == 0
  help fid2niftihdr
  return
end

% set fid extension
if ~isstruct(fidname)
  fidname = setext(fidname,'fid',0);
end

% create an empty header
hdr = cbiCreateNiftiHeader;

% get the qform
[qform44 info] = fid2xform(fidname,verbose,'movepro',movepro,'movepss',movepss);
if isempty(qform44),hdr = [];return;end

% give warning for epi's that have not been processed
if info.isepi && info.compressedFid
  disp(sprintf('(fid2niftihdr) Compressed EPI - needs epibsi processing'));
  hdr = [];
  return
end

% set the qform
hdr = cbiSetNiftiQform(hdr,qform44);

% now set dimensions in header
nDims = length(info.dim);
hdr.dim(1:nDims+1) = [nDims info.dim];

% set voxel dimensions 
hdr.pixdim(2:4) = info.voxsize(1:3);

% set the tr
hdr.pixdim(5) = info.tr*1000;
hdr.xyzt_units = 18;




