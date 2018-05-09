function [fibers,idx] = ea_trk2ftr(trk_in,options)
% transforms .trk to fibers.
% CAVE: fiber information in the trk needs to conform to the MNI space!
%% basis information
template_in = [options.earoot 'templates' filesep 'space' filesep options.prefs.machine.space filesep options.primarytemplate '.nii'];

%% read .trk file
[header,tracks]=ea_trk_read(trk_in);

%% create variables needed
ea_fibformat='1.0';
fourindex=1;
idx=[tracks(:).nPoints]';
idx2 =cumsum(idx);
fibers = NaN(sum(idx),4);

start = 1;
fibercount = numel(idx);
dispercent(0,'Converting trk to fibers.');
for a=1:fibercount
    dispercent(a/fibercount);
    stop = idx2(a);
    fibers(start:stop,1:3)=tracks(a).matrix;
    fibers(start:stop,4)=a;
    start = stop+1;
end
dispercent(100,'end');
fibers=double(fibers);

fibers = fibers';

%% transform fibers to template origin
nii=ea_load_nii(template_in);
nii.mat(1,1)=1;
nii.mat(2,2)=1;
nii.mat(3,3)=1;

fibers_origin = fibers + [nii.mat(1,4) nii.mat(2,4) nii.mat(3,4) 0]';

fibers_final=fibers_origin;
fibers = fibers_final';


function  dispercent(varargin)
percent=round(varargin{1}*100);
if nargin==2
    if strcmp(varargin{2},'end')
        fprintf('\n')
        fprintf('\n')        
        fprintf('\n')        
    else
        fprintf(1,[varargin{2},':     ']);
    end
else
    fprintf(1,[repmat('\b',1,(length(num2str(percent))+1)),'%d','%%'],percent);
end
