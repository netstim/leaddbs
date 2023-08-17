function icn_color=ea_get_icn(varargin)
% This simple function reads in the image files for the Icons of the
% UItoolbar of eAuto-DBS.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn
switch varargin{1}
    case 'captions'
        icn_color=ones(16);
        icn_color(round([18,26,27,34,43,50.54,55,56,59,66,70,72,75,76,77,82,83,84,86,87,89,91,92,93,102,103,...
            114,115,116,119,124,125,126,130,132,135,139,146,147,148,151,154,155,156,...
            162,167,183,178,179,180,183,184,185,187,188,189]))=0;
        icn_color=repmat(icn_color,[1,1,3]);
        
    case 'regions'
        if nargin < 2
            icn_color=rand(16,16,3);
        else
            icn_color=rand(varargin{2},varargin{2},3);
        end
            
    case 'atlas'
        if nargin < 3
            icn_color=zeros(16,16,3);
        else
            icn_color=zeros(varargin{3},varargin{3},3);
        end
        col=varargin{2};
        if ischar(col) % none
            col=[1,1,1];
        end
        icn_color(:,:,1)=col(1);
        icn_color(:,:,2)=col(2);
        icn_color(:,:,3)=col(3);
        
    otherwise
        icn_color = imread(fullfile(ea_getearoot,'icons',[varargin{1},'.png']));
        if numel(size(icn_color)) == 2
            icn_color = repmat(icn_color, [1 1 3]);
        end
end