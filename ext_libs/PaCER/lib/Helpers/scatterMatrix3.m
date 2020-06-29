%% scatterSpheres - plot scattered points 
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function h = scatterMatrix3(pointcloud, options, varargin)

    fixedArgin = {}; %#ok<NASGU>
    if(nargin <= 2)
        fixedArgin = {10};
        options = {}
        
    else
        fixedArgin(1) = varargin(1);
        varargin(1) = [];
    end
    
    
    inParse  = inputParser;
    inParse.KeepUnmatched = true;
    inParse.addParameter('numbers',false, @islogical);
    inParse.parse(options{:});


    
    nm = size(pointcloud);
    if(nm(2) > nm(1) && max(nm) > 3)
        pointcloud = pointcloud';
    end
    
    
    h = scatter3(pointcloud(:,1),pointcloud(:,2), pointcloud(:,3),fixedArgin{:}, varargin{:});
    if(inParse.Results.numbers)
        text(pointcloud(:,1),pointcloud(:,2), pointcloud(:,3),cellstr(num2str((1:size(pointcloud,1))')));
    end
end