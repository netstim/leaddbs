function fibers=FT(FA,VectorF,Roi,parameters)
% This function will perform straight forward fiber tracking from the
% whole brain using linear interpolation
%
% fibers=FT(FA,VectorF,Roi,parameters);
%
% FA: A 3D matrix with the fractional anistropy (FA) of a DTI dataset
% VectorF: A 3D matrix with the main fiber direction of each voxel
% Roi: A 3D matrix with the Roi. The fibers traced from the whole brain
%       which go through this Roi are kept.
% parameters: parameters.FiberLengthMax, parameters.FiberLengthMin
%       parameters.DeviationAngleMax, parameters.Step, parameters.Sampling 
%       and parameters.FiberTrackingComputationTreshold 
%
% fibers: a cell array, with the coordinates of all fibers through the Roi
%
% Example:
%  see FT_example.m
%
% This function is written by D.Kroon University of Twente (August 2008)

if(parameters.textdisplay), disp('Start FT function'); pause(0.1); end

if(parameters.textdisplay), disp('Filter Vector Field and FA to remove high frequency noise'); pause(0.1); end
% First filter the data to remove high frequency noise
VectorF(:,:,:,1)=imfilter(VectorF(:,:,:,1),fspecial3('gaussian',[4 4 4]));
VectorF(:,:,:,2)=imfilter(VectorF(:,:,:,2),fspecial3('gaussian',[4 4 4]));
VectorF(:,:,:,3)=imfilter(VectorF(:,:,:,3),fspecial3('gaussian',[4 4 4]));
FA=imfilter(FA,fspecial3('gaussian',[4 4 4]));

% Number of fibers found
fibern=0;

if(parameters.textdisplay), disp('Perform fiber tracking through the whole brain'); end
% Loop through the voxel locations of the volume
for x=2:parameters.Sampling:(size(Roi,1)-1),
    if(parameters.textdisplay),  disp(['fiber tracking process : ' num2str(round(100*x/(size(Roi,1)-1))) '%']); pause(0.1); end
    for y=2:parameters.Sampling:(size(Roi,2)-1),
        for z=2:parameters.Sampling:(size(Roi,3)-1),
            % fiber is used to store the coordinates of one fiber
            fiber=zeros(600,3);
            % First fiber coordinate is the current position
            fiber(1,:)=[x y z];
            
            % initialize variables
            check=true; f_length=1; f_Roi=false; old_gradient=[]; f_angle=0;
            while(check)
                % Linear interpolation of current location  
                xBas=[0 0 0 0 1 1 1 1]+floor(fiber(f_length,1));
                yBas=[0 0 1 1 0 0 1 1]+floor(fiber(f_length,2));
                zBas=[0 1 0 1 0 1 0 1]+floor(fiber(f_length,3));
                xCom=fiber(f_length,1)-floor(fiber(f_length,1)); 
                yCom=fiber(f_length,2)-floor(fiber(f_length,2)); 
                zCom=fiber(f_length,3)-floor(fiber(f_length,3)); 
                % Linear interpolation percentages.
                perc=[(1-xCom) * (1-yCom) * (1-zCom); (1-xCom) * (1-yCom) * zCom; 
                      (1-xCom) * yCom * (1-zCom); (1-xCom) * yCom * zCom; 
                       xCom * (1-yCom) * (1-zCom); xCom * (1-yCom) * zCom; 
                       xCom * yCom * (1-zCom); xCom * yCom * zCom;];

                % The gradient / fiber direction is determined from a voxel
                % neighborhood of 8.
                gradient=[0 0 0]; f_FA=0;
                for i=1:8,
                    gradient=gradient+(squeeze(VectorF(xBas(i),yBas(i),zBas(i),:))*perc(i))'; 
                    f_FA=f_FA+FA(xBas(i),yBas(i),zBas(i),:)*perc(i);
                end
                gradient=gradient./(sqrt(sum(gradient.^2))+eps);
                
                % Set a step in the direction of the gradient.
                fiber(f_length+1,:)=fiber(f_length,:)+parameters.Step*gradient;

                % Determine if this fiber crosses the Roi
                fiber_r=round(fiber(f_length,:));
                f_Roi=f_Roi||Roi(fiber_r(1),fiber_r(2),fiber_r(3));
               
                % Calculate angle between the current and last on the (fiber) position
                if(f_length>1), f_angle=abs(acos(sum(gradient.*old_gradient))); end
                
                % Do the Fiber stop checks:
                % Stop if the fiber becomes to long
                if(f_length>=parameters.FiberLengthMax), check=false; end
                % Stop if the fiber takes a hard turn.
                if(f_angle>parameters.DeviationAngleMax), check=false; end
                % Stop if the anistropy of the voxel is too low
                if(f_FA<parameters.FiberTrackingComputationTreshold), check=false; end
                % Check if the fiber will grow outside of the volume 
                if((sum(fiber(f_length+1,:)>(size(Roi)-1))+sum(fiber(f_length+1,:)<2))>0), check=false; end
                
                % Update the fiber length
                if(check), f_length=f_length+1; end
                
                % Keep the gradient for next step angle change calculation
                old_gradient=gradient;
            end
            % Keep the fiber if it is long enough and crossed the Roi
            if((f_length>parameters.FiberLengthMin)&&f_Roi)
                fibern=fibern+1;
                fibers{fibern}=fiber(1:f_length,:);
            end
        end
    end
end
if(parameters.textdisplay), disp('FT function finished'); pause(0.1); end 

function h = fspecial3(type,siz)
%FSPECIAL3 Create predefined 3-D filters.
%   H = FSPECIAL3(TYPE,SIZE) creates a 3-dimensional filter H of the
%   specified type and size. Possible values for TYPE are:
%
%     'average'    averaging filter
%     'ellipsoid'  ellipsoidal averaging filter
%     'gaussian'   Gaussian lowpass filter
%     'laplacian'  Laplacian operator
%     'log'        Laplacian of Gaussian filter
%
%   The default SIZE is [5 5 5]. If SIZE is a scalar then H is a 3D cubic
%   filter of dimension SIZE^3.
%
%   H = FSPECIAL3('average',SIZE) returns an averaging filter H of size
%   SIZE. SIZE can be a 3-element vector specifying the dimensions in
%   H or a scalar, in which case H is a cubic matrix.
%
%   H = FSPECIAL3('ellipsoid',SIZE) returns an ellipsoidal averaging filter.
%
%   H = FSPECIAL3('gaussian',SIZE) returns a centered Gaussian lowpass
%   filter of size SIZE with standard deviations defined as SIZE/2/2.354 so
%   that FWHM equals half filter size (http://en.wikipedia.org/wiki/FWHM).
%
%   H = FSPECIAL3('laplacian') returns a 3-by-3-by-3 filter approximating
%   the shape of the three-dimensional Laplacian operator. REMARK: the
%   shape of the Laplacian cannot be adjusted. An infinite number of 3D
%   Laplacian could be defined. If you know any simple formulation allowing
%   one to control the shape, please contact me.
%
%   H = FSPECIAL3('log',SIZE) returns a rotationally symmetric Laplacian of
%   Gaussian filter of size SIZE with standard deviation defined as
%   SIZE/2/2.354.
%
%
%   Class Support
%   -------------
%   H is of class double.
%
%   Example
%   -------
%      I = single(rand(80,40,20));
%      h = fspecial3('gaussian',[9 5 3]); 
%      Inew = imfilter(I,h);
%       
%   See also IMFILTER, FSPECIAL.
%
%   -- Damien Garcia -- 2007/08

type = lower(type);

if nargin==1
        siz = 5;
end

if numel(siz)==1
    siz = round(repmat(siz,1,3));
elseif numel(siz)~=3
    error('Number of elements in SIZ must be 1 or 3')
else
    siz = round(siz(:)');
end

switch type
    
    case 'average'
        h = ones(siz)/prod(siz);
        
    case 'gaussian'
        sig = siz/2/2.354;
        siz   = (siz-1)/2;
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
        h = h/sum(h(:));
        
    case 'ellipsoid'
        R = siz/2;
        h = ones(siz);
        siz = (siz-1)/2;
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        I = (x.*x/R(1)^2+y.*y/R(2)^2+z.*z/R(3)^2)>1;
        h(I) = 0;
        h = convn(h,ones(2,2,2),'same');
        h = h/sum(h(:));
        
    case 'laplacian'
        h = zeros(3,3,3);
        h(:,:,1) = [0 3 0;3 10 3;0 3 0];
        h(:,:,3) = h(:,:,1);
        h(:,:,2) = [3 10 3;10 -96 10;3 10 3];
        
    case 'log'
        sig = siz/2/2.354;
        siz   = (siz-1)/2;
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
        h = h/sum(h(:));
        arg = (x.*x/sig(1)^4 + y.*y/sig(2)^4 + z.*z/sig(3)^4 - ...
            (1/sig(1)^2 + 1/sig(2)^2 + 1/sig(3)^2));
        h = arg.*h;
        h = h-sum(h(:))/prod(2*siz+1);
        
    otherwise
        error('Unknown filter type.')
        
end
        
        
        