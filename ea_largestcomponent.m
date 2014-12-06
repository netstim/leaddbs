function slicebw=ea_largestcomponent(slicebw)

stats=conncomp(slicebw);



if stats.NumObjects>1
    maxlen=0;
    for i=1:stats.NumObjects
    if length(stats.PixelIdxList{i})>maxlen
        maxlen=length(stats.PixelIdxList{i});
        biggestobj=i;
    end
    end
    
    slicebw(:)=0;
    slicebw(stats.PixelIdxList{biggestobj})=1;
end
    


function cc = conncomp(A)
% conncomp         Find connected components in 2D grayscale or label image
%   Four-connected neighborhood is used.
%   
%   CC = CONNCOMP(A) returns the connected components CC found in A. 
%   A is a two dimensional array (image) of sixe I-by-J. Elements in A should
%   be positive, if A(i,j)==0 it is handled as background and not assigned to a
%   component, and cc.Labels(i,j) = 0. class(A) should be uint8, uint16 or double.
%
%   CC is a structure with the following fields:
% 
%      Connectivity   is always 4, two-dimensional four-connected neighborhood
% 
%      ImageSize      Size of A, i.e. [I, J] = size(A).
% 
%      NumObjects     Number of connected components (objects) in A.
% 
%      PixelIdxList   1-by-NumObjects cell array where the kth element
%                     in the cell array is a vector containing the linear
%                     indices of the pixels in the kth object.
%     
%      NumPixels      1-by-NumObjects array where the kth element is the
%                     number of pixels in the kth object.
%     
%      LenBorder      1-by-NumObjects array where the kth element is the
%                     length of the border (to other components) for the kth object 
%
%      OutBorder      1-by-NumObjects array where the kth element is the
%                     length of the border (to background) for the kth object 
%
%      Labels         I-by-J label image
%
%   Example 1
%   ---------
%   A = [1,2,2,2,2,2; 1,1,2,2,2,2; 1,1,1,2,3,2; 1,1,1,3,3,3; 1,1,1,3,3,3];
%   cc = conncomp(A)
%   A([1,2,29,30])=0;
%   cc = conncomp(A)
% 
%   See also bwconncomp  (in Matlab images toolbox)

%   or grayconncomp (an older more complicated version made by K. Skretting, 
%   and which is usually faster). 

%----------------------------------------------------------------------
% Copyright (c) 2012.  Karl Skretting.  All rights reserved.
% University of Stavanger (Stavanger University), Signal Processing Group
% Mail:  karl.skretting@uis.no   Homepage:  http://www.ux.uis.no/~karlsk/
% 
% HISTORY:  dd.mm.yyyy
% Ver. 1.0  11.09.2012  Made function 
% Ver. 1.1  18.09.2012  Added OutBorder field
% Ver. 2.0  28.09.2012  Included conncomp_mex (to speed up the program)
%----------------------------------------------------------------------

mfile = 'conncomp';
[I,J] = size(A); IJ = I*J;

if (exist('conncomp_mex','file') == 3)  % the mex file exist
    
    [nObj, pixList, idxNewObj, inBord, outBord] = conncomp_mex(A, [I,J]);
    % these are with c-indexing (starting with 0), 
    % so this is adjusted and returned
    cc = struct('Connectivity', 4, ...
                'ImageSize', [I,J], ...
                'NumObjects', double(nObj), ... 
                'PixelIdxList', {cell(nObj,1)}, ...     
                'NumPixels', zeros(nObj,1), ...  
                'LenBorder', double(inBord(1:nObj)), ...
                'OutBorder', double(outBord(1:nObj)), ...
                'Labels', zeros(I,J) );
    
    i1 = 1;            
    for i = 1:nObj
        i2 = idxNewObj(i)+1;
        cc.PixelIdxList{i} = double(pixList(i1:i2)+1);
        cc.NumPixels(i) = i2-i1+1;
        i1 = i2 + 1;
        cc.Labels(cc.PixelIdxList{i}) = i;
    end
    
else

cc = struct('Connectivity', 4, ...
            'ImageSize', [I,J], ...
            'NumObjects', 0, ...           % fill in this in the end
            'PixelIdxList', [], ...        % fill in this in the end
            'NumPixels', zeros(IJ,1), ...  % but these as processing goes on
            'LenBorder', zeros(IJ,1), ...
            'OutBorder', zeros(IJ,1), ...
            'Labels', zeros(I,J) );

% an alternative to function n4 = n4list(n,I,J)
% use more memory, and does not increase speed
% n = repmat((1:I)',1,J) + repmat(I*(0:(J-1)),I,1);
% n41 = [zeros(I,1),n(:,2:J)-I];
% n42 = [zeros(1,J); n(2:I,:)-1];
% n43 = [n(1:(I-1),:)+1; zeros(1,J)];
% n44 = [n(:,1:(J-1))+I, zeros(I,1)];
% clear n

% stack is implemented as a simple array here, not a class
lenStack = IJ; % 4000;  % IJ is large enough!
incStack = 2000;  
stack = zeros(IJ,1);  
idxStack = 0;     
maxStack = 0;
idxPlist = zeros(IJ,1);    % idxPlist(k) is index in plist for last pixel in component k
plist = zeros(IJ,1);
count = 0;  % number of labeled pixels
p1 = 0;     % first pixel in previous (or current) component
           
k = 0;     % current component
while 1   % (p <= IJ)
    
    % find next (first) pixel (in new component) to process, 
    p = p1 + 1;
    if (p > IJ); break; end;     % all are processed
    while ((cc.Labels(p) > 0) || (A(p) == 0))
        p = p + 1;
        if (p > IJ); break; end;     % all are processed
    end
    if (p > IJ); break; end;     % all are processed
    
    % p is current pixel, here first pixel in new component
    p1 = p;
    pv = A(p);    % current pixel value
    k = k + 1;    % current component
    cc.Labels(p) = k;
    
    while 1  % process this component pixel by pixel
        count = count + 1;    
        plist(count) = p;
        cc.NumPixels(k) = cc.NumPixels(k)+1;
        % disp(['   count=',int2str(count),', p=',int2str(p),...
        %       ', pv=',int2str(pv),', k=',int2str(k)]);
        
        % add to borders and push neighbors in this component to the stack
        for n4 = n4list(p,I,J)
        % for n4 = [n41(p),n42(p),n43(p),n44(p)]  % is not faster
            if n4  % > 0
                if (A(n4) == pv) && (cc.Labels(n4) == 0)
                    % push n4 to stack
                    idxStack = idxStack + 1;
                    maxStack = max(maxStack, idxStack);
                    if (idxStack > lenStack)
                        stack = [stack; zeros(incStack,1)]; %#ok<AGROW>
                        lenStack = numel(stack);
                    end
                    stack(idxStack) = n4;
                    cc.Labels(n4) = k;   % label it right away
                end
                if (A(n4) ~= pv)
                    if A(n4) %    ~= 0
                        cc.LenBorder(k) = cc.LenBorder(k)+1;
                    else
                        cc.OutBorder(k) = cc.OutBorder(k)+1;
                    end
                end
            else  % (n4 == 0)
                cc.OutBorder(k) = cc.OutBorder(k)+1;
            end
        end
        
        if (idxStack == 0); break;  end
        
        % next pixel
        p = stack(idxStack);
        idxStack = idxStack - 1;
        
    end   % loop for this component
    
    % finishing the component
    idxPlist(k) = count;    % last pixel for component k in plist
    
end

cc.NumObjects = k;
cc.NumPixels = cc.NumPixels(1:k);
cc.LenBorder = cc.LenBorder(1:k);
cc.OutBorder = cc.OutBorder(1:k);
cc.PixelIdxList = cell(k,1);
n = 1;
for k=1:cc.NumObjects
    cc.PixelIdxList{k} = plist(n:idxPlist(k));
    %  cc.PixelIdxList{k} = sort(plist(n:idxPlist(k)));  % ordered perhaps ??
    n = idxPlist(k)+1;
end

disp([mfile,': (m-file) found ',int2str(cc.NumObjects),' components in A (',...
    int2str(I),'-by-',int2str(J),'), count = ',int2str(count), ...
    ', maxStack = ',int2str(maxStack)]);

end

return

function n4 = n4list(n,I,J)
n4 = [n-I, n-1, n+1, n+I];
if (n <= I); n4(1) = 0; end;
if (rem(n-1,I) == 0); n4(2) = 0; end;
if (rem(n,I) == 0); n4(3) = 0; end;
if (n > (J-1)*I); n4(4) = 0; end;
% n4 = nonzeros(n4);
return