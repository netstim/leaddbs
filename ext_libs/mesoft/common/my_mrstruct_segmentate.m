function ret = my_mrstruct_segmentate(mrStruct,queryStr, mmin, mmax, setNum, point, retMR)
%queryStr   string which contains the name of the segmentation algorithmus
%           wanted to use.
%           possible strings:   
%           'treshold'          simple treshold value segmentation
%           'region growing2D'    region growing algorithm
%           'region growing3D'    region growing algorithm
%mmin   minimum treshold value
%mmax   maximum treshold value
%setNum number which is set. 
%point  starting point for region growing. 


%%%%% init and error check
%return value init

if ~exist('retMR') || ~isequal(mrstruct_query(mrStruct, 'sizeAy'), size(retMR))
    retMR= [];
end

ret = [];

dims = mrstruct_query(mrStruct, 'dimensions');
switch (dims)
    case 2
        dims = mrstruct_query(mrStruct, 'sizeAy');
        dim_x = dims(1);
        dim_y = dims(2);
        if isempty(retMR)
            mrSegStruct = zeros(dim_x, dim_y);
        else
            mrSegStruct = retMR;
        end
    case 3 
        dims = mrstruct_query(mrStruct, 'sizeAy');
        dim_x = dims(1);
        dim_y = dims(2);
        dim_z = dims(3);
        if isempty(retMR)
            mrSegStruct = zeros(dim_x, dim_y, dim_z);
        else
            mrSegStruct = retMR;
        end
    otherwise
        warning('mrstruct has two many dimensions!'); return;
end

%Check input arguments
if (nargin<5), warning('insufficient input arguments'); return; end
if (~isstruct(mrStruct) ||  ~ischar(queryStr) ), warning('input argument have wrong type'); return; end

if strcmp(queryStr,'region growing')
   if (nargin<6), warning('insufficient input arguments for region growing segmentation'); return; end
end

%treshold segmentation
if strcmp(queryStr,'treshold')
    temp_points = (mrStruct.dataAy >= mmin ) & (mrStruct.dataAy <= mmax);
    mrSegStruct(temp_points) = setNum;
    ret = mrSegStruct;    
end

%region growing segmentation

neighbourhoodVc= ...
    [[-1,  0,  0]; ...
     [ 1,  0,  0]; ...
     [ 0, -1,  0]; ...
     [ 0,  1,  0]; ...
     [ 0,  0, -1]; ...
     [ 0,  0,  1]];
    



if strcmp(queryStr,'region growing2D')
    dims = mrstruct_query(mrStruct, 'dimensions');
    if (dims ~= 2),  warning('wrong dimensions for 2D segmentation'); return; end
    
    dims = mrstruct_query(mrStruct, 'sizeAy');
    dim_x = dims(1);
    dim_y = dims(2);
    workmap = zeros(dim_x*dim_y, 2);
   
    counter1 = 1;
    counter2 = 2;
    workmap(1, 1) = point(1);
    workmap(1, 2) = point(2);
    
    while (counter1 ~= counter2) 
        for i=-1:1 
            for j=-1:1
                    my_x =  min(max(1, workmap(counter1, 1) + i), dim_x) ;           
                    my_y =  min(max(1, workmap(counter1, 2) + j), dim_y);
                    if ((mrStruct.dataAy(my_x, my_y)>= mmin) ...
                            && (mrStruct.dataAy(my_x, my_y)<= mmax) ...
                            && (mrSegStruct(my_x, my_y) ~= setNum) ) 
                        mrSegStruct(my_x, my_y) = setNum;
                        workmap(counter2, 1) = my_x;
                        workmap(counter2, 2) = my_y;
                        counter2 =counter2 + 1;
                    end
                
            end
        end
        counter1=counter1+1;
    end
    ret = mrSegStruct;  
end

%region growing segmentation
if strcmp(queryStr,'region growing3D')
    dims = mrstruct_query(mrStruct, 'dimensions');
    if (dims ~= 3),  warning('wrong dimensions for 3D segmentation'); return; end
    
    
    if (mrStruct.dataAy(point(1), point(2), point(3))>= mmin) && (mrStruct.dataAy(point(1), point(2), point(3))<= mmax) 
        mrSegStruct(point(1), point(2), point(3)) = setNum;
    else % trivial case. If startvoxel is not set, do nothing
        return;
    end
    
    dims = mrstruct_query(mrStruct, 'sizeAy');
    dim_x = dims(1);
    dim_y = dims(2);
    dim_z = dims(3);
    workmap = zeros(dim_x*dim_y*dim_z, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    tmp= (mrStruct.dataAy >= mmin) & (mrStruct.dataAy <= mmax);
    mask= zeros(dim_x + 2, dim_y + 2, dim_z + 2);
    mask(2:(end - 1), 2:(end - 1), 2:(end - 1))= double(tmp);

    counter= 1;
    workmap(1, :) = reshape(point, [1, length(point)]) + 1; %push

    neighborNo= size(neighbourhoodVc, 1);
    mrSegStruct(point(1), point(2), point(3))= setNum;
    while (counter > 0)
        curPos= workmap(counter, :); counter= counter - 1;  %pop
        for i= 1:neighborNo%:-1:1
            myPos= curPos + neighbourhoodVc(i, :);
            if mask(myPos(1), myPos(2), myPos(3)) == 1
                counter= counter + 1; workmap(counter, :)= myPos; %push
                mask(myPos(1), myPos(2), myPos(3))= 2;
            end
        end
%        if mod(counter, 2000) == 0
%            disp(sprintf('work in progress ... %d max= %d', counter, dim_x*dim_y*dim_z))
%        end
    end
    
    tmp=  mask(2:(end - 1), 2:(end - 1), 2:(end - 1));
    idx= tmp == 2;
    mrSegStruct(idx)= setNum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     counter1 = 1;
%     counter2 = 2;
%     workmap(1, :) = reshape(point, [1, length(point)]);
%     neighborNo= size(neighbourhoodVc, 1);
%     while (counter1 ~= counter2) 
%         for i=1:neighborNo
%             myPos =  min(dims, max([1 1 1], workmap(counter1, :) + neighbourhoodVc(i, :)));           
%             if (mrStruct.dataAy(myPos(1), myPos(2), myPos(3))>= mmin) ...
%                     & (mrStruct.dataAy(myPos(1), myPos(2), myPos(3))<= mmax) ...
%                     & (mrSegStruct(myPos(1), myPos(2), myPos(3)) ~= setNum)
%                 mrSegStruct(myPos(1), myPos(2), myPos(3)) = setNum;
%                 workmap(counter2, :) = myPos;
%                 counter2 =counter2 + 1;
%             end
%         end
%         counter1=counter1+1;
%         if mod(counter1, 2000) == 0
%             disp(sprintf('work in progress ... %d/%d max= %d', counter1, counter2, dim_x*dim_y*dim_z))
%         end
%     end
    
    ret = mrSegStruct;  
end