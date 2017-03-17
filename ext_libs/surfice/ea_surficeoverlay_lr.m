function ea_surficeoverlay_lr(heatmap,threshs,sides,colorbar)
% heatmap can contain wildcards

if ~exist('threshs','var')
    autothresh=1;
else
    if isempty(threshs)
        autothresh=1;
    else
        autothresh=0;
    end
end

if ~exist('colorbar','var')
    colorbar='false';
end
% make sure colorbar is a string
if ~ischar(colorbar)
    if colorbar
        colorbar='true';
    else
        colorbar='false';
    end
end

if ~exist('sides','var')
    sides=1:2;
end

pth=fileparts(heatmap);
hms=dir(heatmap);
for hm=1:length(hms)
   fis{hm}=fullfile(pth,hms(hm).name); 
end
if ~autothresh
    if ~(size(threshs,1)==size(fis,1))
        threshs=repmat(threshs,size(fis,1),1);
    end
end

 sidest={'r','l'};

basedir=[ea_getearoot,'ext_libs',filesep,'surfice',filesep];
if ismac
    surfice = [basedir,'mac',filesep, 'surfice.app',filesep,'Contents',filesep,'MacOS',filesep,'surfice'];
elseif isunix
    ea_libs_helper([basedir, 'linux']);
    surfice = [basedir, 'linux',filesep,'surfice'];
elseif ispc
    surfice = [basedir, 'win',filesep,'surfice.exe'];
end

% get thresholds
if autothresh
    for fi=1:length(fis)
        
        
        
        nii=ea_load_nii(fis{fi});
        nii.img(nii.img==0)=nan;
        thisstd=ea_nanstd(nii.img(:));
        if min(nii.img(:))<0 % has negative values, calculate 4 vals
            threshs(fi,1)=thisstd;
            threshs(fi,2)=1.7*thisstd;
            threshs(fi,3)=-thisstd;
            threshs(fi,4)=-1.7*thisstd;
        else % only positive vals
            threshs(fi,1)=thisstd;
            threshs(fi,2)=1.7*thisstd;
        end
        
        
    end
end
for side=sides

    mesh = [ea_space,'surf_',sidest{side},'.stl'];
    if ~exist(mesh,'file')
        ea_gensurfice_temps;
    end
    script=['BEGIN;',...
        ' RESETDEFAULTS;'];
    for fi=1:length(fis)
            [pth,fn]=fileparts(fis{fi});
    expfn_medial=fullfile(pth,[fn,'_',sidest{side},'_med.png']);
    expfn_lateral=fullfile(pth,[fn,'_',sidest{side},'_lat.png']);
        script=[script,...
            ' MESHLOAD(''',mesh,''');',...
            ' OVERLAYLOAD(''',fis{fi},''');',...
            ' OVERLAYCOLORNAME(1, ''Red-Yellow'');',...
            ' OVERLAYMINMAX(1,',num2str(threshs(fi,1)),',',num2str(threshs(fi,2)),');']; 
        if size(threshs,2)>2
            script=[script,...
                ' OVERLAYLOAD(''',fis{fi},''');',...
                ' OVERLAYCOLORNAME(2, ''Blue-Green'');',...
                ' OVERLAYMINMAX(2,',num2str(threshs(fi,3)),',',num2str(threshs(fi,4)),');'];
        end
        script=[script,...
            ' COLORBARVISIBLE(',colorbar,');',...
            ' AZIMUTHELEVATION(',num2str(90+(180*side-1)),', 0);',...
            ' SAVEBMP(''',expfn_medial,''');',...
            ' AZIMUTHELEVATION(',num2str(270+(180*side-1)),', 0);',...
            ' SAVEBMP(''',expfn_lateral,''');'];
    end
    script=[script,...
        ' QUIT',...
        ' END.'];
    cmd=[surfice,' -S "',script,'"'];
    system(cmd);
    
    [im,map,transp]=imread(expfn_medial);
    [im,transp]=crop_img(im,transp);
    imwrite(im,expfn_medial,'Alpha',transp);
    
        [im,map,transp]=imread(expfn_lateral);
    [im,transp]=crop_img(im,transp);
    imwrite(im,expfn_lateral,'Alpha',transp);
    
end








function [img2,transp] = crop_img(img,transp)
% 
% Crop image by removing edges with homogeneous intensity
% 
% 
%USAGE
%-----
% img2 = crop_img(img)
% img2 = crop_img(img,border)
% 
% 
%INPUT
%-----
% - IMG: MxNxC matrix, where MxN is the size of the image and C is the
%   number of color layers
% - BORDER: maximum number of pixels at the borders (default: 0)
% 
% 
%OUPUT
%-----
% - IMG2: cropped image
% 
% 
%EXAMPLE
%-------
% >> img  = imread('my_pic.png');
% >> img2 = crop_img(img,0);
% >> imwrite(img2,'my_cropped_pic.png')
% 

% Guilherme Coco Beltramini (guicoco@gmail.com)
% 2013-May-29, 12:29 pm


% Input
%==========================================================================

    border = 0;



% Initialize
%==========================================================================
[MM,NN,CC] = size(img);
edge_col   = zeros(2,CC); % image edges (columns)
edge_row   = edge_col;    % image edges (rows)


% Find the edges
%==========================================================================
for cc=1:CC % loop for the colors
    
    
    % Top left corner
    %================
    
    % Find the background
    %--------------------
    img_bg = img(:,:,cc) == img(1,1,cc);
    
    % Columns
    %--------
    cols = sum(img_bg,1);
    if cols(1)==MM % first column is background
        tmp = find(diff(cols),1,'first'); % background width on the left
        if ~isempty(tmp)
            edge_col(1,cc) = tmp + 1 - border;
        else % no background
            edge_col(1,cc) = 1;
        end
    else % no background
        edge_col(1,cc) = 1;
    end
    
    % Rows
    %-----
    rows = sum(img_bg,2);
    if rows(1)==NN % first row is background
        tmp = find(diff(rows),1,'first'); % background height at the top
        if ~isempty(tmp)
            edge_row(1,cc) = tmp + 1 - border;
        else % no background
            edge_row(1,cc) = 1;
        end
    else % no background
        edge_row(1,cc) = 1;
    end
    
    
    % Bottom right corner
    %====================
    
    % Find the background
    %--------------------
    img_bg = img(:,:,cc) == img(MM,NN,cc);
    
    % Columns
    %--------
    cols = sum(img_bg,1);
    if cols(end)==MM % last column is background
        tmp = find(diff(cols),1,'last'); % background width on the right
        if ~isempty(tmp)
            edge_col(2,cc) = tmp + border;
        else % no background
            edge_col(2,cc) = NN;
        end
    else % no background
        edge_col(2,cc) = NN;
    end
    
    % Rows
    %-----
    rows = sum(img_bg,2);
    if rows(end)==NN % last row is background
        tmp = find(diff(rows),1,'last'); % background height at the bottom
        if ~isempty(tmp)
            edge_row(2,cc) = tmp + border;
        else % no background
            edge_row(2,cc) = MM;
        end
    else % no background
        edge_row(2,cc) = MM;
    end
    
    
    % Identify homogeneous color layers
    %==================================
    if edge_col(1,cc)==1 && edge_col(2,cc)==NN && ...
            edge_row(1,cc)==1 && edge_row(2,cc)==MM && ...
            ~any(any(diff(img(:,:,cc),1)))
        edge_col(:,cc) = [NN;1];
        edge_row(:,cc) = [MM;1]; % => ignore layer
    end
    
    
end


% Indices of the edges
%==========================================================================

% Columns
%--------
tmp      = min(edge_col(1,:),[],2);
edge_col = [tmp ; max(edge_col(2,:),[],2)];
if edge_col(1)<1
   edge_col(1) = 1;
end
if edge_col(2)>NN
   edge_col(2) = NN;
end

% Rows
%-----
tmp      = min(edge_row(1,:),[],2);
edge_row = [tmp ; max(edge_row(2,:),[],2)];
if edge_row(1)<1
   edge_row(1) = 1;
end
if edge_row(2)>MM
   edge_row(2) = MM;
end


% Crop the edges
%==========================================================================
img2 = zeros(edge_row(2)-edge_row(1)+1,...
    edge_col(2)-edge_col(1)+1,CC,class(img));
for cc=1:CC % loop for the colors
    img2(:,:,cc) = img(edge_row(1):edge_row(2),edge_col(1):edge_col(2),cc);
end
transp=transp(edge_row(1):edge_row(2),edge_col(1):edge_col(2));


    


