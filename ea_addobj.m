function ea_addobj(ht,obj,resultfig,type,options)

addht=getappdata(resultfig,'addht');
if isempty(addht)
    addht=uitoolbar(resultfig);
end

switch type
    case 'tract'
         % open dialog
        [FileName,PathName]=uigetfile('*.mat','Choose Fibertract to add to scene...',[options.root,options.patientname,filesep]);
        addobj=[PathName,FileName];
        
        
        fs=load(addobj);
            fn = fieldnames(fs);
            
            eval(sprintf('thisset = fs.%s;',fn{1}));
            if size(thisset,2)~=3;
                thisset=thisset';
            end
            fibmax=length(thisset);
            dispercent(0,'Plotting fibers');
            for fib=1:fibmax
                dispercent(fib/fibmax);
                if size(thisset{fib},1)>size(thisset{fib},2)
                    thisset{fib}=thisset{fib}';
                end
                thisset{fib}(4,:)=detcolor(thisset{fib}); % add coloring information to the 4th column.
                
                for dim=1:4
                    thisfib(dim,:)=double(interp1q([1:size(thisset{fib},2)]',thisset{fib}(dim,:)',[1:0.1:size(thisset{fib},2)]')');
                end
                addobjr(fib)=surface([thisfib(1,:);thisfib(1,:)],...
                    [thisfib(2,:);thisfib(2,:)],...
                    [thisfib(3,:);thisfib(3,:)],...
                    [thisfib(4,:);thisfib(4,:)],'facecol','no','edgecol','interp','linew',1.5);
                clear thisfib
                
            end
            dispercent(100,'end');
            
            set(addobjr(:),'EdgeAlpha',0.05);
            
            
            
            
            addbutn=uitoggletool(addht,'CData',ea_get_icn('fibers',options),'TooltipString',FileName,'OnCallback',{@atlasvisible,addobjr},'OffCallback',{@atlasinvisible,addobjr},'State','on');
       
        
    case 'roi' % atlas
        
        % open dialog
        [FileName,PathName]=uigetfile({'*.nii';'*.nii.gz'},'Choose .nii image to add to scene...',[options.root,options.patientname,filesep]);
        addobj=[PathName,FileName];
        
        if ~FileName % user pressed cancel.
           return 
        end
        
        % load nifti
        
        
        nii=load_nii_proxy(addobj,options);
                nii.img=round(nii.img);

        [xx,yy,zz]=ind2sub(size(nii.img),find(nii.img>0)); % find 3D-points that have correct value.
            
            
            if ~isempty(xx)
                
                XYZ=[xx,yy,zz]; % concatenate points to one matrix.
                
                XYZ=map_coords_proxy(XYZ,nii); % map to mm-space
         
                
            end
            
            
            bb=[0,0,0;size(nii.img)];
            
            bb=map_coords_proxy(bb,nii);
            gv=cell(3,1);
            for dim=1:3
                gv{dim}=linspace(bb(1,dim),bb(2,dim),size(nii.img,dim));
            end
            [X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});
            if options.prefs.hullsmooth
                nii.img = smooth3(nii.img,'gaussian',options.prefs.hullsmooth);
            end
            fv=isosurface(X,Y,Z,permute(nii.img,[2,1,3]),max(nii.img(:))/2);
            
            if ischar(options.prefs.hullsimplify)
                
                % get to 700 faces
                simplify=700/length(fv.faces);
                fv=reducepatch(fv,simplify);
                
            else
                if options.prefs.hullsimplify<1 && options.prefs.hullsimplify>0
                    
                    fv=reducepatch(fv,options.prefs.hullsimplify);
                elseif options.prefs.hullsimplify>1
                    simplify=options.prefs.hullsimplify/length(fv.faces);
                    fv=reducepatch(fv,simplify);
                end
            end
            
            
            
            
            % set cdata
            
         atlasc=rand*64;
         jetlist=jet;
         
                cdat=abs(repmat(atlasc,length(fv.vertices),1) ... % C-Data for surface
                +randn(length(fv.vertices),1)*2)';
            
            
            
            % show atlas.
            set(0,'CurrentFigure',resultfig); 
        addobjr=patch(fv,'CData',cdat,'FaceColor',[0.8 0.8 1.0],'facealpha',0.7,'EdgeColor','none','facelighting','phong');
                ea_spec_atlas(addobjr,'',jetlist,1);

              % add toggle button:
              
        addbutn=uitoggletool(addht,'CData',ea_get_icn('atlas',ind2rgb(round(atlasc),jetlist)),'TooltipString',FileName,'OnCallback',{@atlasvisible,addobjr},'OffCallback',{@atlasinvisible,addobjr},'State','on');



                
                
end

storeinfigure(resultfig,addobjr,FileName,type); % store rendering in figure.
axis fill


       setappdata(resultfig,'addht',addht);
 
        
        
        
        
       
function storeinfigure(resultfig,obj,name,type)
AL=getappdata(resultfig,'AL');
% AL has fields of
% ROI
% ROINAMES
% FTS
% FTSNAMES
% MENU
% GUI
% store info..

if isempty(AL) % initialize AL
    AL.FTS=cell(0);
    AL.ROI=cell(0);
    AL.FTSNAMES=cell(0);
    AL.ROINAMES=cell(0);
    AL.MENU=struct;
    AL.GUI=cell(0);
end


switch type
    case 'tract'
        AL.FTS{end+1}=obj;
        AL.FTSNAMES{end+1}=name;
    case 'roi'
        AL.ROI{end+1}=obj;
        AL.ROINAMES{end+1}=name;
end
% store in figure.
setappdata(resultfig,'AL',AL);

% build menu for toggling
try
delete(AL.MENU.MAINMENU)
end
if length(AL.FTS) % only build fibertracking menu if there is at least one fiberset.
AL.MENU.MAINMENU=uimenu(resultfig,'Label','Fibertracking');
for ft=1:length(AL.FTS)
    AL.MENU.FTMENU(ft) = uimenu(AL.MENU.MAINMENU,'Label',AL.FTSNAMES{ft});
    for roi=1:length(AL.ROI)
       AL.MENU.ROIMENU(ft,roi)=uimenu(AL.MENU.FTMENU(ft),'Label',AL.ROINAMES{roi}); 
    end
end
end
axis fill

setappdata(resultfig,'AL',AL);
    



function atlasvisible(hobj,ev,atls)
set(atls, 'Visible', 'on');
%disp([atls,'visible clicked']);

function atlasinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');
%disp([atls,'invisible clicked']);



function coords=map_coords_proxy(XYZ,V)

XYZ=[XYZ';ones(1,size(XYZ,1))];

coords=V.mat*XYZ;
coords=coords(1:3,:)';


function nii=load_nii_proxy(fname,options)

if strcmp(fname(end-2:end),'.gz')
    wasgzip=1;
    gunzip(fname);
    fname=fname(1:end-3);
else
    wasgzip=0;
end
try
nii=spm_vol(fname);

nii.img=spm_read_vols(nii);
catch
    
end

nii.hdr.dime.pixdim=nii.mat(logical(eye(4)));
if ~all(abs(nii.hdr.dime.pixdim(1:3))<=1)
    reslice_nii(fname,fname,[0.5,0.5,0.5],3);
    
nii=spm_vol(fname);
nii.img=spm_read_vols(nii);
end
if wasgzip
    delete(fname); % since gunzip makes a copy of the zipped file.
end

function indcol=detcolor(mat) % determine color based on traversing direction.

xyz=abs(diff(mat,1,2));
rgb=xyz/max(xyz(:));

rgb=[rgb,rgb(:,end)];
rgbim=zeros(1,size(rgb,2),3);
rgbim(1,:,:)=rgb';
indcol=double(rgb2ind(rgbim,jet));