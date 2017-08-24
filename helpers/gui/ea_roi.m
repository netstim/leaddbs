classdef ea_roi < handle
    % ROI class to plot niftis on lead dbs resultfig / 3D Matlab figures
    properties (SetObservable)
        niftiFilename % original nifti filename
        nii % nifti loaded
        threshold % threshold to visualize
        color % color of patch
        alpha=0.7 % alpha of patch
        fv % faces and vertices of patch
        sfv % smoothed version
        visible='on' % turn on/off
        smooth % smooth by FWHM
        hullsimplify % simplify hull
        max % maxvalue in nifti - 0.1 of the way
        min % minvalue in nifti + 0.1 of the way
        controlH % handle to color / threshold control figure
        plotFigureH % handle of figure on which to plot
        patchH % handle of patch
        toggleH % toggle handle
        htH % handle for toggle toolbar
    end
    methods
        function obj=ea_roi(niftiFilename,pobj) % generator function
            if exist('niftiFilename','var') && ~isempty(niftiFilename)
                obj.niftiFilename=niftiFilename;
            end
            
           
            obj.plotFigureH=gcf;
            
            if exist('pobj','var') && ~isempty(pobj)
                try
                    obj.plotFigureH=pobj.plotFigureH;
                end
            end
            obj.htH=getappdata(obj.plotFigureH,'addht');
            if isempty(obj.htH) % first ROI
                obj.htH=uitoolbar(obj.plotFigureH);
                setappdata(obj.plotFigureH,'addht',obj.htH);
            end
            
           set(0,'CurrentFigure',obj.plotFigureH); 
            % set cdata
            if exist('pobj','var') && ~isempty(pobj)
                try
                    obj.color=pobj.color;
                catch
                    obj.color=uisetcolor;
                end
            else
                obj.color=uisetcolor;
            end
            
            
            % load nifti
            obj.nii=ea_load_nii(obj.niftiFilename);
            obj.nii.img(obj.nii.img==0)=nan;
            obj.nii.img=obj.nii.img-nanmin(obj.nii.img(:)); % set min to zero
            obj.nii.img(isnan(obj.nii.img))=0;
            if ~all(abs(obj.nii.voxsize)<=1)
                ea_reslice_nii(obj.niftiFilename,obj.niftiFilename,[0.5,0.5,0.5],0,[],3);
                obj.nii=ea_load_nii(obj.niftiFilename);
            end
            options.prefs=ea_prefs;
            obj.max=ea_nanmax(obj.nii.img(~(obj.nii.img==0)));
            obj.min=ea_nanmin(obj.nii.img(~(obj.nii.img==0)));
            maxmindiff=obj.max-obj.min;
            obj.max=obj.max-0.1*maxmindiff;
            obj.min=obj.min+0.1*maxmindiff;
            
            obj.threshold=obj.max-0.5*maxmindiff;
            obj.smooth=options.prefs.hullsmooth;
            obj.hullsimplify=options.prefs.hullsimplify;
            set(0,'CurrentFigure',obj.plotFigureH);
            obj.patchH=patch;
                        
            %menu=uicontextmenu;
            %entry1=uimenu(menu,'Label','ROI properties...','Callback',{@ea_editroi,obj});
            obj.toggleH=uitoggletool;
            

            % Get the underlying java object using findobj
            jtoggle = findjobj(obj.toggleH);
            
            % Specify a callback to be triggered on any mouse release event
            %set(jtoggle, 'MouseReleasedCallback', @(s,e)rightcallback(h,e))
            
            set(jtoggle, 'MouseReleasedCallback', {@rightcallback,obj})            
            update_roi(obj);
            addlistener(obj,'visible','PostSet',...
                @ea_roi.changeevent);
            addlistener(obj,'color','PostSet',...
                @ea_roi.changeevent);
            addlistener(obj,'threshold','PostSet',...
                @ea_roi.changeevent);
            addlistener(obj,'smooth','PostSet',...
                @ea_roi.changeevent);
            addlistener(obj,'hullsimplify','PostSet',...
                @ea_roi.changeevent);
            addlistener(obj,'alpha','PostSet',...
                @ea_roi.changeevent);
            
            if exist('pobj','var') && isfield(pobj,'openedit') && pobj.openedit
                    ea_editroi([],[],obj)
            end
            
        end
        
        function changeevent(~,event)
            update_roi(event.AffectedObject,event.Source.Name);
        end

        function obj=update_roi(obj,evtnm) % update ROI
            if ~exist('evtnm','var')
                evtnm='all';
            end
            if ismember(evtnm,{'all','threshold','smooth','hullsimplify'}) % need to recalc fv here:
                [xx,yy,zz]=ind2sub(size(obj.nii.img),find(obj.nii.img>0)); %(mean(obj.nii.img(obj.nii.img~=0))/3))); % find 3D-points that have correct value.
                if ~isempty(xx)
                    XYZ=[xx,yy,zz]; % concatenate points to one matrix.
                    XYZ=map_coords_proxy(XYZ,obj.nii); % map to mm-space
                end
                bb=[0,0,0;size(obj.nii.img)];
                bb=map_coords_proxy(bb,obj.nii);
                gv=cell(3,1);
                for dim=1:3
                    gv{dim}=linspace(bb(1,dim),bb(2,dim),size(obj.nii.img,dim));
                end
                [X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});
                %if obj.smooth
                %    obj.nii.simg = smooth3(obj.nii.img,'gaussian',obj.smooth);
                %else
                %    obj.nii.simg=obj.nii.img;
                %end
                
                obj.fv=isosurface(X,Y,Z,permute(obj.nii.img,[2,1,3]),obj.threshold);
                fvc=isocaps(X,Y,Z,permute(obj.nii.img,[2,1,3]),obj.threshold);
                obj.fv.faces=[obj.fv.faces;fvc.faces+size(obj.fv.vertices,1)];
                obj.fv.vertices=[obj.fv.vertices;fvc.vertices];
                
                if obj.smooth
                    obj.sfv=ea_smoothpatch(obj.fv,1,obj.smooth);
                else
                    obj.sfv=obj.fv;
                end
                
                if ischar(obj.hullsimplify)
                    
                    % get to 700 faces
                    simplify=700/length(obj.sfv.faces);
                    obj.sfv=reducepatch(obj.sfv,simplify);
                    
                else
                    if obj.hullsimplify<1 && obj.hullsimplify>0
                        
                        obj.sfv=reducepatch(obj.sfv,obj.hullsimplify);
                    elseif obj.hullsimplify>1
                        simplify=obj.hullsimplify/length(obj.fv.faces);
                        obj.sfv=reducepatch(obj.sfv,simplify);
                    end
                end
            end
            
            %?atlasc=59; %rand*64;
            jetlist=jet;
            
            co=ones(1,1,3);
            co(1,1,:)=obj.color;
            atlasc=double(rgb2ind(co,jetlist));
            
            cdat=abs(repmat(atlasc,length(obj.sfv.vertices),1) ... % C-Data for surface
                +randn(length(obj.sfv.vertices),1)*2)';
            
            % show atlas.
            set(0,'CurrentFigure',obj.plotFigureH);
            set(obj.patchH,...
                {'Faces','Vertices','CData','FaceColor','FaceAlpha','EdgeColor','FaceLighting','Visible'},...
                {obj.sfv.faces,obj.sfv.vertices,cdat,obj.color,obj.alpha,'none','phong',obj.visible});
            
            % add toggle button:
            set(obj.toggleH,...
                {'Parent','CData','TooltipString','OnCallback','OffCallback','State'},...
                {obj.htH,ea_get_icn('atlas',obj.color),stripext(obj.niftiFilename),{@ea_roivisible,'on',obj},{@ea_roivisible,'off',obj},obj.visible});
        end
    end
    
end



function rightcallback(src, evnt,obj)
    if evnt.getButton() == 3
        ea_editroi(src,evnt,obj)
    end
end

function ea_editroi(Hobj,evt,obj)
            obj.controlH=ea_roicontrol(obj);

end

function ea_roivisible(Hobj,evt,onoff,obj)

%try % if not yet initialized wont work
    obj.visible=onoff;
%end
end
function coords=map_coords_proxy(XYZ,V)

XYZ=[XYZ';ones(1,size(XYZ,1))];

coords=V.mat*XYZ;
coords=coords(1:3,:)';
end

function fn=stripext(fn)
[~,fn]=fileparts(fn);
end