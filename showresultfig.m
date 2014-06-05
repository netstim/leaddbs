function [resultfig,coords_mm]=showresultfig(coords_mm,realcoords,fitline,cornii,patientname,options)



%% plot coords

if options.verbose>1; resultfig=figure('name',[patientname,': Auto-Manual-Coherence']); end
if options.verbose>1; pltcrds=plot3(coords_mm(:,1),coords_mm(:,2),coords_mm(:,3),'o','MarkerFaceColor',[0.5,0.5,0.5],'MarkerSize',15); end

if options.verbose>1; hold on; end
try
    if ~isempty(realcoords)
    if options.verbose>1; realcoords_plot=plot3(realcoords(:,1),realcoords(:,2),realcoords(:,3),'*y'); end
    end
end

%% plot lines

for side=options.sides
    
    try
        
        if ~isempty(fitline{side})
        fitline{side}=map_coords(fitline{side}', [options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii'])';
        end
        if options.verbose>1; trajectory_plot=plot3(fitline{side}(:,1),fitline{side}(:,2),fitline{side}(:,3),'color',[0.5,0.5,0.5],'linew',1.5); end
        
    end
end

planecnt=1;
%% plot slices in x and y planes

if options.verbose>1;
    for doxx=0:1
        for side=options.sides
          try  
            sample_width=20;
            meanfitline=genhd_inside(fitline{side});
            clear imat
            %% sample plane left and right from meanfitline
            V=spm_vol(cornii);
            cnt=1;
            for xx=-sample_width:0.1:sample_width
                
                %    map from world to voxel coords (special case)
                %    [XYZ_mm XYZ_vx] = map_coords(XYZ_mm, img); % (XYZ_mm unaltered)
                if doxx % sample in x direction

                    [dummy,thisfitline_vx]=map_coords([meanfitline+[repmat(xx,length(meanfitline),1),zeros(length(meanfitline),2)]]',cornii);
                else % sample in y direction

                    [dummy,thisfitline_vx]=map_coords([meanfitline+[zeros(length(meanfitline),1),repmat(xx,length(meanfitline),1),zeros(length(meanfitline),1)]]',cornii);
                end
                thisfitline_vx=thisfitline_vx';
                    imat(:,cnt)=spm_sample_vol(V,double(thisfitline_vx(:,1)),double(thisfitline_vx(:,2)),double(thisfitline_vx(:,3)),1);

                cnt=cnt+1;
            end
            
            colormap gray
            
            
            if doxx % span surface in x direction
                spanvector=[sample_width,0,0];
            else % span surface in y direction
                spanvector=[0,sample_width,0];
            end
            
            boundingbox=[meanfitline(1,:)-spanvector;...
                meanfitline(1,:)+spanvector;...
                meanfitline(end,:)-spanvector;...
                meanfitline(end,:)+spanvector];
            
            
            xx=[boundingbox(1,1),boundingbox(2,1);boundingbox(3,1),boundingbox(4,1)];
            yy=[boundingbox(1,2),boundingbox(2,2);boundingbox(3,2),boundingbox(4,2)];
            zz=[boundingbox(1,3),boundingbox(2,3);boundingbox(3,3),boundingbox(4,3)];
            
            
            alphamap=imat;
            alphamap(:)=0.7;
            
            
            %imat=repmat(imat,[1,1,3]);
           
%            [meshx,meshy]=meshgrid(
            
            planes(planecnt)=surface('XData',xx,'YData',yy,'ZData',zz,'CData',imat,'alphadata',alphamap,'FaceAlpha', 'texturemap','FaceColor','texturemap','EdgeColor','none','alphadatamapping','none');
            planecnt=planecnt+1;
        end
        end
    end
end

%% Manual height correction here:
if options.manualheightcorrection
    set(resultfig,'Position',[10 400 700 500])
    view(0,0);
    
    disp('Manual correction: Press arrows to adjust, space to end adjustment.');
    
    while 1
        
        pause
        %if k % button press, not mouse click.
        commnd=get (gcf, 'CurrentKey');
        if strcmp(commnd,'space')
            break
        elseif strcmp(commnd,'x') || strcmp(commnd,'a')
            view(0,0);
        elseif strcmp(commnd,'p')
            view(180,0);
        elseif strcmp(commnd,'y') || strcmp(commnd,'l')
            view(90,0);
        elseif strcmp(commnd,'r')
            view(270,0);
            
        else
            coords_mm=correctcoords(coords_mm,fitline,side,commnd);
        end
        %end
        set(pltcrds,'XData',coords_mm(:,1),'YData',coords_mm(:,2),'ZData',coords_mm(:,3))
        refreshdata(pltcrds,'caller')
        drawnow
        
    end
    
    
    
    
    disp('Manual correction done.');

end

try
sample_cuboid_qm(coords_mm(1:4,:),fitline{1},patientname,options,'r')
end
try
sample_cuboid_qm(coords_mm(5:8,:),fitline{2},patientname,options,'l')
end

try
delete(realcoords_plot);
end
delete(planes);

set(gcf,'Renderer','OpenGL')

if options.showatlases
    
    atlassurfs=showatlas(resultfig);
    
end




if options.showfibres

    fib_plots=showfibres(resultfig,coords_mm,options);
    
end

axis equal
axis vis3d
axis off
set(gcf,'color','w');
axis tight

function hdfitline=genhd_inside(fitline)

resolution=20;

hdfitline(:,1)=interp1q([1:length(fitline)]',fitline(:,1),[1:1/resolution:length(fitline)]');
hdfitline(:,2)=interp1q([1:length(fitline)]',fitline(:,2),[1:1/resolution:length(fitline)]');
hdfitline(:,3)=interp1q([1:length(fitline)]',fitline(:,3),[1:1/resolution:length(fitline)]');