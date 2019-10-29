function ea_setatlascline(handles)


contour=getappdata(handles.checkstructures,'contour');

linefiducial=getappdata(handles.checkstructures,'linefiducial');
fiducialview=getappdata(handles.checkstructures,'fiducialview');
bbs=getappdata(handles.checkstructures,'bbs'); % bounding boxes of views
thisbb=bbs(fiducialview);
%linefiducial(:,1)=thisbb.imgdim(1)-linefiducial(:,1);
contour=contour{fiducialview};

[planedim,onedim,secdim]=ea_getdims(fiducialview,1);

drawline=linefiducial(:,[onedim,secdim]);

% map contour to image coords:

expvx=zeros(size(contour'));

expvx(:,1)=(contour(1,:)-smallestentry([thisbb.mm{onedim}(1),thisbb.mm{onedim}(end)]))...
    /abs(diff([thisbb.mm{onedim}(1),thisbb.mm{onedim}(end)]))... % scale from 0 to 1
    *thisbb.imgdim(1)... % scale from 0 to fov size in vx
    ; % shift to correct bb

expvx(:,2)=(contour(2,:)-smallestentry([thisbb.mm{secdim}(1),thisbb.mm{secdim}(end)]))...
    /abs(diff([thisbb.mm{secdim}(1),thisbb.mm{secdim}(end)]))... % scale from 0 to 1
    *thisbb.imgdim(1)... % scale from 0 to fov size in vx
    ; % shift to correct bb
hold on
expvx(:,2)=thisbb.imgdim(2)-expvx(:,2);

%[idx]=knnsearch(expvx,drawline);
idx = ea_near_line(expvx,drawline);

for pt=round(linspace(1,size(drawline,1),size(drawline,1)/7))
arrhandles{pt}=ea_plot_arrow(drawline(pt,1),drawline(pt,2),expvx(idx(pt),1),expvx(idx(pt),2),'linewidth',2,'headwidth',0.25,'headheight',0.33,'color',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5],'edgecolor',[0.5,0.5,0.5]);
end


clinefiducial=linefiducial;
clinefiducial(:,onedim)=expvx(idx,1);
clinefiducial(:,secdim)=expvx(idx,2);
setappdata(handles.checkstructures,'clinefiducial',clinefiducial); % also store corrected linefiducial.
setappdata(handles.checkstructures,'arrhandles',arrhandles);


function v=smallestentry(ay)
ay=sort(ay,'ascend');
v=ay(1);