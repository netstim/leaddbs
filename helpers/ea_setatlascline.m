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


drawline_cum_length = cumsum(sqrt(sum(diff(drawline).^2,2))); % calculate length that each point adds to the total line

minarrow = min(5,size(drawline,1)); % minimum number of arrows

drawline_sample = linspace(drawline_cum_length(1), drawline_cum_length(end), max(minarrow,round(drawline_cum_length(end)/12))); % get equally distant positions
pt = knnsearch(drawline_cum_length, drawline_sample'); % get drawline index
expvx_sample = idx(round(linspace(1,length(idx),length(drawline_sample)))); % sample atlas line with same ammount of pts

for i = 1:length(drawline_sample) 
    arrhandles{i}=ea_plot_arrow(drawline(pt(i),1),drawline(pt(i),2),expvx(expvx_sample(i),1),expvx(expvx_sample(i),2),'linewidth',2,'headwidth',0.25,'headheight',0.33,'color',[0.2 0.8 0.8],'facecolor',[0.2 0.8 0.8],'edgecolor',[0.2 0.8 0.8]);
end


clinefiducial=linefiducial;
clinefiducial(:,onedim)=expvx(idx,1);
clinefiducial(:,secdim)=expvx(idx,2);
setappdata(handles.checkstructures,'clinefiducial',clinefiducial); % also store corrected linefiducial.
setappdata(handles.checkstructures,'arrhandles',arrhandles);


function v=smallestentry(ay)
ay=sort(ay,'ascend');
v=ay(1);