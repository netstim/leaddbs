function Dlist=ea_ftr2aniso(fibers,smri,ten)



[idx,~]=cellfun(@size,fibers');

fibers=cell2mat(fibers);
idxv=zeros(size(fibers,1),1);
lid=1; cnt=1;
for id=idx
    idxv(lid:lid+id-1)=cnt;
    lid=lid+id;
    cnt=cnt+1;
end


% convert wm indices to mm format
[xx,yy,zz]=ind2sub(size(smri.white),find(smri.white));
wmpnts=[xx,yy,zz,ones(size(xx,1),1)]';
wmpnts=smri.hdr.mat*wmpnts;
wmpnts=wmpnts(1:3,:)';

resolution=abs(smri.hdr.mat(logical(eye(4))));
resolution=mean(resolution(1:3));

% figure
% 
% plot3(wmpnts(:,1),wmpnts(:,2),wmpnts(:,3),'r.');
% hold on
% plot3(fibers(:,1),fibers(:,2),fibers(:,3),'b.');

[idx]=rangesearch(fibers,wmpnts,3);
dimensionality=length(idx);
Dlist=zeros(dimensionality,6);
ea_dispercent(0,'Generating tensor estimates from fibers');
for pix=1:dimensionality
    ea_dispercent(pix/dimensionality);
   thispixani=idx{pix}; 
   
   traverses1=abs(fibers(thispixani-1,:)-fibers(thispixani,:));
   traverses2=abs(fibers(thispixani+1,:)-fibers(thispixani,:));
   traverses1(sum(traverses1,2)>10,:)=nan;
   traverses2(sum(traverses2,2)>10,:)=nan;
   
   traverses=cat(3,traverses1,traverses2);
   traverses=ea_nanmean(traverses,3);
   if ~isempty(traverses)
       if size(traverses,1)>2
           [~,S,V]=svd(traverses);
       else
           [~,S,V]=svd(repmat(traverses,3,1));
       end
       D=V'*S(1:3,1:3)*V;
       Dlist(pix,:)=D([1,5,9,2,6,3]);
       
   else
       Dlist(pix,:)=[1,1,1,0,0,0];
   end
end
ea_dispercent(1,'end');


function [thispixani,thisD]=removeneighbors(thispixani)

[thispixani]=sort(thispixani);
thispixani=thispixani(diff(thispixani)>1);
