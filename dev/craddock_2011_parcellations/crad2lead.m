fniis=dir('*.nii');

for fnii=1:length(fniis)
   nii=ea_load_nii(fniis(fnii).name);
   nii.img=round(nii.img);
         of=fopen([fn,'.txt'],'w');
   for c=1:max(nii.img(:));
      [xx,yy,zz]=ind2sub(size(nii.img),find(nii.img==c));
      XYZ=[xx,yy,zz,ones(length(xx),1)]';
      XYZmm=[nii.mat*XYZ]';
      XYZmm=XYZmm(:,1:3);
      centr=mean(XYZmm,1);
      [~,fn]=fileparts(fniis(fnii).name);
      try
fprintf(of,'%d %s \n',c,[num2str(centr(1)),'_',num2str(centr(2)),'_',num2str(centr(3))]);
      catch
          keyboard
      end
   end
   fclose(of);

end