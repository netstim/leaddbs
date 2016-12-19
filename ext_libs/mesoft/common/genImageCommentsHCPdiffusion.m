function genFileTableHCP(folder)


b = load('bvals');
d = load('bvecs');

fi = fopen(fullfile(folder,'imagecomments.txt'),'w');
for k = 1:length(b),
    fprintf(fi,'%s,%i b=%f %f %f %f \n',fullfile(folder,'data.nii'),k,b(k),d(1,k),d(2,k),d(3,k));
end;
    
fclose(fi);

