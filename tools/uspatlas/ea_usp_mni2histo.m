function ea_usp_mni2histo(uspath)

[infile,inpath]=uigetfile('*.nii','Select file to be warped into Histology space...',fullfile(uspath,'histology','masks'));

options=ea_getptopts(fullfile(uspath,'derivatives','leaddbs','sub-USPatlas'));

ea_mkdir(fullfile(uspath,'histology','masks',ea_stripext(infile)));

ea_apply_normalization_tofile(options, fullfile(inpath,infile), ...
    fullfile(uspath,'histology','masks',ea_stripext(infile),[ea_stripext(infile),'_histo.nii']), ...
    1, 1, ...
    fullfile(uspath,'histology','histo_volume.nii.gz'));

outfolder=fullfile(uspath,'histology','masks',ea_stripext(infile));
nii=ea_load_nii(fullfile(outfolder,[ea_stripext(infile),'_histo.nii']));
gzip(fullfile(outfolder,[ea_stripext(infile),'_histo.nii']));
ea_delete(fullfile(outfolder,[ea_stripext(infile),'_histo.nii']));

ea_mkdir(fullfile(outfolder,'blockface'));
ea_mkdir(fullfile(outfolder,'brightfield'));
ea_mkdir(fullfile(outfolder,'darkfield'));


ea_dispercent(0,'Exporting slices');
fullslice=zeros(4912,7360);
xfrom=850:6000;
for slice=1:size(nii.img,3)
    thisslice=nii.img(:,:,slice);
    fullslice(:,xfrom)=thisslice;
    fullslice=fullslice(end:-1:1,:);
    clear thisslice
    imwrite(fullslice,fullfile(outfolder,[sprintf('%03.f',slice),'.png']));
    if any(fullslice(:)) % also export overlays

        h=figure('Visible','off');
        [c,cont] = contour(mean(fullslice,3),2);
        contP = get(cont,'Parent');
        X = contP.XLim;
        Y = contP.YLim;
        hold on
        c=getframe;
        close(h);
        c=single(c.cdata)./255;
               
        bf=imread(fullfile(uspath,'histology','brightfield',[sprintf('%03.f',slice),'.png']));
        cs=imresize(c,[size(bf,1),size(bf,2)]);
        bf=uint8(single(bf).*cs);
        imwrite(bf,fullfile(outfolder,'brightfield',[sprintf('%03.f',slice),'.png']));

        df=imread(fullfile(uspath,'histology','darkfield',[sprintf('%03.f',slice),'.png']));
        cs=imresize(c,[size(df,1),size(df,2)]);
        df=uint8(single(df).*cs);
        imwrite(df,fullfile(outfolder,'darkfield',[sprintf('%03.f',slice),'.png']));

        blf=imread(fullfile(uspath,'histology','blockface',[sprintf('%03.f',slice),'.png']));
        cs=imresize(c,[size(blf,1),size(blf,2)]);
        blf=uint8(single(blf).*cs);
        imwrite(blf,fullfile(outfolder,'blockface',[sprintf('%03.f',slice),'.png']));

    end
    ea_dispercent(slice/size(nii.img,3));
end
ea_dispercent(1,'end');

