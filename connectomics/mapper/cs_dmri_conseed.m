function cs_dmri_conseed(dfold,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask,space)
%# ea_load_nii

switch cmd
    case 'seed'
        map=ea_load_nii(space);
        cfile=dfold;
        load(cfile,'fibers');
        mapsz=size(map.img);

        seedfiles=sfile;
        tree=KDTreeSearcher(fibers(:,1:3));
        for s=1:length(seedfiles)
            map.img(:)=0;
            Vseed=ea_load_nii(seedfiles{s});
            
            if isempty(varargin{4});
                maxdist=abs(Vseed.mat(1))/2;
            else
                maxdist=varargin{4};
            end
            Vseed.img(isnan(Vseed.img))=0;
            
            ixs=find(Vseed.img);
            % subtract nan values from these
            
            ixvals=Vseed.img(ixs);
            [xx,yy,zz]=ind2sub(size(Vseed.img),ixs);
            XYZvx=[xx,yy,zz,ones(length(xx),1)]';
            clear ixs
            XYZmm=Vseed.mat*XYZvx;
            XYZmm=XYZmm(1:3,:)';
            clear Vseed
            
            ids=rangesearch(tree,XYZmm,maxdist,'distance','chebychev');
            % select fibers for each ix
            for ix=1:length(ixvals)
                % assign fibers on map with this weighted value.
                fibnos=unique(fibers(ids{ix},4));
                
                allfibcs=fibers(ismember(fibers(:,4),fibnos),1:3);
                allfibcs=round(map.mat\[allfibcs,ones(size(allfibcs,1),1)]');
                allfibcs(:,logical(sum(allfibcs<1,1)))=[];
                topaint=sub2ind(mapsz,allfibcs(1,:),allfibcs(2,:),allfibcs(3,:));
                map.img(topaint)=map.img(topaint)+ixvals(ix);
            end
            
            [pth,fn]=fileparts(seedfiles{s});
            
            map.fname=fullfile(outputfolder,[fn,'_struc',cmd,'.nii']);
            map.dt=[16,0];
            spm_write_vol(map,map.img);
        end
        
    otherwise
        warning('Structural connectivity only supported for seeds');
end
