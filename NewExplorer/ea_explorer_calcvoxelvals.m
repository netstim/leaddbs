function [AllX,space]=ea_explorer_calcvoxelvals(vatlist,obj)
templateresolution = obj.resolution;
if size(vatlist,2)>1
    sidesuffices={'_r','_l'};
else
    sidesuffices={''};
end
disp('Need to export Efields in proper format, this may take a while');
outdir=[fileparts(obj.leadgroup),filesep,'sweetspots',filesep,obj.ID,filesep];
mkdir(outdir)

disp(['Creating Voxel Template: EFthreshold = ' num2str(obj.calcthreshold) ' V/m; Resolution = ' num2str(templateresolution) ' mm.']);
mins = repmat({[]},1,size(vatlist,2));
maxs = repmat({[]},1,size(vatlist,2));
for side=1:size(vatlist,2)
    for vat=1:size(vatlist,1)
        if ~isempty(vatlist{vat,side})
            xyz=[];
            nii=ea_load_nii(vatlist{vat,side});
            nii.img=nii.img>obj.calcthreshold;
            [xyz(1,:),xyz(2,:),xyz(3,:)]=ind2sub(size(nii.img),find(nii.img));
            xyz(4,:)=1;
            xyz=nii.mat*xyz;
            mins{side} = vertcat(mins{side},min(xyz(1:3,:),[],2)');
            maxs{side} = vertcat(maxs{side},max(xyz(1:3,:),[],2)');
        end
    end
    mins{side}=min(mins{side});
    maxs{side}=max(maxs{side});
end

for side=1:size(vatlist,2)
    templatesize=[];
    templatecenter=[];
    for dim=1:3
        templatesize(dim) = numel(mins{side}(dim):templateresolution:maxs{side}(dim))+1;
        templatecenter(dim) = (mins{side}(dim)+maxs{side}(dim))./2;
    end
    templatecenter = round(templatecenter/templateresolution)*templateresolution;
    ea_createTemplateSpace(templatecenter,templatesize,repmat(templateresolution,1,3),outdir,['template',sidesuffices{side},'.nii']);
    space{side}=ea_load_nii([outdir,'template',sidesuffices{side},'.nii']);
end

% now conform each VTA to space
AllX=cell(1,size(vatlist,2));
for side=1:size(vatlist,2)
    AllX{side}=single(nan(prod(space{side}.dim),size(vatlist,1)));
    for vat=1:size(vatlist,1)
        copyfile(vatlist{vat,side},[outdir,'tmp_efield.nii']);
        ea_conformspaceto([outdir,'template',sidesuffices{side},'.nii'],...
            [outdir,'tmp_efield.nii'],0);
        nii=ea_load_nii([outdir,'tmp_efield.nii']);
        % [~,vatnametmp,~]=fileparts(vatlist{vat,side});
        % copyfile([outdir,'tmp_efield.nii'],[outdir,vatnametmp,'.nii'])
        clear vatnametmp
        AllX{side}(:,vat)=single(nii.img(:));
    end
end
ea_delete([outdir,'template.nii']);
ea_delete([outdir,'template_l.nii']);
ea_delete([outdir,'teamplate_r.nii']);
ea_delete([outdir,'tmp_efield.nii']);
ea_delete(outdir);
