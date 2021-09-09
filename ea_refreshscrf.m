function ispresent=ea_refreshscrf(options,handles,directory)

standardslice=ea_loadrefineslice(directory,options,0);
[refineslice,ispresent]=ea_loadrefineslice(directory,options,1);
set(0,'CurrentFigure',handles.scrf);


handles.scrf.CurrentAxes=handles.standardax;
imshow(standardslice);
handles.scrf.CurrentAxes=handles.scfax;
imshow(refineslice);

% calculate and display transform matrix:
if exist([directory,'scrf',filesep,'scrf_instore.mat'],'file')
    mat=ea_getscrfmat(directory);
    handles.affmatrix.String=sprintf('% 0.2f  % 0.2f  % 0.2f  % 0.2f  \n% 0.2f  % 0.2f  % 0.2f  % 0.2f  \n% 0.2f  % 0.2f  % 0.2f  % 0.2f  \n% 0.2f  % 0.2f  % 0.2f  % 0.2f  ',mat');
    save([directory,'scrf',filesep,'scrf_instore_converted.mat'],'mat');
end

function [slice,ispresent]=ea_loadrefineslice(directory,options,refine)

switch refine
    case 1
        refstr='scrf';
    case 0
        refstr='standard';
end
if isfield(options,'init') && options.init
    if ~exist([directory,'scrf',filesep,refstr,'.png'],'file')
        ea_createrefineslice(directory,options,refine);
    end
else
    ea_createrefineslice(directory,options,refine);
end
try
    slice=imread([directory,'scrf',filesep,refstr,'.png']);
    ispresent=1;
catch
    slice=imread([ea_getearoot,'helpers',filesep,'gui',filesep,'scrf_msg.png']);
    ispresent=0;
end


function ea_createrefineslice(directory,options,refine)

switch refine
    case 1
        scrf='scrf';
    case 0
        scrf='';
end

ea_createbbfiles(directory); % needs to unfortunately be done each time since coregistration may have changed.
ea_createmovim(directory,options);
ea_gencheckregfigs(options, 'brainshift');


function ea_createbbfiles(directory)
[options.root,options.patientname]=fileparts(fileparts(directory));
options.root=[options.root,filesep];
options.earoot=ea_getearoot;
options.prefs=ea_prefs(options.patientname);
options=ea_assignpretra(options);

if ~exist([directory,'scrf',filesep,options.prefs.prenii_unnormalized],'file')
    if ~exist([directory,'scrf'],'dir')
        mkdir([directory,'scrf']);
    end
    to{1}=[directory,'scrf',filesep,'bb.nii'];
    from{1}=[ea_space,'bb.nii'];
    try
        ea_apply_normalization_tofile(options,from,to,[options.root,options.patientname,filesep],1,1);
    catch
        ea_error('Please perform normalization first.');
    end
    ea_crop_nii([directory,'scrf',filesep,'bb.nii']);
    ea_reslice_nii([directory,'scrf',filesep,'bb.nii'],[directory,'scrf',filesep,'bb.nii'],[0.4,0.4,0.4]);
    % do put in primary anat file - needs to be done only once.
    fis={options.prefs.prenii_unnormalized};
    copyfile([directory,fis{1}],[directory,'scrf',filesep,fis{1}])
    ea_conformspaceto([directory,'scrf',filesep,'bb.nii'],[directory,'scrf',filesep,fis{1}]);
    % cleanup:
    delete([directory,'scrf',filesep,'bb.nii']);
end
% apply tonemapping if needed
if strcmp(options.prefs.scrf.tonemap,'tp_')
    if ~exist([directory,options.prefs.scrf.tonemap,options.prefs.ctnii_coregistered],'file') && exist([directory,options.prefs.ctnii_coregistered],'file')
        ea_tonemapct(options, 'native');
    end
end
fis={options.prefs.tranii_unnormalized,options.prefs.cornii_unnormalized,options.prefs.sagnii_unnormalized,[options.prefs.scrf.tonemap,options.prefs.ctnii_coregistered]};
for fi=1:length(fis)
    if exist([directory,fis{fi}],'file')
        copyfile([directory,fis{fi}],[directory,'scrf',filesep,fis{fi}])
        ea_conformspaceto([directory,'scrf',filesep,options.prefs.prenii_unnormalized],[directory,'scrf',filesep,fis{fi}],1);
    end
end


function otherfiles=ea_createmovim(directory,options)
switch options.modality
    case 1
        otherfiles={[directory,'scrf',filesep,options.prefs.tranii_unnormalized],...
            [directory,'scrf',filesep,options.prefs.cornii_unnormalized],...
            [directory,'scrf',filesep,options.prefs.sagnii_unnormalized]};
        cnt=1;
        for ofi=1:length(otherfiles)
            if exist(otherfiles{ofi},'file')
                nii=ea_load_nii(otherfiles{ofi});
                delete(otherfiles{ofi});
                nii.img(abs(nii.img)<0.1)=nan;
                if ~exist('AllX','var')
                    AllX=nii.img;
                else
                    AllX(:,:,:,cnt)=nii.img;
                end
                cnt=cnt+1;
            end
        end
        if ~exist('AllX','var')
        ea_error('Something went wrong. Please make sure that you chose the right modality (MR vs. CT) and there are pre- and postoperative acquisitions in the patient directory.');
        end 
        nii.img=ea_nanmean(AllX,4);
        clear AllX
        nii.fname=[directory,'scrf',filesep,'movim.nii'];
        nii.img(isnan(nii.img))=0;
        %nii.img(~(nii.img==0))=zscore(nii.img(~(nii.img==0)));

        ea_write_nii(nii);
    case 2
        otherfiles={[directory,'scrf',filesep,options.prefs.scrf.tonemap,options.prefs.ctnii_coregistered]};
        copyfile([directory,'scrf',filesep,options.prefs.scrf.tonemap,options.prefs.ctnii_coregistered],[directory,'scrf',filesep,'movim.nii']);
end
