function ea_cvshowvatfmri(resultfig,pX,directory,filesare,handles,pV,selectedparc,mod,options)
%mV=pV; % duplicate labeling handle
%mX=pX; % duplicate labeling data
stims=get(handles.vatseed,'String');
stim=stims{get(handles.vatseed,'Value')};

% check out which vats to use
[usevat,dimensionality,currentseed,sides]=ea_checkvatselection(handles);
if ~dimensionality
    return
end

pX=round(pX);

mod = strrep(mod, 'Patient''s fMRI - ', '');
options.prefs.rest = [mod, '.nii'];

if ~exist([directory,'stimulations',filesep,ea_nt(options),stim,filesep,'vat_',mod,'.mat'],'file')
    ea_warp_vat('rest', options, handles);
    vat_tc=ea_extract_timecourses_vat(options,handles,usevat,dimensionality);
    save([directory,'stimulations',filesep,ea_nt(options),stim,filesep,'vat_',mod,'.mat'],'vat_tc');
else
    load([directory,'stimulations',filesep,ea_nt(options),stim,filesep,'vat_',mod,'.mat']);
end

parcs=get(handles.labelpopup,'String');
tc=load([directory,'connectomics',filesep,parcs{get(handles.labelpopup,'Value')},filesep,ea_stripext(options.prefs.rest),'_tc']);
fn=fieldnames(tc);
tc=eval(['tc.',fn{1},';']);
tc=[vat_tc,tc];

timedim=size(tc,1);
tiwindow=get(handles.timewindow,'String');
tiframe=get(handles.timeframe,'String');

if strcmp(tiwindow,'all') || strcmp(tiframe,'all')
    % use whole CM
    cm=corrcoef(tc);
else
    tiframe=str2double(tiframe);         tiwindow=str2double(tiwindow);
    % check if selected time window is possible:
    if (tiframe+tiwindow)>timedim || tiframe<1 % end is reached
        set(handles.timeframe,'String','1'); tiframe=1; % reset timeframe to 1
        if tiwindow>size(tc,1)
            set(handles.timewindow,'String','1'); tiwindow=1;
        end
    end
    cm=corrcoef(tc(tiframe:tiframe+tiwindow,:)); % actual correlation

    if get(handles.timecircle,'Value')
        % make a step to next timeframe (prepare next iteration).
        if (tiframe+tiwindow+1)>timedim
            set(handles.timeframe,'String','1')
        else
            set(handles.timeframe,'String',num2str(tiframe+1))
        end
    end
end

for iside=1:length(options.sides)
    side=options.sides(iside);
    seedcon=cm(side,:);
    seedcon=seedcon(3:end);
    thresh=get(handles.vatthresh,'String');

    if strcmp(thresh,'auto')
        thresh=nanmean(seedcon)+1*0.5*nanstd(seedcon);
    else
        thresh=str2double(thresh);
    end

    tseedcon=seedcon;
    tseedcon(tseedcon<thresh)=0;
    tseedcon(currentseed)=0;
    pX(pX==0)=nan;
    mX=pX;
    for cs=1:length(tseedcon) % assign each voxel of the corresponding cluster with the entries in tseedcon. Fixme, this should be doable wo forloop..
        mX(ismember(round(pX),cs))=tseedcon(cs);
    end

    Vvat=spm_vol([directory,'stimulations',filesep,ea_nt(options),stim,filesep,'vat_',usevat{side},'.nii,1']);
    Xvat=spm_read_vols(Vvat);
    vatseedsurf{side}=ea_showseedpatch(resultfig,Vvat,Xvat,options);

    %sX=ismember(round(pX),currentseed);
    set(0,'CurrentFigure',resultfig)
    set(handles.vatthreshis,'String',num2str(thresh));
    vatsurf{side}=ea_showconnectivitypatch(resultfig,pV,mX,thresh,[],[],[],1,0);
end

setappdata(resultfig,'vatsurf',vatsurf);
setappdata(resultfig,'vatseedsurf',vatseedsurf);
