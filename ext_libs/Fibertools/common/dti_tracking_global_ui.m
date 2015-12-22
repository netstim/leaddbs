% interface program for batch program
%
% Author: Marco Reisert
% PC 2010

function out = dti_tracking_global_ui(P)

h = findobj(0,'tag','fiberGT_main');
if ~isempty(h)
    delete(h);
end;
 
if isfield(P,'temp') %% GT accum
       
    fiberGT_tool('loadFTR',P.filenameFTR{1});
    ftr = fiberGT_tool('createEFTR',P.numits,P.numsamps,P.temp);
    [path,nam,ext] = fileparts(P.filenameFTR{1});
    [path_out,nam_out,ext_out] = fileparts(P.fname);
    if isempty(path_out),
        path_out = path;
    end;
    ftrname = fullfile(path_out,[nam_out ext_out]);
    ftrstruct_write(ftr,ftrname);
    
     out.ftr = {ftrname};
    
else


    if isfield(P.fname,'filenameHARDI'),
        fiberGT_tool('loadHARDI',P.fname.filenameHARDI{1});
        [path,nam,ext] = fileparts(P.fname.filenameHARDI{1});
        if length(nam) > 6 && strcmp(nam(end-5:end), '_HARDI')
            nam = nam(1:end-6);
        end
    elseif isfield(P.fname,'filenameNII'),
        fiberGT_tool('loadHARDI',P.fname.filenameNII);    
        [path,nam,ext] = fileparts(P.fname.filenameNII{1});
    else
        fiberGT_tool('loadDTD',P.fname.filenameDTD{1},1000);
        [path,nam,ext] = fileparts(P.fname.filenameDTD{1});
        if length(nam) > 4 && strcmp(nam(end-3:end), '_DTD')
            nam = nam(1:end-4);
        end
    end

    h = findobj(0,'tag','fiberGT_main');


    if isfield(P.newfile,'auto')
        out.ftr = {fullfile(path, [nam,'_FTR.mat'])};
        out.fd = {fullfile(path, [nam,'_FTR_fd.mat'])};
        out.epd = {fullfile(path, [nam,'_FTR_epd.mat'])};
    else
        [p n e] = fileparts(P.newfile.out.fname);
        out.ftr = {fullfile(P.newfile.out.dir{1}, P.newfile.out.fname)};
        out.fd = {fullfile(P.newfile.out.dir{1}, fullfile(p,[n '_fd' e]))};    
        out.epd = {fullfile(P.newfile.out.dir{1}, fullfile(p,[n '_epd' e]))};    
    end
    fiberGT_tool('setFTRname',out.ftr{1});



    if isfield(P.trackingarea,'maskstruct')
        name = P.trackingarea.maskstruct.roiid;
        fiberGT_tool('loadMaskStruct',P.trackingarea.maskstruct.filenameMASK{1},P.trackingarea.maskstruct.roiid);    
    elseif isfield(P.trackingarea,'masknii')
        name =  P.trackingarea.masknii.filenameMASKnii;
        threshold = P.trackingarea.masknii.thresholdMask;
        fiberGT_tool('loadMask',name{1},threshold);
    else
        fiberGT_tool('estimateMask');
    end;



    if P.parameters == 0,
        fiberGT_tool('suggestSparse');
    else
        fiberGT_tool('suggestDense');
    end;

    if isfield(P.para_weight,'custom_para_weight'),
        fiberGT_tool('setparam','weight',P.para_weight.custom_para_weight);
    end;

    if isfield(P.para_other,'custom_para_other'),
        fiberGT_tool('setparam','startTemp',P.para_other.custom_para_other(1));
        fiberGT_tool('setparam','stopTemp',P.para_other.custom_para_other(2));
        fiberGT_tool('setparam','steps',P.para_other.custom_para_other(3));
        fiberGT_tool('setparam','itnum',P.para_other.custom_para_other(4));
        fiberGT_tool('setparam','width',P.para_other.custom_para_other(5));
        fiberGT_tool('setparam','length',P.para_other.custom_para_other(6));
        fiberGT_tool('setparam','chemPot2',P.para_other.custom_para_other(7));
        fiberGT_tool('setparam','bfac',P.para_other.custom_para_other(8));
        fiberGT_tool('setparam','connlike',P.para_other.custom_para_other(9));
    end;




    set(findobj(h,'tag','fiberGT_editfiblength'),'string',[num2str(P.minlen) ';' num2str(P.maxlen)]);


    fiberGT_tool('start');
    fiberGT_tool('saveFD');

end


h = findobj(0,'tag','fiberGT_main');
if ~isempty(h)
    delete(h);
end;



