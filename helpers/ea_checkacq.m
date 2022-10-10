function acqTag = ea_checkacq(pixdim)
% check image acq tag (iso, ax, cor or sag)
% Input can be the pixdim or the path of the image

if ischar(pixdim)
    header = ea_fslhd(pixdim);
    pixdim = [header.pixdim1, header.pixdim2, header.pixdim3];
end

if range(pixdim) < 0.05
    acqTag = 'iso';
    return;
end

[C,~, ic] = unique(pixdim);
if numel(C) == 1
    acqTag = 'iso';
else
    if numel(C) == 2
        count = accumarray(ic, 1);
        flag = find(pixdim == C(count==1));
    else
        multi = [pixdim(2)*pixdim(3), pixdim(1)*pixdim(3), pixdim(1)*pixdim(2)];
        flag = find(multi == min(multi));
    end

    switch flag
        case 1
            acqTag = 'sag';
        case 2
            acqTag = 'cor';
        case 3
            acqTag = 'ax';
    end
end
