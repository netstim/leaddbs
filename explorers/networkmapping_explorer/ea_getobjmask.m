function usemask=ea_getobjmask(obj,im)
maskn=ea_mask2maskn(obj);
switch maskn
    case 'equation'
        eq=obj.cvmask(5:end); % first letters will be =Eq:
        % isolate <>~
        im;
        eq=strrep(eq,'X', 'im(:)');
        usemask=eval(eq);
    case 'nifti'
        masknii=ea_load_nii(obj.cvmask(5:end));
        usemask=logical(masknii.img(:));
    otherwise
        usemask=ea_getmask(maskn);
end
end

function maskn=ea_mask2maskn(obj)
switch obj.cvmask
    case 'Gray Matter'
        switch obj.outputspace
            case '222'
                maskn='gray';
            case '111'
                maskn='gray_hd';
            case '555'
                maskn='gray_5';
        end
    case 'Brain'
        switch obj.outputspace
            case '222'
                maskn='brain';
            case '111'
                maskn='brain_hd';
            case '555'
                maskn='brain_5';
        end
    case 'Cortex & Cerebellum'
        switch obj.outputspace
            case '222'
                maskn='cortexcb';
            case '111'
                maskn='cortexcb_hd';
            case '555'
                ea_error('Cortex & Cerebellum Mask not supported for 0.5 mm resolution space');
        end
    case 'Cortex'
        switch obj.outputspace
            case '222'
                maskn='cortex';
            case '111'
                maskn='cortex_hd';
            case '555'
                maskn='cortex_5';
        end
    case 'Cerebellum'
        switch obj.outputspace
            case '222'
                maskn='cb';
            case '111'
                maskn='cb_hd';
            case '555'
                ea_error('Cerebellum Mask not supported for 0.5 mm resolution space');
        end
    otherwise
        if strcmp(obj.cvmask(1:4),'=Eq:') % equation
            maskn='equation';
        elseif strcmp(obj.cvmask(1:4),'=Im:') % path to nifti file
            maskn='nifti';
        end
end
end