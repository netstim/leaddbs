function usemask = ea_unified_getobjmask(obj,im)
    % Determine mask type string
    maskn = ea_mask2maskn(obj);

    switch maskn
        case 'equation'
            eq = obj.cvmask(5:end); % remove =Eq:
            eq = strrep(eq, 'X', 'im(:)');
            usemask = eval(eq);

        case 'nifti'
            masknii = ea_load_nii(obj.cvmask(5:end));
            usemask = logical(masknii.img(:));

        otherwise
            usemask = ea_getmask(maskn);
    end
end

function maskn = ea_mask2maskn(obj)
    % Decide which resolution to use
    res = obj.calcsettings.functionalresolution;  % or structuralresolution

    % Map resolution to suffix
    switch res
        case '2 mm'
            suffix = '';
        case '1 mm'
            suffix = '_hd';
        case '0.5 mm'
            suffix = '_5';
        otherwise
            error('Unsupported resolution: %s', res);
    end

    % Map cvmask to base name
    switch obj.cvmask
        case 'Gray Matter'
            maskn = ['gray' suffix];
        case 'Brain'
            maskn = ['brain' suffix];
        case 'Cortex & Cerebellum'
            if strcmp(res,'0.5 mm')
                ea_error('Cortex & Cerebellum not supported for 0.5 mm');
            else
                maskn = ['cortexcb' suffix];
            end
        case 'Cortex'
            maskn = ['cortex' suffix];
        case 'Cerebellum'
            if strcmp(res,'0.5 mm')
                ea_error('Cerebellum not supported for 0.5 mm');
            else
                maskn = ['cb' suffix];
            end
        otherwise
            if startsWith(obj.cvmask,'=Eq:')
                maskn = 'equation';
            elseif startsWith(obj.cvmask,'=Im:')
                maskn = 'nifti';
            else
                error('Unknown cvmask: %s', obj.cvmask);
            end
    end
end