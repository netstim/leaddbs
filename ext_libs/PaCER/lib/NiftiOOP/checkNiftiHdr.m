%% checkNifti - basic consistency checking of nifti headers according to the nifti-1 standard
%
% Andreas Husch
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
%
% 2018
%
% mail@andreashusch.de, andreas.husch@uni.lu
function checkNiftiHdr(nii)
NIFTI_DOC = ' <a href="https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html">https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html</a>';

if(nii.hdr.hist.qform_code == 0)
    warning(['checkNifti: Nifti qform_code = 0, this is discouraged by the Nifti-1 standard. '...
        'Please check carefully and refer to ' NIFTI_DOC ' section 5. METHOD 1. Consider '...
        'fixing the nifti header of your file before continuing.']);
end

if(nii.hdr.hist.sform_code == nii.hdr.hist.qform_code)
    h = nii.hdr.hist;
    h.pixdim = nii.hdr.dime.pixdim;
    
    if(~all(all(nii.qform - nii.sform < eps(single(1)))))
        warning(['checkNifti: qform_code == sform_code, however the transformation defined in the qform '...
            'differes from the sform! This might indicate a serious flaw in the nifti header '...
            'and lead to unexpected results as different tools/algorithms might deal differently '...
            'with this situation. Fix the nifti header of your file before continuing.']);
    end
end
end

