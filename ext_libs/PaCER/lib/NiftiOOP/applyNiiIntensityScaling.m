function nii = applyNiiIntensityScaling(nii)
% Code extracte from Jimmy Shen's (jimmy@rotman-baycrest.on.ca) xform_nii.m
% method, which is part of the nifti toolbox
%
% Andreas Husch, 2019

    %  if scl_slope field is nonzero, then each voxel value in the
       %  dataset should be scaled as: y = scl_slope * x + scl_inter
       %  I bring it here because hdr will be modified by change_hdr.
       %
       if nii.hdr.dime.scl_slope ~= 0 & ...
        ismember(nii.hdr.dime.datatype, [2,4,8,16,64,256,512,768]) & ...
        (nii.hdr.dime.scl_slope ~= 1 | nii.hdr.dime.scl_inter ~= 0)

          nii.img = ...
        nii.hdr.dime.scl_slope * double(nii.img) + nii.hdr.dime.scl_inter;

          if nii.hdr.dime.datatype == 64

             nii.hdr.dime.datatype = 64;
             nii.hdr.dime.bitpix = 64;
          else
             nii.img = single(nii.img);

             nii.hdr.dime.datatype = 16;
             nii.hdr.dime.bitpix = 32;
          end

          nii.hdr.dime.glmax = max(double(nii.img(:)));
          nii.hdr.dime.glmin = min(double(nii.img(:)));

          %  set scale to non-use, because it is applied in xform_nii
          %
          nii.hdr.dime.scl_slope = 0;

       end

       %  However, the scaling is to be ignored if datatype is DT_RGB24.

       %  If datatype is a complex type, then the scaling is to be applied
       %  to both the real and imaginary parts.
       %
       if nii.hdr.dime.scl_slope ~= 0 & ...
        ismember(nii.hdr.dime.datatype, [32,1792])

          nii.img = ...
        nii.hdr.dime.scl_slope * double(nii.img) + nii.hdr.dime.scl_inter;

          if nii.hdr.dime.datatype == 32
             nii.img = single(nii.img);
          end

          nii.hdr.dime.glmax = max(double(nii.img(:)));
          nii.hdr.dime.glmin = min(double(nii.img(:)));

          %  set scale to non-use, because it is applied in xform_nii
          %
          nii.hdr.dime.scl_slope = 0;

       end
end