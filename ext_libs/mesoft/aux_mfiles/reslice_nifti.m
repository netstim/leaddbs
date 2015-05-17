function vol_res = reslice(vol,hdr,hdrref,oversamp)



T1 = getMat(hdr);
T2 = getMat(hdrref);


SC = diag([[1 1 1]/oversamp 1]);
T2 = T2*SC;

if isfield(hdrref,'dataAy');
    destdim = size(hdrref.dataAy)*oversamp;
else
    destdim = hdrref.dime.dim(2:4)*oversamp;
end;

[X Y Z] = ndgrid(1:destdim(1),1:destdim(2),1:destdim(3));
NR = [X(:)-1 Y(:)-1 Z(:)-1 X(:)*0+1]*(inv(T1)*T2)' +1;
vol_res = reshape(interp3(double(vol),NR(:,2),NR(:,1),NR(:,3)),destdim);


function mat = getMat(hdr);

if isfield(hdr,'dataAy') % a mrstruct
    Q = diag([-1 -1 1 1]); % different world coordinates
    mat = Q*hdr.edges;
         
else % a nifti
    mat = [hdr.hist.srow_x ; hdr.hist.srow_y ; hdr.hist.srow_z ];
    mat(4,4) = 1;
end;

