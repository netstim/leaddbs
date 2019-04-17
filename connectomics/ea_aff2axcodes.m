function axcodes = ea_aff2axcodes(affine, tol, labels)
% Calculcate axis direction codes from affine matrix
% Can be used to set the 'voxel_order' in track file (*.trk) header
%
% tol: {None, float}, optional
% 	  Threshold below which SVD values of the affine are considered zero.
%     If 'tol' is None, and 'S' is an array with singular values for 'affine',
%     and 'eps' is the floating-point relative accuracy, then 'tol' set to
%     'max(S) * max(size(RS)) * eps'
%
% labels : optional, None or Nx2 char array
%     In each row of the char array are the labels for beginning and end of
%     output axis.  That is, if the first row in 'ornt' is [1, 1], and the
%     first row of labels is ('back', 'front') then the first returned axis
%     code will be 'front'.  If the first row in 'ornt' had been [1, -1] then
%     the first returned value would have been 'back'.  If None, equivalent
%     to ['L', 'R'; 'P', 'A'; 'I', 'S'] - that is - RAS axes.
%

%% Orientation of input axes in terms of output axes for 'affine'
% extract the underlying rotation, zoom, shear matrix
RZS = affine(1:size(affine,1)-1,1:size(affine,2)-1);
zooms = sqrt(sum(RZS.*RZS));
% Zooms can be zero, in which case all elements in the column are zero, and
% we can leave them as they are
zooms(zooms==0) = 1;
RS = RZS./repmat(zooms,3,1);
% Transform below is polar decomposition, returning the closest shearless
% matrix R to RS
[P, S, Qs] = svd(RS);
% Threshold the singular values to determine the rank.
if nargin < 2
    tol = max(S) * max(size(RS)) * eps;
end
R = P(:,diag(S)'>tol) * Qs(diag(S)'>tol,:)';
% The matrix R is such that R*R' is projection onto the columns of P(:,keep)
% and R'*R is projection onto the rows of Qs(keep,:)'. R (== R*eye(p)) gives
% rotation of the unit input vectors to output coordinates. Therefore, the
% row index of abs max R(:,N), is the output axis changing most as input
% axis N changes.  In case there are ties, we choose the axes iteratively,
% removing used axes from consideration as we go.
ornt = nan(size(affine,2)-1,2);
for in_ax=1:size(affine,2)-1
    col = R(:,in_ax);
    if ~ea_allclose(col ,0)
        out_ax = find(abs(col) == max(abs(col)));
        ornt(in_ax,1) = out_ax;
        assert(col(out_ax) ~= 0);
        if col(out_ax) < 0
            ornt(in_ax,2) = -1;
        else
            ornt(in_ax,2) = 1;
        end
        % remove the identified axis from further consideration, by zeroing
        % out the corresponding row in R
        R(out_ax, :) = 0;
    end
end

%% Convert orientation 'ornt' to labels for axis directions
if nargin < 3
    labels = ['L', 'R'; 'P', 'A'; 'I', 'S'];
end
axcodes = [];
for i=1:size(ornt,1)
   axno = ornt(i,1);
   direction = ornt(i,2);
   if isnan(axno)
       axcodes = [axcodes nan];
       continue
   end
   if direction == -1
       axcodes = [axcodes,labels(axno,1)];
   elseif direction == 1
       axcodes = [axcodes,labels(axno,2)];
   end
end
