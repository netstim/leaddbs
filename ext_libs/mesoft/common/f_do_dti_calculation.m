%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%% function  do_dti_calculation
%
%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D_trace, FA, RA, EigenVal, eigVect, M_error, sqErr] = f_do_dti_calculation(DE_direction, b_factor, singleSliceDWI_set, sizeIm, tresholdAbs )

D_trace= []; FA= []; RA= [];  EigenVal= []; eigVect= []; M_error= []; sqErr= [];
DTE_mode  = DE_direction ; % transfer DE grad directions read from dicom comment field
%% now time to created the b-matrix
DE = DTE_mode;%squeeze(DTE_mode(1,:,:)); % !! NB now use the data for the 1st slice, if motion correction is on has to be fixed ika 050413
B_tensor =zeros(3,3,length(DE));
DE_size = size(DE);
for i = 2:DE_size(1)
    tmp=DE(i,:)'*DE(i,:);
    if sum(sum(abs(tmp) )) >0, % ika 041118
        B_tensor(:,:,i)= (b_factor/trace(tmp))*tmp;end
end
nsteps =length(B_tensor);
X1 = [  reshape(B_tensor(1,1,:), 1, nsteps , 1)]; %1st diag element
X2 = [  reshape(B_tensor(2,2,:), 1, nsteps , 1)]; %2nd diag element
X3 = [  reshape(B_tensor(3,3,:), 1, nsteps , 1)]; %3rd diag element
X4 = [  reshape(B_tensor(1,2,:), 1, nsteps, 1)];
X5 = [  reshape(B_tensor(1,3,:), 1, nsteps, 1)];
X6 = [  reshape(B_tensor(2,3,:), 1, nsteps, 1)];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   do dti calculation for the paticular slice
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mask =zeros(sizeIm(1), sizeIm(2));
treshold = mean(mean(singleSliceDWI_set(:,:,1 )));  % NB now maks will depend from the slice - for low slices that can be bad
D_trace= mask;
D_iso= mask;

for i=1:sizeIm(1) % vertical
    for j=1:sizeIm(2) % horizontal
        if singleSliceDWI_set(i,j,1) > tresholdAbs %treshold_factor *treshold
            mask(i,j)=1;
        end
    end
end

% reserviere speicher
D_trace= zeros(sizeIm(1:2));
eigVect=zeros(sizeIm(1), sizeIm(2), 3,3);
difComp=zeros(sizeIm(1), sizeIm(2), 3,3);
M_error= zeros(sizeIm(1:2));
sqErr= zeros(sizeIm(1:2));


Y1 =singleSliceDWI_set;
for x_coor =1: sizeIm(2)
    for y_coor =1: sizeIm(1)
        Y  = [];	% clear
        if (mask(y_coor, x_coor)) ~= 0;
            tempY =(Y1(  y_coor, x_coor, :));
            itemp =find(tempY <= 0);
            if(length(itemp) >=1),
                %     'warning! signal is zero or negative',tempY(itemp),  tempY(itemp)= 20.0222; x_coor, y_coor,%  itemp,
                tempY(itemp)= 20.0222;
            end
            y=squeeze( sum( sum( log(tempY ), 2 ),1 ) );  % 2:5 - from 2 to 5th image
            Y= [Y y'];
            anisoX=[ones(size(X1')), X1', X2', X3',2*X4', 2*X5', 2*X6'];
            aniso_a=anisoX\(Y');
            YY = anisoX*aniso_a;
            MaxErr = max(abs(YY- Y'));

            d_tens= [aniso_a(2), aniso_a(5) ,aniso_a(6) ;aniso_a(5) ,aniso_a(3) ,aniso_a(7); aniso_a(6) ,aniso_a(7) ,aniso_a(4) ] ; % check off_diaf elem order!!!!
            [vv, dd] = eig(-d_tens) ;
            %'vv -eigen-vector, dd - eigen-value';
            D_trace(y_coor, x_coor ) = sum(sum(dd)) * 1/3; % ika 2k0121 added 1/3 factor
            eigVect(y_coor, x_coor , :,:) =vv(:,:);
            difComp(y_coor, x_coor , :,:) =dd(:,:);
            % D_raw(y_coor, x_coor,:,: ) = -d_tens;  % ika9041130
            M_error(y_coor, x_coor ) = MaxErr;
            sqErr(y_coor, x_coor ) = mean((YY- Y').^2);
        end %if
    end % for
end % for

EigenVal= sum(difComp, 3);
dd_1= EigenVal(:,:,1).^2 + EigenVal(:,:,2).^2 + EigenVal(:,:,3).^2;
l_mean= (EigenVal(:,:,1) + EigenVal(:,:,2) + EigenVal(:,:,3))/3;
dd_3= (EigenVal(:,:,1)-l_mean).^2 +(EigenVal(:,:,2)-l_mean).^2 +(EigenVal(:,:,3)-l_mean).^2 ;
i6= find(dd_3 ~= 0);
i7= find(l_mean ~= 0);
FA = zeros(size(dd_3));
RA = FA;
FA(i6)=sqrt(3/2)* sqrt(dd_3(i6)./dd_1(i6));
RA(i7)=sqrt(1/3)* sqrt(dd_3(i7))./l_mean(i7);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%       DONE DTI calculation
%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% end of file %%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
