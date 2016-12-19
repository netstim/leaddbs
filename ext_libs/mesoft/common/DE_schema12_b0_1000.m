% DE_schema12_b0_1000
% ika 030124
% script to calculate DTI fron Siem DTI images
% last update 14.03.2007
%
% example : DE_schema12_b0_1000('J.Rick-030124-1.5T','E:\zDataTemp\Rick_Jochen_20030124', DTI, 49, 40)
% explanation: Res =DE_schema12_b0_12(FileName,OutputDir, DTI, N_slices, tresholdAbs, DE_schema, b_factor)
% if b_factor is undefinedm, assumed to be 1000
% if b_factor = BMatrix with size 3x3xNumberOfDESteps, this BMatrix will be used for calculations 
%
%   Writen by:
%   ------------------------------------------------------------
%   Dr.Kamil Il'yasov
%   Section of Medical Physics, Dept. of Diagnostic Radiology
%   University Medical Center of Freiburg
%   Hugstetter Str. 55, D-79106 Freiburg, Germany
%   Tel.: +49 761 270 7392
%   Fax : +49 761 270 3831
%   E_mail: Kamil.Ilyasov@uniklinik-freiburg.de
%   ------------------------------------------------------------

%   Corrected by Susanne Schnell, March 2007
%
%   Changes by ika:
%   ika 031117 added option for exact b-matrix calculation  
%   040726 - correction for b-matrix - to keep identical with fiberviewer_tool
%   ika041121 added Liambda1, 2,3 in save DTD struct
%

function Res =DE_schema12_b0_1000(FileName,OutputDir, DTI, N_slices, tresholdAbs, DE, b_factor)

% if (nargin < 6) 
%     DE_schema = 'DE12';
%     b_factor= 1000;
% end

if (nargin < 7) 
    b_factor= 1000;
end;

if (nargin < 5) 
    'file DE_schema12_b0_12: not enouth input arguments'
end;

treshold_factor = .25 ; %1.0
%tresholdAbs = 30; % absolute value

Max_slice = N_slices;
Roi_size = 0; % ADC will be calculated in square ROI's of size (Roi_size+1)^2

%%%NB!!!!!%%%
DE(:,2) = -DE(:,2); % ika 030114  % 31.10.03 - why I have done that????
%%%NB!!!!!%%%

DE_size = size(DE);
b_factorSize = size(b_factor)
if (length(b_factor)>1 )
    if (length(b_factorSize)==3 && b_factorSize(1)==3 &&  b_factorSize(2)==3) 
        B_tensor = b_factor;
        %%%NB!!!!!%%%
        B_tensor(2,1,:) = -B_tensor(2,1,:);
        B_tensor(1,2,:) = -B_tensor(1,2,:);
        B_tensor(2,3,:) = -B_tensor(2,3,:);
        B_tensor(3,2,:) = -B_tensor(3,3,:);
        %%%NB!!!!!%%%
    elseif(length(b_factorSize)==2 && b_factorSize(1)==1 ) 
        if b_factor(1)==0
            B_tensor(:,:,1)= zeros(3,3);
            b_start = 2;
        else
            b_start = 1;
        end
        i_count = b_start;
        for j= b_start:length(b_factor)
            for i = 2:DE_size(1)
                tmp =DE(i,:)'*DE(i,:);
                B_tensor(:,:,i_count)= (b_factor(j)/trace(tmp))*tmp;
                i_count = i_count +1;
            end
        end   
    else
        'b-factor has unsupported format'
        ret = 0;
        return
    end
else
    for i = 2:DE_size(1)
        tmp =DE(i,:)'*DE(i,:);
        B_tensor(:,:,i)= (b_factor/trace(tmp))*tmp;
    end
end;

nsteps =length(B_tensor);

X1 = [  reshape(B_tensor(1,1,:), 1, nsteps , 1)]; %1st diag element
X2 = [  reshape(B_tensor(2,2,:), 1, nsteps , 1)]; %2nd diag element
X3 = [  reshape(B_tensor(3,3,:), 1, nsteps , 1)]; %3rd diag element
X4 = [  reshape(B_tensor(1,2,:), 1, nsteps, 1)];
X5 = [  reshape(B_tensor(1,3,:), 1, nsteps, 1)];
X6 = [  reshape(B_tensor(2,3,:), 1, nsteps, 1)];

% ++++++++++++++++++++++++++++++++++++=====

imageSize = size(DTI.dataAy);

for Slice=  1:Max_slice    
    B= squeeze(DTI.dataAy(:,:,Slice,:));
    mask =zeros(imageSize(1), imageSize(2));
    treshold = mean(mean(B(:,:,1 )));
    D_trace= mask;
    D_iso= mask;
    
    for i=1:imageSize(1) % vertical
        for j=1:imageSize(2) % horizontal
            if B(i,j,1) > tresholdAbs %treshold_factor *treshold
                mask(i,j)=1;
            end
        end
    end
    
    Y1 =squeeze(DTI.dataAy(:,:,Slice,:  ));
    
    for x_coor =1: imageSize(2)
        for y_coor =1: imageSize(1)
            Y  = [];	% clear      
            if (mask(y_coor, x_coor)) == 0;
                D_trace(y_coor, x_coor ) = 0;
                eigVect(y_coor, x_coor , :,:) =zeros(3,3);
                difComp(y_coor, x_coor , :,:) =zeros(3,3);
                D_raw(y_coor, x_coor , :,:) =zeros(3,3);
                M_error(y_coor, x_coor ) = 0;
            else
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
                D_raw(y_coor, x_coor,:,: ) = -d_tens;
                M_error(y_coor, x_coor ) = MaxErr;
            end %if
        end % for
    end % for
    
    D_trace1(:,:,Slice ) = D_trace;
    eigVect1(:,: ,Slice, :,:) =eigVect;
    D_raw1(:,:,Slice, :,: ) = D_raw;
    Max_error(:,:,Slice  ) = M_error;
    
    EigenVal= sum(difComp, 3); 
    EigenVal1(:,:,Slice,:,:)= sum(difComp, 3); %ika2k1206   
    
end %Slice

% **************************save structures  **************************
b0_image_struc= mrstruct_init('volume', double(DTI.dataAy(:,:,1:Max_slice)), DTI);

eigenVal_struc= mrstruct_init('series3D', squeeze(EigenVal1), DTI);

eigenVec_struc= mrstruct_init('series3DEchos', eigVect1, DTI);

error_struc= mrstruct_init('volume', Max_error, DTI);

text = ['save(''', fullfile(OutputDir, FileName), '_DTD.mat'', ''b0_image_struc'', ''eigenVal_struc'', ''eigenVec_struc'', ''error_struc'');'];
eval(text);
Res ='DTD structure is calculated';
