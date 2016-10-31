function [quatern, qfac, pixdim]=cbiHomogeneousToQuaternion( R )
% Calculates quaternion parameters from a 4x4 matrix
% [quatern, qfac, pixdim]=cbiHomogeneousToQuaternion( M44 )
% 
% - OR - (using a 3x3 rotation matrix)
% [quatern, qfac, pixdim]=cbiHomogeneousToQuaternion( R )
% 
% Modified from nifti_io.c

% check args
    
  if (nargin~=1)
    error('Wrong number of input arguments!');
  end

  if (size(R)==[4,4] | size(R)==[3,3])
    R=R(1:3,1:3);
  else
    error('Input argument must be a homogeneous 4x4 matrix or 3x3 rotation matrix');
  end
  
  if (nargout~=3)
    error('Wrong number of output arguments!');
  end

  [quatern,qfac, pixdim]=calcQuatFromRotMat(R);

  return
  
function [quatern,qfac,pixdim]=calcQuatFromRotMat(R)
  % BEGIN calc
  % assume R set
  
  %  /* compute lengths of each column; these determine grid spacings  */
 

  xd = sqrt( sum(R(:,1).^2) );
  yd = sqrt( sum(R(:,2).^2) );
  zd = sqrt( sum(R(:,3).^2) );
 
  
% $$$    /* if a column length is zero, patch the trouble */
  
  if ( xd == 0 )
    R(1,1) = 1 ; R(2,1) = 0; R(3,1) = 0 ; xd = 1 ;
  end
  if ( yd == 0 )
    R(2,2) = 1 ; R(1,2) = 0; R(3,2) = 0 ; yd = 1 ; 
  end
  if ( zd == 0 )
    R(3,3) = 1 ; R(1,3) = 0; R(2,3) = 0 ; zd = 1 ; 
  end
  pixdim=[xd yd zd];
% $$$    /* normalize the columns */
    
  R(:,1)=R(:,1)./xd;
  R(:,2)=R(:,2)./yd;
  R(:,3)=R(:,3)./zd;
  
  R=polarDecomp33(R);
  zd=det(R);
  
  if (zd>0)
    qfac=1;
  else
    qfac=-1;
    R(:,3)=-R(:,3);
  end
  
  a=R(1,1)+R(2,2)+R(3,3)+1;
  if (a>0.5)
    a = 0.5 * sqrt(a) ;
    b = 0.25 * (R(3,2)-R(2,3)) / a ;
    c = 0.25 * (R(1,3)-R(3,1)) / a ;
    d = 0.25 * (R(2,1)-R(1,2)) / a ;
  else
    xd = 1.0 + R(1,1) - (R(2,2)+R(3,3)) ;
    yd = 1.0 + R(2,2) - (R(1,1)+R(3,3)) ;
    zd = 1.0 + R(3,3) - (R(1,1)+R(2,2)) ;
    if( xd > 1.0 )
      b = 0.5 * sqrt(xd) ;
      c = 0.25* (R(1,2)+R(2,1)) / b ;
      d = 0.25* (R(1,3)+R(3,1)) / b ;
      a = 0.25* (R(3,2)-R(2,3)) / b ;
    elseif ( yd > 1.0 )
      c = 0.5 * sqrt(yd) ;
      b = 0.25* (R(1,2)+R(2,1)) / c ;
      d = 0.25* (R(2,3)+R(3,2)) / c ;
      a = 0.25* (R(1,3)-R(3,1)) / c ;
    else
      d = 0.5 * sqrt(zd) ;
      b = 0.25* (R(1,3)+R(3,1)) / d ;
      c = 0.25* (R(2,3)+R(3,2)) / d ;
      a = 0.25* (R(2,1)-R(1,2)) / d ;
    end
    if( a < 0.0 )
      b=-b; c=-c; d=-d; a=-a; 
    end
    
  end
  quatern=[b c d];
  
  return;

  
  
function OM=polarDecomp33( M );
  
   X = M;

   gam = det(X) ;
   while ( gam == 0.0 )
     gam = 0.00001 * ( 0.001 + norm(X,'inf'));
     X(1,1)=X(1,1)+gam;
     X(2,2)=X(2,2)+gam;
     X(3,3)=X(3,3)+gam;
     gam = det(X);
   end
   dif=1;
   k=0;
   while (1)
     Y = inv(X);
     if( dif > 0.3 )
       alp = sqrt( norm(X,'inf') * norm(X','inf'));
       bet = sqrt( norm(Y,'inf') * norm(Y','inf'));
       gam = sqrt( bet / alp ) ;
       gmi = 1.0 / gam ;
     else
       gam = 1.0;
       gmi = 1.0 ;
     end

     Z(1,1) = 0.5 * ( gam*X(1,1) + gmi*Y(1,1) ) ;     
     Z(1,2) = 0.5 * ( gam*X(1,2) + gmi*Y(2,1) ) ;
     Z(1,3) = 0.5 * ( gam*X(1,3) + gmi*Y(3,1) ) ;
     Z(2,1) = 0.5 * ( gam*X(2,1) + gmi*Y(1,2) ) ;
     Z(2,2) = 0.5 * ( gam*X(2,2) + gmi*Y(2,2) ) ;
     Z(2,3) = 0.5 * ( gam*X(2,3) + gmi*Y(3,2) ) ;
     Z(3,1) = 0.5 * ( gam*X(3,1) + gmi*Y(1,3) ) ;
     Z(3,2) = 0.5 * ( gam*X(3,2) + gmi*Y(2,3) ) ;
     Z(3,3) = 0.5 * ( gam*X(3,3) + gmi*Y(3,3) ) ;

     dif = abs(Z(1,1)-X(1,1))+abs(Z(1,2)-X(1,2)) ...
	   +abs(Z(1,3)-X(1,3))+abs(Z(2,1)-X(2,1)) ...
	   +abs(Z(2,2)-X(2,2))+abs(Z(2,3)-X(2,3)) ...	   
	   +abs(Z(3,1)-X(3,1))+abs(Z(3,2)-X(3,2)) ...
	   +abs(Z(3,3)-X(3,3));

     k = k+1 ;
     if ( k > 100 | dif < 3.e-6 ) 
       break ;
     end
     X = Z ;
   end

   OM=Z;
   return


