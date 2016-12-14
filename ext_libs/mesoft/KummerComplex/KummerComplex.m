function [chg]=KummerComplex(a,b,z)

%
%KUMMERCOMPLEX   Confluent hypergeometric function 1F1 
% (Kummer function).
%
% KUMMERCOMPLEX(a,b,z) is the confluent hypergeometric 
% function 1F1 (Kummer function) for complex 
% parameters a, b and complex variable z.
%
% Example: 
%    KUMMERCOMPLEX(1+0.5i,2-3.1i,4+2i) 
%     equals   0.33300865268261 - 0.02369621687656i



%
% In general case the program calculates the sum of 
% convergent series defining the function until the next 
% term becomes too small (in comparison with the sum of all
% previous terms). The case of large abs(z) is considered 
% separately (e.g., see 13.5.1 in Abramowitz, Stegun 
% "Handbook of special functions", 1964). Some simple cases
% of integer parameter values are considered separately as 
% well.
%
% The function controls the loss of precision and makes 
% check for insufficient number of members in the series. 
% It prints warning if there are any problems. Otherwise, 
% if everything is ok,  the results seem to coincide with 
% Matematica 4.1 with 10-digit precision.
%
% This function is largely based at "Fortran library of 
% special functions" which was converted to Matlab. 
% Unfortunatey, the library can compute confluent 
% hypergeometric function only for real values of a and b. 
% So this file may be considered as its generalization 
% for complex a and b.
%
% This function also requires cgama.m file which computes 
% Gamma function for complex variables. This file was taken
% from just the same "Fortran library" and insignificantly 
% modified.
%



%% Some initialization
precision = 15;


%% Special cases

if (imag(b)==0 && real(b)<=0 && real(b) == fix(real(b)) )   % b==-n   n=1,2,3,..
    
    if ( imag(a)==0 && real(a)<=0 && real(a) == fix(real(a)) && abs(a)<abs(b) ) % a==-m; m=1,2,..
        m=fix(-a);
        cr=complex(1,0);
        chg=complex(1,0);

        cMax = abs(cr);
        
        for  k=1:m;
            cr=cr.*(a+k-1)./k./(b+k-1).*z;
            chg=chg+cr;
            
            cMax=max( cMax , max(abs(cr),abs(chg)) );            
        end;  
        
        precision = 15-fix(Log10(cMax/abs(chg)));       

    elseif ( imag(a)==0 && real(a)<=0 && real(a) == fix(real(a)) && abs(a)==abs(b) ); % a==b;
        format compact;
        '!!!Confluent hypergeometric function is indeterminate!!!'  %, a,b,z
        chg='Error';
        return;
      
    else         
        chg=complex(NaN,NaN);
    end 

elseif (a==0 || z==0)
    chg=1;
    
elseif (a == -1);
    chg=1-z./b;

elseif (a == b);
    chg=exp(z);

elseif ( (a-b) == 1);
    chg=(1+z./b).*exp(z);    
    
elseif ( a==1 && b==2 )
    chg=(exp(z)-1)./z;

    % finite number of elements in a row
elseif ( imag(a)==0 && real(a)<0 && real(a)==fix(real(a)) )    
    m=fix(-a);
    cr=complex(1,0);
    chg=complex(1,0);

    cMax = abs(cr);

    
    for  k=1:m;
        cr=cr.*(a+k-1)./k./(b+k-1).*z;
        chg=chg+cr;
        
        cMax=max( cMax , max(abs(cr),abs(chg)) );            
    end;  
    
    precision = 15-fix(Log10(cMax/abs(chg)));       

elseif  ( abs(z)>10.*abs(a) && abs(z)>10.*abs(b) )      %Abramowitz Stegun 13.5.1.
%%%%%%%%%%
    [g1_real,g1_imag]=cgama(real(a),imag(a),1);  
    g1=complex(g1_real,g1_imag);
    [g2_real,g2_imag]=cgama(real(b),imag(b),1);  
    g2=complex(g2_real,g2_imag);
    ba=b-a;
    [g3_real,g3_imag]=cgama(real(ba),imag(ba),1);  
    g3=complex(g3_real,g3_imag);

    cs1=complex(1.0d0,0.0d0);
    cs2=complex(1.0d0,0.0d0);
    cr1=complex(1.0d0,0.0d0);
    cr2=complex(1.0d0,0.0d0);

    c1Max = abs(cr1);
    c2Max = abs(cr2);
            
    for  j=1:500;
        cr1=-cr1.*(a+j-1.0d0).*(a-b+j)./(z.*j);  
        cr2=cr2.*(b-a+j-1.0d0).*(j-a)./(z.*j);
        cs1=cs1+cr1;
        cs2=cs2+cr2;

        c1Max=max(c1Max,max(abs(cr1),abs(cs1)));
        c2Max=max(c2Max,max(abs(cr2),abs(cs2)));
        
        if ( abs(cr1)/abs(cs1) < 1.d-15 && abs(cr2)/abs(cs2) < 1.d-15 )
            %j
            break; 
        end;
                
        if (j==500) 
            ['Got to the  ' num2str(j) '  limit in the series of confluent hypergeometric function!']
            cs1 = 'Error';
            cs2 = 'Error';            
            chg = 'Error';            
            return;
        end;
                
    end;  
    
    precision = 15-fix(  Log10(  max( c1Max/abs(cs1) , c2Max/abs(cs2) )  )  );
   
    x=real(z);
    y=imag(z);

    if(x == 0.0 && y >= 0.0);
           phi=0.5.*pi;
                
    elseif(x == 0.0 && y <= 0.0);
            phi=-0.5.*pi;
                
    else
            phi=atan(y./x);
    end;
            
    if(phi > -0.5*pi && phi < 1.5*pi)
            ns=1; 
    end;
            
    if(phi > -1.5*pi && phi <= -0.5*pi)
            ns=-1; 
    end;

    cfac=exp(ns.*complex(0,1).*pi.*a);
            
    if(y == 0)
            cfac=cos(pi.*a); 
    end;

    chg1=g2./g3.*z.^(-a).*cfac.*cs1;
    chg2=g2./g1.*exp(z).*z.^(a-b).*cs2;
    chg=chg1+chg2;
             
%%%%%%%%%%
        
    
    
else    % General case
            chg=complex(1.0,0.0);
            crg=complex(1.0,0.0);
            cgMax=abs(crg);
        
            for  j=1:500;
                crg=crg.*(a+j-1.0d0)./(j.*(b+j-1.0d0)).*z;   % Abramowitz Stegun 13.1.2 
                chg=chg+crg;
                
                cgMax = max(cgMax,max(abs(crg),abs(chg)));
                                
                if(abs(crg)/abs(chg)< 1.d-15)
                    %j
                    break; 
                end;
                
                if (j==500) 
                       ['Got to the  ' num2str(j) '  limit in the series of confluent hypergeometric function!'] 
                       chg = 'Error';
                       return;
                end;
                
            end;
            
            precision = 15-fix(Log10(cgMax/abs(chg)));
   
end;    
    
if (precision<=0)
    precision=0;
    chg='Error';
end;

if (precision<10) 
    ['!!! Warning!!! Only about  ' num2str(precision) '  first digits are correct!!!']
end;

return;