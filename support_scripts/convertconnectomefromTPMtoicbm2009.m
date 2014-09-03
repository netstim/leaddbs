function convertconnectomefromTPMtoicbm2009(connfile)




tmat=[    1.0380   -0.0055    0.0021   -0.9238
    0.0080    1.0394   -0.0045    1.0272
   -0.0057   -0.0101    1.0162    0.1113
         0         0         0    1.0000];


load([connfile])

% convert coords_mm



for fib=1:length(gibbsconnectome)
   gibbsconnectome{fib}=tmat*[gibbsconnectome{fib},ones(size(gibbsconnectome{fib},1),1)]';
    gibbsconnectome{fib}=gibbsconnectome{fib}(1:3,:)';
end

save([connfile],'gibbsconnectome','-v7.3');