function hdtrajectory=ea_genhdtrajectory(trajectory,options)

hdtrajectory(:,1)=interp1q([1:length(trajectory)]',trajectory(:,1),[1:1/options.zresolution:length(trajectory)]');
hdtrajectory(:,2)=interp1q([1:length(trajectory)]',trajectory(:,2),[1:1/options.zresolution:length(trajectory)]');
hdtrajectory(:,3)=interp1q([1:length(trajectory)]',trajectory(:,3),[1:1/options.zresolution:length(trajectory)]');