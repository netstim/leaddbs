function coords=ea_findcoords(goodz,trajectory,trajvector,dist,correction,options)

trajectory=ea_genhdtrajectory(trajectory,options);

onecoord=trajectory(abs(trajectory(:,3)-(goodz+correction(3)))<0.001,:);

%% add correction term
orth=null(trajvector);
orth=[orth,trajvector'];
orthx=orth(:,1);
orthy=orth(:,2);

try
    onecoord=[onecoord(1)+orthx(1)*correction(1)+orthy(1)*correction(2); ...
              onecoord(2)+orthx(2)*correction(1)+orthy(2)*correction(2); ...
              onecoord(3)+orthx(3)*correction(1)+orthy(3)*correction(2)]';
catch
    return
end

coords(1,:)=onecoord;

normtrajvector=trajvector./norm(trajvector);

for electrode=2:4
    coords(electrode,:)=coords(1,:)-normtrajvector.*((electrode-1)*dist);
end
