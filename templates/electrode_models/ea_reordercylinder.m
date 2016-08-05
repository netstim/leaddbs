function b=ea_reordercylinder(a,shift)
% highly specialized function.. this assumes zip-like ordering of vertices

if ~exist('shift','var')
    shift=1;
end

   dict=[];
   for d=shift:-1:1
      dict=[dict,d:d:length(a.vertices)];  
   end
   
   b=a;
   b.vertices=a.vertices(dict,:);
   for dim=1:size(b.faces,2)
       [~, ix]=ismember(b.faces(:,dim),dict);
       b.faces(:,dim)=ix;
   end
   
   
