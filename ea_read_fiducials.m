function coordsmm=ea_read_fiducials(path,options)
% Small function to read fiducials (coordinates of electrode contacts) in
% Slicer3 format.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
%
% Andreas Horn

fid=fopen(path);

C=textscan(fid,'%s %f %f %f %f %f','commentStyle', '#','delimiter', ',');

coordsmm=zeros(8,3);
for contact=1:2*options.elspec.numContacts
    coordsmm(contact,1)=C{2}(contact);
    coordsmm(contact,2)=C{3}(contact);
    coordsmm(contact,3)=C{4}(contact);
end

