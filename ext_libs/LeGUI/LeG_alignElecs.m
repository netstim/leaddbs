function lead_sm = LeG_alignElecs(lead,d)
%"lead" is xyz coordinates (patient space mm) of contacts of a lead with
%deepest first (ground truth). "d" is intercontact spacing of contacts in
%mm. Spacing is adjusted between contacts to match "d". Next contact is
%placed "d" away from last in direction of average vector of next two
%contacts.

lead_sm = zeros(size(lead));
lead_sm(1,:) = lead(1,:);       % first contact - ground truth
for j = 1:size(lead,1)-1
    if j < length(lead)-1
        a = lead_sm(j,:);
        b = (lead(j+1,:)+lead(j+2,:))./2;
        v = (b-a)/norm(b-a);    % mean vector of next 2 contacts
    else                        % last contact condition
        a = lead_sm(j,:);
        b = lead(j+1,:);
        v = (b-a)/norm(b-a);    % vector to last contact
    end
    n = a + d*v;                % new coords intercontact dist away
    lead_sm(j+1,:) = n;
end