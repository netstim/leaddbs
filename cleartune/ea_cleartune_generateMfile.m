function S = ea_cleartune_generateMfile(rightsets,leftsets,S,volcur)

rightsets(isnan(rightsets)) = 0;
leftsets(isnan(leftsets)) = 0;

% Generates M. file out of stimulation data
S.label = 'setting';

% Stim OFF
fields1 = fieldnames(S);
for source = 2:9
    fields2 = fieldnames(S.(fields1{source}));
    for contact = 1:8
        S.(fields1{source}).(fields2{contact}).perc = 0;
        S.(fields1{source}).(fields2{contact}).pol = 0;
        S.(fields1{source}).(fields2{contact}).imp = 1;
    end
    S.(fields1{source}).amp = 0;
    S.(fields1{source}).va = volcur;             % mA
    S.(fields1{source}).case.perc = 100;
    S.(fields1{source}).case.pol = 2;
    S.(fields1{source}).case.imp = 1;
end
S.activecontacts{1,1} = [0,0,0,0,0,0,0,0];
S.activecontacts{1,2} = [0,0,0,0,0,0,0,0];
S.amplitude{1,1} = [0,0,0,0];
S.amplitude{1,2} = [0,0,0,0];

% Stim ON
S.Rs1.amp = rightsets(1);
S.Ls1.amp = leftsets(1);
S.active = [1,1];
S.amplitude{1,1} = [rightsets(1),0,0,0];
S.amplitude{1,2} = [leftsets(1),0,0,0];
S.activecontacts{1,1} = abs(rightsets(2:end))>0;
S.activecontacts{1,2} = abs(leftsets(2:end))>0;

% rsum = sum(rightsets(2:end));
% lsum = sum(leftsets(2:end));
% S.Rs1.case.perc = abs(-rsum);
% S.Ls1.case.perc = abs(-lsum);
% if rsum<0
%     S.Rs1.case.pol = 2;
% elseif rsum>0
%     S.Rs1.case.pol = 1;
% else
%     S.Rs1.case.pol = 0;
% end
% 
% if lsum<0
%     S.Ls1.case.pol = 2;
% elseif lsum>0
%     S.Ls1.case.pol = 1;
% else
%     S.Ls1.case.pol = 0;
% end

% %for monopolar cases, pol is always 2 and case is always positive
% S.Rs1.case.perc = 100;
% S.Rs1.case.pol = 2;
% S.Ls1.case.perc = 100;
% S.Rs1.case.pol = 2;

% compute case perc and polarity from amp and other contacts
S.Rs1.case.perc = -1.0 * sum(rightsets(2:end));
if S.Rs1.case.perc > 0.0
    S.Rs1.case.pol = 2;  % anode
elseif S.Rs1.case.perc < 0.0
    S.Rs1.case.pol = 1;  % cathode
else
    S.Rs1.case.pol = 0;
end

S.Ls1.case.perc = -1.0 * sum(leftsets(2:end));
if S.Ls1.case.perc > 0.0
    S.Ls1.case.pol = 2;  % anode
elseif S.Ls1.case.perc < 0.0
    S.Ls1.case.pol = 1;  % cathode
else
    S.Ls1.case.pol = 0;
end



for cont = 1:8
    kr = ['k' num2str(cont-1)];
    S.Rs1.(kr).perc = abs(rightsets(cont+1));

    if rightsets(cont+1) > 0
        S.Rs1.(kr).pol = 2;
    elseif rightsets(cont+1) < 0
        S.Rs1.(kr).pol = 1;
    else
        S.Rs1.(kr).pol = 0;
    end

    kl = ['k' num2str(cont+7)];
    S.Ls1.(kl).perc = abs(leftsets(cont+1));

    if leftsets(cont+1) > 0
        S.Ls1.(kl).pol = 2;
    elseif leftsets(cont+1) < 0
        S.Ls1.(kl).pol = 1;
    else
        S.Ls1.(kl).pol = 0;
    end
end
