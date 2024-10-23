function ea_gentargetreport(M)

thresh = inputdlg('Please enter threshold in mm:','Enter threshold...',1,{'0.5'});

if isempty(thresh)
    return
end
thresh=str2double(thresh{1});

conInd = arrayfun(@num2str, 1:max(M.S.numContacts), 'Uni', 0);
rnames = [strcat('k', conInd, 'R'), strcat('k', conInd, 'L')];

for thr=0:1
    for target=1:length(M.ui.volumeintersections)
        if thr
            rf(target,thr+1)=figure('color','w','Numbertitle','off','name',['Electrode centers residing in ',M.vilist{M.ui.volumeintersections(target)}]);
        else
            rf(target,thr+1)=figure('color','w','Numbertitle','off','name',['Distances of electrode centers to nearest voxel in ',M.vilist{M.ui.volumeintersections(target)}]);
        end

        % Here we suppose electrodes of the select patients have the same
        % number of contacts
        numContacts = size(M.stats(1).ea_stats.conmat{1},1);
        distances=zeros(numContacts*2,length(M.ui.listselect));

        for pt=1:length(M.ui.listselect)
            try
                distances(:,pt)=[M.stats(M.ui.listselect(pt)).ea_stats.conmat{1}(:,M.ui.volumeintersections(target));... % right side
                    M.stats(M.ui.listselect(pt)).ea_stats.conmat{2}(:,M.ui.volumeintersections(target))];
            catch ME
                if strcmp(ME.identifier, 'MATLAB:subsassigndimmismatch')
                    ea_error('Number of contacts not consistent across patients!')
                elseif strcmp(ME.identifier, 'MATLAB:nonExistentField')
                    ea_error('Please calculate DBS stats for all patients first!');
                else
                    ea_error(ME.message);
                end
            end
        end

        if thr
            distances=distances<thresh;
        end

        R{target,thr+1}=distances;

        for xx=1:size(R{target,thr+1},1)
              for yy=1:size(R{target,thr+1},2)
                  side=ceil((xx/numContacts/2)*2);
                  con=xx+(1-side)*size(R{target,thr+1},1)/2;
                  cstring='FFFFFF';
                  try
                      if logical(M.S(yy).activecontacts{side}(con)) && M.ui.hlactivecontcheck
                          cstring='FF9999';
                      end
                  end
                  C{xx,yy}=['<html><table border=0 width=400 bgcolor=#',cstring,'><TR><TD>',num2str(R{target,thr+1}(xx,yy)),'</TR> </table></html>'];
              end
        end

        cnames=M.patient.list(M.ui.listselect');
        [~,cnames]=cellfun(@fileparts,cnames,'UniformOutput',0);

        if M.ui.hlactivecontcheck
            t(target,thr+1)=uitable(rf(target,thr+1),'Data',C,'ColumnName',cnames,'RowName',rnames);
        else
            t(target,thr+1)=uitable(rf(target,thr+1),'Data',R{target,thr+1},'ColumnName',cnames,'RowName',rnames);
        end
        figdims=get(rf(target,thr+1),'Position');
        textend=get(t(target,thr+1),'Extent');
        set(rf(target,thr+1),'Position',[figdims(1:2),textend(3:4)]);
        figdims=get(rf(target,thr+1),'Position');
        set(t(target,thr+1),'Position',[0,0,figdims(3),figdims(4)])
    end
end

assignin('base','R',R);
