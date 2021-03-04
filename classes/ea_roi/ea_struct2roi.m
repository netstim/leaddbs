function roi=ea_struct2roi(s,resultfig,ht)

s.plotFigureH=resultfig;
s.htH=ht;
roi=ea_roi;
props = fieldnames(s);
for p = 1:numel(props)
    if ~ismember(props{p},{'controlH','plotFigureH','patchH','toggleH','htH'})
        roi.(props{p})=s.(props{p});
    end
end
