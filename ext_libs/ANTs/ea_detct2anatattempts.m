function cnt=ea_detct2anatattempts(volumedir)
% determines how many coregistration attempts have been made so far..
cnt=0;
while 1
    if exist([volumedir,'ct2anat',num2str(cnt+1),'.txt'],'file')
        cnt=cnt+1;
    else
        break
    end
end
cnt=cnt+1;