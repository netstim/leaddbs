function absolute_pos = relativePos2absolutePos( relative_pos,parent_position,units )
%RELATIVEPOS2ABSOLUTEPOS returns absolute object position from it's parent
%position and it's position relative to it's parent

if strcmp(units,'normalized')
    scale=parent_position(3:4);
    absolute_pos(1:2)=  relative_pos(1:2).*scale+parent_position(1:2);
    absolute_pos(3:4)=  relative_pos(3:4).*scale;
else
    absolute_pos(1:2)=  relative_pos(1:2)+parent_position(1:2);
    absolute_pos(3:4)=  relative_pos(3:4);
end   

end

