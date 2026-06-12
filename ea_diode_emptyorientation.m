function angles = ea_diode_emptyorientation(numelectrodes)
angles = struct;
for i = 1:numelectrodes
    %%
    angles.postop(i).yaw = NaN;
    angles.postop(i).pitch = NaN;
    angles.postop(i).polar = NaN;
    angles.postop(i).roll = NaN;
    angles.postop(i).orientation = NaN;
    %%
    angles.native(i).yaw = NaN;
    angles.native(i).pitch = NaN;
    angles.native(i).polar = NaN;
    angles.native(i).roll = NaN;
    angles.native(i).orientation = NaN;
    %%
    angles.mni(i).yaw = NaN;
    angles.mni(i).pitch = NaN;
    angles.mni(i).polar = NaN;
    angles.mni(i).roll = NaN;
    angles.mni(i).orientation = NaN;
    %%
    angles.acpc(i).yaw = NaN;
    angles.acpc(i).pitch = NaN;
    angles.acpc(i).polar = NaN;
    angles.acpc(i).roll = NaN;
    angles.acpc(i).orientation = NaN;
    %%
    angles.orientationmethod{i,1} = 'none';
end

end