function tests = test_ea_atlasfiber_cache
% Cached render meshes must never cross atlas versions or display settings.
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
repo=fileparts(fileparts(mfilename('fullpath')));
oldPath=path;
addpath(fullfile(repo,'helpers'));
testCase.addTeardown(@() path(oldPath));
end

function setup(testCase)
assumeTrue(testCase,usejava('jvm'),'Persistent cache keys require Java.');
folder=tempname;
mkdir(folder);
testCase.addTeardown(@() rmdir(folder,'s'));
index=fullfile(folder,'atlas_index.mat');
fid=fopen(index,'w');
fwrite(fid,uint8(1:10));
fclose(fid);
testCase.TestData.index=index;
testCase.TestData.folder=fullfile(folder,'cache');
end

function testRoundTripPreservesEveryMeshValue(testCase)
[index,folder,meshes]=fixture(testCase);
[missing,cachefile,signature]=ea_atlasfiber_cache('load',index,0.3,[],folder);
verifyEmpty(testCase,missing);
verifyFalse(testCase,isfile(cachefile));
ea_atlasfiber_cache('save',index,0.3,meshes,folder,signature);
loaded=ea_atlasfiber_cache('load',index,0.3,[],folder);
verifyEqual(testCase,loaded,meshes);
verifyTrue(testCase,isfile(cachefile));
verifyEqual(testCase,numel(dir(fullfile(folder,'*.mat'))),1);
end

function testSettingsAndMatlabVersionMustMatch(testCase)
[index,folder,meshes]=fixture(testCase);
[~,cachefile]=ea_atlasfiber_cache('save',index,0.3,meshes,folder);
verifyEmpty(testCase,ea_atlasfiber_cache('load',index,'auto',[],folder));
verifyEmpty(testCase,ea_atlasfiber_cache('load',index,0.5,[],folder));
verifyEqual(testCase,ea_atlasfiber_cache('load',index,0.3,[],folder),meshes);
cached=load(cachefile);
cached.metadata.matlabVersion='another MATLAB build';
save(cachefile,'-struct','cached');
verifyEmpty(testCase,ea_atlasfiber_cache('load',index,0.3,[],folder));
end

function testSourceSizeAndModificationTimeInvalidateCache(testCase)
[index,folder,meshes]=fixture(testCase);
ea_atlasfiber_cache('save',index,0.3,meshes,folder);
fid=fopen(index,'a');fwrite(fid,uint8(11));fclose(fid);
verifyEmpty(testCase,ea_atlasfiber_cache('load',index,0.3,[],folder));
ea_atlasfiber_cache('save',index,0.3,meshes,folder);
verifyEqual(testCase,ea_atlasfiber_cache('load',index,0.3,[],folder),meshes);
file=java.io.File(index);
verifyTrue(testCase,file.setLastModified(file.lastModified()+1));
verifyEmpty(testCase,ea_atlasfiber_cache('load',index,0.3,[],folder));
end

function testCorruptAndMalformedCachesFallBackToMiss(testCase)
[index,folder,meshes]=fixture(testCase);
[~,cachefile]=ea_atlasfiber_cache('save',index,0.3,meshes,folder);
fid=fopen(cachefile,'w');fwrite(fid,'invalid MAT-file');fclose(fid);
verifyEmpty(testCase,ea_atlasfiber_cache('load',index,0.3,[],folder));
ea_atlasfiber_cache('save',index,0.3,meshes,folder);
cached=load(cachefile);
cached.meshes={struct('vertices',ones(4,2),'faces',[1 2 3])};
save(cachefile,'-struct','cached');
verifyEmpty(testCase,ea_atlasfiber_cache('load',index,0.3,[],folder));
for invalid=1:4
    cached.meshes=meshes;
    if invalid==1
        cached.meshes{1,1}.vertices(1)=NaN;
    else
        cached.meshes{1,1}.faces=double(cached.meshes{1,1}.faces);
        if invalid==2
            cached.meshes{1,1}.faces(1)=size(cached.meshes{1,1}.vertices,1)+1;
        elseif invalid==3
            cached.meshes{1,1}.faces(1)=1.5;
        else
            cached.meshes{1,1}.faces(1)=Inf;
        end
    end
    save(cachefile,'-struct','cached');
    verifyEmpty(testCase,ea_atlasfiber_cache('load',index,0.3,[],folder));
end
save(cachefile,'meshes'); % Old/incomplete files must not emit load warnings.
verifyWarningFree(testCase,@() ea_atlasfiber_cache('load',index,0.3,[],folder));
verifyEmpty(testCase,ea_atlasfiber_cache('load',index,0.3,[],folder));
end

function testUnavailableCacheFolderDoesNotBlockRendering(testCase)
[index,folder,meshes]=fixture(testCase);
fid=fopen(folder,'w');fwrite(fid,'a file occupies this folder path');fclose(fid);
ea_atlasfiber_cache('save',index,0.3,meshes,folder);
verifyEmpty(testCase,ea_atlasfiber_cache('load',index,0.3,[],folder));
end

function testChangedSourceDuringPreparationIsNotSaved(testCase)
[index,folder,meshes]=fixture(testCase);
[~,cachefile,signature]=ea_atlasfiber_cache('load',index,0.3,[],folder);
fid=fopen(index,'a');fwrite(fid,uint8(11));fclose(fid);
ea_atlasfiber_cache('save',index,0.3,meshes,folder,signature);
verifyFalse(testCase,isfile(cachefile));
verifyEmpty(testCase,ea_atlasfiber_cache('load',index,0.3,[],folder));
end

function [index,folder,meshes]=fixture(testCase)
index=testCase.TestData.index;
folder=testCase.TestData.folder;
mesh=struct('vertices',single([0 0 0;1 0 0;0 1 0;0 0 1]),...
    'faces',uint32([1 2 3;1 2 4;1 3 4;2 3 4]));
meshes=cell(3,2);
meshes{1,1}=mesh;
meshes{2,1}=struct('vertices',zeros(0,3),'faces',zeros(0,3));
mesh.vertices=double(mesh.vertices)+[5 -3 1];
meshes{3,2}=mesh;
end
