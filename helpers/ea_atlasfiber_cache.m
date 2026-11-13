function [meshes,cachefile,signature] = ea_atlasfiber_cache(action,indexFile,hullsimplify,meshes,cacheFolder,expectedSignature)
% Reuse the exact fiber meshes produced for one atlas index and setting.
% A cache miss or unavailable cache directory is harmless. Source atlas data
% is never changed. Optional cacheFolder isolates tests/custom cache storage.
% Pass the third load output as expectedSignature when saving to reject an
% index that changed while its meshes were being prepared.

if nargin<4
    meshes=[];
end
if nargin<5 || isempty(cacheFolder)
    cacheFolder=fullfile(prefdir,'LeadDBS','atlas-render-cache');
end
if nargin<6
    expectedSignature=[];
end
cachefile='';
signature=[];
action=char(action);
if ~any(strcmp(action,{'load','save'}))
    error('ea_atlasfiber_cache:InvalidAction','Unknown cache action: %s',action);
end
if strcmp(action,'load')
    meshes=[];
end

try
    indexFile=char(java.io.File(char(indexFile)).getCanonicalPath());
    indexHandle=java.io.File(indexFile);
    source=dir(indexFile);
    if isempty(source) || source.isdir
        return
    end
    signature=struct('schema',1,'indexFile',indexFile,'bytes',source.bytes,...
        'datenum',source.datenum,'modifiedMilliseconds',int64(indexHandle.lastModified()),...
        'hullsimplify',hullsimplify,...
        'matlabVersion',version,'platform',computer);
    digest=java.security.MessageDigest.getInstance('SHA-256');
    digest.update(unicode2native(indexFile,'UTF-8'));
    key=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
    cachefile=fullfile(cacheFolder,[key,'.mat']);

    if strcmp(action,'load')
        if ~isfile(cachefile)
            return
        end
        contents=whos('-file',cachefile);
        if ~all(ismember({'metadata','meshes'},{contents.name}))
            return
        end
        cached=load(cachefile,'-mat','metadata','meshes');
        if isfield(cached,'metadata') && isequaln(cached.metadata,signature) && ...
                isfield(cached,'meshes') && validmeshes(cached.meshes,true)
            meshes=cached.meshes;
        end
    elseif validmeshes(meshes)
        if ~isempty(expectedSignature) && ~isequaln(expectedSignature,signature)
            return
        end
        if ~isfolder(cacheFolder)
            [created,~]=mkdir(cacheFolder);
            if ~created
                return
            end
        end
        temporary=[tempname(cacheFolder),'.mat'];
        cleanup=onCleanup(@() removetemporary(temporary)); %#ok<NASGU>
        metadata=signature; %#ok<NASGU>
        % Geometry arrays are already prepared; compression would spend CPU
        % on every cache fill and decompression on every subsequent display.
        save(temporary,'metadata','meshes','-v7.3','-nocompression');
        movefile(temporary,cachefile,'f');
    end
catch
    % Cache corruption, read-only preferences, or unavailable Java must not
    % prevent atlas display: the caller can always compute the meshes again.
    if strcmp(action,'load')
        meshes=[];
    end
end
end


function valid=validmeshes(meshes,checkValues)
if nargin<2
    checkValues=false;
end
valid=iscell(meshes) && ~isempty(meshes);
if ~valid
    return
end
for k=1:numel(meshes)
    mesh=meshes{k};
    if isempty(mesh)
        continue
    end
    if ~isstruct(mesh) || ~isscalar(mesh) || ...
            ~isfield(mesh,'vertices') || ~isfield(mesh,'faces') || ...
            ~isnumeric(mesh.vertices) || ~isreal(mesh.vertices) || ...
            ~ismatrix(mesh.vertices) || ...
            size(mesh.vertices,2)~=3 || ~isnumeric(mesh.faces) || ~isreal(mesh.faces) || ...
            ~ismatrix(mesh.faces) || size(mesh.faces,2)<3
        valid=false;
        return
    end
    if checkValues
        % Bound temporary arrays while scanning large cached meshes. An
        % invalid cache must fall back before its data reaches patch().
        batchSize=1000000;
        for first=1:batchSize:numel(mesh.vertices)
            last=min(first+batchSize-1,numel(mesh.vertices));
            if any(~isfinite(mesh.vertices(first:last)))
                valid=false;
                return
            end
        end
        vertexCount=size(mesh.vertices,1);
        for first=1:batchSize:numel(mesh.faces)
            last=min(first+batchSize-1,numel(mesh.faces));
            indices=mesh.faces(first:last);
            if any(~isfinite(indices) | indices<1 | indices>vertexCount | ...
                    indices~=fix(indices))
                valid=false;
                return
            end
        end
    end
end
end


function removetemporary(filename)
if isfile(filename)
    delete(filename);
end
end
