function mesh = read_mesh(filename)
%READ_MESH Read a Gmsh .msh file (MSH 4.1 ASCII)
%
% Output:
%   mesh.nodes           [Nn x 3] coordinates in meters
%   mesh.cells           struct with volume elements
%   mesh.boundaryFaces   struct with boundary surface elements
%   mesh.physicalNames   containers.Map (key = physical tag, value = name)



    fid = fopen(filename, 'r');
    cleanup = onCleanup(@() fclose(fid));

    lengthScale = 1e-3;     %input mesh in mm -> output in m

    %Initialize Output
    mesh = struct();
    mesh.nodes = [];
    mesh.cells = struct( ...
        'id', [], ...
        'physicalTag', [], ...
        'conn', [] ...
    );
    mesh.boundaryFaces = struct( ...
        'id', [], ...
        'physicalTag', [], ...
        'conn', [] ...
    );

    mesh.physicalNames = containers.Map('KeyType', 'int32', 'ValueType', 'char');

    % Main Loop
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if ~ischar(line)
            break;
        end

        switch line
            case '$PhysicalNames'
                mesh.physicalNames = parsePhysicalNames(fid);
            case '$Entities'
                entityPhys = parseEntities(fid);
            case '$Nodes'
                mesh.nodes = parseNodes(fid, lengthScale);
            case '$Elements'
                [mesh.cells, mesh.boundaryFaces] = parseElements(fid, entityPhys);
            otherwise
                % ignore unknown or empty lines
        end
    end
end


%% LOCAL FUNCTIONS

function physicalNames = parsePhysicalNames(fid)

    physicalNames = containers.Map('KeyType', 'int32', 'ValueType', 'char');
   
    n = sscanf(strtrim(fgetl(fid)), '%d', 1);
    for i = 1:n
        line = strtrim(fgetl(fid));
        tok = regexp(line, '^(\d+)\s+(\d+)\s+"(.*)"$', 'tokens', 'once');
        physicalTag = int32(str2double(tok{2}));
        name = tok{3};
        physicalNames(physicalTag) = name;
    end
end

function entityPhys = parseEntities(fid)
    % entityPhys{dim-1}(entityTag) = [physicalTags...]

    entityPhys = {
        containers.Map('KeyType','int32','ValueType','any')  % dim 2
        containers.Map('KeyType','int32','ValueType','any')  % dim 3
    };

    counts = sscanf(strtrim(fgetl(fid)), '%d');
    nPoints   = counts(1);
    nCurves   = counts(2);
    nSurfaces = counts(3);
    nVolumes  = counts(4);

    % Points
    for i = 1:nPoints
        fgetl(fid);
    end

    % Curves
    for i = 1:nCurves
        fgetl(fid);
    end

    % Surfaces
    for i = 1:nSurfaces
        vals = sscanf(strtrim(fgetl(fid)), '%f').';
        tag = int32(vals(1));
        nPhys = vals(8);
        if nPhys > 1
            error('An entity has multiple physicalNames');
        elseif nPhys==0
            error('An entity has no physicalName')
        end
        phys = int32(vals(9));
        entityPhys{1}(tag) = phys;
    end

    % Volumes
    for i = 1:nVolumes
        vals = sscanf(strtrim(fgetl(fid)), '%f').';
        tag = int32(vals(1));
        nPhys = vals(8);
        if nPhys > 1
            error('An entity has multiple physicalNames');
        elseif nPhys==0
            error('An entity has no physicalName')
        end
        phys = int32(vals(9));
        entityPhys{2}(tag) = phys;
    end
end

function nodes = parseNodes(fid, lengthScale)
    hdr = sscanf(strtrim(fgetl(fid)), '%d');

    numEntityBlocks = hdr(1);
    numNodes        = hdr(2);
    maxNodeTag      = hdr(4);

    nodes = nan(maxNodeTag, 3);

    countedNodes = 0;

    for b = 1:numEntityBlocks
        bh = sscanf(strtrim(fgetl(fid)), '%d');

        numNodesInBlock = bh(4);

        nodeTags = zeros(numNodesInBlock,1);
        for i = 1:numNodesInBlock
            nodeTags(i) = sscanf(strtrim(fgetl(fid)), '%d', 1);
        end

        for i = 1:numNodesInBlock
            vals = sscanf(strtrim(fgetl(fid)), '%f').';
            xyz = vals(1:3) * lengthScale;   % mm -> m
            nodes(nodeTags(i), :) = xyz;
        end

        countedNodes = countedNodes + numNodesInBlock;
    end

    if countedNodes ~= numNodes
        warning('Node count mismatch');
    end
end

function [cells, boundaryFaces] = parseElements(fid, entityPhys)
    hdr = sscanf(strtrim(fgetl(fid)), '%d');

    numPhysicalBlocks = hdr(1);
    numElements = hdr(2);

    cell_id = [];
    cell_phys = [];
    cell_conn = {};

    face_id = [];
    face_phys = [];
    face_conn = {};

    countedElems = 0;

    for b = 1:numPhysicalBlocks
        bh = sscanf(strtrim(fgetl(fid)), '%d');

        entityDim        = bh(1);
        entityTag      = int32(bh(2));
        physicalTag = entityPhys{entityDim-1}(entityTag);
        numElemsInBlock  = bh(4);

        for i = 1:numElemsInBlock
            vals = sscanf(strtrim(fgetl(fid)), '%d').';
            elemTag = vals(1);
            conn = vals(2:end);

            if entityDim == 3
                cell_id(end+1,1) = elemTag; 
                cell_phys(end+1,1) = physicalTag; 
                cell_conn{end+1,1} = conn;

            elseif entityDim == 2
                face_id(end+1,1) = elemTag; 
                face_phys(end+1,1) = physicalTag;
                face_conn{end+1,1} = conn;
            end
        end

        countedElems = countedElems + numElemsInBlock;
    end

    if countedElems ~= numElements
        warning('Element count mismatch');
    end

    cells = struct();
    cells.id = cell_id;
    cells.physicalTag = cell_phys;
    cells.conn = cell_conn;

    boundaryFaces = struct();
    boundaryFaces.id = face_id;
    boundaryFaces.physicalTag = face_phys;
    boundaryFaces.conn = face_conn;
end





