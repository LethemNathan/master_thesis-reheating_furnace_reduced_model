function marker = read_marker_nonParticipatingMedia(filename)
% READ_BLOCKED_MARKER_MSH Read a Gmsh MSH 4.1 ASCII marker file
% Stores only PhysicalNames + 2D/3D Entities with bbox and physical tags.
% Output bboxes are converted from mm to m.



    fid = fopen(filename, 'r');
    cleanup = onCleanup(@() fclose(fid));

    lengthScale = 1e-3;     %input mesh in mm -> output in m

    %Initialize Output
    marker = struct();
    marker.physicalNames = containers.Map('KeyType','int32','ValueType','char');
    marker.surfaces = struct(...
        'physicalTag',{},...
        'bbox',{});
    marker.volumes  = struct(...
        'physicalTag',{},...
        'bbox',{});

    % Main Loop
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if ~ischar(line)
            break
        end

        switch line
            case '$PhysicalNames'
                marker.physicalNames = parsePhysicalNames(fid);

            case '$Entities'
                [marker.surfaces, marker.volumes] = parseEntities(fid, lengthScale);
            otherwise
            % ignore unknown or empty lines
        end
    end
end


%% LOCAL FUNCTIONS

function physicalNames = parsePhysicalNames(fid)
    physicalNames = containers.Map('KeyType','int32','ValueType','char');
    n = sscanf(strtrim(fgetl(fid)), '%d', 1);

    for i = 1:n
        line = strtrim(fgetl(fid));
        tok = regexp(line, '^(\d+)\s+(\d+)\s+"(.*)"$', 'tokens', 'once');
        ptag = int32(str2double(tok{2}));
        pname = tok{3};
        physicalNames(ptag) = pname;
    end
end


function [surfaces, volumes] = parseEntities(fid, lengthScale)
    counts = sscanf(strtrim(fgetl(fid)), '%d');

    nSurfaces = counts(3);
    nVolumes  = counts(4);

    surfaces = struct(...
        'physicalTag',{},...
        'bbox',{});
    volumes  = struct(...
        'physicalTag',{},...
        'bbox',{});

    % Ignore points and curves
    for i = 1:(sum(counts)-nSurfaces-nVolumes)
        fgetl(fid);
    end

    % surfaces
    for i = 1:nSurfaces
        vals = sscanf(strtrim(fgetl(fid)), '%f').';
        bbox = vals(2:7) * lengthScale;   % mm -> m
        nPhys = vals(8);
        
        if nPhys > 0
            ptag = int32(vals(9));
            surfaces(end+1).physicalTag = ptag;
            surfaces(end).bbox = bbox;
        end
    end

    % volumes
    for i = 1:nVolumes
        vals = sscanf(strtrim(fgetl(fid)), '%f').';
        bbox = vals(2:7) * lengthScale;   % mm -> m
        nPhys = vals(8);

        if nPhys == 0
            error('Not all volumes are assigned: some volumes have no physical tag.');
        else
            ptag = int32(vals(9));
        end

        volumes(end+1).physicalTag = ptag;
        volumes(end).bbox = bbox;
    end
end