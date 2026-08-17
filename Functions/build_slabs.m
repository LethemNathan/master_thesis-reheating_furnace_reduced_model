function slabs = build_slabs(slabsBbox)

    ids = cell2mat(keys(slabsBbox));
    ids = sort(ids);

    slabs = struct('id', {}, 'bbox', {}, 'surfaces', {});

    for i = 1:numel(ids)
        id = ids(i);
        bbox = slabsBbox(id);

        slabs(i) = build_slabs_from_bbox(id, bbox);
    end
end


function slabs = build_slabs_from_bbox(id, bbox)

    xmin = bbox(1); ymin = bbox(2); zmin = bbox(3);
    xmax = bbox(4); ymax = bbox(5); zmax = bbox(6);

    % Map des physical names (UNIQUEMENT pour les surfaces)
    surfacePhysicalNames = containers.Map('KeyType', 'int32', 'ValueType', 'char');

 
    surfaceData = struct('faceTag', {}, 'bbox', {});

    surfaceDefs = {
        1, 'symmetry', [xmin ymin zmin xmax ymin zmax];  % ymin
        5, 'side',     [xmin ymax zmin xmax ymax zmax];  % ymax
        4, 'rear',     [xmax ymin zmin xmax ymax zmax];  % xmax
        6, 'front',    [xmin ymin zmin xmin ymax zmax];  % xmin
        8, 'top',      [xmin ymin zmax xmax ymax zmax];  % zmax
        2, 'bottom',   [xmin ymin zmin xmax ymax zmin];  % zmin
    };

    for k = 1:size(surfaceDefs, 1)
        yy = surfaceDefs{k, 1};
        faceName = surfaceDefs{k, 2};
        faceBbox = surfaceDefs{k, 3};

        tag = makeTag(id, yy);
        physName = sprintf('wall_slab-%d-%s', id, faceName);

        surfaceData(k).faceTag = tag;
        surfaceData(k).bbox = faceBbox;

        surfacePhysicalNames(tag) = physName;
    end

    % Structure surfaces
    surfaces = struct();
    surfaces.physicalNames = surfacePhysicalNames;
    surfaces.faceTag = [surfaceData.faceTag]';
    surfaces.bbox = {surfaceData.bbox}';

    % Structure slab (SANS physicalNames)
    slabs = struct();
    slabs.id = id;
    slabs.bbox = bbox;
    slabs.surfaces = surfaces;
end


function tag = makeTag(id, yy)
    tag = int32(110000 + double(id)*100 + yy);
end