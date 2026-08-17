function [geom, mesh] = build_topology(mesh, marker)
%BUILD_TOPOLOGY Build face-based topology and geometry from mesh
%
% Inputs:
%   mesh.nodes         [Nn x 3]
%   mesh.cells         struct from read_mesh.m
%   mesh.boundaryFaces struct from read_mesh.m
%   marker             marker .msh parsed by read_blocked_marker_msh.m
%                      pass [] if not used
%
% Outputs:
%   geom.faceArea      [Nf x 1]
%   geom.faceNormal    [Nf x 3] oriented owner->neighbour for internal faces,
%                      outward for boundary faces
%   geom.owner         [Nf x 1]
%   geom.neighbour     [Nf x 1] = 0 for boundary face
%   geom.cellVolume    [Nc x 1]
%   geom.cellCenter    [Nc x 3]
%
% Additional fields:
%   geom.faceCenter    [Nf x 3]
%   geom.faceNodes     cell(Nf,1)
%   geom.faceType      [Nf x 1] 2=tri, 3=quad
%   geom.boundaryTag   [Nf x 1] boundary physical tag, 0 otherwise
%
% Notes:
%   - If marker is provided, blocked-off volumes are projected onto cells:
%       * mesh.cells.physicalTag is updated
%       * mesh.cells.entityTag is updated
%       * mesh.cells.isBlocked is created

    nodes = mesh.nodes;
    Nc = numel(mesh.cells.id);

    % Cell geometric properties
    cellCenter = zeros(Nc, 3);
    cellVolume = zeros(Nc, 1);

    for c = 1:Nc
        conn = mesh.cells.conn{c};
        pts = nodes(conn, :);

        [Vc, Cc] = polyhedronVolumeCentroid(pts);
        cellVolume(c) = Vc;
        cellCenter(c, :) = Cc;
    end


    % Relabel mesh from marker volumes
    geomTmp = struct();
    geomTmp.cellCenter = cellCenter;
    geomTmp.cellVolume = cellVolume;
    mesh = relabel_blocked_cells(mesh, geomTmp, marker);


    % Build all faces from volume cells
    faceMap = containers.Map('KeyType', 'char', 'ValueType', 'int32');

    faceNodes = {};
    owner = [];
    neighbour = [];

    for c = 1:Nc
        conn = mesh.cells.conn{c};
        localFaces = {
                [4 3 2 1]
                [1 5 8 4]
                [5 6 7 8]
                [2 3 7 6]
                [8 7 3 4]
                [1 2 6 5]
            };

        for lf = 1:numel(localFaces)
            fn = conn(localFaces{lf});   % oriented according to local cell convention
            key = faceKey(fn);           % canonical key for matching faces

            if ~isKey(faceMap, key)
                f = int32(numel(faceNodes) + 1);
                faceMap(key) = f;

                faceNodes{f,1} = fn(:).';
                owner(f,1) = c;
                neighbour(f,1) = 0;

            else
                f = faceMap(key);
                if neighbour(f) ~= 0
                    error('A face is shared by more than two cells');
                end
                neighbour(f) = c;
            end
        end
    end

    Nf = numel(faceNodes);


    % Compute face center, area, normal

    faceCenter = zeros(Nf, 3);
    faceArea = zeros(Nf, 1);
    faceNormal = zeros(Nf, 3);

    for f = 1:Nf
        fn = faceNodes{f};
        pts = nodes(fn, :);

        [Af, Cf, nf] = polygonAreaCenterNormal(pts);

        co = cellCenter(owner(f), :);

        if neighbour(f) ~= 0
            cn = cellCenter(neighbour(f), :);
            d = cn - co;

            if dot(nf, d) < 0
                %nf = -nf;
                error('convention de signe incompatible');
            end
        else
            % Boundary face: outward from owner
            d = Cf - co;
            if dot(nf, d) < 0
                %nf = -nf;
                error('convention de signe incompatible');
            end
        end

        faceArea(f) = Af;
        faceCenter(f, :) = Cf;
        faceNormal(f, :) = nf;
    end

    % Face tags:
    %    - Face category (0 internal, 1 interface/boundary)
    %    - FaceTag (physicalTag of the face when faceCategory == 1)

    faceTag = zeros(Nf, 1, 'int32');
    faceCategory = zeros(Nf, 1, 'uint8');

    Nbf = numel(mesh.boundaryFaces.id);
    for bf = 1:Nbf
        fn  = mesh.boundaryFaces.conn{bf};
        key = faceKey(fn);

        if ~isKey(faceMap, key)
            error('Boundary face #%d (element id %d) does not match any face built from volume cells.', ...
                bf, mesh.boundaryFaces.id(bf));
        end

        f = faceMap(key);

        % A boundary face must have no neighbour (only one owner cell)
        if neighbour(f) ~= 0
            error('Face #%d is tagged as a boundary face but is shared by two cells (owner=%d, neighbour=%d).', ...
                f, owner(f), neighbour(f));
        end

        faceCategory(f) = 1;
        faceTag(f)      = mesh.boundaryFaces.physicalTag(bf);
    end

    % Consistency check:
    for f = 1:Nf
        if neighbour(f) == 0 && faceCategory(f) ~= 1
            error('Warning: a face is detected as a boundary but has not been assigned a boundaryTag (face index %d, owner cell %d).', ...
                f, owner(f));
        end
    end


    % Loop over all faces to detect interfaces involving blocked-off cells
    % (non-participating medium) 
    for f = 1:Nf
        ownerC = owner(f);
        neighC = neighbour(f);

        ownerBlocked = mesh.cells.isBlocked(ownerC);

        if neighC == 0
            % Domain boundary face (only an owner cell exists)
            isInterface = true;
            if ownerBlocked
                blockedCell = ownerC;
            else
                blockedCell = 0;   % boundary not adjacent to a blocked cell
            end
        else
            neighBlocked = mesh.cells.isBlocked(neighC);
            if xor(ownerBlocked, neighBlocked)
                % Exactly one side is blocked -> interface face
                isInterface = true;
                if ownerBlocked
                    blockedCell = ownerC;
                else
                    blockedCell = neighC;
                end
            else
                % Both blocked or both non-blocked -> not an interface
                isInterface = false;
                blockedCell = 0;
            end
        end

        if isInterface
            faceCategory(f) = 1;

            if blockedCell ~= 0
                % Retrieve the physical tag of the blocked cell and the
                % corresponding physical name from the mesh
                blockedPhysTag  = mesh.cells.physicalTag(blockedCell);
                blockedPhysName = mesh.physicalNames(blockedPhysTag);

                newFaceTag = getInterfaceFaceTag(blockedPhysName, faceCenter(f, :), faceNormal(f,:), marker);

                faceTag(f) = newFaceTag;
            end
        end
    end



    % Output
    geom = struct();
    geom.faceArea = faceArea;
    geom.faceNormal = faceNormal;
    geom.owner = owner;
    geom.neighbour = neighbour;
    geom.cellVolume = cellVolume;
    geom.cellCenter = cellCenter;
    geom.faceCenter = faceCenter;
    geom.faceNodes = faceNodes;
    geom.faceTag = faceTag;
    geom.faceCategory = faceCategory;
    geom.isLoadFace = false(Nf,1);
    geom.loadTag = zeros(Nf,1);
end


%% LOCAL FUNCTIONS


function mesh = relabel_blocked_cells(mesh, geom, marker, tol)
%RELABEL_BLOCKED_CELLS
% Update mesh.cells.physicalTag if cell center
% lies inside a nonParticipatingMedium volume.
% Also merges marker physical names into mesh.physicalNames.

    if nargin < 4
        tol = 1e-6;
    end

    Nc = numel(mesh.cells.id);

    % Merge marker physical names into mesh.physicalNames
    markerKeys = marker.physicalNames.keys;
    for i = 1:numel(markerKeys)
        k = markerKeys{i};
        if ~isKey(mesh.physicalNames, k)
            mesh.physicalNames(k) = marker.physicalNames(k);
        else       
            error('Keys already exist.');
        end
    end

    mesh.cells.isBlocked = false(Nc,1);
    mesh.cells.isLoad = false(Nc,1);
    mesh.cells.loadTag = zeros(Nc,1);
    mesh.cells.zoneId = zeros(Nc,1);

    keysList = mesh.physicalNames.keys;
    blockedTags = []; 
    for i = 1:numel(keysList)
        k = keysList{i};
        name = mesh.physicalNames(k);
    
        if startsWith(name, 'nonParticipatingMedium')
            blockedTags(end+1) = k;
        end
    end

    tags = [marker.volumes.physicalTag];

    isBlockedVolume = ismember(tags, blockedTags);
    blockedVolumes = marker.volumes(isBlockedVolume);

    for c = 1:Nc
        xc = geom.cellCenter(c,:);

        for k = 1:numel(blockedVolumes)
            bb = blockedVolumes(k).bbox;

            if pointInBBox(xc, bb, tol)
                mesh.cells.physicalTag(c) = blockedVolumes(k).physicalTag;
                mesh.cells.isBlocked(c)   = true;
                break
            end
        end
    end
end


function inside = pointInBBox(x, bb, tol)
    inside = (x(1) >= bb(1)-tol) && (x(1) <= bb(4)+tol) && ...
             (x(2) >= bb(2)-tol) && (x(2) <= bb(5)+tol) && ...
             (x(3) >= bb(3)-tol) && (x(3) <= bb(6)+tol);
end


function key = faceKey(nodeIds)
% Canonical key for a face, independent of orientation
    s = sort(nodeIds(:).');
    key = sprintf('%d_', s);
end


function [A, C, n] = polygonAreaCenterNormal(pts)
% Area, centroid, unit normal

    p1 = pts(1,:);
    pa = pts(2,:) - p1;
    pb = pts(3,:) - p1;

    triN = 0.5*cross(pa, pb);
    A = 2*norm(triN);
    n = triN/norm(triN);
    C=0.5*(p1+pts(3,:));
end


function [V, C] = polyhedronVolumeCentroid(pts)
% Volume and centroid of a convex polyhedron from its vertices
% hull triangulation

    x0 = mean(pts, 1);
    K = convhulln(pts);

    V = 0.0;
    C = zeros(1,3);

    for k = 1:size(K,1)
        a = pts(K(k,1), :);
        b = pts(K(k,2), :);
        c = pts(K(k,3), :);

        vol = abs(dot(a - x0, cross(b - x0, c - x0))) / 6.0;
        ctr = (x0 + a + b + c) / 4.0;

        V = V + vol;
        C = C + vol * ctr;
    end

    if V <= 0
        error('Degenerate cell with zero volume detected.');
    end

    C = C / V;
end


function physTag = getInterfaceFaceTag(blockedPhysName, faceCenter, faceNormal, marker)

    blockedPrefix = 'nonParticipatingMedium_';
    if ~startsWith(blockedPhysName, blockedPrefix)
        error(['Blocked cell physical name "%s" does not follow the ', ...
               'expected convention "%s<BLOCK_NAME>".'], ...
               blockedPhysName, blockedPrefix);
    end
    blockName = extractAfter(blockedPhysName, blockedPrefix);

    interfacePrefix = ['interface_', blockName, '_'];

    bboxTol  = 0.089;   % m, tolerance for bbox inclusion
    axisTol  = 1e-6;   % unitless, tolerance for axis alignment (1 - |dot|)
    degenTol = 1e-9;   % m, tolerance for detecting a degenerate bbox dim.

    % Determine the axis of the face normal
    faceAxis = normalAxis(faceNormal, axisTol);
    if faceAxis == 0
        error(['Face normal (%.6g, %.6g, %.6g) is not aligned with any ', ...
               'coordinate axis (X, Y or Z). Non-axis-aligned interfaces ', ...
               'are not supported.'], ...
               faceNormal(1), faceNormal(2), faceNormal(3));
    end

    % Search among marker interface surfaces of the current block
    physTag = int32(0);
    found   = false;

    for i = 1:numel(marker.surfaces)
        sPhysTag  = marker.surfaces(i).physicalTag;
        sPhysName = marker.physicalNames(sPhysTag);

        % Keep only interface surfaces attached to the current block
        if ~startsWith(sPhysName, interfacePrefix)
            continue;
        end

        bbox = marker.surfaces(i).bbox;   % [xmin ymin zmin xmax ymax zmax]

        % Determine the normal axis of this marker surface from its bbox
        surfAxis = markerSurfaceAxis(bbox, degenTol);
        if surfAxis == 0
            error(['Marker surface "%s" (physicalTag %d) is not aligned ', ...
                   'with a coordinate plane (no single degenerate bbox ', ...
                   'dimension).'], sPhysName, sPhysTag);
        end

        % Filter out marker surfaces not parallel to the face plane
        if surfAxis ~= faceAxis
            continue;
        end

        % Bounding box inclusion test on the face centre
        inside = faceCenter(1) >= bbox(1) - bboxTol && faceCenter(1) <= bbox(4) + bboxTol && ...
                 faceCenter(2) >= bbox(2) - bboxTol && faceCenter(2) <= bbox(5) + bboxTol && ...
                 faceCenter(3) >= bbox(3) - bboxTol && faceCenter(3) <= bbox(6) + bboxTol;

        if inside
            physTag = int32(sPhysTag);
            found   = true;
            break
        end
    end

    if ~found
        error(['No marker interface surface matches face centre ', ...
               '(%.6g, %.6g, %.6g) with normal aligned with axis %d ', ...
               'for block "%s" (expected physical name prefix "%s").'], ...
               faceCenter(1), faceCenter(2), faceCenter(3), faceAxis, ...
               blockName, interfacePrefix);
    end
end


function ax = normalAxis(n, tol)
%NORMALAXIS Return 1, 2, or 3 if n is parallel to X, Y or Z (up to sign),
% else 0. n need not be unit; it is normalised internally.
    nn = n(:).' / norm(n);
    if     abs(abs(nn(1)) - 1) < tol, ax = 1;
    elseif abs(abs(nn(2)) - 1) < tol, ax = 2;
    elseif abs(abs(nn(3)) - 1) < tol, ax = 3;
    else,                              ax = 0;
    end
end


function ax = markerSurfaceAxis(bbox, tol)
%MARKERSURFACEAXIS Return 1, 2, or 3 corresponding to the degenerate
% dimension of an axis-aligned planar surface bbox, else 0.
%   bbox = [xmin ymin zmin xmax ymax zmax]
%   ax = 1 if xmin == xmax (surface in YZ plane, normal along X)
%   ax = 2 if ymin == ymax (surface in XZ plane, normal along Y)
%   ax = 3 if zmin == zmax (surface in XY plane, normal along Z)
    dx = abs(bbox(4) - bbox(1));
    dy = abs(bbox(5) - bbox(2));
    dz = abs(bbox(6) - bbox(3));
    if     dx < tol && dy >= tol && dz >= tol, ax = 1;
    elseif dy < tol && dx >= tol && dz >= tol, ax = 2;
    elseif dz < tol && dx >= tol && dy >= tol, ax = 3;
    else,                                       ax = 0;
    end
end
