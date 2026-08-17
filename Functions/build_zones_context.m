

function zonesCtx = build_zones_context(zonesDefinition, mesh, geom)
%BUILD_ZONES_CONTEXT Build per-zone context for the FVM energy equation.
%
% In addition to the bbox-based topology (zoneCenter, zoneNeighbor, zoneFace),
% this function classifies every face of geom according to the zone(s) it
% borders, with a sign convention indicating the flux direction relative
% to the zone (geom.faceNormal is oriented owner -> neighbour) :
%
%   * -f  if the face normal points OUT of the zone (the zone is on the
%         OWNER side of f).
%   * +f  if the face normal points INTO the zone (the zone is on the
%         NEIGHBOUR side of f).
%
% Three new fields are populated on top of the existing topology:
%
%   zonesCtx.zoneFaces : Nz x 1 cell array. zoneFaces{k} is the int32
%       column vector of signed face indices for the zone at position k
%       (i.e. with tag zonesCtx.id(k)).
%
%   zonesCtx.byFaceTag : Nz x 1 cell array. byFaceTag{k} is a struct
%       with two parallel fields
%         .faceTag (1 x K1_k int32, unique values of geom.faceTag
%                   bordering zone at position k)
%         .faceIdx (1 x K1_k cell of int32 vectors of signed face
%                   indices grouped by faceTag, for that zone only).
%       Populated by:
%         - external boundary faces (case A) bordering the zone,
%         - fluid|blocked obstacle interfaces (case B2 with isBlocked).
%
%   zonesCtx.byLoadTag : Nz x 1 cell array. byLoadTag{k} is a struct
%       with two parallel fields
%         .loadTag (1 x K2_k int32, unique values of geom.loadTag(f)
%                   for slab-fluid interface faces bordering zone k,
%                   i.e. subzone-specific tags from
%                   slabs(i).surfaces.faceTag rather than the bare
%                   slab id ; format = 110000 + slabId*100 + yy with
%                   yy in {1,2,4,5,6,8} for symmetry / bottom / rear /
%                   side / front / top).
%         .faceIdx (1 x K2_k cell of int32 vectors).
%       Populated by fluid|load obstacle interfaces (case B2 with
%       isLoad), which have geom.faceTag == 0 but geom.loadTag != 0.
%
% Inputs
%   zonesDefinition : struct from the user, with .numZones and
%                     .volumes(j).bbox / .volumes(j).zoneTag.
%   mesh            : with .cells.zoneId, .cells.isBlocked, .cells.isLoad,
%                     .cells.loadTag.
%   geom            : with .owner, .neighbour, .faceTag.

    Nz = zonesDefinition.numZones;

    %% --- Existing bbox-based topology ---
    zoneCenter = zeros(Nz, 3);
    for j = 1:Nz
        zoneCenter(j, :) = bboxCenter(zonesDefinition.volumes(j).bbox);
    end

    bboxes              = vertcat(zonesDefinition.volumes.bbox);
    zoneNeighbor        = findBBoxNeighborsMatrix(bboxes);   % default tol
    [zoneFace, faceAxis] = buildFaceCellsFromNeighbors(zoneNeighbor);

    %% --- Per-zoneFace face center (Nf x 3) ---
    %  The face center is derived from the cells' bboxes plus the axis info
    %  recorded by buildFaceCellsFromNeighbors. Convention reminder:
    %  column 1 -> column 2 follows -x (axis 1), -y (axis 2), +z (axis 3).
    Nf_z = size(zoneFace, 1);
    zoneFaceCenter = zeros(Nf_z, 3);
    for f = 1:Nf_z
        zoneFaceCenter(f, :) = computeBboxFaceCenter( ...
            zoneFace(f, 1), zoneFace(f, 2), faceAxis(f), bboxes);
    end

    %% --- zoneFaceByDir (Nz x 6) ---
    %  For each zone, store the row index in zoneFace of the face on each
    %  of the 6 directions (E/W/N/S/U/D), or 0 if no face. Direction is
    %  identified from (zoneFaceCenter - zoneCenter): the dominant
    %  component in absolute value gives the axis, its sign gives +/-.


    cdir  = loadTopologyCodes();
    DIR_E = cdir.DIR_E;   % +y
    DIR_W = cdir.DIR_W;   % -y
    DIR_N = cdir.DIR_N;   % +z
    DIR_S = cdir.DIR_S;   % -z
    DIR_U = cdir.DIR_U;   % +x
    DIR_D = cdir.DIR_D;   % -x

    zoneFaceByDir = zeros(Nz, 6, 'int32');
    for f = 1:Nf_z
        twoCells = zoneFace(f, :);
        for kk = 1:2
            cellIdx = twoCells(kk);
            if cellIdx == 0
                continue
            end
            delta = zoneFaceCenter(f, :) - zoneCenter(cellIdx, :);
            d     = dirFromDelta(delta, DIR_E, DIR_W, DIR_N, DIR_S, DIR_U, DIR_D);
            zoneFaceByDir(cellIdx, d) = int32(f);
        end
    end

    %% --- zoneTag -> position 1..Nz (vector lookup, fast) ---
    zoneTagsList = int32([zonesDefinition.volumes.zoneTag]');
    if any(zoneTagsList <= 0)
        error(['build_zones_context: zone tags must be strictly positive ', ...
               'integers (got %d).'], min(zoneTagsList));
    end
    maxTag = double(max(zoneTagsList));
    tagToPos = zeros(maxTag, 1, 'int32');
    tagToPos(zoneTagsList) = int32(1:Nz);

    %% --- faceBC (Nz x 6) ---
    %  Default: 1 (radiative) on every boundary direction (no zone neighbour).
    %  Override to 2 (inlet) on the lateral / roof burner zones.
    faceBC = zeros(Nz, 6, 'int32');
    faceBC(zoneNeighbor == 0) = 1;

    if isfield(zonesDefinition, 'lateralBurnerLoc')
        faceBC = override_burner_bc(faceBC, ...
            zonesDefinition.lateralBurnerLoc, [DIR_E, DIR_W], 'lateral', ...
            tagToPos, maxTag);
    end
    if isfield(zonesDefinition, 'roofBurnerLoc')
        faceBC = override_burner_bc(faceBC, ...
            zonesDefinition.roofBurnerLoc, [DIR_N, DIR_S], 'roof', ...
            tagToPos, maxTag);
    end
    if isfield(zonesDefinition, 'inletLoc')
        faceBC = override_burner_bc(faceBC, ...
            zonesDefinition.inletLoc, [DIR_U, DIR_D], 'roof', ...
            tagToPos, maxTag);
    end

    

    %% --- Walk over all faces and classify ---
    Nf          = numel(geom.owner);
    owner       = geom.owner;
    neigh       = geom.neighbour;
    faceTag     = geom.faceTag;
    faceLoadTag = geom.loadTag;             % per-face subzone tag (from slabs.surfaces.faceTag)
    zoneId      = mesh.cells.zoneId;
    isBlock     = mesh.cells.isBlocked;
    isLoad      = mesh.cells.isLoad;

    % Pre-allocated flat accumulators (over-allocated, trimmed at the end).
    zoneAccPos  = zeros(2*Nf, 1, 'int32');   % zone position 1..Nz
    zoneAccSign = zeros(2*Nf, 1, 'int32');   % signed face idx
    nZA = 0;

    ftAccKey  = zeros(Nf, 1, 'int32');
    ftAccSign = zeros(Nf, 1, 'int32');
    ftAccZone = zeros(Nf, 1, 'int32');   % zone position 1..Nz
    nFT = 0;

    ltAccKey  = zeros(Nf, 1, 'int32');
    ltAccSign = zeros(Nf, 1, 'int32');
    ltAccZone = zeros(Nf, 1, 'int32');
    nLT = 0;

    for f = 1:Nf
        o  = owner(f);
        n  = neigh(f);
        zO = int32(zoneId(o));

        if n == 0
            % --- Case A : external boundary face ---
            if zO ~= 0
                signedF = -int32(f);
                zPos    = tagToPos(zO);

                nZA = nZA + 1;
                zoneAccPos(nZA)  = zPos;
                zoneAccSign(nZA) = signedF;

                nFT = nFT + 1;
                ftAccKey(nFT)  = int32(faceTag(f));
                ftAccSign(nFT) = signedF;
                ftAccZone(nFT) = zPos;
            end
            continue
        end

        % --- Case B : internal face ---
        zN = int32(zoneId(n));

        if zO == zN
            continue
        end

        if zO ~= 0 && zN ~= 0
            % B1 : zone | zone interface (no entry in by*Tag)
            sFowner = -int32(f);
            sFneigh = +int32(f);

            nZA = nZA + 1;
            zoneAccPos(nZA)  = tagToPos(zO);
            zoneAccSign(nZA) = sFowner;

            nZA = nZA + 1;
            zoneAccPos(nZA)  = tagToPos(zN);
            zoneAccSign(nZA) = sFneigh;
        else
            % B2 : zone | obstacle (one of zO, zN is 0)
            if zO ~= 0
                zoneTagB2    = zO;
                obstacleCell = n;
                signedF      = -int32(f);    % zone on owner side
            else
                zoneTagB2    = zN;
                obstacleCell = o;
                signedF      = +int32(f);    % zone on neighbour side
            end

            zPos = tagToPos(zoneTagB2);

            nZA = nZA + 1;
            zoneAccPos(nZA)  = zPos;
            zoneAccSign(nZA) = signedF;

            if isBlock(obstacleCell)
                nFT = nFT + 1;
                ftAccKey(nFT)  = int32(faceTag(f));
                ftAccSign(nFT) = signedF;
                ftAccZone(nFT) = zPos;
            elseif isLoad(obstacleCell)
                nLT = nLT + 1;
                ltAccKey(nLT)  = int32(faceLoadTag(f));   % subzone-specific tag
                ltAccSign(nLT) = signedF;
                ltAccZone(nLT) = zPos;
            else
                error(['build_zones_context: face %d has zoneId == 0 on ', ...
                       'one side, but the corresponding cell %d is neither ', ...
                       'isBlocked nor isLoad (inconsistent topology).'], ...
                      f, obstacleCell);
            end
        end
    end

    % Trim
    zoneAccPos  = zoneAccPos(1:nZA);
    zoneAccSign = zoneAccSign(1:nZA);
    ftAccKey    = ftAccKey(1:nFT);
    ftAccSign   = ftAccSign(1:nFT);
    ftAccZone   = ftAccZone(1:nFT);
    ltAccKey    = ltAccKey(1:nLT);
    ltAccSign   = ltAccSign(1:nLT);
    ltAccZone   = ltAccZone(1:nLT);

    %% --- Group into per-zone cell array ---
    zoneFaces = cell(Nz, 1);
    for k = 1:Nz
        zoneFaces{k} = zoneAccSign(zoneAccPos == k);
    end

    %% --- Per-zone byFaceTag (Nz x 1 cell of parallel-arrays struct) ---
    byFaceTag = cell(Nz, 1);
    for k = 1:Nz
        mask         = (ftAccZone == k);
        byFaceTag{k} = group_by_key(ftAccKey(mask), ftAccSign(mask), 'faceTag');
    end

    %% --- Per-zone byLoadTag (Nz x 1 cell of parallel-arrays struct) ---
    byLoadTag = cell(Nz, 1);
    for k = 1:Nz
        mask         = (ltAccZone == k);
        byLoadTag{k} = group_by_key(ltAccKey(mask), ltAccSign(mask), 'loadTag');
    end

    %% --- Output ---
    zonesCtx = struct();
    zonesCtx.id               = [zonesDefinition.volumes.zoneTag]';
    zonesCtx.zoneCenter       = zoneCenter;
    zonesCtx.zoneNeighbor     = zoneNeighbor;
    zonesCtx.zoneFace         = zoneFace;
    zonesCtx.zoneFaceCenter   = zoneFaceCenter;
    zonesCtx.zoneFaceByDir    = zoneFaceByDir;
    zonesCtx.faceBC           = faceBC;
    zonesCtx.globaToLocalFace = 0;
    zonesCtx.zoneFaces        = zoneFaces;
    zonesCtx.byFaceTag        = byFaceTag;
    zonesCtx.byLoadTag        = byLoadTag;
end


%% LOCAL FUNCTIONS

function out = group_by_key(keyVec, signVec, keyFieldName)
%GROUP_BY_KEY Build a struct with parallel fields (keyFieldName, faceIdx)
% from flat (key, signed face) accumulators.
    if isempty(keyVec)
        out = struct(keyFieldName, int32([]), 'faceIdx', {{}});
        return
    end
    [uK, ~, bin] = unique(keyVec);
    out = struct();
    out.(keyFieldName) = uK(:).';                       % 1 x K int32
    out.faceIdx        = cell(1, numel(uK));            % 1 x K cell
    for i = 1:numel(uK)
        out.faceIdx{i} = signVec(bin == i);
    end
end


function center = bboxCenter(bbox)
% bboxCenter computes the center of a rectangular parallelepiped bbox.
    if numel(bbox) ~= 6
        error('bbox must contain 6 values: [xmin ymin zmin xmax ymax zmax].');
    end
    center = [(bbox(1) + bbox(4))/2, ...
              (bbox(2) + bbox(5))/2, ...
              (bbox(3) + bbox(6))/2];
end


function cellNeighbor = findBBoxNeighborsMatrix(bboxes, tol)
% findBBoxNeighborsMatrix finds the 6-direction neighbor of each bbox.
    if nargin < 2
        tol = 1e-9;
    end

    c = loadTopologyCodes();
    DIR_E = c.DIR_E;   % +y
    DIR_W = c.DIR_W;   % -y
    DIR_N = c.DIR_N;   % +z
    DIR_S = c.DIR_S;   % -z
    DIR_U = c.DIR_U;   % +x
    DIR_D = c.DIR_D;   % -x

    Nc = size(bboxes,1);
    cellNeighbor = zeros(Nc,6);

    for i = 1:Nc
        bb_i = bboxes(i,:);
        for j = 1:Nc
            if i == j
                continue
            end
            bb_j = bboxes(j,:);

            if abs(bb_i(5) - bb_j(2)) < tol && overlapXZ(bb_i, bb_j, tol)
                cellNeighbor(i,DIR_E) = j;
            elseif abs(bb_i(2) - bb_j(5)) < tol && overlapXZ(bb_i, bb_j, tol)
                cellNeighbor(i,DIR_W) = j;
            elseif abs(bb_i(6) - bb_j(3)) < tol && overlapXY(bb_i, bb_j, tol)
                cellNeighbor(i,DIR_N) = j;
            elseif abs(bb_i(3) - bb_j(6)) < tol && overlapXY(bb_i, bb_j, tol)
                cellNeighbor(i,DIR_S) = j;
            elseif abs(bb_i(4) - bb_j(1)) < tol && overlapYZ(bb_i, bb_j, tol)
                cellNeighbor(i,DIR_U) = j;
            elseif abs(bb_i(1) - bb_j(4)) < tol && overlapYZ(bb_i, bb_j, tol)
                cellNeighbor(i,DIR_D) = j;
            end
        end
    end
end


function tf = overlapXY(a,b,tol)
    tf = intervalsOverlap(a(1),a(4),b(1),b(4),tol) && ...
         intervalsOverlap(a(2),a(5),b(2),b(5),tol);
end

function tf = overlapXZ(a,b,tol)
    tf = intervalsOverlap(a(1),a(4),b(1),b(4),tol) && ...
         intervalsOverlap(a(3),a(6),b(3),b(6),tol);
end

function tf = overlapYZ(a,b,tol)
    tf = intervalsOverlap(a(2),a(5),b(2),b(5),tol) && ...
         intervalsOverlap(a(3),a(6),b(3),b(6),tol);
end

function tf = intervalsOverlap(a1,a2,b1,b2,tol)
    tf = min(a2,b2) - max(a1,b1) > tol;
end


function [faceCells, faceAxis] = buildFaceCellsFromNeighbors(cellNeighbor)
% buildFaceCellsFromNeighbors creates a [Nf x 2] face-cell adjacency table.
% Also returns faceAxis [Nf x 1] with values 1, 2, or 3 for the axis
% normal to each face (1 = x, 2 = y, 3 = z).
    DIR_E = 1;   % +y
    DIR_W = 2;   % -y
    DIR_N = 3;   % +z
    DIR_S = 4;   % -z
    DIR_U = 5;   % +x
    DIR_D = 6;   % -x

    Nc = size(cellNeighbor,1);
    faceCells = [];
    faceAxis  = zeros(0, 1);

    for i = 1:Nc

        % Faces normal to y.
        % DIR_E adds every +y face from cell i (internal AND boundary).
        % DIR_W adds only boundary -y faces. This way every internal
        % y-face is added exactly once (by the cell that sees the face
        % on its +y side), and there is no need for an i < j guard.
        j = cellNeighbor(i,DIR_E);
        if j > 0
            faceCells(end+1,:) = [j i];   %#ok<AGROW>
            faceAxis(end+1,1)  = 2;       %#ok<AGROW>
        elseif j == 0
            faceCells(end+1,:) = [0 i];   %#ok<AGROW>
            faceAxis(end+1,1)  = 2;       %#ok<AGROW>
        end

        j = cellNeighbor(i,DIR_W);
        if j == 0
            faceCells(end+1,:) = [i 0];   %#ok<AGROW>
            faceAxis(end+1,1)  = 2;       %#ok<AGROW>
        end

        % Faces normal to z (same logic with DIR_N taking every +z face).
        j = cellNeighbor(i,DIR_N);
        if j > 0
            faceCells(end+1,:) = [i j];   %#ok<AGROW>
            faceAxis(end+1,1)  = 3;       %#ok<AGROW>
        elseif j == 0
            faceCells(end+1,:) = [i 0];   %#ok<AGROW>
            faceAxis(end+1,1)  = 3;       %#ok<AGROW>
        end

        j = cellNeighbor(i,DIR_S);
        if j == 0
            faceCells(end+1,:) = [0 i];   %#ok<AGROW>
            faceAxis(end+1,1)  = 3;       %#ok<AGROW>
        end

        % Faces normal to x (same logic with DIR_U taking every +x face).
        j = cellNeighbor(i,DIR_U);
        if j > 0
            faceCells(end+1,:) = [j i];   %#ok<AGROW>
            faceAxis(end+1,1)  = 1;       %#ok<AGROW>
        elseif j == 0
            faceCells(end+1,:) = [0 i];   %#ok<AGROW>
            faceAxis(end+1,1)  = 1;       %#ok<AGROW>
        end

        j = cellNeighbor(i,DIR_D);
        if j == 0
            faceCells(end+1,:) = [i 0];   %#ok<AGROW>
            faceAxis(end+1,1)  = 1;       %#ok<AGROW>
        end
    end
end


function center = computeBboxFaceCenter(c1, c2, ax, bboxes)
%COMPUTEBBOXFACECENTER Center of the bbox-face f = [c1 c2] normal to axis ax.
%
% Convention recap (as in buildFaceCellsFromNeighbors):
%   column 1 -> column 2 follows -x (axis 1), -y (axis 2), +z (axis 3).
% Equivalently, column 1 cell is on +x / +y / -z side of column 2 cell.

    % --- Coordinate of the face plane along the normal axis ---
    if c1 > 0 && c2 > 0
        bb1 = bboxes(c1, :);
        if ax == 1
            cAx = bb1(1);            % col1 is +x of col2 -> shared coord = bb1.xmin
        elseif ax == 2
            cAx = bb1(2);            % col1 is +y of col2 -> shared coord = bb1.ymin
        else
            cAx = bb1(6);            % col1 is -z of col2 -> shared coord = bb1.zmax
        end
        bbA = bboxes(c1, :);
        bbB = bboxes(c2, :);
    else
        if c1 > 0
            % column 2 is missing : missing cell would be on -x/-y/+z of c1
            bb = bboxes(c1, :);
            isCol1 = true;
        else
            % column 1 is missing : missing cell would be on +x/+y/-z of c2
            bb = bboxes(c2, :);
            isCol1 = false;
        end
        if ax == 1
            if isCol1, cAx = bb(1); else, cAx = bb(4); end   % -x or +x face
        elseif ax == 2
            if isCol1, cAx = bb(2); else, cAx = bb(5); end   % -y or +y face
        else
            if isCol1, cAx = bb(6); else, cAx = bb(3); end   % +z or -z face
        end
        bbA = bb;
        bbB = bb;
    end

    % --- Center on the two transverse axes : midpoint of overlap (= bbox
    %     midpoint when only one cell exists, since bbA == bbB). ---
    center      = zeros(1, 3);
    center(ax)  = cAx;
    otherAxes   = [1 2 3];
    otherAxes(ax) = [];
    for k = otherAxes
        lo = max(bbA(k),   bbB(k));
        hi = min(bbA(k+3), bbB(k+3));
        center(k) = 0.5 * (lo + hi);
    end
end


function d = dirFromDelta(delta, DIR_E, DIR_W, DIR_N, DIR_S, DIR_U, DIR_D)
%DIRFROMDELTA Pick the topology direction (1..6) matching the dominant
% component of (faceCenter - cellCenter).
    [~, ax] = max(abs(delta));
    switch ax
        case 1   % x
            if delta(1) > 0, d = DIR_U; else, d = DIR_D; end
        case 2   % y
            if delta(2) > 0, d = DIR_E; else, d = DIR_W; end
        case 3   % z
            if delta(3) > 0, d = DIR_N; else, d = DIR_S; end
    end
end


function faceBC = override_burner_bc(faceBC, zoneTagList, dirs, label, ...
                                     tagToPos, maxTag)
%OVERRIDE_BURNER_BC Set faceBC to 2 (inlet) on the requested directions
% for every zone whose tag is in zoneTagList. Only entries currently at 1
% (radiative boundary) are overwritten; if none of the requested
% directions on a zone is a boundary, a warning is emitted.
    for kk = 1:numel(zoneTagList)
        zt = int32(zoneTagList(kk));
        if zt <= 0 || zt > maxTag || tagToPos(zt) == 0
            error(['build_zones_context: %s burner zone tag %d is not ', ...
                   'a valid zone.'], label, zt);
        end
        i = double(tagToPos(zt));
        anyApplied = false;
        for dd = dirs
            if faceBC(i, dd) == 1
                faceBC(i, dd) = 2;
                anyApplied = true;
            end
        end
        if ~anyApplied
            warning(['build_zones_context: %s burner zone %d (pos %d) ', ...
                     'has no boundary face along the requested ', ...
                     'directions; faceBC unchanged.'], label, zt, i);
        end
    end
end
