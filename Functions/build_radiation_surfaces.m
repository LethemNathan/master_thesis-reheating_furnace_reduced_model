function radSurf = build_radiation_surfaces(zonesCtx, mesh, geom)
%BUILD_RADIATION_SURFACES Enumerate the radiative elements of the zonal model.
%
% Builds the list of "elements" between which the exchange areas will be
% computed. Three kinds:
%
%   * GAS elements : one per furnace zone (the volume zones of
%     zonesCtx.id). Emitting / absorbing gas volume with uniform temperature.
%
%   * WALL / OBSTACLE surfaces : within each zone, every physical-tag group
%     of zonesCtx.byFaceTag{k} (external walls, heaters / blocked obstacles,
%     skids, ...) forms ONE individual diffuse-gray surface.
%
%   * SLAB (load) surfaces : defined INDEPENDENTLY of the gas subzones.
%     Each slab is cut into 3 equidistant slices along its LENGTH (y axis ;
%     the width along x is NOT subdivided). Its faces are classified by
%     orientation into bottom / top / lateral (lateral = front + rear +
%     side), the symmetry face being ignored. This yields up to
%     3 slices x 3 orientation classes = 9 surfaces per slab.
%
% Zone-zone interfaces are internal fluid faces (transparent to radiation)
% and are NOT surface elements ; radiation crosses them through the medium
% and is captured by the gas-gas exchange areas.
%
% Element global numbering: all GAS elements first (order of zonesCtx.id),
% then WALL/OBSTACLE surfaces (per zone, byFaceTag order), then SLAB
% surfaces (per slab, orientation, slice). This ordering is fixed and
% reused by compute_direct_exchange_areas and the validation / coupling.
%
% Slab face-tag decoding (geom.loadTag): loadTag = 110000 + slabId*100 + yy,
% with yy in {1,2,4,5,6,8} = symmetry / bottom / rear / side / front / top.
%
% Inputs
%   zonesCtx : from build_zones_context ; uses .id, .byFaceTag{k}.
%   mesh     : with .cells.zoneId (fluid cells carry their zone tag).
%   geom     : with .faceArea (Nf x 1), .CellVolume (Nc x 1),
%              .faceCenter (Nf x 3), .loadTag (Nf x 1).
%
% Output struct radSurf
%   .nElem, .nGas, .nSurf
%   .elem : 1 x nElem struct array, each with
%             .idx      global element index
%             .type     'gas' | 'surf'
%             .zoneTag  owning zone tag (0 for slab surfaces: subzone-independent)
%             .tag      zoneTag (gas) | faceTag (wall) | composite (slab)
%             .source   'gas' | 'faceTag' | 'load'
%             .faces    (Mf x 1 int32) abs geom face indices  (surf ; [] gas)
%             .cells    (Mc x 1 int32) cell indices           (gas ; [] surf)
%             .area     total surface area  [m^2]  (surf ; 0 gas)
%             .volume   total gas volume    [m^3]  (gas ; 0 surf)
%             .slabId   slab id (slab surfaces ; [] otherwise)
%             .orient   'bottom'|'top'|'lateral' (slab surfaces ; '' otherwise)
%             .bin      slice index 1..3 (slab surfaces ; [] otherwise)
%   .gasElem, .surfElem : global indices of each kind
%   .gasOfZone : containers.Map(double zoneTag -> gas element index)
%   .zoneTags  : (nGas x 1) zone tags aligned with gasElem
%
% NOTE ON AXES: the slab length is assumed along y (geom.faceCenter(:,2))
% and the width along x. If your convention is swapped, change the column
% index used for the slicing (marked below).

    LENGTH_AXIS   = 2;      % y = slab length (change to 1 if x is the length)
    N_SLICE       = 3;
    SKIP_FACETAGS = 14;     % symmetry plane (radBcType 4): not a radiative surface

    faceArea   = double(geom.faceArea);
    cellVolume = double(geom.cellVolume);
    faceCenter = double(geom.faceCenter);
    loadTagVec = double(geom.loadTag);
    zoneId     = mesh.cells.zoneId;

    zoneTags = double(zonesCtx.id(:));
    Nz       = numel(zoneTags);

    tmpl = struct('idx', [], 'type', '', 'zoneTag', [], 'tag', [], ...
                  'source', '', 'faces', int32([]), 'cells', int32([]), ...
                  'area', 0, 'volume', 0, 'slabId', [], 'orient', '', 'bin', []);

    elem = repmat(tmpl, 1, 0);

    %% --- 1) GAS elements (one per zone) ---
    gasOfZone = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for k = 1:Nz
        ztag  = zoneTags(k);
        cells = int32(find(zoneId == ztag));

        s         = tmpl;
        s.idx     = numel(elem) + 1;
        s.type    = 'gas';
        s.zoneTag = ztag;
        s.tag     = ztag;
        s.source  = 'gas';
        s.cells   = cells;
        s.volume  = sum(cellVolume(cells));

        elem(end+1)     = s;             %#ok<AGROW>
        gasOfZone(ztag) = s.idx;
    end
    nGas = numel(elem);

    %% --- 2) WALL / OBSTACLE surfaces (byFaceTag, per zone) ---
    for k = 1:Nz
        ztag = zoneTags(k);
        bt   = zonesCtx.byFaceTag{k};
        if isempty(bt) || ~isfield(bt, 'faceTag')
            continue
        end
        for j = 1:numel(bt.faceTag)
            if ismember(double(bt.faceTag(j)), SKIP_FACETAGS)
                continue    % skip symmetry planes (not radiative surfaces)
            end
            faces = int32(abs(double(bt.faceIdx{j}(:))));
            if isempty(faces), continue; end
            s         = tmpl;
            s.idx     = numel(elem) + 1;
            s.type    = 'surf';
            s.zoneTag = ztag;
            s.tag     = double(bt.faceTag(j));
            s.source  = 'faceTag';
            s.faces   = faces;
            s.area    = sum(faceArea(faces));
            elem(end+1) = s;             %#ok<AGROW>
        end
    end

    %% --- 3) SLAB surfaces (subzone-independent, 9 per slab) ---
    loadFaces = find(loadTagVec ~= 0);
    base      = loadTagVec(loadFaces) - 110000;
    slabIdF   = floor(base / 100);
    yyF       = base - slabIdF * 100;

    % Keep only bottom(2), top(8), lateral(4,5,6) ; drop symmetry(1) & other.
    keep      = ismember(yyF, [2 8 4 5 6]);
    loadFaces = loadFaces(keep);
    slabIdF   = slabIdF(keep);
    yyF       = yyF(keep);

    orient          = strings(numel(yyF), 1);
    orient(yyF==2)  = "bottom";
    orient(yyF==8)  = "top";
    orient(yyF==4)  = "rear";
    orient(yyF==5)  = "side";
    orient(yyF==6)  = "front";

    yFace   = faceCenter(loadFaces, LENGTH_AXIS);
    classes = ["bottom", "top", "rear", "side", "front"];
    slabIds = unique(slabIdF);

    for si = 1:numel(slabIds)
        sid    = slabIds(si);
        inSlab = (slabIdF == sid);
        ys     = yFace(inSlab);
        ymin   = min(ys);
        ymax   = max(ys);
        if ymax <= ymin
            edges = [ymin, ymin, ymin, ymin + eps];   % degenerate: all in bin 1
        else
            edges = linspace(ymin, ymax, N_SLICE + 1);
        end

        for ci = 1:numel(classes)
            cls = classes(ci);
            for b = 1:N_SLICE
                if b < N_SLICE
                    inBin = (yFace >= edges(b)) & (yFace < edges(b+1));
                else
                    inBin = (yFace >= edges(b)) & (yFace <= edges(b+1));
                end
                sel   = inSlab & (orient == cls) & inBin;
                faces = int32(loadFaces(sel));
                if isempty(faces), continue; end

                s         = tmpl;
                s.idx     = numel(elem) + 1;
                s.type    = 'surf';
                s.zoneTag = 0;                       % independent of subzones
                s.tag     = sid*1000 + ci*10 + b;    % composite slab-surface tag
                s.source  = 'load';
                s.faces   = faces;
                s.area    = sum(faceArea(faces));
                s.slabId  = sid;
                s.orient  = char(cls);
                s.bin     = b;
                elem(end+1) = s;         %#ok<AGROW>
            end
        end
    end

    nElem = numel(elem);
    nSurf = nElem - nGas;

    %% --- Pack ---
    radSurf           = struct();
    radSurf.nElem     = nElem;
    radSurf.nGas      = nGas;
    radSurf.nSurf     = nSurf;
    radSurf.elem      = elem;
    radSurf.gasElem   = (1:nGas).';
    radSurf.surfElem  = (nGas+1:nElem).';
    radSurf.gasOfZone = gasOfZone;
    radSurf.zoneTags  = zoneTags;

    nSlab = numel(slabIds);
    fprintf(['build_radiation_surfaces: %d elements (%d gas zones, ', ...
             '%d surfaces ; %d slabs -> up to %d slab surfaces).\n'], ...
            nElem, nGas, nSurf, nSlab, 15*nSlab);
end
