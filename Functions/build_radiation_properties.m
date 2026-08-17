function radProps = build_radiation_properties(geom, mesh, fields, ...
                                                Twall_zone_table, ...
                                                kappaA_value, physicalConst, ...
                                                radBcType2Tags, radBcType3Tags)
%BUILD_RADIATION_PROPERTIES  (optimised, drop-in)
%
% Same output as the original, but the STATIC part (wall classification,
% radBcType, and the per-wall-face temperature-source mapping) is purely
% geometric and is computed ONCE and cached (persistent). Each call then only
% updates the temperatures (Ib_cell, Ib_wall) from fields.T and
% Twall_zone_table, fully vectorised.
%
% ASSUMPTION: the KEYS of Twall_zone_table are stable across the coupling loop
% (only the values change). Call  clear build_radiation_properties  to force a
% rebuild if the mesh, the wall tags, or the set of table keys change.

    persistent C

    Nc    = numel(mesh.cells.id);
    Nf    = numel(geom.faceArea);
    sigma = physicalConst.sigma;

    if nargin < 7 || isempty(radBcType2Tags), radBcType2Tags = []; end
    if nargin < 8 || isempty(radBcType3Tags), radBcType3Tags = []; end

    %% ---- Build / refresh the static cache if the geometry changed ----
    sig = [Nf, Nc, double(sum(geom.faceArea)), ...
           numel(radBcType2Tags), numel(radBcType3Tags), ...
           sum(double(radBcType2Tags(:))), sum(double(radBcType3Tags(:)))];
    if isempty(C) || ~isequal(C.sig, sig)
        C = build_static_cache(geom, mesh, Twall_zone_table, ...
                               radBcType2Tags, radBcType3Tags, sig);
    end

    %% ---- Per-call temperature update (vectorised) ----
    Ib_cell = (sigma/pi) * fields.T(:).^4;

    Tw_face = zeros(Nf, 1);
    Tw_face(C.isWall) = 1300;                       % default on every wall

    % (a) slab (load) faces -> surface temperature T_s (from the slab-conduction
    %     surface node, fields.T_surf) when available, else the adjacent
    %     load-cell temperature (fallback: NaN entry or no T_surf field).
    if any(C.slabMask)
        sm      = C.slabMask;
        Tw_slab = fields.T(C.loadCellOfFace(sm));        % fallback (cell temp)
        if isfield(fields, 'T_surf') && ~isempty(fields.T_surf)
            ts  = fields.T_surf(sm);
            use = ~isnan(ts);
            Tw_slab(use) = ts(use);
        end
        Tw_face(sm) = Tw_slab;
    end

    % (b) type-1 refractory faces -> Twall_zone_table lookup via precomputed slots
    if C.nSlot > 0
        slotVals = nan(C.nSlot, 1);
        for s = 1:C.nSlot
            k = C.slotKey{s};
            if isKey(Twall_zone_table, k)
                slotVals(s) = Twall_zone_table(k);
            end
        end
        rf   = C.refracIdx;
        slot = C.slotOfRefrac;                      % >=1 resolved, 0 = missing
        Tw_r = 1300 * ones(numel(rf), 1);
        ok   = slot > 0;
        vals = slotVals(max(slot,1));
        use  = ok & ~isnan(vals);
        Tw_r(use) = vals(use);
        Tw_face(rf) = Tw_r;
    end

    % (c) faceTag-15 faces refined from the (static) Fluent box mean
    if any(C.fluentMask)
        Tw_face(C.fluentMask) = C.TwFluent(C.fluentMask);
    end

    Ib_wall = (sigma/pi) * Tw_face.^4;
    Ib_wall(~C.isWall) = 0;

    %% ---- Pack output ----
    radProps           = struct();

    radProps.kappa_a   = kappaA_value * ones(Nc, 1);

    % unitDigit = mod(mesh.cells.zoneId, 10);
    % radProps.kappa_a(unitDigit == 1) = 0.2209;
    % radProps.kappa_a(unitDigit == 2) = 0.2055;

    radProps.sigma_s   = 0;
    radProps.eps_wall  = C.eps_wall;
    radProps.isWall    = C.isWall;
    radProps.radBcType = C.radBcType;
    radProps.Tw_face   = Tw_face;
    radProps.Ib_cell   = Ib_cell;
    radProps.Ib_wall   = Ib_wall;
end


%% ======================= STATIC CACHE (built once) =======================
function C = build_static_cache(geom, mesh, Twall_zone_table, ...
                                radBcType2Tags, radBcType3Tags, sig)

    Nf = numel(geom.faceArea);

    o       = double(geom.owner(:));
    n       = double(geom.neighbour(:));
    faceTag = double(geom.faceTag(:));
    cat1    = geom.faceCategory(:) == 1;
    isLoadF = logical(geom.isLoadFace(:));

    isLoadCell    = logical(mesh.cells.isLoad(:));
    isBlockedCell = logical(mesh.cells.isBlocked(:));
    zoneId        = double(mesh.cells.zoneId(:));

    %% --- Vectorised wall classification (verified vs the original loop) ---
    loadO  = isLoadCell(o);   blockO = isBlockedCell(o);
    hasN   = n > 0;
    loadN  = false(Nf,1);     blockN = false(Nf,1);
    loadN(hasN)  = isLoadCell(n(hasN));
    blockN(hasN) = isBlockedCell(n(hasN));

    bnd   = ~hasN;
    condB = (~isLoadF) & cat1 & bnd  & (~loadO) & (~blockO);
    condO = (~isLoadF) & cat1 & hasN & blockO & (~blockN) & (~loadN);
    condN = (~isLoadF) & cat1 & hasN & blockN & (~blockO) & (~loadO);
    isWall = isLoadF | condB | condO | condN;

    radBcType = zeros(Nf, 1, 'int32');
    radBcType(isWall) = 1;
    if ~isempty(radBcType3Tags)
        radBcType(isWall & ismember(faceTag, radBcType3Tags)) = 3;
    end
    if ~isempty(radBcType2Tags)   % type-2 precedence over type-3 (as in elseif)
        radBcType(isWall & ismember(faceTag, radBcType2Tags)) = 2;
    end
    radBcType(isWall & (faceTag == 14)) = 4;   % type-4 highest precedence

    eps_wall = zeros(Nf, 1);
    eps_wall(isWall) = 1;

    %% --- Slab (load) faces: adjacent load-cell index ---
    slabMask = isWall & isLoadF & (radBcType == 1);
    loadCellOfFace = zeros(Nf, 1);
    if any(slabMask)
        oo = o(slabMask);  nn = n(slabMask);
        useO = isLoadCell(oo);
        lc = nn;  lc(useO) = oo(useO);
        loadCellOfFace(slabMask) = lc;
    end

    %% --- Type-1 refractory faces: fluid cell, zone tag, table slot ---
    refracMask = isWall & (~isLoadF) & (radBcType == 1);
    refracIdx  = find(refracMask);

    oo = o(refracIdx);  nn = n(refracIdx);
    oFluid = ~isLoadCell(oo) & ~isBlockedCell(oo);
    hN     = nn > 0;
    nFluid = false(size(nn));
    nFluid(hN) = ~isLoadCell(nn(hN)) & ~isBlockedCell(nn(hN));
    fc = zeros(size(oo));
    fc(oFluid) = oo(oFluid);
    useN = (~oFluid) & nFluid;
    fc(useN) = nn(useN);

    zoneTagR = zeros(size(fc));
    good = fc > 0;
    zoneTagR(good) = zoneId(fc(good));
    faceTagR = faceTag(refracIdx);

    % Resolve each refractory face to a table slot ONCE (keys are stable).
    slotOfRefrac = zeros(numel(refracIdx), 1);
    slotKey  = {};
    keyCanon = containers.Map('KeyType','char','ValueType','double');
    missing  = zeros(0, 2);
    for r = 1:numel(refracIdx)
        zt = zoneTagR(r);  ft = faceTagR(r);
        if fc(r) == 0 || zt == 0
            continue                     % no fluid cell / zone -> stays 1300
        end
        pairStr = sprintf('%d_%d', zt, ft);
        if isKey(Twall_zone_table, pairStr)
            keyObj = pairStr;  canon = ['P:' pairStr];
        elseif isKey(Twall_zone_table, zt)
            keyObj = zt;       canon = sprintf('Z:%d', zt);
        else
            if ft ~= 0 && ft ~= 14, missing(end+1,:) = [zt ft]; end %#ok<AGROW>
            continue                     % missing -> stays 1300
        end
        if isKey(keyCanon, canon)
            slotOfRefrac(r) = keyCanon(canon);
        else
            slotKey{end+1,1} = keyObj;                 %#ok<AGROW>
            keyCanon(canon)  = numel(slotKey);
            slotOfRefrac(r)  = numel(slotKey);
        end
    end
    nSlot = numel(slotKey);

    if ~isempty(missing)
        missing = unique(missing, 'rows');
        warning(['build_radiation_properties: missing Twall values for ', ...
                 'some (zoneTag, faceTag) pairs. Missing pairs: %s'], ...
                 mat2str(missing));
    end

    %% --- Fluent faceTag-15 refinement (static box means) ---
    fluentMask = false(Nf, 1);
    TwFluent   = zeros(Nf, 1);
    fluentFile = 'input/inputFluidFurnaceTempData';
    if isfile(fluentFile)
        fd = readmatrix(fluentFile, 'FileType', 'text');
        fY = fd(:,3);  fZ = fd(:,4);  fT = fd(:,5);
        dy = 0.0901;  dz = 0.0673;  tol = 1e-8;
        idx15 = refracIdx(faceTagR == 15);
        for q = 1:numel(idx15)
            f   = idx15(q);
            fcC = geom.faceCenter(f,:);
            in  = fY >= fcC(2)-dy/2-tol & fY <= fcC(2)+dy/2+tol & ...
                  fZ >= fcC(3)-dz/2-tol & fZ <= fcC(3)+dz/2+tol;
            if any(in)
                TwFluent(f)   = mean(fT(in));
                fluentMask(f) = true;
            end
        end
    end

    C = struct('sig', sig, 'isWall', isWall, 'radBcType', radBcType, ...
        'eps_wall', eps_wall, 'slabMask', slabMask, ...
        'loadCellOfFace', loadCellOfFace, 'refracIdx', refracIdx, ...
        'slotOfRefrac', slotOfRefrac, 'slotKey', {slotKey}, 'nSlot', nSlot, ...
        'fluentMask', fluentMask, 'TwFluent', TwFluent);
end
