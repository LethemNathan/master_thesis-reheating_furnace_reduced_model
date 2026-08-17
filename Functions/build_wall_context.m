function wallCtx = build_wall_context(wallDataFurnace, zonesCtx, geom, ...
                                      matProps, targetTags)
%BUILD_WALL_CONTEXT Build per (zone, wall-tag) 1D multilayer conduction context.
%
% For every furnace zone that borders one of the requested wall tags, this
% function precomputes the geometry- and material-only data needed to solve
% a 1D steady (later unsteady) multilayer conduction problem through the
% wall thickness. Everything temperature-dependent (k(T), sigma*T^4, ...) is
% handled later by solve_wall_conduction ; NOTHING here depends on T.
%
% Granularity: one 1D column per (zoneTag, faceTag) couple. All the mesh
% faces carrying that faceTag inside that zone are lumped into a single
% column ; the incident radiation and the imposed fluxes are area-averaged
% over them (done in the solver/driver) and the resulting interior-surface
% temperature is applied back to all of them. This matches the (zone, tag)
% granularity of Twall_zone_table.
%
% Axis convention: x = 0 at the INTERIOR (furnace) surface, x = Ltot at the
% AMBIENT surface. layers(1) is the interior layer, layers(end) the ambient
% (outermost) one.
%
% FVM layout: each layer l is discretised into
%       nCell(l) = max(1, round(nodesPerMeter(l) * thickness(l)))
% uniform cells of size dx(l) = thickness(l) / nCell(l). Cells are numbered
% 1..N from the interior surface outward. The face conductance used later by
% the solver is the exact series resistance between two adjacent cell
% centres, R = (dx_i/2)/k_i + (dx_{i+1}/2)/k_{i+1}, which on a uniform mesh
% reduces to harmonic_mean(k_i, k_{i+1})/dx and handles layer interfaces
% (different material AND different dx) correctly.
%
% Inputs
%   wallDataFurnace : struct array (see wallDataFurnace_Generation.m) with
%                     .tag, .wallName, .nLayers, .hAmb, .Tinf, .TextRad and
%                     .layers(l).materialName / .thickness / .nodesPerMeter.
%   zonesCtx        : from build_zones_context ; uses .id (Nz x 1 zone tags)
%                     and .byFaceTag{k}.faceTag / .faceIdx.
%   geom            : uses .faceArea (Nf x 1).
%   matProps        : uses .extWall.<sanitizedMaterialName>.k / .cp / .rho.
%   targetTags      : vector of wall face tags to model
%                     (e.g. int32([18 16 19 865 866])).
%
% Output
%   wallCtx : 1 x Ncol struct array, one element per (zone, tag) column:
%     .zoneTag    scalar          furnace zone tag
%     .faceTag    scalar          wall face tag
%     .wallName   char            human-readable name
%     .meshFaces  (Mf x 1 int32)  absolute geom face indices lumped here
%     .area       scalar          total wall area = sum(faceArea)
%     .nLayers    scalar
%     .layerName  (1 x nL cell)   sanitized material names, interior->ambient
%     .N          scalar          total number of FVM cells across thickness
%     .dx         (N x 1)         cell thickness
%     .xc         (N x 1)         cell-centre coordinate from interior surface
%     .cellLayer  (N x 1 int32)   layer index (1..nL) of each cell
%     .kC         (N x nkC)       k(T) polynomial coeffs per cell (ascending)
%     .cp         (N x 1)         cp per cell (constant per material for now)
%     .rho        (N x 1)         rho per cell
%     .dInt       scalar          interior surface -> first cell-centre dist
%     .dExt       scalar          last cell-centre -> exterior surface dist
%     .dFace      (N-1 x 1)       inter cell-centre distances
%     .hAmb, .Tinf, .TextRad      ambient-side BC parameters

    targetTags = int32(targetTags(:).');
    wtags      = int32([wallDataFurnace.tag]);

    % Empty 1x0 struct array with the full field set (so wallCtx(end+1)=s works).
    wallCtx = struct('zoneTag', {}, 'faceTag', {}, 'wallName', {}, ...
                     'meshFaces', {}, 'area', {}, 'nLayers', {}, ...
                     'layerName', {}, 'N', {}, 'dx', {}, 'xc', {}, ...
                     'cellLayer', {}, 'kC', {}, 'cp', {}, 'rho', {}, ...
                     'dInt', {}, 'dExt', {}, 'dFace', {}, ...
                     'hAmb', {}, 'Tinf', {}, 'TextRad', {});

    Nz = numel(zonesCtx.id);

    for k = 1:Nz
        bt = zonesCtx.byFaceTag{k};
        if isempty(bt) || ~isfield(bt, 'faceTag') || isempty(bt.faceTag)
            continue
        end
        zoneTag = double(zonesCtx.id(k));

        for j = 1:numel(bt.faceTag)
            tg = int32(bt.faceTag(j));
            if ~ismember(tg, targetTags)
                continue
            end

            wi = find(wtags == tg, 1);
            if isempty(wi)
                error(['build_wall_context: no wallDataFurnace entry for ', ...
                       'wall tag %d (zone %d).'], tg, zoneTag);
            end
            wentry = wallDataFurnace(wi);

            faceSigned = bt.faceIdx{j};
            meshFaces  = int32(abs(double(faceSigned(:))));
            area       = sum(double(geom.faceArea(meshFaces)));

            s = build_one_column(wentry, matProps);
            s.zoneTag   = zoneTag;
            s.faceTag   = double(tg);
            s.meshFaces = meshFaces;
            s.area      = area;

            wallCtx(end+1) = s; %#ok<AGROW>
        end
    end

    if isempty(wallCtx)
        warning(['build_wall_context: no (zone, tag) column built. Check ', ...
                 'that targetTags actually border some zone in ', ...
                 'zonesCtx.byFaceTag.']);
    end
end


%% LOCAL FUNCTIONS

function s = build_one_column(wentry, matProps)
%BUILD_ONE_COLUMN Geometry + material discretisation for a single wall.

    nL = wentry.nLayers;

    nCell = zeros(nL, 1);
    dxL   = zeros(nL, 1);
    for l = 1:nL
        t = wentry.layers(l).thickness;
        if ~isfinite(t) || t <= 0
            error(['build_wall_context: wall tag %d, layer %d has invalid ', ...
                   'thickness (%g). This tag cannot be modelled.'], ...
                  wentry.tag, l, t);
        end
        npm      = wentry.layers(l).nodesPerMeter;
        nCell(l) = max(1, round(npm * t));
        dxL(l)   = t / nCell(l);
    end

    N         = sum(nCell);
    dx        = zeros(N, 1);
    cellLayer = zeros(N, 1, 'int32');
    p = 0;
    for l = 1:nL
        idx            = p + (1:nCell(l));
        dx(idx)        = dxL(l);
        cellLayer(idx) = int32(l);
        p              = p + nCell(l);
    end

    xc    = cumsum(dx) - dx / 2;      % cell centres from interior surface
    dInt  = dx(1) / 2;
    dExt  = dx(N) / 2;
    dFace = dx(1:end-1) / 2 + dx(2:end) / 2;

    % --- Material properties per layer, then broadcast to cells ---
    layerName = cell(1, nL);
    kByLayer  = cell(nL, 1);
    cpByLayer = zeros(nL, 1);
    rhoByLayer = zeros(nL, 1);
    maxK = 0;
    for l = 1:nL
        nm = strrep(wentry.layers(l).materialName, '-', '_');
        layerName{l} = nm;
        if ~isfield(matProps.extWall, nm)
            error(['build_wall_context: material "%s" (wall tag %d, layer ', ...
                   '%d) not found in matProps.extWall.'], nm, wentry.tag, l);
        end
        mp           = matProps.extWall.(nm);
        kByLayer{l}  = double(mp.k(:).');       % ascending-order coeffs
        maxK         = max(maxK, numel(kByLayer{l}));
        cpByLayer(l) = double(mp.cp(1));        % constant per material
        rhoByLayer(l) = double(mp.rho(1));
    end

    kC  = zeros(N, maxK);
    cp  = zeros(N, 1);
    rho = zeros(N, 1);
    for l = 1:nL
        mask          = (cellLayer == l);
        cc            = kByLayer{l};
        kC(mask, 1:numel(cc)) = repmat(cc, nnz(mask), 1);
        cp(mask)      = cpByLayer(l);
        rho(mask)     = rhoByLayer(l);
    end

    % --- Pack (zone-specific fields filled by the caller) ---
    s = struct();
    s.zoneTag   = [];
    s.faceTag   = [];
    s.wallName  = wentry.wallName;
    s.meshFaces = [];
    s.area      = [];
    s.nLayers   = nL;
    s.layerName = layerName;
    s.N         = N;
    s.dx        = dx;
    s.xc        = xc;
    s.cellLayer = cellLayer;
    s.kC        = kC;
    s.cp        = cp;
    s.rho       = rho;
    s.dInt      = dInt;
    s.dExt      = dExt;
    s.dFace     = dFace;
    s.hAmb      = wentry.hAmb;
    s.Tinf      = wentry.Tinf;
    s.TextRad   = wentry.TextRad;
end
