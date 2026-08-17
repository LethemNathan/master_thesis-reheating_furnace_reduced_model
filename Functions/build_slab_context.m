function slabCtx = build_slab_context(slab, mesh, geom, bcConfig)
%BUILD_SLAB_CONTEXT Build local data structure for one slab's conduction solver.
%
% Inputs
%   slab     : a single slab returned by build_slab() (fields: id, bbox, surfaces)
%   mesh     : full mesh struct (after relabel_load_cells)
%   geom     : geometry struct (after load_face_tag)
%   bcConfig : config struct with bcConfig.slabs.<subzone> in {"adiabatic",
%              "radiative", "dirichlet"} for each subzone in
%              {symmetry, bottom, side, top, front, rear}.
%
% Output
%   slabCtx : struct with fields
%     .id              slab id
%     .Nc              number of cells in the slab
%     .globalCellIdx   [Nc x 1]   global cell indices in mesh.cells
%     .globalToLocal   [NcGlob x 1] map global -> local (0 if cell not in slab)
%     .cellNeighbor    [Nc x 6]   local index of neighbour cell, 0 if BC face
%     .cellFace        [Nc x 6]   global face index
%     .faceBC          [Nc x 6]   BC type code (codes(): BC_*)
%     .faceSubzone     [Nc x 6]   sub-zone yy code (0 if internal)
%     .cellVolume      [Nc x 1]
%     .cellCenter      [Nc x 3]
%     .faceArea        [Nc x 6]
%     .faceDist        [Nc x 6]   centre-to-centre if internal,
%                                 cell-centre-to-face-centre if BC
%     .surfaceFaces    struct with one field per sub-zone (symmetry, bottom,
%                      rear, side, front, top), each containing
%                        .cellLocal  [Ns x 1]
%                        .dir        [Ns x 1]
%                        .faceGlobal [Ns x 1]
%
% Direction column ordering (codes(): DIR_*):
%   1=E (+y), 2=W (-y), 3=N (+z), 4=S (-z), 5=U (+x), 6=D (-x)

    c = loadTopologyCodes();

    %% Step A: identify slab cells
    inSlab = (mesh.cells.loadTag == slab.id) & mesh.cells.isLoad;
    globalCellIdx = find(inSlab);
    Nc = numel(globalCellIdx);

    if Nc == 0
        error('build_slab_context: no cells found for slab id %d.', slab.id);
    end

    %% Step B: build global -> local map
    NcGlob = numel(mesh.cells.id);
    globalToLocal = zeros(NcGlob, 1, 'int32');
    globalToLocal(globalCellIdx) = int32(1:Nc);

    %% Step C+D: per (cell, direction) topology and geometry
    cellNeighbor = zeros(Nc, 6, 'int32');
    cellFace     = zeros(Nc, 6, 'int32');
    faceBC       = -ones(Nc, 6, 'int32');   % -1 = sentinel "unassigned"
    faceSubzone  = zeros(Nc, 6, 'int32');
    faceArea     = zeros(Nc, 6);
    faceDist     = zeros(Nc, 6);

    Nf = numel(geom.owner);
    for f = 1:Nf
        oG = geom.owner(f);
        nG = geom.neighbour(f);

        oL = double(globalToLocal(oG));
        if nG > 0
            nL = double(globalToLocal(nG));
        else
            nL = 0;
        end

        % Skip faces with no slab cell on either side
        if oL == 0 && nL == 0
            continue;
        end

        % Encode each "side that is in the slab" as a row:
        %   [iLoc, iGlob, otherL, otherG, normalSign]
        sides = zeros(2, 5);
        nSides = 0;
        if oL > 0
            nSides = nSides + 1;
            sides(nSides, :) = [oL, oG, nL, nG, +1];
        end
        if nL > 0
            nSides = nSides + 1;
            sides(nSides, :) = [nL, nG, oL, oG, -1];
        end

        for s = 1:nSides
            iLoc   = sides(s, 1);
            iGlob  = sides(s, 2);
            otherL = sides(s, 3);
            otherG = sides(s, 4);
            nSign  = sides(s, 5);

            n_out = nSign * geom.faceNormal(f, :);
            d = direction_from_normal(n_out, c);

            % Detect double assignment (a cell shouldn't have two faces in
            % the same direction for an axis-aligned hex)
            if cellFace(iLoc, d) ~= 0
                error(['build_slab_context: cell %d (global %d) already has ', ...
                       'face %d in direction %d, trying to add face %d.'], ...
                       iLoc, iGlob, cellFace(iLoc, d), d, f);
            end

            cellFace(iLoc, d) = f;
            faceArea(iLoc, d) = geom.faceArea(f);

            if otherL > 0
                % Internal face within the slab
                cellNeighbor(iLoc, d) = int32(otherL);
                faceBC(iLoc, d)       = c.BC_INTERNAL;
                faceSubzone(iLoc, d)  = 0;
                faceDist(iLoc, d)     = norm( ...
                    geom.cellCenter(otherG, :) - geom.cellCenter(iGlob, :));
            else
                % Surface face: BC depends on sub-zone
                cellNeighbor(iLoc, d) = 0;
                lt = double(geom.loadTag(f));
                yy = mod(lt, 100);

                if yy == 0
                    error(['build_slab_context: surface face %d of cell %d ', ...
                           '(global %d) has loadTag=%d but no sub-zone code.'], ...
                           f, iLoc, iGlob, lt);
                end

                faceSubzone(iLoc, d) = int32(yy);
                faceBC(iLoc, d)      = subzone_to_bc(yy, bcConfig, c);
                faceDist(iLoc, d)    = norm( ...
                    geom.faceCenter(f, :) - geom.cellCenter(iGlob, :));
            end
        end
    end

    %% Step E: sanity checks
    if any(cellFace(:) == 0)
        [iBad, dBad] = find(cellFace == 0, 1);
        error(['build_slab_context: cell %d (global %d) is missing a face ', ...
               'in direction %d. Hex cell expected to have 6 faces.'], ...
               iBad, globalCellIdx(iBad), dBad);
    end
    if any(faceBC(:) == -1)
        error('build_slab_context: some faces have unassigned BC type.');
    end
    for i = 1:Nc
        if numel(unique(cellFace(i, :))) ~= 6
            error(['build_slab_context: cell %d has duplicate face indices ', ...
                   'across directions.'], i);
        end
    end

    %% Step F: aggregate surface faces by sub-zone
    surfaceFaces = struct();
    subzoneList = {
        'symmetry', c.SUBZONE_SYMMETRY;
        'bottom',   c.SUBZONE_BOTTOM;
        'rear',     c.SUBZONE_REAR;
        'side',     c.SUBZONE_SIDE;
        'front',    c.SUBZONE_FRONT;
        'top',      c.SUBZONE_TOP;
    };
    for k = 1:size(subzoneList, 1)
        name = subzoneList{k, 1};
        code = subzoneList{k, 2};
        mask = (faceSubzone == code);
        [iLocs, dirs] = find(mask);
        sub = struct();
        sub.cellLocal  = int32(iLocs);
        sub.dir        = int32(dirs);
        sub.faceGlobal = cellFace(mask);
        surfaceFaces.(name) = sub;
    end

    %% Build output struct
    slabCtx = struct();
    slabCtx.id            = slab.id;
    slabCtx.Nc            = Nc;
    slabCtx.globalCellIdx = int32(globalCellIdx);
    slabCtx.globalToLocal = globalToLocal;
    slabCtx.cellNeighbor  = cellNeighbor;
    slabCtx.cellFace      = cellFace;
    slabCtx.faceBC        = faceBC;
    slabCtx.faceSubzone   = faceSubzone;
    slabCtx.cellVolume    = geom.cellVolume(globalCellIdx);
    slabCtx.cellCenter    = geom.cellCenter(globalCellIdx, :);
    slabCtx.faceArea      = faceArea;
    slabCtx.faceDist      = faceDist;
    slabCtx.surfaceFaces  = surfaceFaces;
end


%% LOCAL FUNCTIONS

function d = direction_from_normal(n_out, c)
%DIRECTION_FROM_NORMAL Map an axis-aligned outward normal to a direction code.
    [maxAbs, ax] = max(abs(n_out));
    if maxAbs < 0.9
        error(['direction_from_normal: outward normal (%.3g, %.3g, %.3g) ', ...
               'is not axis-aligned (max component %.3g < 0.9).'], ...
               n_out(1), n_out(2), n_out(3), maxAbs);
    end
    sgn = n_out(ax) > 0;
    if ax == 1
        if sgn, d = c.DIR_U; else, d = c.DIR_D; end
    elseif ax == 2
        if sgn, d = c.DIR_E; else, d = c.DIR_W; end
    else
        if sgn, d = c.DIR_N; else, d = c.DIR_S; end
    end
end


function bcType = subzone_to_bc(yy, bcConfig, c)
%SUBZONE_TO_BC Resolve a sub-zone code + bcConfig into a BC type code.
    switch yy
        case c.SUBZONE_SYMMETRY, cfg = bcConfig.slabs.symmetry;
        case c.SUBZONE_BOTTOM,   cfg = bcConfig.slabs.bottom;
        case c.SUBZONE_REAR,     cfg = bcConfig.slabs.rear;
        case c.SUBZONE_SIDE,     cfg = bcConfig.slabs.side;
        case c.SUBZONE_FRONT,    cfg = bcConfig.slabs.front;
        case c.SUBZONE_TOP,      cfg = bcConfig.slabs.top;
        otherwise
            error('subzone_to_bc: unknown sub-zone code %d.', yy);
    end

    switch lower(string(cfg))
        case "adiabatic", bcType = c.BC_ADIABATIC;
        case "radiative", bcType = c.BC_RADIATIVE;
        case "dirichlet", bcType = c.BC_DIRICHLET;
        otherwise
            error('subzone_to_bc: unknown BC type "%s" (sub-zone yy=%d).', ...
                   string(cfg), yy);
    end
end
