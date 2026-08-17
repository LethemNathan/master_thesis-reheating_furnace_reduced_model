function slabCtxs = build_slab_rear_coupling(slabCtxs, geom, bcConfig)

    Ns = numel(slabCtxs);
    if Ns < 1
        error('build_slab_rear_coupling: empty slabCtxs.');
    end

    %% Slab 1: imposed Dirichlet
    sc1 = slabCtxs(1);
    rear1 = sc1.surfaceFaces.rear;
    Nrear1 = numel(rear1.cellLocal);

    if Nrear1 == 0
        error('build_slab_rear_coupling: slab %d has no rear face.', sc1.id);
    end

    if ~isfield(bcConfig.slabs, 'firstSlabRearT')
        error(['build_slab_rear_coupling: bcConfig.slabs.firstSlabRearT is ', ...
               'missing.']);
    end

    Timposed = bcConfig.slabs.firstSlabRearT;
    if isscalar(Timposed)
        Tvec = double(Timposed) * ones(Nrear1, 1);
    elseif numel(Timposed) == Nrear1
        Tvec = double(Timposed(:));
    else
        error(['build_slab_rear_coupling: bcConfig.slabs.firstSlabRearT must ', ...
               'be scalar or a vector of length %d (got %d).'], ...
               Nrear1, numel(Timposed));
    end

    rearBC = struct();
    rearBC.type           = "imposed";
    rearBC.cellLocal      = rear1.cellLocal;
    rearBC.dir            = rear1.dir;
    rearBC.faceGlobal     = rear1.faceGlobal;
    rearBC.T              = Tvec;
    rearBC.prevCellGlobal = zeros(Nrear1, 1, 'int32');
    rearBC.prevCellLocal  = zeros(Nrear1, 1, 'int32');
    rearBC.prevSlabIdx    = int32(0);

    slabCtxs(1).rearBC = rearBC;

    %% Slabs 2..N
    for i = 2:Ns
        slabCtxs(i).rearBC = match_rear_to_prev_front( ...
            slabCtxs(i), slabCtxs(i-1), int32(i-1), geom);
    end

    %% Build per-slab (cellLocal, dir) -> rearBC.T index map.
    % Used by assemble_slab_system to find T_face in O(1) for a Dirichlet
    % face. Entry is 0 for any (P, d) that is not a rear-Dirichlet face.
    for i = 1:Ns
        Nc_i      = slabCtxs(i).Nc;
        rb        = slabCtxs(i).rearBC;
        rearBCIdx = zeros(Nc_i, 6, 'int32');
        for k = 1:numel(rb.cellLocal)
            rearBCIdx(rb.cellLocal(k), rb.dir(k)) = int32(k);
        end
        slabCtxs(i).rearBCIdx = rearBCIdx;
    end
end


%% LOCAL FUNCTIONS

function rearBC = match_rear_to_prev_front(scI, scIm1, prevIdx, geom, tol)
%MATCH_REAR_TO_PREV_FRONT For each rear face of slab i, find the matching
% front face of slab i-1 by (y, z) coordinates of face centres.

    if nargin < 5, tol = 1e-6; end

    rearI    = scI.surfaceFaces.rear;
    frontIm1 = scIm1.surfaceFaces.front;

    NrearI    = numel(rearI.cellLocal);
    NfrontIm1 = numel(frontIm1.cellLocal);

    if NrearI == 0
        error('match_rear_to_prev_front: slab %d has no rear face.', scI.id);
    end
    if NfrontIm1 == 0
        error(['match_rear_to_prev_front: slab %d (previous) has no front ', ...
               'face.'], scIm1.id);
    end
    if NrearI ~= NfrontIm1
        error(['match_rear_to_prev_front: slab %d has %d rear faces but ', ...
               'previous slab %d has %d front faces. Mesh mismatch.'], ...
               scI.id, NrearI, scIm1.id, NfrontIm1);
    end

    % (y, z) coordinates of the face centres
    rearYZ  = geom.faceCenter(rearI.faceGlobal,    2:3);   % NrearI x 2
    frontYZ = geom.faceCenter(frontIm1.faceGlobal, 2:3);   % NfrontIm1 x 2

    prevCellLocal  = zeros(NrearI, 1, 'int32');
    prevCellGlobal = zeros(NrearI, 1, 'int32');
    used           = false(NfrontIm1, 1);

    for k = 1:NrearI
        diffs = frontYZ - rearYZ(k, :);
        d2    = sum(diffs.^2, 2);
        [mind2, idx] = min(d2);

        if sqrt(mind2) > tol
            error(['match_rear_to_prev_front: no front face match for rear ', ...
                   'face %d at (y,z)=(%.6g, %.6g) within tol=%g (min dist=%.3e).'], ...
                   rearI.faceGlobal(k), rearYZ(k,1), rearYZ(k,2), tol, sqrt(mind2));
        end
        if used(idx)
            error(['match_rear_to_prev_front: front face %d (slab %d) ', ...
                   'matched twice. Mapping is not bijective.'], ...
                   frontIm1.faceGlobal(idx), scIm1.id);
        end
        used(idx) = true;

        prevCellLocal(k)  = frontIm1.cellLocal(idx);
        prevCellGlobal(k) = scIm1.globalCellIdx(frontIm1.cellLocal(idx));
    end

    rearBC = struct();
    rearBC.type           = "coupled";
    rearBC.cellLocal      = rearI.cellLocal;
    rearBC.dir            = rearI.dir;
    rearBC.faceGlobal     = rearI.faceGlobal;
    rearBC.T              = zeros(NrearI, 1);   % filled at runtime
    rearBC.prevCellGlobal = prevCellGlobal;
    rearBC.prevCellLocal  = prevCellLocal;
    rearBC.prevSlabIdx    = prevIdx;
end
