function [Twall_zone_table, wallState, info] = solve_walls_conduction( ...
        wallCtx, radFields, geom, Twall_zone_table, physicalConst, opts)
%SOLVE_WALLS_CONDUCTION Driver: solve every 1D wall column and feed the
% interior-surface temperatures back into Twall_zone_table.
%
% For each (zone, tag) column in wallCtx this function:
%   1. area-averages the incident radiation radFields.G over the column's
%      mesh faces  ->  G_col = sum(G_f * A_f) / sum(A_f)   [W/m^2] ;
%   2. looks up the imposed interior convective flux q_conv_in and the
%      corrector flux q_corr for that (zoneTag, faceTag) ;
%   3. solves the 1D steady multilayer conduction (solve_wall_conduction),
%      warm-started from the previous outer iteration if available ;
%   4. writes the interior-surface temperature Ts_int into Twall_zone_table
%      at key '<zoneTag>_<faceTag>', so the next build_radiation_properties
%      call uses the PREDICTED wall temperature for these tags.
%
% This is meant to be called once per outer (radiation/conduction) iteration
% in main.m, after compute_G_field has refreshed radFields.G.
%
% Inputs
%   wallCtx          : struct array from build_wall_context.
%   radFields        : with .G (Nf x 1, W/m^2) from compute_G_field.
%   geom             : with .faceArea (Nf x 1).
%   Twall_zone_table : containers.Map, key '<zoneTag>_<faceTag>' -> Tw [K].
%   physicalConst    : with .sigma.
%   opts             : optional struct
%       .qConvMap    interior convective flux source. Either a numeric/
%                    sparse matrix indexed (zoneTag, faceTag) [e.g.
%                    extWallReport.convectiveHeatFlux], or a containers.Map
%                    keyed '<zoneTag>_<faceTag>'. Missing entries -> 0.
%                    SIGN CONVENTION: positive = flux INTO the wall.
%       .qCorrMap    corrector flux, same format. Default: 0 everywhere.
%       .solverParams  passed to solve_wall_conduction (defaults there).
%       .wallState   1 x Ncol struct array from a previous call, used for
%                    warm starts (field .u). Default: cold start.
%       .verbose     print a per-column line. Default false.
%
% Outputs
%   Twall_zone_table : same map, updated in place (also returned).
%   wallState        : 1 x Ncol struct array with per-column results:
%                        .zoneTag .faceTag .Tsi .Tse .u .qIn .qOut
%                        .energyImbalance .iters .converged
%                      Pass it back in via opts.wallState next iteration.
%   info             : struct with
%                        .nCols, .nNotConverged, .maxAbsImbalance,
%                        .maxIters, .Tsi (Ncol x 1), .keys (Ncol x 1 cell)

    if nargin < 6 || isempty(opts), opts = struct(); end
    qConvMap     = getopt(opts, 'qConvMap',     []);
    qCorrMap     = getopt(opts, 'qCorrMap',     []);
    solverParams = getopt(opts, 'solverParams', struct());
    prevState    = getopt(opts, 'wallState',    []);
    verbose      = getopt(opts, 'verbose',      false);

    Ncol = numel(wallCtx);

    % Template so wallState keeps a consistent field set.
    stateTmpl = struct('zoneTag', [], 'faceTag', [], 'Tsi', [], 'Tse', [], ...
                       'u', [], 'qIn', [], 'qOut', [], ...
                       'energyImbalance', [], 'iters', [], 'converged', []);
    wallState = repmat(stateTmpl, 1, Ncol);

    faceArea = double(geom.faceArea);
    G        = double(radFields.G);

    nNotConv    = 0;
    maxAbsImbal = 0;
    maxIters    = 0;
    TsiVec      = zeros(Ncol, 1);
    keyList     = cell(Ncol, 1);

    for c = 1:Ncol
        wc      = wallCtx(c);
        zoneTag = wc.zoneTag;
        faceTag = wc.faceTag;
        faces   = double(wc.meshFaces);

        % --- 1. Area-averaged incident radiation over the column ---
        Af    = faceArea(faces);
        sumAf = sum(Af);
        if sumAf <= 0
            warning(['solve_walls_conduction: zero total area for (zone %d, ', ...
                     'tag %d) ; skipping.'], zoneTag, faceTag);
            continue
        end
        G_col = sum(G(faces) .* Af) / sumAf;

        % --- 2. Imposed interior fluxes for this (zoneTag, faceTag) ---
        drive         = struct();
        drive.G       = G_col;
        drive.qConvIn = lookup_flux(qConvMap, zoneTag, faceTag);
        drive.qCorr   = lookup_flux(qCorrMap, zoneTag, faceTag);

        % --- 3. Warm start from previous outer iteration if compatible ---
        u0 = [];
        if ~isempty(prevState) && c <= numel(prevState) ...
                && isfield(prevState(c), 'u') && ~isempty(prevState(c).u) ...
                && numel(prevState(c).u) == wc.N + 2
            u0 = prevState(c).u;
        end

        % --- 4. Solve the column ---
        [out, cinfo] = solve_wall_conduction(wc, drive, physicalConst, ...
                                             solverParams, u0);

        % --- 5. Write predicted interior-surface temperature back ---
        key = sprintf('%d_%d', zoneTag, faceTag);
        Twall_zone_table(key) = out.Tsi;

        % --- Bookkeeping ---
        wallState(c).zoneTag         = zoneTag;
        wallState(c).faceTag         = faceTag;
        wallState(c).Tsi             = out.Tsi;
        wallState(c).Tse             = out.Tse;
        wallState(c).u               = out.u;
        wallState(c).qIn             = out.qIn;
        wallState(c).qOut            = out.qOut;
        wallState(c).energyImbalance = cinfo.energyImbalance;
        wallState(c).iters           = cinfo.iters;
        wallState(c).converged       = cinfo.converged;

        TsiVec(c)  = out.Tsi;
        keyList{c} = key;

        nNotConv    = nNotConv + double(~cinfo.converged);
        maxAbsImbal = max(maxAbsImbal, abs(cinfo.energyImbalance));
        maxIters    = max(maxIters, cinfo.iters);

        if verbose
            fprintf(['  [wall %2d/%2d] %s : G=%.3e W/m2 qConv=%.3e qCorr=%.3e ', ...
                     '-> Tsi=%.1f K, Tse=%.1f K (imbal=%.2e, it=%d)\n'], ...
                    c, Ncol, key, G_col, drive.qConvIn, drive.qCorr, ...
                    out.Tsi, out.Tse, cinfo.energyImbalance, cinfo.iters);
        end
    end

    info = struct();
    info.nCols           = Ncol;
    info.nNotConverged   = nNotConv;
    info.maxAbsImbalance = maxAbsImbal;
    info.maxIters        = maxIters;
    info.Tsi             = TsiVec;
    info.keys            = keyList;
end


%% LOCAL FUNCTIONS

function v = lookup_flux(src, zoneTag, faceTag)
%LOOKUP_FLUX Read a per-(zoneTag,faceTag) flux from a matrix or a Map.
% Returns 0 if src is empty or the pair is absent. Result in W/m^2.
    v = 0;
    if isempty(src)
        return
    end
    if isa(src, 'containers.Map')
        key = sprintf('%d_%d', zoneTag, faceTag);
        if isKey(src, key)
            v = full(double(src(key)));
        end
        return
    end
    % Numeric / sparse matrix indexed by (zoneTag, faceTag).
    [nr, nc] = size(src);
    if zoneTag >= 1 && zoneTag <= nr && faceTag >= 1 && faceTag <= nc
        v = full(double(src(zoneTag, faceTag)));
    end
end

function v = getopt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end
