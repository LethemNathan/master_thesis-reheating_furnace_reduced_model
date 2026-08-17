function [T_slab, info] = solve_slab_conduction(slabCtx, fields, matProps, mode, dt, solverParams)
%SOLVE_SLAB_CONDUCTION Picard + direct (sparse Cholesky) solver on one slab.
%
% The radiative boundary condition uses an EXPLICIT surface temperature node
% per radiative face, distinct from the first-cell temperature (see
% assemble_slab_system). The linearised FVM system is symmetric positive
% definite, so it is solved with a direct sparse Cholesky factorisation at
% each Picard iteration; no iterative linear solver is used.
%
% WARM START. The cell temperatures start from fields.T (previous outer
% iteration). The surface temperatures start from fields.T_surf (the surfaces
% converged at the previous outer iteration), falling back to the adjacent
% first-cell temperature when T_surf is not available (first outer iteration).
%
% Iteration scheme (one slab, one time step):
%   x_lin = [T_cells ; T_surf]                          (warm start)
%   for k = 1..picardMaxIter
%       [A, b] = assemble_slab_system(..., x_lin, ..., smap)
%       R      = chol(A)                    % A = R'*R  (SPD)
%       x_dir  = R \ (R' \ b)
%       x_new  = alpha*x_dir + (1 - alpha)*x_lin        (under-relaxation)
%       if max|x_new - x_lin| < picardTol, break
%       x_lin  = x_new
%   end
%
% This function does NOT write the result back into fields.T; the caller
% (solve_slabs_conduction) is responsible for that. The converged surface
% temperatures are returned in info.T_surface so that the caller can feed
% them to the radiation solver for the radiative faces.
%
% Inputs
%   slabCtx       : slab context, with rearBC.T already refreshed by
%                   update_slab_dirichlet_bc.
%   fields        : struct
%                     .T      global Nc_glob x 1, cell warm start (via
%                             globalCellIdx).
%                     .G      global Nf_glob x 1, irradiation per face.
%                     .T_surf global Nf_glob x 1, surface warm start per
%                             radiative face (optional; NaN entries ignored).
%                     .T_n    global Nc_glob x 1, T at previous time step
%                             (REQUIRED if mode == "unsteady").
%   matProps      : see assemble_slab_system.
%   mode          : "steady" or "unsteady".
%   dt            : time step (used only when mode == "unsteady").
%   solverParams  : struct with optional fields (defaults applied below)
%                     .picardTol      (default 1e-3, in K)
%                     .picardMaxIter  (default 30)
%                     .relaxAlpha     (default 0.7, in (0, 1])
%                     .verbose        (default false)
%
% Outputs
%   T_slab : Nc x 1 converged CELL temperature on the slab (local).
%   info   : diagnostic struct with fields
%              .slabId, .converged, .picardIters
%              .picardResidual (vector of length picardIters)
%              .T_surface      (Ns x 1 converged surface temperatures)
%              .surfaceMap     (smap: radiative-face -> surface-node mapping)
%              .nSurf          (Ns)

    Nc = slabCtx.Nc;
    g  = slabCtx.globalCellIdx;

    %% Defaults for solverParams
    sp = solverParams;
    if ~isfield(sp, 'picardTol'),     sp.picardTol     = 1e-3;  end
    if ~isfield(sp, 'picardMaxIter'), sp.picardMaxIter = 30;    end
    if ~isfield(sp, 'relaxAlpha'),    sp.relaxAlpha    = 0.7;   end
    if ~isfield(sp, 'verbose'),       sp.verbose       = false; end
    if sp.relaxAlpha <= 0 || sp.relaxAlpha > 1
        error('solve_slab_conduction: relaxAlpha must be in (0, 1].');
    end

    %% Surface map: one surface node per radiative face
    smap = build_surface_map(slabCtx);
    Ns   = smap.count;

    %% Initial guess and T_n
    if ~isfield(fields, 'T') || isempty(fields.T)
        error('solve_slab_conduction: fields.T is missing or empty.');
    end
    T_cell0 = double(fields.T(g));               % Nc x 1

    % --- Surface warm start ---
    % Default: adjacent (first-cell) temperature. Then override with the
    % previous outer-iteration surface temperature fields.T_surf (per global
    % face) wherever it is available (not NaN).
    T_surf0 = T_cell0(smap.P);                   % Ns x 1  (fallback)
    if Ns > 0 && isfield(fields, 'T_surf') && ~isempty(fields.T_surf)
        linPd = sub2ind(size(slabCtx.cellFace), smap.P, smap.d);
        fg    = double(slabCtx.cellFace(linPd));         % global face of each node
        ts    = nan(Ns, 1);
        ok    = fg > 0 & fg <= numel(fields.T_surf);
        ts(ok) = fields.T_surf(fg(ok));
        use   = ~isnan(ts);
        T_surf0(use) = ts(use);
    end

    x_lin = [T_cell0; T_surf0];                  % N x 1 (augmented)

    isUnsteady = strcmpi(string(mode), "unsteady");
    if isUnsteady
        if ~isfield(fields, 'T_n') || isempty(fields.T_n)
            error(['solve_slab_conduction: fields.T_n is missing for ', ...
                   'unsteady mode.']);
        end
        T_n = double(fields.T_n(g));             % Nc x 1 (cells only)
    else
        T_n = [];
    end

    %% Diagnostic accumulators
    picardRes = zeros(sp.picardMaxIter, 1);
    converged = false;
    nDone     = 0;

    %% Picard loop
    for k = 1:sp.picardMaxIter
        % --- 1. Assemble linearised augmented FVM system at x_lin ---
        [A, b] = assemble_slab_system(slabCtx, x_lin, T_n, fields, ...
                                      matProps, mode, dt, smap);

        %% 2-3. Cholesky with fill-reducing reordering (AMD) + solve
        [R, pchol, s] = chol(A, 'vector');        % R'*R = A(s,s), s = AMD permutation
        if pchol == 0
            x_dir = zeros(size(b));
            x_dir(s) = R \ (R' \ b(s));           % solve on the reordered system
        else
            warning(['solve_slab_conduction: chol failed (slab %d, picard iter %d): ', ...
                     'A not SPD (flag=%d). Falling back to backslash.'], ...
                     slabCtx.id, k, pchol);
            x_dir = A \ b;
        end

        % --- 3. Under-relaxation ---
        alpha = sp.relaxAlpha;
        x_new = alpha * x_dir + (1 - alpha) * x_lin;

        % --- 4. Picard convergence test (cells and surfaces) ---
        deltaT = max(abs(x_new - x_lin));
        picardRes(k) = deltaT;
        nDone        = k;

        if sp.verbose
            fprintf('  [slab %2d] Picard %2d: ||dT||inf = %.3e K\n', ...
                    slabCtx.id, k, deltaT);
        end

        x_lin = x_new;
        if deltaT < sp.picardTol
            converged = true;
            break;
        end
    end

    if ~converged
        warning(['solve_slab_conduction: slab %d did not reach picardTol ', ...
                 '%.2e in %d iterations (last ||dT||inf = %.3e K).'], ...
                slabCtx.id, sp.picardTol, sp.picardMaxIter, picardRes(nDone));
    end

    %% Output
    T_slab = x_lin(1:Nc);

    info = struct();
    info.slabId         = slabCtx.id;
    info.converged      = converged;
    info.picardIters    = nDone;
    info.picardResidual = picardRes(1:nDone);
    info.T_surface      = x_lin(Nc+1:end);       % Ns x 1 surface temperatures
    info.surfaceMap     = smap;
    info.nSurf          = Ns;
end


%% LOCAL FUNCTIONS

function smap = build_surface_map(slabCtx)
%BUILD_SURFACE_MAP One surface node per radiative boundary face.
%
% Returns
%   smap.count      : Ns, number of radiative faces on the slab.
%   smap.surfOfFace : Nc x 6, local surface index (1..Ns) of each radiative
%                     face, 0 for non-radiative faces.
%   smap.P, smap.d  : Ns x 1, the (cell, direction) of each surface node.

    c  = loadTopologyCodes();
    Nc = slabCtx.Nc;

    surfOfFace = zeros(Nc, 6);
    Plist = zeros(6 * Nc, 1);
    dlist = zeros(6 * Nc, 1);
    s = 0;
    for P = 1:Nc
        for d = 1:6
            if slabCtx.faceBC(P, d) == c.BC_RADIATIVE
                s = s + 1;
                surfOfFace(P, d) = s;
                Plist(s) = P;
                dlist(s) = d;
            end
        end
    end

    smap.count      = s;
    smap.surfOfFace = surfOfFace;
    smap.P          = Plist(1:s);
    smap.d          = dlist(1:s);
end
