function [I, info, radProps] = solve_radiation(geom, mesh, radProps, radCtx, ...
                                                D, ang, q_imposed, ...
                                                physicalConst, solverParams)
%SOLVE_RADIATION Outer Picard loop driving the sweep-based FVM radiation
% solver on the FINE gmsh mesh.
%
% At each iteration :
%   1. update_wall_intensities  -> refresh radProps.Ib_wall on walls of
%                                  radBcType in {2, 3} using the I from
%                                  the previous iteration.
%   2. For m = 1..M :
%        sweep_radiation_direction(..., I)  -> solve I(:, m) in place
%        (Gauss-Seidel for the symmetry BC ; once a direction is swept,
%         its values are visible to subsequent symmetry lookups in the
%         same iteration).
%   3. residual = ||I - I_old|| / ||I||   (L2, relative)
%
% The iteration stops when residual < picardTol or picardMaxIter is
% reached.
%
% Note : with epsilon_w = 1 and sigma_s = 0, the only coupling between
% directions comes from
%   - walls of radBcType in {2, 3} (Ib_w depends on the full incoming
%     radiation),
%   - walls of radBcType = 4 (symmetry, I^m at the face = I^{m_refl} of
%     the adjacent active cell).
% If you don't have any of these BCs, the loop converges in 2 iterations
% (one to populate I, one to confirm).
%
% Inputs
%   geom          : standard topology.
%   mesh          : with cells.isLoad, cells.isBlocked, cells.id.
%   radProps      : built by build_radiation_properties. Will be returned
%                   updated (Ib_wall, Tw_face on radBcType 2/3 faces).
%   radCtx        : with cellFaces, cellFaceSign (from build_rad_context).
%   D             : Nf x M directional weights.
%   ang           : output of build_angles.
%   q_imposed     : scalar W/m^2 for radBcType == 3 (positive = wall ->
%                   gas, entering the domain).
%   physicalConst : with .sigma.
%   solverParams  : optional struct, with fields
%                     .picardTol      (default 1e-4 on relative L2 dI)
%                     .picardMaxIter  (default 50)
%                     .verbose        (default true)
%                     .I0             (default zeros) initial guess
%                                     (Nc x M)
%
% Outputs
%   I        : Nc x M, converged intensity field.
%   info     : struct with .iter, .residuals (1 x iter), .converged.
%   radProps : same struct as input with Ib_wall / Tw_face refreshed on
%              radBcType 2 and 3 walls at convergence.

    %% --- Defaults ---
    picardTol     = 1e-4;
    picardMaxIter = 50;
    verbose       = true;
    I0            = [];

    if nargin >= 9 && ~isempty(solverParams)
        if isfield(solverParams, 'picardTol'),     picardTol     = solverParams.picardTol; end
        if isfield(solverParams, 'picardMaxIter'), picardMaxIter = solverParams.picardMaxIter; end
        if isfield(solverParams, 'verbose'),       verbose       = solverParams.verbose; end
        if isfield(solverParams, 'I0'),            I0            = solverParams.I0; end
    end

    Nc = numel(mesh.cells.id);
    M  = ang.M;

    %% --- Initial guess ---
    if isempty(I0)
        I = zeros(Nc, M);
    else
        if ~isequal(size(I0), [Nc, M])
            error(['solve_radiation: I0 must be Nc x M = %d x %d ', ...
                   '(got %d x %d).'], Nc, M, size(I0, 1), size(I0, 2));
        end
        I = I0;
    end

    info           = struct();
    info.iter      = 0;
    info.residuals = zeros(picardMaxIter, 1);
    info.converged = false;

    residual = inf;

    %% --- Outer Picard loop ---
    for iter = 1:picardMaxIter
        I_old = I;

        % 1. Refresh wall blackbody intensities on radBcType {2, 3}
        radProps = update_wall_intensities(geom, mesh, radProps, D, ang, ...
                                            I_old, q_imposed, physicalConst);

        % 2. Sweep every direction (Gauss-Seidel for the symmetry BC)
        for m = 1:M
            I(:, m) = sweep_radiation_direction(geom, mesh, radProps, ...
                                                 radCtx, D, ang, m, I);
        end

        % 3. Residual
        normI    = norm(I(:));
        residual = norm(I(:) - I_old(:)) / max(normI, eps);
        info.residuals(iter) = residual;
        info.iter            = iter;

        if verbose
            fprintf('  Rad iter %3d : ||dI|| / ||I|| = %.3e\n', iter, residual);
        end

        if residual < picardTol && iter > 1
            info.converged = true;
            break
        end
    end

    info.residuals = info.residuals(1:info.iter);

    if ~info.converged
        warning(['solve_radiation: did not converge in %d iterations ', ...
                 '(last residual %.3e).'], picardMaxIter, residual);
    end
end
