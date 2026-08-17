function [out, info] = solve_wall_conduction(wc, drive, physicalConst, ...
                                             solverParams, u0)
%SOLVE_WALL_CONDUCTION Steady 1D multilayer conduction on a single wall column.
%
% Solves the temperature field through one wall (one (zone, tag) column built
% by build_wall_context), given the interior-side drive (incident radiation
% and imposed convective / corrector fluxes) and the ambient-side BC stored
% in the column.
%
% Unknowns  u = [ Ts_int ; T_1 ; ... ; T_N ; Ts_ext ]   (length N+2)
%   Ts_int : interior (furnace) surface temperature  [K]  <- the quantity
%            fed back into the radiation solver.
%   T_i    : cell-centre temperatures across the thickness [K].
%   Ts_ext : ambient surface temperature [K].
%
% Physics (all fluxes per unit area, W/m^2), x = 0 interior -> x = Ltot ambient:
%
%   Interior surface energy balance (net flux entering the wall):
%       G - sigma*Ts_int^4 + q_conv_in + q_corr  =  k_1/dInt * (Ts_int - T_1)
%     with eps = 1 (blackbody emission of the wall face). G is the incident
%     radiation, q_conv_in the imposed interior convective flux, q_corr the
%     corrector flux (second convective flux for calibration).
%
%   Interior cells: steady conduction, no volumetric source
%       Gf_{i-1}(T_{i-1}-T_i) - Gf_i(T_i-T_{i+1}) = 0
%     with face conductance (series resistance, exact at layer interfaces)
%       Gf_i = 1 / ( dx_i/(2 k_i) + dx_{i+1}/(2 k_{i+1}) ).
%
%   Exterior surface energy balance (net flux leaving the wall):
%       k_N/dExt * (T_N - Ts_ext) = h_amb*(Ts_ext - Tinf)
%                                   + sigma*(Ts_ext^4 - TextRad^4)   (eps = 1)
%
% Summing every equation telescopes to the global balance
%       q_in(Ts_int) = q_out(Ts_ext),
% so energy is conserved to solver tolerance at convergence.
%
% Nonlinearity: k(T) (temperature-dependent conductivity) and the two
% sigma*T^4 surface terms. Solved by Picard iteration: conductances are
% frozen at the current T, the two radiative terms are Newton-linearised
% (sigma*T^4 ~ 4*sigma*T0^3*T - 3*sigma*T0^4), the resulting TRIDIAGONAL
% linear system is solved, under-relaxed, and repeated. Same pattern and
% conventions as solve_slab_conduction.
%
% Inputs
%   wc            : one element of the wallCtx array (build_wall_context).
%   drive         : struct with scalar fields
%                     .G        area-averaged incident radiation  [W/m^2]
%                     .qConvIn  imposed interior convective flux  [W/m^2]
%                     .qCorr    corrector flux                    [W/m^2]
%                   Missing fields default to 0.
%   physicalConst : with .sigma (Stefan-Boltzmann).
%   solverParams  : optional struct
%                     .picardTol     (default 1e-6, in K)
%                     .picardMaxIter (default 200)
%                     .relaxAlpha    (default 0.9, in (0,1])
%                     .verbose       (default false)
%   u0            : optional warm start (N+2 x 1). If omitted, a radiative-
%                   equilibrium guess is used.
%
% Outputs
%   out : struct with
%           .Tsi     interior surface temperature [K]  (= wall Tw for rad)
%           .Tse     exterior surface temperature [K]
%           .Tcells  (N x 1) cell temperatures [K]
%           .u       (N+2 x 1) full solution vector
%           .qIn     net interior flux entering the wall [W/m^2]
%           .qOut    net exterior flux leaving the wall  [W/m^2]
%           .qCond   conductive flux interior surface -> cell 1 [W/m^2]
%   info : struct with .converged, .iters, .dTmax (last update),
%          .energyImbalance (= qIn - qOut, W/m^2).

    sigma = physicalConst.sigma;
    N     = wc.N;

    if N < 2
        error('solve_wall_conduction: column needs N >= 2 cells (got %d).', N);
    end

    %% --- Drive defaults ---
    if nargin < 2 || isempty(drive), drive = struct(); end
    G       = getfielddef(drive, 'G',       0);
    qConvIn = getfielddef(drive, 'qConvIn', 0);
    qCorr   = getfielddef(drive, 'qCorr',   0);
    qSrcInt = G + qConvIn + qCorr;    % T-independent part of interior drive

    %% --- Solver params defaults ---
    if nargin < 4 || isempty(solverParams), solverParams = struct(); end
    sp = solverParams;
    if ~isfield(sp, 'picardTol'),     sp.picardTol     = 1e-6; end
    if ~isfield(sp, 'picardMaxIter'), sp.picardMaxIter = 200;  end
    if ~isfield(sp, 'relaxAlpha'),    sp.relaxAlpha     = 0.9; end
    if ~isfield(sp, 'verbose'),       sp.verbose        = false; end
    if sp.relaxAlpha <= 0 || sp.relaxAlpha > 1
        error('solve_wall_conduction: relaxAlpha must be in (0,1].');
    end

    %% --- Geometry / material shortcuts ---
    dx   = wc.dx(:);
    dInt = wc.dInt;
    dExt = wc.dExt;
    kC   = wc.kC;
    h    = wc.hAmb;
    Tinf = wc.Tinf;
    Trad = wc.TextRad;

    M = N + 2;                 % total unknowns
    iSi = 1;                   % index of Ts_int
    iC  = (2:N+1).';           % indices of T_1..T_N
    iSe = N + 2;               % index of Ts_ext

    %% --- Initial guess ---
    if nargin >= 5 && ~isempty(u0)
        u = double(u0(:));
        if numel(u) ~= M
            error('solve_wall_conduction: u0 must have length N+2 = %d.', M);
        end
    else
        % Radiative-equilibrium interior estimate, linear ramp to ambient.
        Tsi0 = max((max(G, 0) / sigma)^0.25, Tinf);
        if ~isfinite(Tsi0) || Tsi0 <= 0, Tsi0 = 1000; end
        u        = zeros(M, 1);
        u(iSi)   = Tsi0;
        u(iC)    = linspace(Tsi0, Tinf + 50, N).';
        u(iSe)   = Tinf + 50;
    end

    %% --- Picard loop ---
    converged = false;
    dTmax     = NaN;
    it        = 0;

    for it = 1:sp.picardMaxIter
        Told = u;

        Tcell  = u(iC);
        Tsi    = u(iSi);
        Tse    = u(iSe);

        % Conductivity per cell and face conductances (frozen this iteration).
        kcell = eval_k(kC, Tcell);                 % N x 1, > 0
        Gint  = kcell(1) / dInt;                    % surface -> cell 1
        Gext  = kcell(N) / dExt;                    % cell N -> surface
        if N >= 2
            Gf = 1 ./ ( dx(1:N-1) ./ (2*kcell(1:N-1)) ...
                      + dx(2:N)   ./ (2*kcell(2:N)) ); % (N-1) x 1
        else
            Gf = [];
        end

        % Newton linearisation of the two sigma*T^4 surface terms.
        aInt = 4 * sigma * Tsi^3;   cInt = 3 * sigma * Tsi^4;
        aExt = 4 * sigma * Tse^3;   cExt = 3 * sigma * Tse^4;

        % Assemble tridiagonal system A u = b (row ordering matches unknowns).
        I = zeros(3*M, 1); J = I; V = I; nz = 0;
        b = zeros(M, 1);

        % Row 1 : interior surface balance
        %   (aInt + Gint)*Tsi - Gint*T_1 = qSrcInt + cInt
        [I,J,V,nz] = addv(I,J,V,nz, iSi, iSi, aInt + Gint);
        [I,J,V,nz] = addv(I,J,V,nz, iSi, iC(1), -Gint);
        b(iSi) = qSrcInt + cInt;

        % Cell rows
        for i = 1:N
            r = iC(i);
            if i == 1
                % -Gint*Tsi + (Gint+Gf_1)*T_1 - Gf_1*T_2 = 0
                [I,J,V,nz] = addv(I,J,V,nz, r, iSi, -Gint);
                if N == 1
                    % (handled by N>=2 guard, kept for completeness)
                else
                    [I,J,V,nz] = addv(I,J,V,nz, r, r,      Gint + Gf(1));
                    [I,J,V,nz] = addv(I,J,V,nz, r, iC(2), -Gf(1));
                end
            elseif i == N
                % -Gf_{N-1}*T_{N-1} + (Gf_{N-1}+Gext)*T_N - Gext*Tse = 0
                [I,J,V,nz] = addv(I,J,V,nz, r, iC(N-1), -Gf(N-1));
                [I,J,V,nz] = addv(I,J,V,nz, r, r,        Gf(N-1) + Gext);
                [I,J,V,nz] = addv(I,J,V,nz, r, iSe,     -Gext);
            else
                % -Gf_{i-1}*T_{i-1} + (Gf_{i-1}+Gf_i)*T_i - Gf_i*T_{i+1} = 0
                [I,J,V,nz] = addv(I,J,V,nz, r, iC(i-1), -Gf(i-1));
                [I,J,V,nz] = addv(I,J,V,nz, r, r,        Gf(i-1) + Gf(i));
                [I,J,V,nz] = addv(I,J,V,nz, r, iC(i+1), -Gf(i));
            end
        end

        % Row N+2 : exterior surface balance
        %   -Gext*T_N + (Gext + h + aExt)*Tse = h*Tinf + cExt + sigma*Trad^4
        [I,J,V,nz] = addv(I,J,V,nz, iSe, iC(N), -Gext);
        [I,J,V,nz] = addv(I,J,V,nz, iSe, iSe,    Gext + h + aExt);
        b(iSe) = h * Tinf + cExt + sigma * Trad^4;

        A  = sparse(I(1:nz), J(1:nz), V(1:nz), M, M);
        us = A \ b;

        % Under-relaxation
        u     = sp.relaxAlpha * us + (1 - sp.relaxAlpha) * Told;
        dTmax = max(abs(u - Told));

        if sp.verbose
            fprintf(['  [wall z%d t%d] Picard %3d: ||dT||inf = %.3e K ', ...
                     '(Tsi=%.1f, Tse=%.1f)\n'], wc.zoneTag, wc.faceTag, ...
                    it, dTmax, u(iSi), u(iSe));
        end

        if dTmax < sp.picardTol
            converged = true;
            break
        end
    end

    if ~converged
        warning(['solve_wall_conduction: (zone %d, tag %d) did not reach ', ...
                 'picardTol %.2e in %d iters (last ||dT||inf = %.3e K).'], ...
                wc.zoneTag, wc.faceTag, sp.picardTol, sp.picardMaxIter, dTmax);
    end

    %% --- Outputs and diagnostics (evaluated at the converged state) ---
    Tsi   = u(iSi);
    Tse   = u(iSe);
    Tcell = u(iC);
    kcell = eval_k(kC, Tcell);
    Gint  = kcell(1) / dInt;

    qIn   = qSrcInt - sigma * Tsi^4;                       % net interior flux in
    qOut  = h * (Tse - Tinf) + sigma * (Tse^4 - Trad^4);   % net exterior flux out
    qCond = Gint * (Tsi - Tcell(1));

    out = struct();
    out.Tsi    = Tsi;
    out.Tse    = Tse;
    out.Tcells = Tcell;
    out.u      = u;
    out.qIn    = qIn;
    out.qOut   = qOut;
    out.qCond  = qCond;

    info = struct();
    info.converged      = converged;
    info.iters          = it;
    info.dTmax          = dTmax;
    info.energyImbalance = qIn - qOut;
end


%% LOCAL FUNCTIONS

function kc = eval_k(kC, T)
%EVAL_K Vectorised Horner evaluation of per-cell k(T). kC is N x p ascending.
    p  = size(kC, 2);
    kc = kC(:, p);
    for j = p-1:-1:1
        kc = kc .* T + kC(:, j);
    end
end

function [I,J,V,nz] = addv(I,J,V,nz, r, c, v)
%ADDV Append one (row, col, val) triplet to the sparse accumulators.
    nz    = nz + 1;
    I(nz) = r;  J(nz) = c;  V(nz) = v;
end

function v = getfielddef(s, f, d)
%GETFIELDDEF s.(f) if present and non-empty, else default d.
    if isfield(s, f) && ~isempty(s.(f))
        v = double(s.(f));
    else
        v = d;
    end
end
