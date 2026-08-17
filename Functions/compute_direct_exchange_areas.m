function [S, info] = compute_direct_exchange_areas(radSurf, geom, mesh, ...
        radPropsBase, radCtx, D, ang, physicalConst, opts)
%COMPUTE_DIRECT_EXCHANGE_AREAS Direct exchange areas by influence coefficients.
%
% Exploits the LINEARITY of the gray, non-scattering DOM (eps = 1 walls,
% fixed kappa): the map {emissive powers} -> {absorbed powers} is linear, so
% emitting one element at a time with unit emissive power and reading what
% every element absorbs directly gives the direct-exchange-area matrix.
%
% For each source element j:
%   * build a radProps where ONLY j emits (unit emissive power), all walls
%     are forced BLACK & COLD absorbers (radBcType = 1, Ib_wall = 0),
%     symmetry faces (radBcType = 4) kept as symmetry, gas cold (Ib_cell = 0)
%     except the source zone ;
%   * solve the transport once with the existing solve_radiation ;
%   * integrate the absorbed power on every element i:
%       surface i : sum_f  G_face(f) * A_f          (black -> absorbed = incident)
%       gas i     : sum_c  kappa(c) * G_cell(c) * V_c
%   * store it in column j.
%
% Unit emissive power convention: a source emits with E = 1 W/m^2, i.e.
% Ib = E/pi = 1/pi on its faces (surface) or cells (gas). The absorbed power
% at i is then numerically equal to the direct exchange area between i and j
% (m^2): S(i,j) = s_i s_j / g_i s_j / g_i g_j depending on the types.
%
% Element ordering and membership come entirely from radSurf
% (build_radiation_surfaces): gas elements first, then surfaces.
%
% Inputs
%   radSurf       : from build_radiation_surfaces.
%   geom, mesh    : standard.
%   radPropsBase  : a build_radiation_properties output (any current one).
%                   Used only for its geometry-consistent classification
%                   (isWall, radBcType, eps_wall, kappa_a). Emissions are
%                   overwritten here.
%   radCtx        : from build_rad_context.
%   D, ang        : directional weights / angular quadrature.
%   physicalConst : with .sigma.
%   opts          : optional struct
%       .solverParams : passed to solve_radiation (default verbose off,
%                       picardMaxIter 10).
%       .qImposed     : scalar for radBcType 3 (irrelevant here, default 0).
%       .symmetrize   : average S with S' to enforce reciprocity (default true).
%       .verbose      : progress printout (default true).
%
% Outputs
%   S    : nElem x nElem direct exchange areas [m^2]. S(i,j) = exchange area
%          between elements i and j (symmetric if opts.symmetrize).
%   info : struct with
%          .emitted      (nElem x 1) theoretical emission per element
%                        (A_j for surfaces, 4*sum(kappa*V) for gas zones)
%          .colSum       (nElem x 1) sum_i S_raw(i,j)  (before symmetrise)
%          .consRatio    (nElem x 1) colSum ./ emitted  (should be ~1)
%          .maxConsErr   max |consRatio - 1|
%          .reciprocityErr  max|S_raw - S_raw'| / max(S_raw)  (before sym)
%          .elapsed      seconds.

    if nargin < 9 || isempty(opts), opts = struct(); end
    sp         = getf(opts, 'solverParams', struct('verbose', false, ...
                                                   'picardMaxIter', 10));
    if ~isfield(sp, 'verbose'),       sp.verbose = false;      end
    if ~isfield(sp, 'picardMaxIter'), sp.picardMaxIter = 10;   end
    qImposed   = getf(opts, 'qImposed',   0);
    doSym      = getf(opts, 'symmetrize', true);
    verbose    = getf(opts, 'verbose',    true);
    printEvery = getf(opts, 'printEvery', 5);   % progress line every N sources

    Nc = numel(mesh.cells.id);
    Nf = numel(geom.faceArea);

    faceArea   = double(geom.faceArea);
    cellVolume = double(geom.cellVolume);

    nElem = radSurf.nElem;
    elem  = radSurf.elem;

    %% --- Extraction template: black cold absorbers, symmetry preserved ---
    rp0 = radPropsBase;
    rp0.sigma_s = 0;
    wallMask = logical(rp0.isWall);
    symMask  = wallMask & (rp0.radBcType == 4);
    rp0.radBcType(wallMask & ~symMask) = int32(1);   % force black absorbing
    rp0.eps_wall(wallMask) = 1;
    rp0.Ib_cell = zeros(Nc, 1);
    rp0.Ib_wall = zeros(Nf, 1);
    rp0.Tw_face = zeros(Nf, 1);

    invPi = 1 / pi;

    %% --- Theoretical emission per element (for conservation check) ---
    emitted = zeros(nElem, 1);
    for i = 1:nElem
        e = elem(i);
        if strcmp(e.type, 'surf')
            emitted(i) = e.area;                                % A_i
        else
            c = double(e.cells);
            emitted(i) = 4 * sum(rp0.kappa_a(c) .* cellVolume(c)); % 4*sum(kV)
        end
    end

    %% --- Loop over sources ---
    S = zeros(nElem, nElem);
    tStart = tic;

    for j = 1:nElem
        e  = elem(j);
        rp = rp0;
        if strcmp(e.type, 'gas')
            rp.Ib_cell(double(e.cells)) = invPi;      % E_g = 1
        else
            rp.Ib_wall(double(e.faces)) = invPi;      % E   = 1
        end

        I = solve_radiation(geom, mesh, rp, radCtx, D, ang, qImposed, ...
                            physicalConst, sp);
        [Gf, Gc] = compute_G_field(I, geom, mesh, rp, D, ang);

        % Absorbed power on every element -> column j
        col = zeros(nElem, 1);
        for i = 1:nElem
            ei = elem(i);
            if strcmp(ei.type, 'surf')
                f      = double(ei.faces);
                col(i) = sum(Gf(f) .* faceArea(f));
            else
                c      = double(ei.cells);
                col(i) = sum(rp0.kappa_a(c) .* Gc(c) .* cellVolume(c));
            end
        end
        S(:, j) = col;

        if verbose && (mod(j, printEvery) == 0 || j == nElem || j == 1)
            elapsedNow = toc(tStart);
            rate       = elapsedNow / j;              % s per source
            etaSec     = rate * (nElem - j);          % remaining
            fprintf(['  exch. areas: %4d / %-4d (%5.1f%%) | %-4s z%-4g | ', ...
                     'elapsed %6.1fs | ETA %6.1fs (%.2fs/src)\n'], ...
                    j, nElem, 100*j/nElem, e.type, e.zoneTag, ...
                    elapsedNow, etaSec, rate);
        end
    end

    elapsed = toc(tStart);

    %% --- Diagnostics (on the raw, pre-symmetrised matrix) ---
    colSum    = sum(S, 1).';
    consRatio = colSum ./ max(emitted, eps);
    reciErr   = max(abs(S - S.'), [], 'all') / max(max(S(:)), eps);

    if doSym
        S = 0.5 * (S + S.');
    end

    info = struct();
    info.emitted        = emitted;
    info.colSum         = colSum;
    info.consRatio      = consRatio;
    info.maxConsErr     = max(abs(consRatio - 1));
    info.reciprocityErr = reciErr;
    info.elapsed        = elapsed;

    if verbose
        fprintf(['compute_direct_exchange_areas: %d sources in %.1f s | ', ...
                 'max conservation error = %.2e | reciprocity error = %.2e\n'], ...
                nElem, elapsed, info.maxConsErr, info.reciprocityErr);
    end
end


%% LOCAL FUNCTIONS
function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
