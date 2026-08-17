function [radFields, exchInfo] = solve_radiation_exchange(S, radSurf, ...
        geom, mesh, radProps, physicalConst, radFields, opts)
%SOLVE_RADIATION_EXCHANGE Reduced radiation "solve" via exchange areas,
% with MIXED boundary conditions (imposed T and imposed flux).
%
% Replaces the full DOM (solve_radiation + compute_G_field) inside the outer
% coupling loop by a single small linear system + one matrix-vector product.
%
% Element classification (from the real radBcType carried by radProps):
%   * T-imposed elements  : gas zones (imposed T), wall surfaces (predicted
%     T from wall conduction), slab surfaces (predicted T). Their emissive
%     power E is KNOWN.
%   * flux-imposed elements : heaters / skids (radBcType 2 or 3). A NET
%     radiative flux is imposed over the WHOLE surface (global, not per
%     mesh face). Their E (temperature) is UNKNOWN.
%
% Black-surface energy balance (eps = 1), net flux leaving surface i into
% the enclosure:
%       Q_i = A_i E_i - sum_j S(i,j) E_j.
%   - T-imposed  : E_i known  -> Q_i is an output.
%   - flux-imposed: Q_i known -> solve for E_i from
%       (diag(A_F) - S_FF) E_F = Q_F + S_FK E_K.
%
% Symmetry is already embedded in S (extracted with radBcType 4 faces kept
% reflecting), so no special treatment here.
%
% After E is fully known, the incident (absorbed) power on every surface is
% sum_j S(i,j) E_j ; broadcast as a uniform flux density G_i = absorbed_i/A_i
% onto that surface's mesh faces to rebuild radFields.G for the wall and
% slab conduction solvers.
%
% Inputs
%   S             : nElem x nElem exchange areas.
%   radSurf       : element list (build_radiation_surfaces).
%   geom, mesh    : standard (faceArea, CellVolume).
%   radProps      : build_radiation_properties output at the CURRENT
%                   temperatures (real radBcType ; Ib_wall/Ib_cell carry the
%                   T-imposed emission ; flux-surface Ib is ignored here).
%   physicalConst : with .sigma.
%   radFields     : struct ; its .G field (Nf x 1) is (re)built and returned.
%   opts          : optional
%       .qImposedType3 : net flux density on radBcType 3 surfaces [W/m^2],
%                        positive = wall -> gas (default 0).
%       .qImposedType2 : net flux density on radBcType 2 surfaces [W/m^2]
%                        (default 0). >>> confirm convention vs
%                        update_wall_intensities.m <<<
%       .verbose       : default false.
%
% Outputs
%   radFields : with .G (Nf x 1) rebuilt (uniform per surface element ; 0 on
%               symmetry / internal faces).
%   exchInfo  : struct with
%       .E           (nElem x 1) emissive powers (solved)
%       .absorbed    (nElem x 1) incident absorbed power [W]
%       .Gsurf       (nSurf x 1) incident flux density per surface [W/m^2]
%       .isFlux      (nElem x 1 logical)
%       .Tflux       temperatures of the flux-imposed surfaces [K]
%       .fluxElem    global indices of flux-imposed surfaces
%       .gasAbsorbed (nGas x 1) absorbed power in each gas zone [W]
%                    (net gas radiative source available for later coupling)

    if nargin < 8 || isempty(opts), opts = struct(); end
    qT3     = getf(opts, 'qImposedType3', 0);
    qT2     = getf(opts, 'qImposedType2', 0);
    verbose = getf(opts, 'verbose', false);

    sigma      = physicalConst.sigma;
    faceArea   = double(geom.faceArea);
    cellVolume = double(geom.cellVolume);
    Nf         = numel(faceArea);

    elem  = radSurf.elem;
    nElem = radSurf.nElem;

    %% --- Element emissive powers (known) + flux classification ---
    E      = zeros(nElem, 1);
    isFlux = false(nElem, 1);
    Qimp   = zeros(nElem, 1);      % imposed NET flux (W), wall -> gas positive
    areaEl = zeros(nElem, 1);

    for i = 1:nElem
        e = elem(i);
        if strcmp(e.type, 'gas')
            c = double(e.cells);
            w = radProps.kappa_a(c) .* cellVolume(c);
            E(i) = sum(w .* (pi * radProps.Ib_cell(c))) / max(sum(w), eps);
        else
            f          = double(e.faces);
            areaEl(i)  = e.area;
            rt         = double(radProps.radBcType(f(1)));   % uniform per element
            if any(double(radProps.radBcType(f)) ~= rt)
                warning(['solve_radiation_exchange: element %d has mixed ', ...
                         'radBcType ; using the first face (%d).'], i, rt);
            end
            if rt == 3 || rt == 2
                isFlux(i) = true;
                qdens     = (rt == 3) * qT3 + (rt == 2) * qT2;
                Qimp(i)   = qdens * e.area;                 % global over surface
            else
                w    = faceArea(f);
                E(i) = sum(w .* (pi * radProps.Ib_wall(f))) / sum(w);
            end
        end
    end

    %% --- Solve the mixed-BC closure for the unknown E of flux surfaces ---
    F = find(isFlux);
    K = find(~isFlux);
    if ~isempty(F)
        S_FF = S(F, F);
        S_FK = S(F, K);
        Mff  = diag(areaEl(F)) - S_FF;
        rhs  = Qimp(F) + S_FK * E(K);
        E(F) = Mff \ rhs;
    end

    %% --- Incident absorbed power and per-surface flux density ---
    absorbed = S * E;                         % nElem x 1  [W]

    radFields.G = zeros(Nf, 1);
    sIdx  = radSurf.surfElem;
    Gsurf = zeros(numel(sIdx), 1);
    for kk = 1:numel(sIdx)
        i        = sIdx(kk);
        e        = elem(i);
        Gi       = absorbed(i) / max(e.area, eps);
        Gsurf(kk)= Gi;
        radFields.G(double(e.faces)) = Gi;    % uniform broadcast on the group
    end

    %% --- Pack diagnostics ---
    Tflux = (max(E(F), 0) / sigma) .^ 0.25;

    exchInfo             = struct();
    exchInfo.E           = E;
    exchInfo.absorbed    = absorbed;
    exchInfo.Gsurf       = Gsurf;
    exchInfo.isFlux      = isFlux;
    exchInfo.Tflux       = Tflux;
    exchInfo.fluxElem    = F;
    exchInfo.gasAbsorbed = absorbed(radSurf.gasElem);

    if verbose
        fprintf(['solve_radiation_exchange: %d flux surfaces solved ', ...
                 '(T in [%.0f, %.0f] K) | incident G in [%.2e, %.2e] W/m2\n'], ...
                numel(F), min([Tflux; inf]), max([Tflux; -inf]), ...
                min(Gsurf), max(Gsurf));
    end
end


%% LOCAL FUNCTIONS
function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
