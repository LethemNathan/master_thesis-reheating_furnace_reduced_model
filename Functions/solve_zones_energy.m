function [T_sub, info] = solve_zones_energy(zEnCtx, gasProps, exchInfo, ...
        T_ez0, src, physicalConst, opts)
%SOLVE_ZONES_ENERGY Steady gas energy balance on the (merged) energy zones.
%
% Solves, for each energy zone I, the nonlinear steady balance (W):
%
%   advection_I(T) + Qburner_I + Qreactant_I + Qinlet_I - Qconv_I
%       + absorbed_I - 4*kappaV_I*sigma*T_I^4 - mdotOut_I*h(T_I) = 0
%
%   * Qreactant_I : sensible enthalpy carried in by the burner air/fuel mass
%     at its inlet temperature, mdot_burner*h(T_in). Q_burner being only the
%     (efficiency-adjusted) combustion release, the zone must additionally
%     heat this cool injected mass up to T_zone before it leaves with the
%     gas -> requires enforceMassConservation = true for consistency.
%
%   * advection_I : upwind enthalpy transport on the aggregated inter-EZ
%     flows,  sum over edges  (+/-) mdot * h(T_upwind) ;
%   * Qburner_I, Qinlet_I : imposed heat inputs [W] (burners, inlet hot gas) ;
%   * Qconv_I : imposed convective loss gas->walls [W] (from zEnCtx) ;
%   * absorbed_I : gas radiative absorption [W], LAGGED from the last
%     radiation solve (exchInfo.gasAbsorbed) ;
%   * 4*kappaV_I*sigma*T_I^4 : gas emission at the current (unknown) T ;
%   * mdotOut_I*h(T_I) : enthalpy leaving the domain at the upper inlet.
%
% Solved by Newton on the energy-zone temperatures (small system), with the
% analytic Jacobian (emission ~T^4, outflow ~h(T), upwind advection ~cp).
% Radiative absorption is frozen within the Newton solve and refreshed by
% the outer coupling loop.
%
% Mass balance. The per-EZ mass residual  B*mdot - s  (s = burner + inlet -
% outlet mass sources) is always computed and returned. If
% opts.enforceMassConservation is true, the inter-EZ flows are projected
% (minimum-adjustment least squares) onto the conservation constraint
% before the energy solve, so you can measure the impact of your CFD
% post-processing mass errors.
%
% Input vector orderings (default zone lists ; each overridable via
% src.<name>Zones). All boundary mass flows given as POSITIVE magnitudes:
%   src.Qburner  -> [111 211 311 312 322 332]  [W]
%   src.Qinlet   -> [121 131]                  [W]   (111 is an outflow, not an inflow)
%   src.mdotOut  -> [111 112 122 132]          [kg/s] (leaves the domain)
%   src.mdotBurnerLat  -> [111 211 311]        [kg/s] (mass bal. + reactant enthalpy)
%   src.mdotBurnerRoof -> [312 322 332]        [kg/s] (mass bal. + reactant enthalpy)
%   src.mdotInlet      -> [121 131]            [kg/s] (mass balance)
%   src.Tin_burnerLat  -> scalar [K] burner air/fuel inlet T (default 578.83)
%   src.Tin_burnerRoof -> scalar [K] roof burner inlet T     (default 578.83)
%
% Inputs
%   zEnCtx        : build_zone_energy_context output.
%   gasProps      : build_gas_properties output (.cp, .h handles).
%   exchInfo      : solve_radiation_exchange output (.gasAbsorbed, per gas
%                   element in radSurf order).
%   T_ez0         : (nEZ x 1) warm-start energy-zone temperatures [K].
%   src           : struct of source vectors (see orderings above).
%   physicalConst : with .sigma.
%   opts          : optional
%       .newtonTol      (default 1e-6, K)
%       .newtonMaxIter  (default 50)
%       .relax          (default 1.0, under-relaxation of the Newton step)
%       .Tmin, .Tmax    (default 250, 3000) clamp
%       .enforceMassConservation (default false)
%       .verbose        (default false)
%
% Outputs
%   T_sub : (nSub x 1) temperature per SUB-zone (energy-zone T broadcast to
%           its members), aligned with zEnCtx.subTags -> use to update
%           radFields.T on the fluid cells.
%   info  : struct with .T_ez, .converged, .iters, .dTmax,
%           .massResidual (nEZ), .maxAbsMassResidual, .flowMdotUsed,
%           .terms (per-EZ W breakdown: adv, burner, inlet, conv, radNet, out)

    sigma = physicalConst.sigma;

    if nargin < 7 || isempty(opts), opts = struct(); end
    tol      = getf(opts,'newtonTol',1e-6);
    maxIter  = getf(opts,'newtonMaxIter',50);
    relax    = getf(opts,'relax',1.0);
    Tmin     = getf(opts,'Tmin',250);
    Tmax     = getf(opts,'Tmax',3000);
    enforce  = getf(opts,'enforceMassConservation',false);
    verbose  = getf(opts,'verbose',false);

    nEZ = zEnCtx.nEZ;

    %% --- Map source vectors to energy zones ---
    % Zone lists are configurable via src.*Zones (defaults below). A zone may
    % appear both as an outflow (mdotOut) and elsewhere; slots are independent.
    zB   = getf(src,'QburnerZones',        [111 211 311 312 322 332]);
    zQi  = getf(src,'QinletZones',         [121 131]);
    zOut = getf(src,'mdotOutZones',        [111 112 122 132]);
    zmBL = getf(src,'mdotBurnerLatZones',  [111 211 311]);
    zmBR = getf(src,'mdotBurnerRoofZones', [312 322 332]);
    zmI  = getf(src,'mdotInletZones',      [121 131]);

    Qburner_ez = scatterEZ(zEnCtx, getf(src,'Qburner', zeros(numel(zB),1)),   zB);
    Qinlet_ez  = scatterEZ(zEnCtx, getf(src,'Qinlet',  zeros(numel(zQi),1)),  zQi);
    mdotOut_ez = scatterEZ(zEnCtx, getf(src,'mdotOut', zeros(numel(zOut),1)), zOut);

    mBurnLat_ez  = scatterEZ(zEnCtx, getf(src,'mdotBurnerLat',  zeros(numel(zmBL),1)), zmBL);
    mBurnRoof_ez = scatterEZ(zEnCtx, getf(src,'mdotBurnerRoof', zeros(numel(zmBR),1)), zmBR);
    mInlet_ez    = scatterEZ(zEnCtx, getf(src,'mdotInlet',      zeros(numel(zmI),1)),  zmI);

    %% --- Reactant sensible enthalpy carried IN by the burner mass ---
    % Q_burner is the (efficiency-adjusted) combustion heat release. The
    % injected air/fuel mass ALSO enters at its inlet temperature T_in and,
    % once heated to T_zone, leaves with the gas (mass conservation). Its
    % inlet sensible enthalpy  mdot_burner * h(T_in)  is therefore a real
    % energy inflow that must be added ; the zone then "pays" to heat this
    % cool mass up to T_zone. Requires enforceMassConservation = true so the
    % burner mass is actually advected out at T_zone.
    TinLat  = getf(src,'Tin_burnerLat',  578.83);   % K, editable
    TinRoof = getf(src,'Tin_burnerRoof', 578.83);   % K, editable
    mLat_v  = getf(src,'mdotBurnerLat',  zeros(numel(zmBL),1));
    mRoof_v = getf(src,'mdotBurnerRoof', zeros(numel(zmBR),1));
    qReacLat  = mLat_v(:)  * gasProps.hReac(TinLat);
    qReacRoof = mRoof_v(:) * gasProps.hReac(TinRoof);
    qReac_ez  = scatterEZ(zEnCtx, qReacLat, zmBL) + scatterEZ(zEnCtx, qReacRoof, zmBR);

    Qconv_ez = zEnCtx.Qconv_ez;

    %% --- Lagged gas radiative absorption per energy zone [W] ---
    absorbed_sub = exchInfo.gasAbsorbed(zEnCtx.gasElemOfSub);   % nSub
    absorbed_ez  = accumarray(zEnCtx.ezOfSub, absorbed_sub, [nEZ 1]);

    %% --- Mass balance (diagnostic ; optional correction) ---
    flowFrom = zEnCtx.flowFrom;  flowTo = zEnCtx.flowTo;  mdot = zEnCtx.flowMdot;
    nFlow    = numel(mdot);

    B = zeros(nEZ, nFlow);              % incidence: +1 leaves from, -1 enters to
    for e = 1:nFlow
        B(flowFrom(e), e) = B(flowFrom(e), e) + 1;
        B(flowTo(e),   e) = B(flowTo(e),   e) - 1;
    end
    s_mass = (mBurnLat_ez + mBurnRoof_ez + mInlet_ez) - mdotOut_ez;  % net injected

    massResidual = B*mdot - s_mass;     % should be ~0 if mass conserved

    if enforce && nFlow > 0
        % minimum-adjustment projection: mdot <- mdot - B'*pinv(B*B')*(B*mdot - s)
        mdot = mdot - B.' * (pinv(B*B.') * (B*mdot - s_mass));
        massResidual = B*mdot - s_mass;
    end

    %% --- Newton iteration on energy-zone temperatures ---
    T = T_ez0(:);
    cp = gasProps.cp;  h = gasProps.h;

    converged = false;  dTmax = NaN;  it = 0;
    for it = 1:maxIter
        R = zeros(nEZ,1);
        J = zeros(nEZ,nEZ);

        % Local terms: burner + inlet - conv + absorbed - emission - outflow
        emit = 4 * zEnCtx.kappaV_ez .* sigma .* T.^4;
        R = R + Qburner_ez + qReac_ez + Qinlet_ez - Qconv_ez + absorbed_ez - emit ...
              - mdotOut_ez .* h(T);
        for I = 1:nEZ
            J(I,I) = J(I,I) - 16*zEnCtx.kappaV_ez(I)*sigma*T(I)^3 ...
                            - mdotOut_ez(I)*cp(T(I));
        end

        % Advection (upwind = from-zone, since mdot >= 0)
        for e = 1:nFlow
            a = flowFrom(e);  b = flowTo(e);  m = mdot(e);
            ha = h(T(a));  cpa = cp(T(a));
            R(a) = R(a) - m*ha;      % leaves a
            R(b) = R(b) + m*ha;      % enters b
            J(a,a) = J(a,a) - m*cpa;
            J(b,a) = J(b,a) + m*cpa;
        end

        dT = -J\R;
        T  = T + relax*dT;
        T  = min(max(T, Tmin), Tmax);
        dTmax = max(abs(dT));

        if verbose
            fprintf('  [zones-energy] Newton %2d: ||dT||inf = %.3e K\n', it, dTmax);
        end
        if dTmax < tol
            converged = true;
            break
        end
    end

    if ~converged
        warning(['solve_zones_energy: Newton did not reach tol %.1e in %d ', ...
                 'iters (last ||dT||inf = %.3e K).'], tol, maxIter, dTmax);
    end

    %% --- Broadcast to sub-zones ---
    T_sub = T(zEnCtx.ezOfSub);

    %% --- Per-term breakdown (at converged T) for diagnostics ---
    adv = zeros(nEZ,1);
    for e = 1:nFlow
        a = flowFrom(e); b = flowTo(e); m = mdot(e); ha = h(T(a));
        adv(a) = adv(a) - m*ha;  adv(b) = adv(b) + m*ha;
    end
    terms = struct('adv',adv,'burner',Qburner_ez,'reactant',qReac_ez, ...
                   'inlet',Qinlet_ez, 'conv',-Qconv_ez, ...
                   'radNet', absorbed_ez - 4*zEnCtx.kappaV_ez.*sigma.*T.^4, ...
                   'out', -mdotOut_ez.*h(T));

    info = struct();
    info.T_ez               = T;
    info.converged          = converged;
    info.iters              = it;
    info.dTmax              = dTmax;
    info.massResidual       = massResidual;
    info.maxAbsMassResidual = max(abs(massResidual));
    info.flowMdotUsed       = mdot;
    info.terms              = terms;

    if verbose
        fprintf(['solve_zones_energy: %d EZ, T in [%.1f, %.1f] K | ', ...
                 'max|massRes|=%.3e kg/s (enforce=%d)\n'], ...
                nEZ, min(T), max(T), info.maxAbsMassResidual, enforce);
    end
end


%% LOCAL FUNCTIONS
function v_ez = scatterEZ(zEnCtx, values, zoneList)
    values = double(values(:));
    v_sub  = zeros(zEnCtx.nSub, 1);
    for k = 1:numel(zoneList)
        p = find(zEnCtx.subTags == zoneList(k), 1);
        if isempty(p)
            warning('solve_zones_energy: zone %d not found; value ignored.', zoneList(k));
            continue
        end
        v_sub(p) = v_sub(p) + values(k);
    end
    v_ez = accumarray(zEnCtx.ezOfSub, v_sub, [zEnCtx.nEZ 1]);
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
