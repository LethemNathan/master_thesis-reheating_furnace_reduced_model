function zEnCtx = build_zone_energy_context(zonesCtx, radSurf, geom, ...
        flowData, convHeatFluxMap, kappa, opts)
%BUILD_ZONE_ENERGY_CONTEXT Static context for the gas-zone energy balance.
%
% Prepares everything time-INDEPENDENT needed by solve_zones_energy:
%   * the (possibly merged) "energy zones" and the sub-zone <-> energy-zone
%     maps (hybrid model) ;
%   * per energy-zone volume and kappa*V (for the 4*kappa*V*sigma*T^4 gas
%     emission) ;
%   * the inter-energy-zone mass-flow connectivity, aggregated from the
%     internal-surface report (flows internal to a merged group cancel ;
%     flows across a group boundary are SUMMED) ;
%   * the boundary mass flows (report rows touching a domain boundary) ;
%   * the imposed convective sink gas->walls per energy-zone (W) ;
%   * the mapping sub-zone -> gas element index (into S / exchInfo.gasAbsorbed).
%
% Energy zones (hybrid). By default every sub-zone is its own energy zone.
% opts.mergeGroups lets you MERGE sub-zones that must share a single
% temperature in the energy equation only (radiation stays per sub-zone):
%   opts.mergeGroups = { [112 122 132], [212 222 232], ... }
% Each cell is a set of sub-zone tags fused into one energy zone.
%
% Mass-flow mapping. flowData is an N x 2 array [faceID, mdot] where faceID
% is the ROW index r into zonesCtx.zoneFace and mdot is the (already
% sign-corrected) mass flow in kg/s. zoneFace(r,:) = [c1 c2] are positions
% into zonesCtx.id, and mdot > 0 means flow from id(c1) to id(c2) (i.e.
% toward -x/-y/+z). A zero in c1 or c2 marks a domain-boundary face.
%
% Convective sink. For every wall surface element of radSurf (source
% 'faceTag'), the imposed convective flux gas->wall is
% convHeatFluxMap(zoneTag, faceTag) [W/m^2] times its area ; summed per
% energy zone. Slab (load) surfaces are excluded (no gas<->slab convection).
%
% Inputs
%   zonesCtx        : .id (Nz x 1 tags), .zoneFace (Nf_z x 2 positions).
%   radSurf         : build_radiation_surfaces output (.elem, .gasOfZone,
%                     .gasElem, .zoneTags).
%   geom            : (kept for interface area if needed later ; not required).
%   flowData        : N x 2 [faceID(row of zoneFace), mdot(kg/s)].
%   convHeatFluxMap : sparse/num matrix indexed (zoneTag, faceTag), W/m^2,
%                     positive = into wall (out of gas). Missing -> 0.
%   kappa           : gas absorption coefficient (scalar, matProps.fluid.kappa).
%   opts            : optional
%       .mergeGroups : cell array of sub-zone tag vectors to merge. Default {}.
%
% Output struct zEnCtx
%   .subTags     (nSub x 1) sub-zone tags (== zonesCtx.id)
%   .nSub
%   .nEZ
%   .subToEZ     containers.Map(subTag -> EZ index)
%   .ezOfSub     (nSub x 1) EZ index per sub-zone (aligned with subTags)
%   .ezMembers   (nEZ x 1 cell) sub-zone tags in each EZ
%   .ezMemberPos (nEZ x 1 cell) positions (into subTags) of members
%   .V_sub, .V_ez        volumes [m^3]
%   .kappaV_sub, .kappaV_ez   kappa*V [m^2]
%   .gasElemOfSub (nSub x 1) gas element index (into exchInfo.gasAbsorbed)
%   .Qconv_ez    (nEZ x 1) imposed convective sink gas->walls [W]
%   .flowFrom, .flowTo (nFlow x 1) EZ indices ; .flowMdot (nFlow x 1) kg/s
%                (aggregated, positive = From -> To)
%   .bndZoneEZ, .bndMdot   boundary-face flows (EZ index, kg/s) if any
%
% NOTE: burner / inlet energy inputs and domain-outflow rates are NOT here ;
% they are passed at solve time (they may change / be tuned).

    if nargin < 7 || isempty(opts), opts = struct(); end
    mergeGroups = getf(opts, 'mergeGroups', {});

    % Accept flowData as a table (readInternalSurfaceReport surfaceReport)
    % or a numeric array. Expect columns [faceID, mdot].
    if istable(flowData)
        flowData = table2array(flowData);
    end

    subTags = double(zonesCtx.id(:));
    nSub    = numel(subTags);
    tagPos  = containers.Map(num2cell(subTags), num2cell(1:nSub));  % tag->pos

    %% --- 1) Energy-zone assignment (merge map) ---
    ezOfSub = (1:nSub).';        % start: each sub-zone its own EZ
    for g = 1:numel(mergeGroups)
        grp = double(mergeGroups{g});
        pos = zeros(numel(grp),1);
        for k = 1:numel(grp)
            if ~isKey(tagPos, grp(k))
                error('build_zone_energy_context: merge tag %d is not a zone.', grp(k));
            end
            pos(k) = tagPos(grp(k));
        end
        ezOfSub(pos) = min(ezOfSub(pos));   % collapse to a common label
    end
    % Re-index EZ labels to a contiguous 1..nEZ
    [uEZ, ~, newIdx] = unique(ezOfSub, 'stable');
    ezOfSub = newIdx;
    nEZ     = numel(uEZ);

    subToEZ = containers.Map(num2cell(subTags), num2cell(ezOfSub));

    ezMembers   = cell(nEZ,1);
    ezMemberPos = cell(nEZ,1);
    for e = 1:nEZ
        ezMemberPos{e} = find(ezOfSub == e);
        ezMembers{e}   = subTags(ezMemberPos{e});
    end

    %% --- 2) Volumes, kappa*V, gas element mapping ---
    V_sub        = zeros(nSub,1);
    gasElemOfSub = zeros(nSub,1);
    for p = 1:nSub
        ztag            = subTags(p);
        ge              = radSurf.gasOfZone(ztag);   % global elem idx (== gas order)
        gasElemOfSub(p) = ge;
        V_sub(p)        = radSurf.elem(ge).volume;
    end
    kappaV_sub = kappa * V_sub;

    V_ez      = accumarray(ezOfSub, V_sub,      [nEZ 1]);
    kappaV_ez = accumarray(ezOfSub, kappaV_sub, [nEZ 1]);

    %% --- 3) Imposed convective sink gas->walls, per energy zone ---
    Qconv_sub = zeros(nSub,1);
    for i = 1:radSurf.nElem
        e = radSurf.elem(i);
        if ~strcmp(e.type,'surf') || ~strcmp(e.source,'faceTag')
            continue          % only wall/heater/skid surfaces ; skip gas & slabs
        end
        z  = e.zoneTag;
        ft = e.tag;
        q  = lookup_flux(convHeatFluxMap, z, ft);   % W/m^2, into wall
        if q == 0, continue; end
        if isKey(tagPos, z)
            Qconv_sub(tagPos(z)) = Qconv_sub(tagPos(z)) + q * e.area;
        end
    end
    Qconv_ez = accumarray(ezOfSub, Qconv_sub, [nEZ 1]);

    %% --- 4) Mass-flow connectivity (aggregate per EZ pair) ---
    zoneFace = zonesCtx.zoneFace;
    nZF      = size(zoneFace,1);

    faceID = double(flowData(:,1));
    mdot   = double(flowData(:,2));

    % Aggregate internal EZ-pair flows into a map keyed (minEZ,maxEZ),
    % storing the net signed flow oriented low->high EZ index.
    pairKey = containers.Map('KeyType','char','ValueType','double');

    bndZoneEZ = [];
    bndMdot   = [];

    for k = 1:numel(faceID)
        r = faceID(k);
        if r < 1 || r > nZF
            warning('build_zone_energy_context: faceID %d out of zoneFace range; skipped.', r);
            continue
        end
        c1 = zoneFace(r,1);
        c2 = zoneFace(r,2);
        m  = mdot(k);                 % >0 : id(c1) -> id(c2)

        if c1 > 0 && c2 > 0
            eA = ezOfSub(c1);         % from-side (positive flow)
            eB = ezOfSub(c2);         % to-side
            if eA == eB
                continue              % internal to a merged group : cancels
            end
            lo = min(eA,eB); hi = max(eA,eB);
            s  = m; if eA > eB, s = -m; end   % orient low->high
            key = sprintf('%d_%d', lo, hi);
            if isKey(pairKey,key), pairKey(key) = pairKey(key) + s;
            else,                  pairKey(key) = s; end
        else
            % boundary face : keep the nonzero-zone EZ and the flow.
            if c1 > 0, ez = ezOfSub(c1); else, ez = ezOfSub(c2); end
            bndZoneEZ(end+1,1) = ez;   %#ok<AGROW>
            bndMdot(end+1,1)   = m;    %#ok<AGROW>
        end
    end

    keys_   = pairKey.keys;
    nFlow   = numel(keys_);
    flowFrom = zeros(nFlow,1); flowTo = zeros(nFlow,1); flowMdot = zeros(nFlow,1);
    for k = 1:nFlow
        parts = sscanf(keys_{k}, '%d_%d');
        lo = parts(1); hi = parts(2);
        m  = pairKey(keys_{k});
        if m >= 0
            flowFrom(k) = lo; flowTo(k) = hi; flowMdot(k) = m;
        else
            flowFrom(k) = hi; flowTo(k) = lo; flowMdot(k) = -m;
        end
    end

    %% --- Pack ---
    zEnCtx = struct();
    zEnCtx.subTags     = subTags;
    zEnCtx.nSub        = nSub;
    zEnCtx.nEZ         = nEZ;
    zEnCtx.subToEZ     = subToEZ;
    zEnCtx.ezOfSub     = ezOfSub;
    zEnCtx.ezMembers   = ezMembers;
    zEnCtx.ezMemberPos = ezMemberPos;
    zEnCtx.V_sub       = V_sub;
    zEnCtx.V_ez        = V_ez;
    zEnCtx.kappaV_sub  = kappaV_sub;
    zEnCtx.kappaV_ez   = kappaV_ez;
    zEnCtx.gasElemOfSub = gasElemOfSub;
    zEnCtx.Qconv_ez    = Qconv_ez;
    zEnCtx.flowFrom    = flowFrom;
    zEnCtx.flowTo      = flowTo;
    zEnCtx.flowMdot    = flowMdot;
    zEnCtx.bndZoneEZ   = bndZoneEZ;
    zEnCtx.bndMdot     = bndMdot;

    fprintf(['build_zone_energy_context: %d sub-zones -> %d energy zones | ', ...
             '%d inter-EZ flows | %d boundary flows | sum|Qconv|=%.3g W\n'], ...
            nSub, nEZ, nFlow, numel(bndMdot), sum(abs(Qconv_ez)));
end


%% LOCAL FUNCTIONS
function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function v = lookup_flux(src, zoneTag, faceTag)
    v = 0;
    if isempty(src), return; end
    if isa(src,'containers.Map')
        key = sprintf('%d_%d', zoneTag, faceTag);
        if isKey(src,key), v = full(double(src(key))); end
        return
    end
    [nr,nc] = size(src);
    if zoneTag>=1 && zoneTag<=nr && faceTag>=1 && faceTag<=nc
        v = full(double(src(zoneTag, faceTag)));
    end
end
