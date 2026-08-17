function radProps = update_wall_intensities(geom, mesh, radProps, D, ang, ...
                                              I, q_imposed, physicalConst)
%UPDATE_WALL_INTENSITIES Recompute radProps.Ib_wall (and the diagnostic
% radProps.Tw_face) on every wall face with radBcType in {2, 3}, using
% the current cell intensity field I. Wall faces with radBcType in
% {1, 4} are NOT touched (1 already has its imposed Tw, 4 = symmetry is
% handled directly in the sweep via the reflected direction).
%
% Physics
% -------
% For a diffuse-gray wall with epsilon_w = 1, the radiative balance is
%
%       q_net_{wall -> gas}  =  pi * Ib_w  -  q_in
%
% where q_in is the incoming radiative flux on the wall (positive,
% W/m^2) :
%
%       q_in(f) = sum over m : D_oriented(f, m) > 0   of
%                 D_oriented(f, m) * I(fluidCell, m)
%
% with D_oriented(f, m) = sg * D(f, m), sg = +1 if the active fluid
% cell is the owner of f, -1 if it is the neighbour.
%
%   radBcType == 2  (adiabatic)   :  q_net = 0
%                                    -> Ib_w = q_in / pi
%   radBcType == 3  (imposed)     :  q_net = q_imposed   (W/m^2,
%                                    positive = wall emits more
%                                    than it absorbs, i.e. flux
%                                    entering the gas domain)
%                                    -> Ib_w = (q_imposed + q_in) / pi
%
% If Ib_w comes out negative (unphysical, can happen at the very first
% Picard iterations with cold q_in or with an aggressive q_imposed), we
% clamp it to 0 and emit a single aggregated warning at the end.
%
% Inputs
%   geom, mesh    : standard.
%   radProps      : current radiation properties.
%   D             : Nf x M directional weights from
%                   precompute_face_direction_weights.
%   ang           : with .M.
%   I             : Nc x M, current intensity field (typically the one
%                   coming out of the previous outer Picard iteration).
%   q_imposed     : scalar, W/m^2, value for radBcType == 3 faces.
%                   Positive = wall -> gas (entering the domain).
%   physicalConst : with .sigma (Stefan-Boltzmann).
%
% Output
%   radProps : same struct, with .Ib_wall and .Tw_face refreshed on
%              wall faces of type 2 and 3.

    sigma         = physicalConst.sigma;
    isLoadCell    = mesh.cells.isLoad;
    isBlockedCell = mesh.cells.isBlocked;
    owner         = geom.owner;
    neighbour     = geom.neighbour;

    Nf     = numel(radProps.isWall);
    nNegIb = 0;

    for f = 1:Nf
        if ~radProps.isWall(f)
            continue
        end
        bcType = radProps.radBcType(f);
        if bcType ~= 2 && bcType ~= 3
            continue
        end

        % --- Identify the fluid cell adjacent to f, with its sign ---
        o = owner(f);
        n = neighbour(f);
        if n == 0
            if isLoadCell(o) || isBlockedCell(o)
                continue           % wall not adjacent to a fluid cell
            end
            fluidC = o;
            sg     = +1;
        elseif ~isLoadCell(o) && ~isBlockedCell(o)
            fluidC = o;
            sg     = +1;
        elseif ~isLoadCell(n) && ~isBlockedCell(n)
            fluidC = n;
            sg     = -1;
        else
            continue               % both sides inactive : skip
        end

        % --- Incoming flux at the wall (W/m^2) ---
        D_oriented = sg * D(f, :);                 % 1 x M, + outgoing / - incoming
        incoming   = D_oriented > 0;               % directions that hit the wall
        q_in       = sum(D_oriented(incoming) .* I(fluidC, incoming));

        % --- Wall blackbody intensity from the flux balance ---
        if bcType == 2
            Ib_w = q_in / pi;
        else   % bcType == 3
            Ib_w = (q_imposed + q_in) / pi;
        end

        if Ib_w < 0
            nNegIb = nNegIb + 1;
            Ib_w   = 0;            % clamp for stability
        end

        radProps.Ib_wall(f) = Ib_w;
        radProps.Tw_face(f) = (Ib_w * pi / sigma)^(1/4);
    end

    if nNegIb > 0
        warning(['update_wall_intensities: %d wall face(s) produced a ', ...
                 'negative Ib_w and were clamped to 0. Check the value ', ...
                 'of q_imposed vs the local q_in.'], nNegIb);
    end
end
