function [G_face, G_cell] = compute_G_field(I, geom, mesh, radProps, D, ang)
%COMPUTE_G_FIELD Compute the incident radiation field on the FINE mesh,
% both per cell and per wall face.
%
% Cell field (W/m^2, useful as a diagnostic or as a volumetric source
% for an energy equation that includes a radiative source) :
%
%       G_cell(c) = sum_m I(c, m) * w_m
%
% Wall face field (W/m^2 ; this is the quantity consumed by
% solve_slabs_conduction via radFields.G, and by compute_qrad_per_zone
% to evaluate the net flux at each wall) :
%
%       G_face(f) = sum over m : sg*D(f, m) > 0    of
%                   sg * D(f, m) * I(fluidCell, m)
%
% where sg = +1 if the fluid cell adjacent to f is the owner, -1 if it
% is the neighbour. G_face is the INCOMING radiative flux on the wall
% (positive). It equals q_in in update_wall_intensities terminology.
%
% Non-wall faces and faces with radBcType == 4 (symmetry) are left at 0
% in G_face : the slab solver does not consume them, the per-zone qrad
% accumulator skips them, and they have no physical meaning as a wall
% incident flux.
%
% Inputs
%   I        : Nc x M, converged intensity field.
%   geom     : with owner, neighbour, faceArea.
%   mesh     : with cells.isLoad, cells.isBlocked.
%   radProps : with isWall, radBcType.
%   D        : Nf x M directional weights.
%   ang      : with .w (M x 1 solid-angle weights).
%
% Outputs
%   G_face : Nf x 1, in W/m^2.
%   G_cell : Nc x 1, in W/m^2.

    Nf = numel(radProps.isWall);

    %% --- Cell field : simple weighted sum over directions ---
    G_cell = I * ang.w(:);            % Nc x 1

    %% --- Wall-face field : same logic as update_wall_intensities ---
    G_face = zeros(Nf, 1);

    isLoadCell    = mesh.cells.isLoad;
    isBlockedCell = mesh.cells.isBlocked;

    for f = 1:Nf
        if ~radProps.isWall(f)
            continue
        end
        if radProps.radBcType(f) == 4
            continue   % symmetry : no physical incident wall flux
        end

        o = geom.owner(f);
        n = geom.neighbour(f);

        if n == 0
            if isLoadCell(o) || isBlockedCell(o)
                continue
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
            continue
        end

        D_oriented = sg * D(f, :);                 % 1 x M
        incoming   = D_oriented > 0;
        G_face(f)  = sum(D_oriented(incoming) .* I(fluidC, incoming));
    end
end
