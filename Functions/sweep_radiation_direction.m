function I_m = sweep_radiation_direction(geom, mesh, radProps, radCtx, ...
                                         D, ang, m, I_prev)
%SWEEP_RADIATION_DIRECTION Solve the FVM radiation equation for ONE
% direction m using the sweep method (no matrix assembly, no linear
% solve). Returns the intensity I^m on every cell.
%
% Cells are processed in upstream -> downstream order along direction
% s_m, determined by sorting (cellCenter . s_m) ascending. For each
% active cell P, the equation
%
%   I^m_P = ( sum_f max(-A_f * D_f, 0) * I^m_{f, upwind}
%             + kappa_a * Ib_cell(P) * dV * dOmega^m )
%         / ( sum_f max(A_f * D_f, 0) + kappa_a * dV * dOmega^m )
%
% is solved directly using already-swept upwind values. Inactive cells
% (slab / blocked) are forced to 0 (blocked-off).
%
% Wall faces dispatch on radProps.radBcType :
%   1 : std radiative   -> I_upwind = Ib_wall(f)
%   2 : adiabatic       -> I_upwind = Ib_wall(f) [pre-updated by
%                                                  update_wall_intensities]
%   3 : imposed q net   -> I_upwind = Ib_wall(f) [pre-updated]
%   4 : symmetry        -> I_upwind = I_prev(c, ang.reflY(m))
%
% Inputs
%   geom     : with .faceArea, .cellVolume, .cellCenter,
%              .owner, .neighbour.
%   mesh     : with .cells.isLoad, .cells.isBlocked, .cells.id.
%   radProps : with .kappa_a, .Ib_cell, .Ib_wall, .isWall, .radBcType.
%   radCtx   : with .cellFaces, .cellFaceSign (from build_rad_context).
%   D        : Nf x M directional weights (from
%              precompute_face_direction_weights).
%   ang      : output of build_angles, with .s (M x 3), .w (M x 1),
%              .reflY (M x 1).
%   m        : scalar direction index in 1..M.
%   I_prev   : Nc x M, intensities from the previous outer Picard
%              iteration (used for symmetry BCs). For the very first
%              iteration just pass zeros(Nc, M).
%
% Output
%   I_m      : Nc x 1, intensity field for direction m at this
%              iteration.

    Nc      = numel(mesh.cells.id);
    s_m     = ang.s(m, :);
    w_m     = ang.w(m);
    m_refl  = ang.reflY(m);

    isLoadCell    = mesh.cells.isLoad;
    isBlockedCell = mesh.cells.isBlocked;

    cellFaces    = radCtx.cellFaces;
    cellFaceSign = radCtx.cellFaceSign;

    faceArea  = geom.faceArea;
    cellVol   = geom.cellVolume;
    owner     = geom.owner;
    neighbour = geom.neighbour;

    kappa_a   = radProps.kappa_a;
    Ib_cell   = radProps.Ib_cell;
    Ib_wall   = radProps.Ib_wall;
    isWall    = radProps.isWall;
    radBcTy   = radProps.radBcType;

    %% 1. Sweep order : sort cells by (cellCenter . s_m) ascending
    proj = geom.cellCenter * s_m.';            % Nc x 1
    [~, sweepOrder] = sort(proj, 'ascend');    % upstream cells first

    %% 2. Cell-by-cell sweep
    I_m = zeros(Nc, 1);

    for idx = 1:Nc
        c = sweepOrder(idx);

        % --- Blocked-off : inactive cells forced to 0 ---
        if isLoadCell(c) || isBlockedCell(c)
            I_m(c) = 0;
            continue
        end

        % --- Active cell : assemble local equation ---
        a_P = kappa_a(c) * cellVol(c) * w_m;
        rhs = kappa_a(c) * Ib_cell(c) * cellVol(c) * w_m;

        for k = 1:6
            f = cellFaces(c, k);
            if f == 0
                continue   % defensive ; should not happen for hex cells
            end

            sg  = double(cellFaceSign(c, k));   % +1 owner, -1 neighbour
            D_f = sg * D(f, m);                 % signed for cell c : + out / - in
            A_f = faceArea(f);

            if D_f >= 0
                % --- Outgoing face : own value on the diagonal ---
                a_P = a_P + A_f * D_f;
            else
                % --- Incoming face : need I_upwind at face f ---
                if isWall(f)
                    bcType = radBcTy(f);
                    switch bcType
                        case {1, 2, 3}
                            I_up = Ib_wall(f);
                        case 4
                            % Symmetry plane : reflected direction
                            % evaluated at the SAME active cell c
                            I_up = I_prev(c, m_refl);
                        otherwise
                            error(['sweep_radiation_direction: unknown ', ...
                                   'radBcType %d at face %d.'], bcType, f);
                    end
                else
                    % Internal fluid | fluid face : upwind is the neighbour
                    if sg > 0
                        nc = neighbour(f);
                    else
                        nc = owner(f);
                    end
                    if nc == 0 || isLoadCell(nc) || isBlockedCell(nc)
                        % Should not happen on a non-wall face if
                        % radProps was built consistently. Fall back
                        % to 0 (treats it as a dark wall).
                        I_up = 0;
                    else
                        I_up = I_m(nc);   % upwind already swept
                    end
                end

                rhs = rhs + (-A_f * D_f) * I_up;
            end
        end

        if a_P <= 0
            error(['sweep_radiation_direction: non-positive a_P (%.3e) ', ...
                   'for cell %d, direction %d.'], a_P, c, m);
        end

        I_m(c) = rhs / a_P;
    end
end
