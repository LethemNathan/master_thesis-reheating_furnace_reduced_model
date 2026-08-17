function c = loadTopologyCodes()
%CODES Numeric constants used across the conduction solver.
%
% Returns a struct with three groups of codes:
%   - Direction codes (1..6) used as columns in (Nc x 6) per-direction arrays
%   - BC type codes used in faceBC arrays
%   - Sub-zone codes (yy values from build_slab.m) used in faceSubzone arrays
%
% Convention:
%   U = x sup (rear)    D = x inf (front)
%   E = y sup (side)    W = y inf (symmetry)
%   N = z sup (top)     S = z inf (bottom)

    % Direction codes (column index in Nc x 6 arrays)
    c.DIR_E = 1;   % +y
    c.DIR_W = 2;   % -y
    c.DIR_N = 3;   % +z
    c.DIR_S = 4;   % -z
    c.DIR_U = 5;   % +x
    c.DIR_D = 6;   % -x

    % BC type codes
    c.BC_INTERNAL  = 0;
    c.BC_ADIABATIC = 1;
    c.BC_RADIATIVE = 2;
    c.BC_DIRICHLET = 3;

    % Sub-zone codes (must match yy in build_slab.m)
    c.SUBZONE_SYMMETRY = 1;
    c.SUBZONE_BOTTOM   = 2;
    c.SUBZONE_REAR     = 4;
    c.SUBZONE_SIDE     = 5;
    c.SUBZONE_FRONT    = 6;
    c.SUBZONE_TOP      = 8;
end
