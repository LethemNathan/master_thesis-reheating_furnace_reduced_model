function radCtx = build_rad_context(geom, mesh)
%BUILD_RAD_CONTEXT Build the cell -> faces connectivity needed by the
% sweep-based radiation FVM solver.
%
% This is the "rad" analogue of slabCtx / zonesCtx : a one-time
% precomputation of per-cell topology, currently containing only the
% cell -> faces map (with an owner / neighbour sign per face entry),
% but designed to be extended later (per-direction sweep orders, etc.)
% without changing its public interface.
%
% Inputs
%   geom : with .owner, .neighbour (one entry per face, 1-based cell
%          indices ; 0 means "no cell" -> boundary face on that side).
%   mesh : with .cells.id (used only to get Nc).
%
% Output struct radCtx :
%   .cellFaces    (Nc x 6 int32)
%       Indices of the (exactly 6 for hex cells) faces touching each
%       cell. Order within a row matches the order in which the faces
%       appear in the geom face list ; it has no semantic meaning
%       beyond "these are the faces of this cell".
%
%   .cellFaceSign (Nc x 6 int8)
%       +1 if the cell is the OWNER of the face (face normal points
%       outward from the cell), -1 if the cell is the NEIGHBOUR
%       (face normal points into the cell). Used by the sweep to know
%       the outward direction of each face from the cell's viewpoint :
%
%           D_outward(c, f, m) = cellFaceSign(c, k) * D(f, m)
%
%       where k is the local index of f in cellFaces(c, :).

    Nc = numel(mesh.cells.id);
    Nf = numel(geom.owner);

    MAX_FACES_PER_CELL = 6;     % hex cells in build_topology

    cellFaces    = zeros(Nc, MAX_FACES_PER_CELL, 'int32');
    cellFaceSign = zeros(Nc, MAX_FACES_PER_CELL, 'int8');
    nFacesInCell = zeros(Nc, 1, 'int32');

    for f = 1:Nf
        o = geom.owner(f);
        if o > 0
            nFacesInCell(o) = nFacesInCell(o) + 1;
            k = nFacesInCell(o);
            if k > MAX_FACES_PER_CELL
                error(['build_rad_context: cell %d has more than %d ', ...
                       'faces (got %d). Mesh is not pure hex ?'], ...
                      o, MAX_FACES_PER_CELL, k);
            end
            cellFaces(o, k)    = int32(f);
            cellFaceSign(o, k) = int8(+1);
        end

        n = geom.neighbour(f);
        if n > 0
            nFacesInCell(n) = nFacesInCell(n) + 1;
            k = nFacesInCell(n);
            if k > MAX_FACES_PER_CELL
                error(['build_rad_context: cell %d has more than %d ', ...
                       'faces (got %d). Mesh is not pure hex ?'], ...
                      n, MAX_FACES_PER_CELL, k);
            end
            cellFaces(n, k)    = int32(f);
            cellFaceSign(n, k) = int8(-1);
        end
    end

    %% Sanity : every hex cell should have exactly 6 faces in the list.
    badCells = find(nFacesInCell ~= MAX_FACES_PER_CELL);
    if ~isempty(badCells)
        nShow = min(10, numel(badCells));
        msg = sprintf(['build_rad_context: %d cell(s) do not have ', ...
                       'exactly %d faces. First %d (id, faceCount) :\n'], ...
                      numel(badCells), MAX_FACES_PER_CELL, nShow);
        for k = 1:nShow
            c = badCells(k);
            msg = [msg, sprintf('    %d : %d\n', c, nFacesInCell(c))];  %#ok<AGROW>
        end
        error('%s', msg);
    end

    radCtx              = struct();
    radCtx.cellFaces    = cellFaces;
    radCtx.cellFaceSign = cellFaceSign;
end
