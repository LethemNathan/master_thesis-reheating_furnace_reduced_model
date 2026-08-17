function [mesh, geom]=build_topology_slabs(mesh,geom,slabs)

    Nc = numel(mesh.cells.id);
    Nf = numel(geom.faceCategory);

    mesh.cells.isLoad  = false(Nc,1);
    mesh.cells.loadTag = zeros(Nc,1);

    geom.isLoadFace = false(Nf,1);
    geom.loadTag    = zeros(Nf,1);


    % Relabel mesh from marker volumes
    mesh = relabel_load_cells(mesh, geom, slabs);

    % Face Tag
    Nf=numel(geom.faceCategory);

    geom.isLoadFace = false(Nf,1);
    geom.loadTag = zeros(Nf,1);

    for f = 1:Nf
        ownerC = geom.owner(f);
        neighC = geom.neighbour(f);

        if neighC == 0
            if mesh.cells.isLoad(ownerC)
                isInterface = true;
                blockedCell = ownerC;
            else
                isInterface = false;
                blockedCell=0;
            end
        else
            ownerBlocked = mesh.cells.isLoad(ownerC);
            neighBlocked = mesh.cells.isLoad(neighC);
            if xor(ownerBlocked, neighBlocked)
                isInterface = true;
                if ownerBlocked
                    blockedCell = ownerC;
                else
                    blockedCell = neighC;
                end
            else
                isInterface = false;
                blockedCell=0;
            end
        end

        if isInterface
            geom.isLoadFace(f) = true;
            blockedPhysTag = mesh.cells.loadTag(blockedCell);
            newFaceTag = getInterfaceLoadTag(blockedPhysTag, ...
                                             geom.faceCenter(f,:), ...
                                             geom.faceNormal(f,:), ...
                                             slabs);
            geom.loadTag(f) = newFaceTag;
        end
    end
end



%% LOCAL FUNCTIONS

function mesh = relabel_load_cells(mesh, geom, slabs, tol)

    Nc = numel(mesh.cells.id);

    mesh.cells.isLoad = false(Nc,1);
    mesh.cells.loadTag = zeros(Nc,1);

    if nargin < 4
        tol = 1e-6;
    end

    for j = 1:numel(slabs)
        markerKeys = slabs(j).surfaces.physicalNames.keys;
        for i = 1:numel(markerKeys)
            k = markerKeys{i};
            if ~isKey(mesh.physicalNames, k)
                mesh.physicalNames(k) = slabs(j).surfaces.physicalNames(k);
            else       
                error('Keys already exist.');
            end
        end

        for c = 1:Nc
            xc = geom.cellCenter(c,:);
            bb = slabs(j).bbox;
    
            if pointInBBox(xc, bb, tol)
                mesh.cells.loadTag(c) = slabs(j).id;
                mesh.cells.isLoad(c)   = true;
            end
        end
    end
end


function inside = pointInBBox(x, bb, tol)
    inside = (x(1) >= bb(1)-tol) && (x(1) <= bb(4)+tol) && ...
             (x(2) >= bb(2)-tol) && (x(2) <= bb(5)+tol) && ...
             (x(3) >= bb(3)-tol) && (x(3) <= bb(6)+tol);
end


function physTag = getInterfaceLoadTag(blockedPhysTag, faceCenter, faceNormal, slabs)
    blockIdx = find([slabs.id] == blockedPhysTag, 1);
    if isempty(blockIdx)
        error('No slab with id %d found in slabs array.', blockedPhysTag);
    end
    blockSlabSurface = slabs(blockIdx).surfaces;

    bboxTol  = 0.089;
    axisTol  = 1e-6;
    degenTol = 1e-9;

    faceAxis = normalAxis(faceNormal, axisTol);
    if faceAxis == 0
        error(['Face normal (%.6g, %.6g, %.6g) is not aligned with any ', ...
               'coordinate axis (X, Y or Z). Non-axis-aligned interfaces ', ...
               'are not supported.'], ...
               faceNormal(1), faceNormal(2), faceNormal(3));
    end

    physTag = int32(0);
    found   = false;

    for i = 1:numel(blockSlabSurface.faceTag)
        sPhysTag = blockSlabSurface.faceTag(i);
        bbox     = blockSlabSurface.bbox{i};

        surfAxis = markerSurfaceAxis(bbox, degenTol);
        if surfAxis == 0
            error(['Slab surface (faceTag %d) is not aligned with a ', ...
                   'coordinate plane (no single degenerate bbox dimension).'], ...
                   sPhysTag);
        end

        if surfAxis ~= faceAxis
            continue;
        end

        inside = faceCenter(1) >= bbox(1) - bboxTol && faceCenter(1) <= bbox(4) + bboxTol && ...
                 faceCenter(2) >= bbox(2) - bboxTol && faceCenter(2) <= bbox(5) + bboxTol && ...
                 faceCenter(3) >= bbox(3) - bboxTol && faceCenter(3) <= bbox(6) + bboxTol;

        if inside
            if found
                error(['Ambiguous interface match: face centre ', ...
                       '(%.6g, %.6g, %.6g) is contained in several slab ', ...
                       'surface bboxes for slab id %d.'], ...
                       faceCenter(1), faceCenter(2), faceCenter(3), blockedPhysTag);
            end
            physTag = int32(sPhysTag);
            found   = true;
        end
    end

    if ~found
        error(['No slab interface surface matches face centre ', ...
               '(%.6g, %.6g, %.6g) with normal aligned with axis %d ', ...
               'for slab id %d.'], ...
               faceCenter(1), faceCenter(2), faceCenter(3), faceAxis, ...
               blockedPhysTag);
    end
end


function ax = normalAxis(n, tol)
    nn = n(:).' / norm(n);
    if     abs(abs(nn(1)) - 1) < tol, ax = 1;
    elseif abs(abs(nn(2)) - 1) < tol, ax = 2;
    elseif abs(abs(nn(3)) - 1) < tol, ax = 3;
    else,                              ax = 0;
    end
end

function ax = markerSurfaceAxis(bbox, tol)
    dx = abs(bbox(4) - bbox(1));
    dy = abs(bbox(5) - bbox(2));
    dz = abs(bbox(6) - bbox(3));
    if     dx < tol && dy >= tol && dz >= tol, ax = 1;
    elseif dy < tol && dx >= tol && dz >= tol, ax = 2;
    elseif dz < tol && dx >= tol && dy >= tol, ax = 3;
    else,                                       ax = 0;
    end
end