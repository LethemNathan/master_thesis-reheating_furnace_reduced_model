

function mesh=build_topology_zones(mesh, geom, zonesDefinition)
    
    mesh.cells.zoneId = zeros(numel(mesh.cells.id), 1, 'int32');

    inFluid = (~mesh.cells.isBlocked) & (~mesh.cells.isLoad);
    globalCellIdx = find(inFluid);

    Nc = numel(globalCellIdx);

    for c = 1:Nc
        idx = globalCellIdx(c);
        xc = geom.cellCenter(idx,:);
        isFound = false;

        for i = 1:zonesDefinition.numZones
            bb = zonesDefinition.volumes(i).bbox;
            if pointInBBox(xc, bb)
                mesh.cells.zoneId(idx) = zonesDefinition.volumes(i).zoneTag;
                isFound = true;
                break
            end

        end

        if ~isFound
            error('Fluid cell is not assigned to a zone')
        end
    end
end


function inside = pointInBBox(x, bb, tol)

    if nargin < 3
        tol = 1e-6;
    end

    inside = (x(1) >= bb(1)-tol) && (x(1) <= bb(4)+tol) && ...
             (x(2) >= bb(2)-tol) && (x(2) <= bb(5)+tol) && ...
             (x(3) >= bb(3)-tol) && (x(3) <= bb(6)+tol);
end