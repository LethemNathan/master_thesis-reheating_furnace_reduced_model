function radFields = assign_zone_temperature_to_radFields(radFields, zoneFields, mesh, zonesCtx)
%ASSIGN_ZONE_TEMPERATURE_TO_RADFIELDS Project zone temperatures to mesh cells.
%
% Inputs
%   radFields    : struct containing radiation fields
%   zoneFields.T : Nz x 1 vector of zone temperatures
%   mesh         : mesh struct with mesh.cells.zoneId
%   zonesCtx     : zone context with zonesCtx.id
%
% Output
%   radFields.T  : Nc x 1 cell temperature field

    Nc = numel(mesh.cells.zoneId);


    for c = 1:Nc

        zoneTag = mesh.cells.zoneId(c);

        if zoneTag == 0
            continue
        end

        idx = find(zonesCtx.id == zoneTag, 1);

        if isempty(idx)
            error('Zone tag %d from mesh.cells.zoneId(%d) not found in zonesCtx.id.', ...
                  zoneTag, c);
        end

        radFields.T(c) = zoneFields.T(idx);
    end
end