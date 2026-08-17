function slabCtx = update_slab_dirichlet_bc(slabCtx, fields)
%UPDATE_SLAB_DIRICHLET_BC Refresh the rear-face Dirichlet temperatures of a
% slab from the global temperature field, before solving conduction in this
% slab.
%
% Behaviour depends on slabCtx.rearBC.type:
%   - "imposed"  : no-op. The first slab has its rear T fixed at init by
%                  bcConfig.slabs.firstSlabRearT.
%   - "coupled"  : reads fields.T at the global indices stored in
%                  slabCtx.rearBC.prevCellGlobal (which point to the front
%                  cells of the previous slab) and writes them into
%                  slabCtx.rearBC.T. This is the rear/front recycling:
%                  T_U(slab i) = T_D(slab i-1) at the cell level.
%
% Inputs
%   slabCtx : a slab context with a populated .rearBC sub-struct
%             (built by build_slab_rear_coupling)
%   fields  : struct with .T being a global cell-centred temperature vector
%             of length numel(mesh.cells.id)
%
% Output
%   slabCtx : same struct with slabCtx.rearBC.T updated (when "coupled").

    if ~isfield(slabCtx, 'rearBC')
        error(['update_slab_dirichlet_bc: slabCtx has no .rearBC field. ', ...
               'Did you forget to call build_slab_rear_coupling?']);
    end
    if ~isfield(fields, 'T') || isempty(fields.T)
        error('update_slab_dirichlet_bc: fields.T is missing or empty.');
    end

    rb = slabCtx.rearBC;

    switch rb.type
        case "imposed"
            % Nothing to do; T was set at init and stays.
            return;

        case "coupled"
            idx = rb.prevCellGlobal;
            if any(idx <= 0) || any(idx > numel(fields.T))
                error(['update_slab_dirichlet_bc: prevCellGlobal contains ', ...
                       'invalid indices (out of range of fields.T).']);
            end
            slabCtx.rearBC.T = fields.T(idx);

        otherwise
            error('update_slab_dirichlet_bc: unknown rearBC.type "%s".', ...
                   string(rb.type));
    end
end
