function [fields, slabCtxs, info] = solve_slabs_conduction(slabCtxs, fields, matProps, mode, dt, solverParams)
%SOLVE_SLABS_CONDUCTION Single forward pass over all slabs in id ascending
% order.
%
% The inter-slab coupling is strictly uni-directional:
%   T_rear(slab i) = T_front(slab i-1)
% with slab 1 having an imposed rear T. No feedback exists from slab i+1
% back to slab i, so a single forward pass solves the coupling exactly.
%
% For each slab i = 1..Ns:
%   1. Refresh slab i's rear Dirichlet BC from fields.T
%      (no-op for i=1, reads previous slab's just-solved front for i>=2).
%   2. Solve slab i with solve_slab_conduction (Picard + direct Cholesky).
%   3. Write the slab cell solution into fields.T.
%   4. Write the slab radiative-face surface temperatures into fields.T_surf
%      (read by the radiation model to emit each slab face at its true
%      surface temperature T_s, not the first-cell temperature).
%
% Inputs
%   slabCtxs     : Ns x 1 struct array, sorted by id ascending, each with
%                  .rearBC and .rearBCIdx populated.
%   fields       : struct
%                    .T   global Nc_glob x 1, current temperature field
%                         (used as warm start; updated in place over slab
%                         cells).
%                    .G   global Nf_glob x 1, irradiation per face.
%                    .T_surf global Nf_glob x 1, surface temperature per
%                         radiative face (updated in place over slab faces;
%                         allocated here if missing).
%                    .T_n global Nc_glob x 1, T at the previous time step
%                         (REQUIRED if mode == "unsteady"). NEVER modified.
%   matProps     : material properties (see assemble_slab_system).
%   mode         : "steady" or "unsteady".
%   dt           : time step (used only when mode == "unsteady").
%   solverParams : struct, passed through to solve_slab_conduction.
%
% Outputs
%   fields    : same as input, with fields.T updated in slab cells only and
%               fields.T_surf updated on slab radiative faces only.
%               fields.T_n is NOT modified.
%   slabCtxs  : same as input, with slabCtxs(i).rearBC.T set to the values
%               used at the solve.
%   info      : diagnostic struct
%                 .allConverged   (bool, true iff every slab converged)
%                 .slabInfo       (Ns x 1 cell of info structs from
%                                  solve_slab_conduction)

    Ns = numel(slabCtxs);
    if Ns == 0
        error('solve_slabs_conduction: empty slabCtxs.');
    end

    %% Validate fields
    if ~isfield(fields, 'T') || isempty(fields.T)
        error('solve_slabs_conduction: fields.T is missing or empty.');
    end
    if ~isfield(fields, 'G') || isempty(fields.G)
        error('solve_slabs_conduction: fields.G is missing or empty.');
    end
    isUnsteady = strcmpi(string(mode), "unsteady");
    if isUnsteady && (~isfield(fields, 'T_n') || isempty(fields.T_n))
        error(['solve_slabs_conduction: fields.T_n is missing for ', ...
               'unsteady mode.']);
    end

    %% Ensure the surface-temperature field exists (one value per global face)
    if ~isfield(fields, 'T_surf') || isempty(fields.T_surf)
        fields.T_surf = nan(size(fields.G));
    end

    %% Forward pass over all slabs
    slabInfo = cell(Ns, 1);
    for i = 1:Ns
        % 1. Refresh rear Dirichlet BC from fields.T
        slabCtxs(i) = update_slab_dirichlet_bc(slabCtxs(i), fields);

        % 2. Solve slab i (Picard + direct Cholesky)
        [T_slab, info_i] = solve_slab_conduction(slabCtxs(i), fields, ...
                                                 matProps, mode, dt, ...
                                                 solverParams);

        % 3. Write back the cell temperatures: the next slab's rearBC update
        %    will read these freshly written values.
        fields.T(slabCtxs(i).globalCellIdx) = T_slab;

        % 4. Write back the radiative-face surface temperatures so the
        %    radiation model emits each slab face at its true surface
        %    temperature T_s (rather than the first-cell temperature).
        smap = info_i.surfaceMap;
        for s = 1:info_i.nSurf
            f_glob = double(slabCtxs(i).cellFace(smap.P(s), smap.d(s)));
            fields.T_surf(f_glob) = info_i.T_surface(s);
        end

        slabInfo{i} = info_i;
    end

    %% Aggregate diagnostics
    info = struct();
    info.slabInfo     = slabInfo;
    info.allConverged = all(cellfun(@(x) x.converged, slabInfo));
end
