function [A, b] = assemble_slab_system(slabCtx, T_lin, T_n, fields, matProps, mode, dt)
%ASSEMBLE_SLAB_SYSTEM Build the linearised FVM system A*T = b for one slab.
%
% Cell-centred FVM with centred (CDS) gradient and harmonic mean for k at
% internal faces. Boundary conditions are dispatched by face from
% slabCtx.faceBC. The radiative BC is linearised around T_lin assuming
% T_surface ~ T_P (zero-gradient at the face).
%
% Inputs
%   slabCtx  : slab context (build_slab_context + build_slab_rear_coupling),
%              with .rearBC.T already refreshed by update_slab_dirichlet_bc
%              and .rearBCIdx populated.
%   T_lin    : Nc x 1 vector. Current Picard iterate (local). Used to
%              evaluate k(T), cp(T), rho(T) and to linearise sigma*T^4.
%   T_n      : Nc x 1 vector. Temperature at the previous time step (local).
%              Used only when mode == "unsteady"; pass [] for steady.
%   fields   : struct, must contain .G (Nf_glob x 1, irradiation per global
%              face, W/m^2). G is read at f_glob = slabCtx.cellFace(P, d)
%              for radiative faces.
%   matProps : .slabs.k, .slabs.cp, .slabs.rho   (poly coeffs ASCENDING)
%              .slabs.eps                         (scalar emissivity)
%   mode     : "steady" or "unsteady".
%   dt       : time step (used only when mode == "unsteady").
%
% Output
%   A : Nc x Nc sparse matrix. SPD as soon as at least one BC pins the
%       solution (Dirichlet, radiative or time term).
%   b : Nc x 1 right-hand side.

    c     = loadTopologyCodes();
    sigma = 5.670374419e-8;          

    Nc    = slabCtx.Nc;
    eps_s = matProps.slabs.eps;

    %% Material properties at T_lin (vectorised)
    k_cell   = eval_poly(matProps.slabs.k,   T_lin);     % Nc x 1
    cp_cell  = eval_poly(matProps.slabs.cp,  T_lin);     % Nc x 1
    rho_cell = eval_poly(matProps.slabs.rho, T_lin);     % Nc x 1

    %% Mode handling and basic checks
    isUnsteady = strcmpi(string(mode), "unsteady");
    if isUnsteady
        if isempty(T_n)
            error('assemble_slab_system: T_n is empty but mode == "unsteady".');
        end
        if numel(T_n) ~= Nc
            error(['assemble_slab_system: T_n has %d entries, expected %d ', ...
                   '(slab Nc).'], numel(T_n), Nc);
        end
        if ~(isnumeric(dt) && isscalar(dt) && dt > 0)
            error('assemble_slab_system: dt must be a positive scalar in unsteady mode.');
        end
    end
    if numel(T_lin) ~= Nc
        error(['assemble_slab_system: T_lin has %d entries, expected %d ', ...
               '(slab Nc).'], numel(T_lin), Nc);
    end

    %% Pre-allocate triplet storage
    % Worst case per cell: 1 diagonal + 6 off-diagonals + (1 diag from time
    % term) -> we budget 8*Nc to be safe. sparse() will sum duplicates.
    nMax = 8 * Nc;
    I   = zeros(nMax, 1);
    J   = zeros(nMax, 1);
    V   = zeros(nMax, 1);
    cnt = 0;
    b   = zeros(Nc, 1);

    %% Time term (unsteady only)
    if isUnsteady
        m = rho_cell .* cp_cell .* slabCtx.cellVolume / dt;     % Nc x 1
        idx = (cnt + 1):(cnt + Nc);
        I(idx) = (1:Nc).';
        J(idx) = (1:Nc).';
        V(idx) = m;
        cnt = cnt + Nc;
        b = b + m .* T_n;
    end

    %% Loop over cells x directions
    for P = 1:Nc
        T_P = T_lin(P);
        T3  = T_P * T_P * T_P;
        T4  = T3  * T_P;

        for d = 1:6
            bc  = slabCtx.faceBC(P, d);
            A_f = slabCtx.faceArea(P, d);

            switch bc
                case c.BC_INTERNAL
                    nb  = double(slabCtx.cellNeighbor(P, d));
                    d_f = slabCtx.faceDist(P, d);
                    k_f = harmonic_mean(k_cell(P), k_cell(nb));
                    a_f = k_f * A_f / d_f;

                    % Diagonal contribution: +a_f
                    cnt = cnt + 1;
                    I(cnt) = P;  J(cnt) = P;  V(cnt) = a_f;
                    % Off-diagonal: -a_f
                    cnt = cnt + 1;
                    I(cnt) = P;  J(cnt) = nb; V(cnt) = -a_f;

                case c.BC_ADIABATIC
                    % Zero flux; no contribution to A or b.

                case c.BC_DIRICHLET
                    kIdx = double(slabCtx.rearBCIdx(P, d));
                    if kIdx == 0
                        error(['assemble_slab_system: cell %d direction %d ', ...
                               'is BC_DIRICHLET but has no rearBCIdx entry. ', ...
                               'Did build_slab_rear_coupling run?'], P, d);
                    end
                    T_face = slabCtx.rearBC.T(kIdx);
                    d_Pf   = slabCtx.faceDist(P, d);
                    k_f    = k_cell(P);                  % cell-side conductivity
                    a_f    = k_f * A_f / d_Pf;

                    cnt = cnt + 1;
                    I(cnt) = P;  J(cnt) = P;  V(cnt) = a_f;
                    b(P) = b(P) + a_f * T_face;

                case c.BC_RADIATIVE
                    f_glob = double(slabCtx.cellFace(P, d));
                    if f_glob <= 0 || f_glob > numel(fields.G)
                        error(['assemble_slab_system: face global index %d ', ...
                               'out of range of fields.G (numel=%d).'], ...
                               f_glob, numel(fields.G));
                    end
                    G_f = fields.G(f_glob);

                    % Linearisation: sigma*T_P^4 ~ sigma*(4*T_lin^3*T_P - 3*T_lin^4)
                    a_rad = 4 * eps_s * sigma * A_f * T3;

                    cnt = cnt + 1;
                    I(cnt) = P;  J(cnt) = P;  V(cnt) = a_rad;
                    b(P) = b(P) + eps_s * A_f * (G_f + 3 * sigma * T4);

                otherwise
                    error('assemble_slab_system: unknown BC code %d at cell %d, dir %d.', ...
                           bc, P, d);
            end
        end
    end

    %% Trim and assemble
    I = I(1:cnt);
    J = J(1:cnt);
    V = V(1:cnt);
    A = sparse(I, J, V, Nc, Nc);
end
