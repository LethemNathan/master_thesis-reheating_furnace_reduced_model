function ang = build_angles(Nth, Nph)
%BUILD_ANGLES Build a simple 3D angular quadrature for DOM
%
% Inputs:
%   Nth : number of polar angles theta in [0, pi]
%   Nph : number of azimuthal angles phi in [0, 2*pi]
%
% Output:
%   ang.s           [M x 3]  direction unit vectors (sx, sy, sz)
%   ang.w           [M x 1]  solid-angle weights ; sum(ang.w) = 4*pi
%   ang.theta       [M x 1]  polar angle (cone midpoint)
%   ang.phi         [M x 1]  azimuthal angle (cone midpoint)
%   ang.thetaBounds [M x 2]  [theta_min, theta_max] of each cone
%   ang.phiBounds   [M x 2]  [phi_min, phi_max]   of each cone
%   ang.dth         scalar   theta step = pi / Nth
%   ang.dph         scalar   phi step   = 2*pi / Nph
%   ang.Nth, ang.Nph  scalars (the input discretization)
%   ang.M           scalar   total number of directions = Nth * Nph
%   ang.reflY       [M x 1]  reflection map about the y = 0 plane
%
% Notes:
%   - Midpoint quadrature in theta and phi.
%   - The bounds are kept so that downstream code can integrate
%     (s . n_i) analytically over each cone DeltaOmega^m (used by
%     precompute_face_direction_weights for the radiation FVM).

    dth = pi / Nth;
    dph = 2*pi / Nph;

    M = Nth * Nph;

    s           = zeros(M, 3);
    w           = zeros(M, 1);
    theta       = zeros(M, 1);
    phi         = zeros(M, 1);
    thetaBounds = zeros(M, 2);
    phiBounds   = zeros(M, 2);

    m = 0;
    for it = 1:Nth
        th     = (it - 0.5) * dth;        % midpoint in theta band
        thMin  = (it - 1)   * dth;
        thMax  =  it        * dth;

        sth = sin(th);
        cth = cos(th);

        for ip = 1:Nph
            ph     = (ip - 0.5) * dph;    % midpoint in phi band
            phMin  = (ip - 1)   * dph;
            phMax  =  ip        * dph;

            cph = cos(ph);
            sph = sin(ph);

            m = m + 1;

            s(m,1) = sth * cph;
            s(m,2) = sth * sph;
            s(m,3) = cth;

            w(m) = sth * dth * dph;

            theta(m)         = th;
            phi(m)           = ph;
            thetaBounds(m,:) = [thMin, thMax];
            phiBounds(m,:)   = [phMin, phMax];
        end
    end

    w = w * (4*pi / sum(w));    % normalize so sum(w) = 4*pi exactly

    reflY = build_reflection_map(s, [0 1 0]);

    ang             = struct();
    ang.s           = s;
    ang.w           = w;
    ang.theta       = theta;
    ang.phi         = phi;
    ang.thetaBounds = thetaBounds;
    ang.phiBounds   = phiBounds;
    ang.dth         = dth;
    ang.dph         = dph;
    ang.Nth         = Nth;
    ang.Nph         = Nph;
    ang.M           = M;
    ang.reflY       = reflY;
end


%% LOCAL FUNCTIONS

function reflIdx = build_reflection_map(s, normal)

    tol = 1e-10;

    n = normal(:).';
    nn = norm(n);
    if nn <= 0
        error('build_reflection_map: normal must be non-zero.');
    end
    n = n / nn;

    M = size(s,1);
    reflIdx = zeros(M,1);

    for m = 1:M
        dir  = s(m,:);
        sRef = dir - 2 * dot(dir, n) * n;

        err = vecnorm(s - sRef, 2, 2);
        [emin, j] = min(err);

        if emin > tol
            error(['build_reflection_map: no reflected direction found for m=%d. ' ...
                   'Check that the angular quadrature is symmetric. Min error = %.3e'], ...
                   m, emin);
        end

        reflIdx(m) = j;
    end
end
