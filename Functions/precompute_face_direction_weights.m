function D = precompute_face_direction_weights(geom, ang)
%PRECOMPUTE_FACE_DIRECTION_WEIGHTS Pre-compute the directional weights
% D^m_i used by the radiation FVM (Eq. (11) of Kim 2007), extended to 3D:
%
%       D^m_i = integral over DeltaOmega^m of (s . n_i) dOmega
%
% where n_i is geom.faceNormal(i, :) and DeltaOmega^m is the cone
% [theta^-, theta^+] x [phi^-, phi^+] from build_angles. Since the
% integrand is linear in n_i, we can decompose:
%
%       D^m_i = n_x . Dx^m + n_y . Dy^m + n_z . Dz^m
%
% where (Dx^m, Dy^m, Dz^m) are the three axial integrals (independent
% of the face, depend only on the cone) :
%
%       Dx^m = (sin(phi^+) - sin(phi^-)) . I_theta
%       Dy^m = (cos(phi^-) - cos(phi^+)) . I_theta
%       Dz^m = 0.5 . (phi^+ - phi^-) . (sin^2(theta^+) - sin^2(theta^-))
%
%       I_theta = (theta^+ - theta^-) / 2 - (sin(2 theta^+) - sin(2 theta^-)) / 4
%
% These formulas are EXACT for the midpoint quadrature cones produced
% by build_angles. They work for any face normal (axis-aligned or not).
%
% Inputs
%   geom : with .faceNormal (Nf x 3). Convention : oriented from owner
%          to neighbour (for internal faces) ; outward from the domain
%          (for boundary faces). D is computed with this normal ; the
%          sign for the owner / neighbour viewpoint is handled at the
%          assembly step.
%   ang  : output of build_angles, must have .thetaBounds (M x 2),
%          .phiBounds (M x 2), .M.
%
% Output
%   D : Nf x M matrix of integral values (no face area factor).

    if size(geom.faceNormal, 2) ~= 3
        error('precompute_face_direction_weights: geom.faceNormal must be Nf x 3.');
    end
    if ~isfield(ang, 'thetaBounds') || ~isfield(ang, 'phiBounds')
        error(['precompute_face_direction_weights: ang must contain ', ...
               'thetaBounds and phiBounds (rebuild via the updated ', ...
               'build_angles).']);
    end

    Nf = size(geom.faceNormal, 1);
    M  = ang.M;

    thMin = ang.thetaBounds(:, 1);   % M x 1
    thMax = ang.thetaBounds(:, 2);
    phMin = ang.phiBounds(:, 1);
    phMax = ang.phiBounds(:, 2);

    % Axial integrals (M x 1 each)
    I_theta = 0.5 * (thMax - thMin) ...
            - 0.25 * (sin(2 * thMax) - sin(2 * thMin));

    Dx = (sin(phMax) - sin(phMin))                .* I_theta;
    Dy = (cos(phMin) - cos(phMax))                .* I_theta;
    Dz = 0.5 * (phMax - phMin) .* (sin(thMax).^2 - sin(thMin).^2);

    DAxial = [Dx, Dy, Dz];           % M x 3

    % D(f, m) = geom.faceNormal(f, :) * DAxial(m, :).'  ->  Nf x M
    D = geom.faceNormal * DAxial.';
end
