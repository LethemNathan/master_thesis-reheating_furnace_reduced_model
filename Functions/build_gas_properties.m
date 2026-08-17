function gasProps = build_gas_properties(gasComposition, opts)
%BUILD_GAS_PROPERTIES cp(T) and h(T) for the flue-gas mixture and for the
% burner reactant mixture.
%
% PRODUCT (flue-gas) mixture. From a fixed molar composition [O2 N2 H2O CO2],
% builds the mass-based mixture specific heat cp_mix(T) [J/(kg K)] and enthalpy
% h_mix(T) [J/kg] used by the advective / outflow terms of the zone energy
% balance.
%
% REACTANT mixture (added). The burners inject AIR (O2/N2) and NATURAL GAS
% (CH4/C2H6), not flue gas. The sensible enthalpy carried in by this stream,
% h_reac(T), is therefore built from a separate, mass-based reactant mixture
% and returned as gasProps.hReac. It is used by solve_zones_energy for the
% Q_reactant term; the product mixture (gasProps.h) still handles advection and
% outflow.
%
% Species cp are per-mass polynomials in T (ascending order, 7 coeffs):
%       cp_k(T) = sum_{j=0}^{6} a_{k,j} T^j        [J/(kg K)]
% Mixture cp is the mass-weighted average, itself a degree-6 polynomial.
% Enthalpy is the definite integral from Tref:  h(T) = int_{Tref}^{T} cp dT'.
%
% Inputs
%   gasComposition : [x_O2 x_N2 x_H2O x_CO2] molar fractions of the PRODUCT
%                    mixture (renormalised internally).
%   opts           : optional
%       .Tref         reference temperature for h [K] (default 298.15).
%       .reactantAirMolar [x_O2 x_N2] molar fractions in air
%                     (default [0.2080 0.7920]).
%       .reactantNGMolar  [x_CH4 x_C2H6] molar fractions in natural gas
%                     (default [0.9790 0.0210]).
%       .reactantYair air MASS fraction of the injected stream (default from
%                     the lateral burner flows: 0.4056 air + 0.02156 NG).
%
% Output struct gasProps
%   .cp,  .h        product-mixture cp(T) [J/(kg K)] and h(T) [J/kg] handles
%   .cpReac, .hReac reactant-mixture cp(T) and h(T) handles
%   .cpCoeff, .hCoeff, .hRef              product polynomial data
%   .cpReacCoeff, .hReacCoeff, .hReacRef  reactant polynomial data
%   .Y      (1 x 4) product mass fractions [O2 N2 H2O CO2]
%   .Yreac  (1 x 4) reactant mass fractions [O2 N2 CH4 C2H6]
%   .M, .Mmix, .Tref

    if nargin < 2 || isempty(opts), opts = struct(); end
    Tref = 298.15;
    if isfield(opts,'Tref') && ~isempty(opts.Tref), Tref = opts.Tref; end

    x = double(gasComposition(:)).';
    if numel(x) ~= 4
        error('build_gas_properties: gasComposition must be [O2 N2 H2O CO2].');
    end
    x = x / sum(x);                         % renormalise molar fractions

    %% ================= PRODUCT (flue-gas) mixture =================

    % --- Species per-mass cp polynomials (ascending, 7 coeffs) : O2 N2 H2O CO2 ---
    cpSpecies = [ ...
        975.974, -0.66073,   0.00223857, -2.451160e-06,  1.281630e-09, -3.236540e-13,  3.172820e-17;   % O2
       1116.63,  -0.615175,  0.0015183,  -1.296410e-06,  5.498740e-10, -1.167650e-13,  9.910130e-18;   % N2
       1875.84,  -0.419047,  0.00157042, -9.671970e-07,  2.533150e-10, -2.497410e-14,  0;              % H2O
        452.206,  1.6783,   -0.00140075,  6.435410e-07, -1.534180e-10,  1.475850e-14,  0 ];            % CO2

    % --- Molar masses [g/mol] : O2 N2 H2O CO2 ---
    M = [31.9988, 28.0134, 18.01534, 44.00995];

    % --- Mass fractions of the product mixture ---
    xM   = x .* M;
    Mmix = sum(xM);
    Y    = xM / Mmix;                       % 1 x 4

    % --- Mixture cp polynomial and antiderivative ---
    cpCoeff = Y * cpSpecies;                % 1 x 7 ascending
    hCoeff  = antideriv(cpCoeff);
    hRef    = eval_poly(hCoeff, Tref);

    %% ================= REACTANT (air + natural gas) mixture =================

    % --- Reactant species per-mass cp polynomials (ascending, 7 coeffs) :
    %     O2 N2 CH4 C2H6  (C2H6 has 4 coeffs, zero-padded) ---
    cpReacSpecies = [ ...
        975.974, -0.66073,   0.00223857, -2.451160e-06,  1.281630e-09, -3.236540e-13,  3.172820e-17;   % O2
       1116.63,  -0.615175,  0.0015183,  -1.296410e-06,  5.498740e-10, -1.167650e-13,  9.910130e-18;   % N2
       2490.75,  -4.86764,   0.0189244,  -1.958560e-05,  9.608440e-09, -2.283220e-12,  2.112000e-16;   % CH4
        354.153,  5.44874,  -0.0020021,   2.524810e-07,  0,             0,             0 ];             % C2H6

    % --- Reactant stream composition ---
    % Air and natural gas are given as MOLAR fractions (consistent with the
    % product side) and converted to mass fractions here, using the reactant
    % molar masses.
    Mreac = [31.9988, 28.0134, 16.04246, 30.06904];     % g/mol : O2 N2 CH4 C2H6

    xAir = getf(opts,'reactantAirMolar', [0.2080 0.7920]);  % [O2 N2]    molar frac in air
    xNG  = getf(opts,'reactantNGMolar',  [0.9790 0.0210]);  % [CH4 C2H6] molar frac in natural gas
    xAir = xAir / sum(xAir);
    xNG  = xNG  / sum(xNG);

    % Convert each stream from molar to mass fractions
    yAir = xAir .* Mreac(1:2);  yAir = yAir / sum(yAir);    % [O2 N2]    mass frac in air
    yNG  = xNG  .* Mreac(3:4);  yNG  = yNG  / sum(yNG);     % [CH4 C2H6] mass frac in NG

    % Air / fuel MASS split of the injected stream. Default from the lateral
    % burner mass flows (0.4056 kg/s air + 0.02156 kg/s natural gas). Lateral
    % and screen burners share the same fuel-to-air ratio, so a single reactant
    % mixture applies to both.
    mAirRef = 0.4056;  mNGRef = 0.02156;
    Yair = getf(opts,'reactantYair', mAirRef / (mAirRef + mNGRef));
    Yng  = 1 - Yair;
    Yreac = [Yair*yAir(1), Yair*yAir(2), Yng*yNG(1), Yng*yNG(2)];   % [O2 N2 CH4 C2H6] mass

    % --- Reactant mixture cp polynomial and antiderivative ---
    cpReacCoeff = Yreac * cpReacSpecies;    % 1 x 7 ascending
    hReacCoeff  = antideriv(cpReacCoeff);
    hReacRef    = eval_poly(hReacCoeff, Tref);

    %% ================= Pack output =================
    gasProps             = struct();
    gasProps.cp          = @(T) eval_poly(cpCoeff, T);
    gasProps.h           = @(T) eval_poly(hCoeff, T) - hRef;
    gasProps.cpReac      = @(T) eval_poly(cpReacCoeff, T);
    gasProps.hReac       = @(T) eval_poly(hReacCoeff, T) - hReacRef;
    gasProps.cpCoeff     = cpCoeff;
    gasProps.hCoeff      = hCoeff;
    gasProps.hRef        = hRef;
    gasProps.cpReacCoeff = cpReacCoeff;
    gasProps.hReacCoeff  = hReacCoeff;
    gasProps.hReacRef    = hReacRef;
    gasProps.Y           = Y;
    gasProps.Yreac       = Yreac;
    gasProps.M           = M;
    gasProps.Mmix        = Mmix;
    gasProps.Tref        = Tref;

    fprintf(['build_gas_properties: Mmix=%.2f g/mol | Yprod=[%.3f %.3f %.3f %.3f] | ', ...
             'cp(1500)=%.1f J/kgK\n'], ...
            Mmix, Y(1),Y(2),Y(3),Y(4), gasProps.cp(1500));
    fprintf(['build_gas_properties: Yreac=[%.4f %.4f %.4f %.4f] (O2 N2 CH4 C2H6) | ', ...
             'hReac(578.83)=%.0f J/kg  (h_prod=%.0f)\n'], ...
            Yreac(1),Yreac(2),Yreac(3),Yreac(4), ...
            gasProps.hReac(578.83), gasProps.h(578.83));
end


%% ================= LOCAL FUNCTIONS =================

function hCoeff = antideriv(cpCoeff)
%ANTIDERIV Antiderivative coefficients of a cp polynomial:
%   int( sum a_j T^j ) = sum a_j/(j+1) T^{j+1}
%   -> ascending coeffs [0, a0/1, a1/2, ...] (constant term 0).
    n      = numel(cpCoeff);
    hCoeff = zeros(1, n+1);
    for j = 0:n-1
        hCoeff(j+2) = cpCoeff(j+1) / (j+1);
    end
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
