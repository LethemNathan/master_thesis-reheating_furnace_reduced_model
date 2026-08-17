function val = eval_poly(coeffs, T)
%EVAL_POLY Evaluate a polynomial in T from coefficients given in ASCENDING order.
%
%   val = eval_poly(coeffs, T)
%
% Convention (ascending order):
%   coeffs(1) = a0  (constant)
%   coeffs(2) = a1
%   coeffs(3) = a2
%   ...
%   coeffs(n+1) = an
% so val = a0 + a1*T + a2*T^2 + ... + an*T^n.
%
% A scalar coeffs (or coeffs of length 1) returns a constant a0 broadcast over
% the shape of T (useful for constant rho, etc.).
%
% Inputs
%   coeffs : numeric vector, ascending order
%   T      : array of any shape (typically the slab's local T vector)
%
% Output
%   val    : array with same shape as T
%
% Implementation: Horner's method, vectorised in T.

    if isempty(coeffs)
        error('eval_poly: coeffs is empty.');
    end

    coeffs = double(coeffs(:)).';   % row vector
    n = numel(coeffs);

    % Horner's scheme: start with the highest-order coefficient,
    % multiply by T, add the next coefficient, repeat.
    val = coeffs(n) * ones(size(T));
    for k = (n-1):-1:1
        val = val .* T + coeffs(k);
    end
end
