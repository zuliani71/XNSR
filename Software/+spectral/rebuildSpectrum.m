function X = rebuildSpectrum(H, n, symmetry, dim)
%REBUILDSPECTRUM Rebuild a full spectrum from its non-negative half.
%   X = spectral.rebuildSpectrum(H,N) uses Hermitian symmetry.
%   SYMMETRY can be "hermitian" or "mirror".

arguments
    H {mustBeNumeric}
    n (1,1) double {mustBeInteger, mustBePositive}
    symmetry (1,1) string {mustBeMember(symmetry,["hermitian","mirror"])} = "hermitian"
    dim (1,1) double {mustBeInteger, mustBePositive} = 1
end

expected = floor(n / 2) + 1;
if size(H, dim) ~= expected
    error("spectral:rebuildSpectrum:InvalidHalfLength", ...
        "Dimension %d of H must contain %d values for N=%d.", ...
        dim, expected, n);
end

indices = repmat({':'}, 1, max(ndims(H), dim));
if mod(n, 2) == 0
    indices{dim} = (size(H, dim) - 1):-1:2;
else
    indices{dim} = size(H, dim):-1:2;
end
tail = H(indices{:});
if symmetry == "hermitian"
    tail = conj(tail);
end
X = cat(dim, H, tail);
end
