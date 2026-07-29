% Regression tests for triangular spectral smoothing.

scriptPath = fileparts(mfilename('fullpath'));
addpath(scriptPath);

f = (0.2:0.1:1.0)';
x = f.^2;

% A non-empty window must preserve the historical calculation.
fcRegular = 0.5;
regularMask = f >= 0.4 & f <= 0.6;
expectedRegular = mean(x(regularMask).*triang(nnz(regularMask)));
actualRegular = triangFilter(x,f,fcRegular,20);
assert(abs(actualRegular-expectedRegular) <= 10*eps(expectedRegular));

% A window between FFT bins must fall back to the nearest bin.
fcEmpty = 0.26;
actualEmpty = triangFilter(x,f,fcEmpty,1);
[~,nearestBin] = min(abs(f-fcEmpty));
assert(actualEmpty == x(nearestBin));
assert(isfinite(actualEmpty));

% The fallback must work independently on every matrix column.
xMatrix = [x,2*x,3*x];
actualMatrix = triangFilter(xMatrix,f,fcEmpty,1);
assert(isequal(actualMatrix,xMatrix(nearestBin,:)));
assert(all(isfinite(actualMatrix),'all'));

% Three-dimensional inputs must retain columns and pages.
x3 = cat(3,xMatrix,4*xMatrix);
actual3 = triangFilter(x3,f,[fcEmpty;fcRegular],1);
assert(isequal(size(actual3),[2,3,2]));
assert(all(isfinite(actual3),'all'));
assert(isequal(squeeze(actual3(1,:,:)),squeeze(x3(nearestBin,:,:))));

fprintf('TUTTI I TEST TRIANGFILTER SONO SUPERATI\n');
