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

% Precomputed sparse weights must reproduce the direct calculation.
fcVector = [fcEmpty;fcRegular;0.8];
weights = triangFilterWeights(f,fcVector,20);
assert(issparse(weights));
direct = triangFilter(x3,f,fcVector,20);
precomputed = triangFilter(x3,f,fcVector,20,weights);
assert(max(abs(direct(:)-precomputed(:))) <= ...
    10*eps(max(abs(direct(:)))));

% Compare against the historical per-frequency definition.
manual = zeros(numel(fcVector),size(xMatrix,2));
for centerIndex = 1:numel(fcVector)
    mask = f >= fcVector(centerIndex)*(1-20/100) & ...
        f <= fcVector(centerIndex)*(1+20/100);
    if ~any(mask)
        [~,index] = min(abs(f-fcVector(centerIndex)));
        mask(index) = true;
    end
    manual(centerIndex,:) = ...
        sum(xMatrix(mask,:).*triang(nnz(mask)),1)/nnz(mask);
end
optimized = triangFilter(xMatrix,f,fcVector,20,weights);
assert(max(abs(manual(:)-optimized(:))) <= ...
    10*eps(max(abs(manual(:)))));

fprintf('TUTTI I TEST TRIANGFILTER SONO SUPERATI\n');
