% Regression and performance tests for the optimized Konno-Ohmachi filter.
rng(41);

n = 127;
m = 53;
k = 5;
pages = 3;
f = linspace(0.1,30,n)';
fc = linspace(0.2,20,m)';
b = 40;
X = abs(randn(n,k,pages) + 1i*randn(n,k,pages));

% Independent implementation of the historical REPMAT/ACCUMARRAY method.
expected = legacyKonnoOhmachiFilter(X,f,fc,b);

% New direct call and precomputed-weight call must match the old result.
actual = KonnoOhmachiFilter(X,f,fc,b);
weights = KonnoOhmachiWeights(f,fc,b);
reused = KonnoOhmachiFilter(X,f,fc,b,weights);

tolerance = 5e-13 * max(1,max(abs(expected(:))));
assert(max(abs(actual(:)-expected(:))) < tolerance);
assert(max(abs(reused(:)-expected(:))) < tolerance);
assert(all(abs(sum(weights,1)-1) < 5e-14));
assert(isequal(size(actual),[m,k,pages]));

% Verify vectors, matrices and complex spectra.
complexMatrix = randn(n,k) + 1i*randn(n,k);
complexExpected = legacyKonnoOhmachiFilter(complexMatrix,f,fc,b);
complexActual = KonnoOhmachiFilter(complexMatrix,f,fc,b,weights);
complexTolerance = 5e-13 * max(1,max(abs(complexExpected(:))));
assert(max(abs(complexActual(:)-complexExpected(:))) < complexTolerance);

vectorResult = KonnoOhmachiFilter(X(:,1,1),f,fc,b,weights);
assert(isequal(size(vectorResult),[m,1]));

% Exercise the same page-wise PARFOR pattern used by XN_Cruncher.
parallelResult = zeros(m,k,pages);
parfor pageIndex = 1:pages
    parallelResult(:,:,pageIndex) = KonnoOhmachiFilter( ...
        X(:,:,pageIndex),f,fc,b,weights);
end
assert(max(abs(parallelResult(:)-expected(:))) < tolerance);

% Benchmark both implementations on the same representative array.
nBench = 512;
mBench = 250;
kBench = 12;
pagesBench = 4;
fBench = linspace(0.1,30,nBench)';
fcBench = linspace(0.2,20,mBench)';
XBench = abs(randn(nBench,kBench,pagesBench) + ...
    1i*randn(nBench,kBench,pagesBench));
weightsBench = KonnoOhmachiWeights(fBench,fcBench,b);

tic;
legacyBench = legacyKonnoOhmachiFilter(XBench,fBench,fcBench,b);
legacyTime = toc;

tic;
optimizedBench = KonnoOhmachiFilter( ...
    XBench,fBench,fcBench,b,weightsBench);
optimizedTime = toc;

benchmarkTolerance = 5e-13 * max(1,max(abs(legacyBench(:))));
assert(max(abs(optimizedBench(:)-legacyBench(:))) < benchmarkTolerance);

fprintf('Konno-Ohmachi legacy:     %.6f s\n',legacyTime);
fprintf('Konno-Ohmachi optimized:  %.6f s\n',optimizedTime);
fprintf('Speed-up:                 %.2fx\n',legacyTime/optimizedTime);
disp('TUTTI I TEST KONNO-OHMACHI SONO SUPERATI');

function Y = legacyKonnoOhmachiFilter(X,f,fc,b)
% Reference copy of the original algorithm, retained only for regression.
inputSize = size(X);
n = numel(f);
m = numel(fc);
X = reshape(X,n,[]);
k = size(X,2);

win = KonnoOhmachiSmoothingWindow(f,fc,b);
win = reshape(win,m*n,1);
signal = reshape(repmat(X,m,1),m*k*n,1);
win = repmat(win,k,1);
products = signal.*win;
subs = reshape(repmat(1:m*k,n,1),n*m*k,1);
normalization = accumarray(subs,win);
Y = reshape(accumarray(subs,products)./normalization,m,k);

inputSize(1) = m;
Y = reshape(Y,inputSize);
end
