function weights = triangFilterWeights(f,fc,pcent)
%TRIANGFILTERWEIGHTS Build normalized triangular smoothing weights.
%   WEIGHTS = TRIANGFILTERWEIGHTS(F,FC,PERCENT) returns a sparse
%   numel(F)-by-numel(FC) matrix. Each column represents one output
%   frequency. If a percentage window contains no input bin, the nearest
%   input frequency receives unit weight.

validateattributes(f,{'numeric'}, ...
    {'vector','real','finite','nonempty'},mfilename,'f',1);
validateattributes(fc,{'numeric'}, ...
    {'vector','real','finite','nonempty'},mfilename,'fc',2);
validateattributes(pcent,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'},mfilename,'pcent',3);

f = f(:);
fc = fc(:);
frequencyCount = numel(f);
centerCount = numel(fc);
lowerBound = fc-fc*pcent/100;
upperBound = fc+fc*pcent/100;

rowIndices = cell(centerCount,1);
columnIndices = cell(centerCount,1);
values = cell(centerCount,1);
for centerIndex = 1:centerCount
    bins = find(f >= lowerBound(centerIndex) & ...
        f <= upperBound(centerIndex));
    if isempty(bins)
        [~,bins] = min(abs(f-fc(centerIndex)));
    end
    binCount = numel(bins);
    rowIndices{centerIndex} = bins;
    columnIndices{centerIndex} = ...
        repmat(centerIndex,binCount,1);
    values{centerIndex} = triang(binCount)/binCount;
end

weights = sparse(vertcat(rowIndices{:}), ...
    vertcat(columnIndices{:}),vertcat(values{:}), ...
    frequencyCount,centerCount);
end
