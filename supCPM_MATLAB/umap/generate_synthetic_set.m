function [synData,synLabels] = generate_synthetic_set(data,labels,sz)
%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
WARNING_THRESHOLD = 100;

if nargin < 3
    sz = 1e5;
end

[R,C] = size(data);

subsetIds=unique(labels);
nIds = length(subsetIds);
idxs=false(R,nIds);
synSize = zeros(1, nIds);
means = zeros(nIds,C);
covariances=zeros(C,C,nIds);

for i = 1:nIds
    id = subsetIds(i);
    idxs(:,i)= labels==id;
    nIdxs=nnz(idxs(:,i));
    synSize(i) = round(nIdxs*sz/R);
    means(i,:) = mean(data(idxs(:,i),:));
    covariances(:,:,i) = cov(data(idxs(:,i),:));
end

if any(synSize < WARNING_THRESHOLD)
    warning(['At least one of your classes in the synthetic dataset will have fewer than ' num2str(WARNING_THRESHOLD) ' data points. Proceed with caution']);
end

sz = sum(synSize);
synData = zeros(sz, C);
synLabels = zeros(sz, 1);
partialSizes=cumsum(synSize);

synData(1:partialSizes(1),:) = mvnrnd(means(1,:), covariances(:,:,1),synSize(1));
synLabels(1:partialSizes(1)) = subsetIds(1);

for i = 1:nIds-1
    synData((partialSizes(i)+1):partialSizes(i+1), :) = mvnrnd(means(i+1,:), covariances(:,:,i+1),synSize(i+1));
    synLabels((partialSizes(i)+1):partialSizes(i+1)) = subsetIds(i+1);
end

end

