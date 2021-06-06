function mappedX = mds(X, no_dims,type)
%MDS Run MDS on the data to get a low-dimensional visualization
% 
%   mappedX = mds(X, no_dims)
%
% Run multidimensional scaling on the dataset X to get a two-dimensional 
% visualization. The low-dimensional representation is returned in mappedX.
% It has dimensionality no_dims (default = 2).
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if ~exist('no_dims', 'var')
        no_dims = 2;
    end

    % NOTE: Classical scaling is identical to performing PCA, except the
    % input data is different. Specifying pairwise similarity data is not
    % yet supported by the toolbox.
    sumX = sum(X .^ 2, 2);
    D  = bsxfun(@plus,sumX, bsxfun(@plus,sumX',-2*(X*X')));
    if type ==1
    D = cdist(D, no_dims);
    D = D.^2;
    end
    n = size(D,1);
    D(1:(n+1):end)=0;
    B = bsxfun(@plus, bsxfun(@minus, bsxfun(@minus, D, sum(D,1)/n),sum(D,2)/n), sum(D(:))/(n^2))*(-0.5);
    mappedX = compute_mapping(B,'PCA',no_dims);
    