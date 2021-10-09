function D = distance(X)
    sumX = sum(X .^ 2, 2);
    D =  bsxfun(@plus, sumX, bsxfun(@plus, sumX', -2 * (X * X')));
end