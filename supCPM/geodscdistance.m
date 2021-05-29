function D = geodscdistance(X,k,percent,l)
% This function computes the geodesic distance for multi-components dataset
% X: N x P data matrix 
% k: number of neighbours
% Author: Xiaopeng Zhang
n = size(X,1);
sum_X = sum(X.^2,2);
DD = bsxfun(@plus,sum_X,bsxfun(@plus,sum_X',-2*X*X'));
DD(1:n+1:end) = 0;
DD = real(sqrt(DD));
D = real(find_nn(X,k));
D = full(D);
blocks = components(D)';
unique(blocks)
for i = 1:(max(blocks)-1)
    for j = i+1:max(blocks)
        c1 = find(blocks==i);
        c2 = find(blocks==j);
        Ds = DD(c1,c2);
        [a,b] = find(Ds==min(Ds(:)));
        s = sort(Ds(:));
        quan = s(ceil(length(c1)*length(c2)*percent)+1);
        disp([i,j,quan])
        D(c1(a),c2(b)) = 2*quan;
        D(c2(b),c1(a)) = 2*quan;
    end
end
[D,P]=dijkstra(D,1:n);
end