function D = geodscdistance(X,k)
% This function computes the geodesic distance for multi-components dataset
% X: N x P data matrix 
% k: number of neighbours
% Author: Xiaopeng Zhang
n = size(X,1);
sum_X = sum(X.^2,2);
DD = bsxfun(@plus,sum_X,bsxfun(@plus,sum_X',-2*X*X'));
DD(1:n+1:end) = 0;
DD = sqrt(DD);
 D = real(find_nn(X,k));
 D = full(D);
blocks = components(D)';
for i = 1:(max(blocks)-1)
    for j = i+1:max(blocks)
        c1 = find(blocks==i);
        c2 = find(blocks==j);
        Ds = DD(c1,c2);
        [a,b] = find(Ds==min(Ds(:)));
        D(c1(a),c2(b)) = 2*DD(c1(a),c2(b));
        D(c2(b),c1(a)) = 2*DD(c2(b),c1(a));
    end
end
if n<500
   D=dijkstra(D,1:n);
else
   D = approx_geodesic(D);
end


end