function D =  compute_dist(X,geodesic,k)

 [n,p]=size(X);
 initial_dims = min(500,p);
 X = X - min(X(:));
 X = X / max(X(:));
 X = bsxfun(@minus, X, mean(X, 1));
 if size(X, 2) < size(X, 1)
            C = X' * X;
 else
            C = (1 / size(X, 1)) * (X * X');
 end
 [M, lambda] = eig(C);
 [lambda, ind] = sort(diag(lambda), 'descend');
 M = M(:,ind(1:initial_dims));
 lambda = lambda(1:initial_dims);
 if ~(size(X, 2) < size(X, 1))
      M = bsxfun(@times, X' * M, (1 ./ sqrt(size(X, 1) .* lambda))');
 end
 X = X * M;
    
    %distance matrix
 if geodesic ~=1
    sum_X=sum(X.^2,2);
    D=bsxfun(@plus,sum_X,bsxfun(@plus,sum_X',-2*X*X'));
    D(1:n+1:end)=0;
    D = max(D,0);
 else
 
    D=geodscdistance(X,k);
    D=D.^2;
  end
end