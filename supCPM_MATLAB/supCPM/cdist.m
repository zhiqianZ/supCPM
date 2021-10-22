function Dst = cdist(D,no_dims,compel_force)
% The function computes the Capacity adjusted distance based on the pairwise
% distance in the original metric
% D: Pairwise distance matrix
% no_dims: mapped dimension
% (optional) w (=2 or 3): determines the dimension correction power 
     n = size(D,1);   
     nscales = 7; % number of scales for dimension estimation
     [val, I] = sort(sqrt(D(:)),'ascend');
     [dim1,r] = compute_dim(sqrt(D),nscales,0.01,val);
     if compel_force == 1
         dim1(2) = 1.5*dim1(2);
     end
     max_dim = 10*(no_dims-1);   % set an upper bound on the estimated dimension to avoid divergence when the data dimension is too big
     ind = (dim1<max_dim);
     dim = ones(size(dim1))*max_dim;
     dim(ind) = dim1(ind);
     dim(dim1<no_dims) = no_dims; %the lower bound on the estimated dimension is set to the mapped-to dimension
    
  %   r1 = sqrt(D(I(ceil(n*n/100))));
     r1 = prctile(sort(sqrt(D(:))),1); % a small postive number preventing dividing by 0 when defining the probability
     curr_dim = zeros(n,n);
     r0 = [-1;r];
     for k = 2:(length(r)+1)
      curr_dim(sqrt(D)<=r0(k)&sqrt(D)>r0(k-1)) = dim(k-1);
     end
      Dst = (r1^2 + D).^(curr_dim/no_dims);
      Dst(I) =  make_mono(sqrt(Dst(I)),sqrt(D(I)+r1^2));  % impose monotonicity
      Dst = sqrt(Dst);
end
function d = find_dim(r,dim,a)
% find the estimated dimension at the scale r
i = 1;
while i<= length(r)
    if a <= r(i)
        d = dim(i);
        break;
    else
        i = i+1;
    end
end
%d = dim(find(r>=a,1));
end
function x=make_mono(x,y)
d = max(x(1)-y(1),0);
for i = 2:length(x)
    x(i) = max([x(i),x(i-1),y(i)+d]);
    d = max(d,max(x(i)-y(i),0));
end
x = x.^2;
end

