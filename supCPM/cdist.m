function Dst = cdist(D,no_dims,compel_force)
% The function computes the Capacity adjusted distance based on the pairwise
% distance in the original metric
% D: Pairwise distance matrix
% no_dims: mapped dimension
% (optional) w (=2 or 3): determines the dimension correction power 
     n = size(D,1);   
     nscales = 7; % number of scales for dimension estimation
     [dim1,r] = compute_dim(sqrt(D),nscales,0.01);
     if compel_force == 1
         dim1(2) = 1.5*dim1(2);
     end
     max_dim = 10*(no_dims-1);   % set an upper bound on the estimated dimension to avoid divergence when the data dimension is too big
     ind = (dim1<max_dim);
     dim = ones(size(dim1))*max_dim;
     dim(ind) = dim1(ind);
     dim(dim1<no_dims) = no_dims; %the lower bound on the estimated dimension is set to the mapped-to dimension
    
    % for i=1:nscales-1
    %     dim(i+1)=max(dim(i),dim(i+1));
    % end
     r1 = prctile(sort(sqrt(D(:))),1); % a small postive number preventing dividing by 0 when defining the probability
     [~, I] = sort(D(:),'ascend');
     Dst = zeros(n,n);
     for  i = 1:n
         for j = 1:n
          curr_dim = find_dim(r,dim,sqrt(D(i,j))); % find the estimated dimension at the scale r
          Dst(i,j) = (r1^2+(D(i,j)))^((curr_dim)/(no_dims));   % defined the new distance
         end
     end
 
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
end
function x=make_mono(x,y)
d = max(x(1)-y(1),0);
for i = 2:length(x)
    x(i) = max([x(i),x(i-1),y(i)+d]);
    d = max(d,max(x(i)-y(i),0));
end
x = x.^2;
end

