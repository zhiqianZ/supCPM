function ydata = supCPM(data,label,no_dims,degree,distance,w,ratio,k,factor,change,niter,seed)
% data: N x P data matrix (N is the sample size and P is the feature number)
% no_dims: mapping dimensions
% degree: the degree of freedom of t-distribution
% distance (value of 0 and 1): choise of distance, 0-euclidean; 1- geodesic
% w (value betweeen 2-10): determines the dimension correction power
% ratio (value between 0-1): the trade-off between KL-divergence and supervised term
% k: k nearest neighbors for geodesic distance
% factor: a constant to multiply on distance matrix
% change: the iteration time for the first phase
% niter: total iteration
% seed: random seed for reproduction

% Compute the pairwise distances
n = size(data,1);
   if distance == 0 
     disp('computing the euclidean distance')
     D = compute_dist(data,0); 
   end
   if distance == 1 
      disp('computing the geodesic distance')
      D = compute_dist(data,1,k,0.01,label); 
   end
      D = D - min(D(:));
      D = D / max(D(:));
      D(1:n+1:end)=0;
      D = max(D,0);
  
  % compute the Capacity adjusted distance
  if exist('w')
       D1 = cdist(D,no_dims,w);   
  else
       D1 = cdist(D,no_dims);
  end
      L_tilde = zeros(n);
      for i = unique(label)'
           cl = find(label==i);
           L_tilde(cl,cl)=1/length(cl);
      end
        D1(L_tilde==0) = D1(L_tilde==0)*factor;
  % change the distances to probabilities
     vP_D = 1./((D1.^2./degree).^((1+degree)/2));
     P_D  = reshape(vP_D,size(vP_D));
     P  = P_D./sum(P_D(:));
     P(1:n+1:end) = 0;
     P = max(P ./ sum(P(:)), realmin); 
     P = 0.5 * (P + P'); 
     const = real(sum(P(:).*log(P(:))));
   % minimize KL divergence w.r.t. P_D
   
    P = P*4;
    if exist('seed')
       rng(seed);
    end
    ydata = 0.0001*randn(n,no_dims);
    
    
    %const statement
    momentum = 0.5;
    final_momentum = 0.8;                               % value to which momentum is changed
    mom_switch_iter = 500;                              % iteration at which momentum is changed
    stop_lying_iter = 100;                    % iteration at which lying about P-values is stopped
    if exist('niter')
        max_iter = niter;
    else
        max_iter = 1000;
    end                                               
    epsilon = 1000;                                     
    min_gain = .01;  
    y_incs  = zeros(size(ydata));
    gains = ones(size(ydata));
    
   
    I = eye(n);
    one = ones(n,1);  
    Sb = (L_tilde-one*one'/n);
    Sw = (I - L_tilde);   
        
    % Run the iterations
    for iter=1:max_iter
        sum_ydata = sum(ydata .^ 2, 2);
        num = 1 ./ (1 + bsxfun(@plus, sum_ydata, bsxfun(@plus, sum_ydata', -2 * (ydata * ydata')))); % Student-t distribution
        num(1:n+1:end) = 0;     % set diagonal to zero
        Q = max(num ./ sum(num(:)), realmin);                               % normalize to get probabilities
         
        % gradient 
        L = (P - Q) .* num;
        y_grads = 4 * (diag(sum(L, 1)) - L) * ydata;
        
       if iter >= change
          y = ydata';
          y_grads = (1-ratio) * y_grads + ratio * 2 *(((y*Sw)'+(y*Sb)')/(trace(y*Sw*y')+trace(y*Sb*y')) - (y*Sb)'/trace(y*Sb*y'));
       end
       
        gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...     
              + (gains * .8) .* (sign(y_grads) == sign(y_incs));
        gains(gains < min_gain) = min_gain;
        if iter == change
            epsilon = epsilon/1e4;
        end
        
        y_incs = momentum * y_incs - epsilon * (gains .* y_grads);
        ydata = ydata + y_incs;
        ydata = bsxfun(@minus, ydata, mean(ydata, 1));
        
        % 
        if iter == mom_switch_iter
            momentum = final_momentum;
        end
        if iter == stop_lying_iter
            P = P ./ 4;
        end
        
        
        % print error,keep iteration
        if ~rem(iter, 10) && iter > change
            cost = (1-ratio) * const - (1-ratio) * sum(P(:) .* log(Q(:))) + ratio * (log(trace(y*Sw*y')+trace(y*Sb*y'))-log(trace(y*Sb*y')));
            disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]); 
        end
        if ~rem(iter, 10) && iter <= change
           cost = const - sum(P(:) .* log(Q(:)));
        disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]); 
        end
 
    end
end
