function ydata= supCPM(data,label,no_dims,compel_force,geodist,degree,ratio,k,change,niter,seed,factor)
% data: N x P data matrix
% no_dims: mapping dimension
% (optional)compel_force=1 if user wants to pull clusters a bit apart
%                       =0 if user wants the best preservation of the geometry
% (optional) niter : number of iterations
% (optional) D: the pairwise distance matrix to be used to define probabilities, 
%               if not an input, it will be computed from data

% Compute the pairwise Euclidean(geodesic) distances
 if ~exist('D')
  if geodist ==0 
     D = compute_dist(data,0); 
  else
      disp('computing the geodesic distance')
      D = compute_dist(data,1,k); 
  end
 end
 if ~exist('compel_force')
     compel_force = 0;
 end
 n = size(D,1);
       D = D - min(D(:));
      D = D / max(D(:));
      D(1:n+1:end)=0;
      D = max(D,0);
 % compute the Capacity adjusted distance
 D1 = cdist(D,no_dims,compel_force);

 L_tilde = zeros(n);
 for i = unique(label)'
     cl = find(label==i);
     L_tilde(cl,cl)=1/length(cl);
 end
 D1(L_tilde==0) = D1(L_tilde==0)*factor;
 
  % change the distances to probabilities
     vP_D = 1./((D1.^2).^((1+degree)/2));
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
    epsilon = 500;                                      
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
        num(1:n+1:end) = 0;                                                 % set diagonal to zero
        Q = max(num ./ sum(num(:)), realmin);                               % normalize to get probabilities
        
        % gradient 
        L = (P - Q) .* num;
        y_grads = 4 * (diag(sum(L, 1)) - L) * ydata;
        
        if iter >= change
          y = ydata';
          y_grads = (1-ratio) * y_grads +  ratio * 2 * (trace(y*Sb*y')*Sw*y'-trace(y*Sw*y')*Sb*y')/trace(y*Sb*y')^2;
        end
       
        gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...     
              + (gains * .8) .* (sign(y_grads) == sign(y_incs));
        gains(gains < min_gain) = min_gain;
        
        if iter == change
            epsilon = epsilon/5e3;
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
            cost =  (1-ratio) * const - (1-ratio) * sum(P(:) .* log(Q(:))) + ratio * trace(y*Sw*y')/trace(y*Sb*y');
            disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]); 
     %   figure(1); scatter(ydata(:,1),ydata(:,2),[],label); pause(.2)
        end
        if ~rem(iter, 10) && iter <= change
            cost = const - sum(P(:) .* log(Q(:)));
            disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]);
      %    figure(1); scatter(ydata(:,1),ydata(:,2),[],label); pause(.2)
        end
    end
end

