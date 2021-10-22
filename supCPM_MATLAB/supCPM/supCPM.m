function ydata= supCPM(data,label,no_dims,compel_force,geodist,degree,ratio,k,change,niter,seed,factor)
% data: N x P data matrix
% no_dims: mapping dimension
% (optional)compel_force=1 if user wants to pull clusters a bit apart
%                       =0 if user wants the best preservation of the geometry
% (optional) niter : number of iterations
% (optional) D: the pairwise distance matrix to be used to define probabilities, 
%               if not an input, it will be computed from data

% Compute the pairwise Euclidean(geodesic) distances
 n = size(data,1);
 if min(label)==0
    label =label+1;
 end
 [label,I] = sort(label);
 I_recov=I;
 I_recov(I)=1:n;
 data = data(I,:);
 summ = tabulate(label);
 cluster_size = summ(:,2);
 cluster_num = length(unique(label));
 if ~exist('D')
  if geodist ==0 
      disp('computing the euclidean distance')
     D = compute_dist(data,0); 
  else
      disp('computing the geodesic distance')
      D = compute_dist(data,1,k); 
  end
 end
 if ~exist('compel_force')
     compel_force = 0;
 end
 D = D / max(D(:));
 D(1:n+1:end)=0;
 % compute the Capacity adjusted distance
 D1 = cdist(D,no_dims,compel_force);
 %clear D;
 L_tilde = zeros(n);
 for i = unique(label)'
     cl = find(label==i);
     L_tilde(cl,cl)=1;
 end
 D1(L_tilde==0) = (D1(L_tilde==0)).* factor;
  % change the distances to probabilities
     vP_D = 1./((D1.^2/degree).^((degree+1)/2));
     vP_D(1:n+1:end) = 0;
     P_D  = reshape(vP_D,size(vP_D));
     P = max(P_D/sum(P_D(:)), realmin); 
     P = 0.5 * (P + P'); 
     const = real(sum(P(:).*log(P(:))));
   % minimize KL divergence w.r.t. P_D
    
    P = P*4;
    if exist('seed')
       rng(seed);
    end
    ydata = 0.0001*randn(n,no_dims);
%     if init == 1
%         disp('computing initialization')
%         D_temp = (D'+D)/2;
%         D_temp(1:n+1:end) = 0;
%         ydata = cmdscale(D_temp,2);
%         %ydata = ydata(:,1:2);
%     end
    
%    clearvars cl D D1 data L_tilde P_D vP_D 
    %const statement
    momentum = 0.5;
    final_momentum = 0.8;                               % value to which momentum is changed
    mom_switch_iter = 500;                              % iteration at which momentum is changed
    stop_lying_iter = 200;                    % iteration at which lying about P-values is stopped
    if exist('niter')
        max_iter = niter;
    else
        max_iter = 1000;
    end                                               
    epsilon = 500;                                      
    min_gain = .01;  
    y_incs  = zeros(size(ydata));
    gains = ones(size(ydata));
    
    
    % Run the iterations
   
    for iter=1:max_iter
        [y_grads,Q] = unsup_grad(ydata,P);
        
        if iter >= change
            num = 1./(1+distance(ydata));
            num(1:n+1:end) = 0;          
            Q = max(num ./ sum(num(:)), realmin); 
            L = (P-Q) .* num;
            y_grads = 4*(sum(L, 1))'.*ydata;
            y_grads = y_grads - 4 * L *ydata;
           [grad,TraceSw,TraceSb] = sup_grad(ydata,label,cluster_num,cluster_size);
%            if iter == change
%             L1  = const - sum(P(:) .* log(Q(:)));
%             L2 = TraceSb/TraceSw;
%             ratio = ratio*L1/(L2+ratio*L1);
%            end
           y_grads = (1-ratio) * y_grads +  2 * ratio * grad;
        end
       
        gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...     
              + (gains * .8) .* (sign(y_grads) == sign(y_incs));
        gains(gains < min_gain) = min_gain;
        
         if iter == change
             epsilon = epsilon/1e3;
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
        
        if iter==change+1
            cost =  (1-ratio) * const - (1-ratio) * sum(P(:) .* log(Q(:))) + ratio * TraceSb/TraceSw;
            disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]); 
        end
        % print error,keep iteration
        if  iter >= change && ~rem(iter, 10) 
            cost =  (1-ratio) * const - (1-ratio) * sum(P(:) .* log(Q(:))) + ratio * TraceSb/TraceSw;
            disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]); 
     %   figure(1); scatter(ydata(:,1),ydata(:,2),[],label); pause(.2)
        end
        if iter < change && ~rem(iter, 10) 
            cost = const - sum(P(:) .* log(Q(:)));% + sum((1-P(:)) .* log(1-P(:))) - sum((1-P(:)) .* log(1-Q(:)));
            disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]);
      %    figure(1); scatter(ydata(:,1),ydata(:,2),[],label); pause(.2)
        end
    end
    ydata = ydata(I_recov,:);
end

