function [dim,r] = compute_dim(D,nscales,smallest_scale)
%This function estimates the intrinsic dimension at various scales
% nscales: total number of scales
% smallest_scale: the smallest scale in consideration


%compute the radius at each scale
      factor = (1/smallest_scale)^(1/nscales);
      val = sort((D(:)));
      r = zeros(nscales,1);
      for i = 1:nscales
          r(i) = prctile(val,factor^i/factor^nscales*100);
      end
       r1 = prctile(val,1/factor^nscales*100);
       dim = zeros(nscales,1);
       dim(1) = dim_pow(val,1e-5,r1,r(1));
      for i = 1:nscales-1
          r2 = r(i);                % at each scale, three radii are needed to compute the dimension
          if i>1
              r11 = r(i-1);
          else
              r11 = r1;
          end
         dim(i+1) = dim_pow(val,r11,(r2-r11)/2+r11,r2); % estimate the dimension at the current scale
      end
     for i = 1:nscales
        if i > nscales-2
            dim(i) = max(dim(i-1),dim(i-1)); % at large scale, the data is too sparse for an estimation, use the dimension of the previous scale
        end
     end        
end
     
  