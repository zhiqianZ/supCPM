function [dim,r] = compute_dim(D,nscales,smallest_scale,val)
%This function estimates the intrinsic dimension at various scales
% nscales: total number of scales
% smallest_scale: the smallest scale in consideration


%compute the radius at each scale
      factor = (1/smallest_scale)^(1/nscales);
      r = zeros(nscales,1);
      D5 = [];
      n = size(D,1);
%        for i=1:n
%                 %dist = sort(D(:,i));
%                 dist = mink(D(:,i),10);
%                 D5=[D5; dist(2:10)];             
%        end
        D5 = mink(D,10,2);
        D5 = D5(:,2:10);
       r(1) = prctile(D5(:),100);
       r1 = median(D5(:));
       percens = [];
       for i = 2:nscales
           percens = [percens,factor^i/factor^nscales*100];
       end
       %percens = prctile(val,percens);
       percens = val(ceil(percens*length(val)/100));
       for i = 2:nscales
          r(i) = percens(i-1);
          r(i) = max(r(i-1)*1.1,r(i));
       end
       %r1 = prctile(val,1/factor^nscales*100);
       dim = zeros(nscales,1);
       %dim(1) = dim_pow(val,1e-5,r1,r(1));
       [dim(1),c] =  corr_dim(D,r1,r(1));
        s = zeros(nscales,1);
      

        for j = 1:nscales
           mat = max(0,sign(r(j)-D));
           s(j) = sum(mat(:));
        end
        clear mat;
        for j=1:nscales
            Cr(j) = (1 / (n * (n - 1))) * s(j);
        end
        for j=1:nscales-1
            rho(j) = Cr(j+1)-Cr(j);
        end
        for j=1:nscales-1
            a = rho(j);
            b = r(j);
            c = max(c,1e-6);
            fun = @(x)a-c*x*b^(x-1);
            dim(j+1) = fzero(fun,1);
        end           
end

function [no_dims,c] = corr_dim(D,r1,r2)
  n = size(D, 1);
           
            s1 = 0; s2 = 0;
       
            for i=1:n         
                dist = D(:,i);
                s1 = s1 + length(find(dist < r1));
                s2 = s2 + length(find(dist < r2));
            end
            Cr1 = (1 / (n * (n - 1))) * s1;
            Cr2 = (1 / (n * (n - 1))) * s2;

            % Estimate intrinsic dimensionality
            no_dims = (log(Cr2) - log(Cr1)) / (log(r2) - log(r1));
            c = [r1^no_dims; r2^no_dims]\[Cr1;Cr2];
end      
  