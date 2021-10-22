function dim = dim_pow(val,r1,r2,r3)
% This function estimates the intrinsic dimension
% val: the vector that stores the pairwise distances at a certain scale
% r1-r3: three radii at this scale
% Author: Rongrong Wang

            if nargin <= 2
             r1 = prctile(val,5);
             r2 = prctile(val,10);
             r3 =  prctile(val,12);
            end
          
               Dist = val(val>=r1&val<=r3);
               val = sort(Dist);
                N1 = length(val);
                N = 200;
               
            
               [E1,C1]= int(r1,r2,N,N1,val);
               [E2,C2] = int(r1,r3,N,N1,val);
               alpha1 = C1/log(r2/r1);
               alpha2 = C2/log(r3/r1);
              
               alpha = (C2-C1)/C1;
               beta = alpha*log(r2/r1)/log(r3/r2);
               dim = (beta-1)/(log(r2/r1)-E1); %initialization
            
               for i=1:100 
                   dim = min(1/E2,dim);
                   C0 = -(E1-1/dim)*alpha1;
                   C0_a = -(E2-1/dim)*alpha2;
                   if C0 < 1e10 
                       if C0_a<0
                       C0_a = max(C0_a,C2/2);
                       end
                       dim = log((C2+C0_a)/(C1+C0_a))/log(r3/r2);
                   else
                       if C0<0
                          C0=max(C0,C1/2);
                       end
                       dim = log((C1+C0)/C0)/log(r2/r1);
                   end
               end
               
              
              
end
  

function [E,C] = int(r1,r2,N,N1,val)
delta = (r2-r1)/N;
      j = 1;
      E = 0;
            for i = 1:N
                r = i*delta+r1;
                if j == N1
                    break;
                end
                while val(j)<r
                    
                    j=j+1;
                    if j== N1
                    break;
                    end
                end
                if exist('j_ini')
                    E = E+(j-j_ini)/r*delta;
                else
                E = E+j/r*delta;
                end
                if i == 1
                    j_ini = j;
                end
            end
            C = j-j_ini;
            E = E/C;
            
end
