function [y_grads,Q] = unsup_grad(ydata,P)
           n = size(ydata,1);
        num = 1./(1+distance(ydata));
        num(1:n+1:end) = 0;          
        Q = max(num ./ sum(num(:)), realmin); 
        L = (P - Q) .* num;
        y_grads = 4*(sum(L, 1))'.*ydata;
        y_grads = y_grads - 4 * L *ydata;
         
%           D = distance(ydata);
%           num = 1./(1+D);
%           num(1:n+1:end) = 0;
%           L = P.*num;
%           y_grads = 2*(sum(L,1))'.* ydata - 2*L*ydata;
%           L = ((1-P).*num) ./D;
%           L(1:n+1:end)=0;
%           y_grads = y_grads - 2*(sum(L,1))'.* ydata + 2*L*ydata;
%           Q = max(num, realmin); 
 end