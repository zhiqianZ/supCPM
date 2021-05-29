function [Z U] = SPCA(X, Y, d, param)
%Copyright Barshan, Ghodsi 2009
%Paper: Supervised principal component analysis: Visualization, classification and
%regression on subspaces and submanifolds.

% [Z U] = SPCA(X,Y,d,param)
% Input:
%       X:  explanatory variable (pxn)
%       Y:  response variables (lxn)
%       d:  dimension of effective subspaces
%       param: 
%             param.ktype_y : kernel type of the response variable
%             param.kparam_y : kernel parameter of the response variable

% Output:
%       Z:  dimension reduced data (dxn)
%       U:  orthogonal projection matrix (pxd)

if size(X,2)~=size(Y,2)
    error('X and Y must be the same length')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computing Kernel Function of Labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[l,n] = size(Y);
L = repmat(0,n,n);
%making L full rank for classification
if strcmp(param.ktype_y,'delta_cls')
     L = L+eye(n);
     param.ktype_y = 'delta';
end
for i = 1:n
    for j = 1:n
        L(i,j) = L(i,j) + kernel(param.ktype_y,Y(:,i),Y(:,j),param.kparam_y,[]);
    end
end



H = eye(n)-1/n*(ones(n,n));
[p,n] = size(X);
if n>p
    tmp = X*H*L*H*X';
    [U D] = eigendec(tmp,d,'LM');
else
   [u s v] = svd(L);
   phi_Y = s^.5 * v';
   tmp = phi_Y*H*X'*X*H*phi_Y';
   [V D] = eigendec(tmp,d,'LM');
   U = X*H*phi_Y'*V*inv(diag(D)^.5);
end
Z = U'*X;


