function [cv,knc,rho,fisher,knn]= embedding_quality(X,Y,label,knc_class,knn_class,subsetsize)  
   covar_high = zeros(1,length(unique(label)));
   covar_low = zeros(1,length(unique(label)));
   k = 1;
   for i = unique(label)'
      cl = find(label==i);
      covar_high(k) = trace(cov(Y(cl,:)));
      covar_low(k) = trace(cov(X(cl,:)));
      k = k+1;
   end
   cv = corr(covar_high',covar_low','Type','Spearman');
   cl = unique(label);
   C = size(cl,1);
   mu1 = zeros(C,size(X,2));
   mu2 = zeros(C,size(Y,2));
   for c = 1:C
       mu1(c,:) = mean(X(label==cl(c),:));
       mu2(c,:) = mean(Y(label==cl(c),:));
   end
   
   nbrs1 = knnsearch(mu1,mu1,'k',knc_class+1);
   nbrs2 = knnsearch(mu2,mu2,'k',knc_class+1);
   
   intersection = 0;
   
   for i = 1:C
       intersection = intersection + length(intersect(nbrs1(i,2:end),nbrs2(i,2:end)));
   end
   knc = intersection / C / knc_class;
   
   randinx = randperm(size(X,1));
   subset = randinx(1:subsetsize);
   d1 = pdist(X(subset,:));
   d2 = pdist(Y(subset,:));
   rho =  corr(d1',d2','Type','Spearman');
   
   number = zeros(C,1);
   for c = 1:C
       number(c,:) = length(find(label==c));
   end
   
   n = length(label);
   L_tilde = zeros(n);
   for i = unique(label)'
           cl = find(label==i);
           L_tilde(cl,cl)=1/length(cl);
   end
   I = eye(n);
   one = ones(n,1);  
   Sb = (L_tilde-one*one'/n);
   Sw = (I - L_tilde);   
   fisher = trace(X'*Sw*X)/trace(X'*Sb*X); 

   nbrs1 = knnsearch(X,X,'k',knn_class+1);
   nbrs2 = knnsearch(Y,Y,'k',knn_class+1);
   intersection = 0;
   for i = 1:n
       intersection = intersection + length(intersect(nbrs1(i,:),nbrs2(i,:)));
   end
   knn = intersection / (n*knn_class);
end