function D = approx_geodesic(W)
% fast geodesic distance matrix computation based on Nystrom method
%
 W = max(W,W');
 n = size(W,1);
 G = graph(W);
 select = 50;
 randindex = randperm(n);
 test_set = randindex(end-10:end);
 test_column = zeros(length(test_set),n)';
 for i=1:length(test_set)
     [T,dist] = shortestpathtree(G,test_set(i),'all','OutputForm','cell');
     test_column(:,i) = dist';
 end
 for draw = 1:floor(n/select)
     for i = (1+(draw-1)*select) : select*draw
        [~,dist] = shortestpathtree(G,randindex(i),'all','OutputForm','cell');
         sample_column(:,i) = dist';
     end
     Y = nystrom(sample_column.^2,randindex(1:i),test_set).^0.5;
     if norm(Y-test_column)<3e-2*norm(test_column)
         break;
         draw
     end
     %temp_st = find(min([sample_column';sample_column'])> prctile(min([sample_column';sample_column']),0));
   %  kk = randperm(length(temp_st));
   %  set = [set,temp_st(kk(1))];
 end
 D = nystrom(sample_column.^2,randindex(1:i),1:n).^0.5;
 D(W~=0)=W(W~=0);
end

function Y = nystrom(X,set,test_set)
n = length(set);
[U S V] = svd(X(set,:)');
n0 = round(n*1/2);
lambda = S(n0,n0);
S(S<=lambda)=0;
Y = X*pinv(V*S*U')*X(test_set,:)';
Y = max(0,Y);
end
 