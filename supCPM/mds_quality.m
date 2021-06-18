function [MDS,knc,cpd,D,D2] = mds_quality(X,label,knn_class)
    center_high = zeros(length(unique(label)),size(X,2));
    n = 1;
    for i = unique(label)'
        cl = find(label==i);
        center_high(n,:)  = mean(X(cl,:));
        n = n+1;
    end
    MDS = mds(center_high,2,0);
    D  = pdist2(center_high,center_high);
    D2 = pdist2(MDS,MDS);
   nbrs1 = knnsearch(center_high,center_high,'k',knn_class+1);
   nbrs2 = knnsearch(MDS,MDS,'k',knn_class+1);
   
   intersection = 0;
   C = length(unique(label)');
   for i = 1:C
       intersection = intersection + length(intersect(nbrs1(i,2:end),nbrs2(i,2:end)));
   end
   knc = intersection / C / knn_class;
   
   d1 = pdist(center_high);
   d2 = pdist(MDS);
   cpd =  corr(d1',d2','Type','Spearman');
end
    