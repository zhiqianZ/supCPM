function ydata = supCPM_downsample(alpha,beta,data,label,no_dims,compel_force,geodist,degree,ratio,k,change,niter,seed,factor)
% downsample
d = pdist2(data,data);
[~,I]=sort(d(:),'ascend');
[n,~] = size(data);
I(1:n) = []; % remove the diag
rind = I(ceil(length(I)*alpha));
r = d(rind);

merge = zeros(n,1);
indtable= zeros(n,2);
k=1;
cl_label = unique(label);
large_cl = [];
for i = cl_label
   if(length(find(label==i))>n*beta)
      large_cl= [large_cl,i];
   end
end
for i = large_cl
    cl = find(label==i);
    for j = 1:length(cl)
        if merge(cl(j))==0
            smallind = find(d(cl(j),cl)<r);
            smallind = intersect(smallind,find(merge(cl)==0));
            if length(smallind)>1
                merge(cl(smallind))=k;
                indtable(k,:) = [cl(j),k];
                k = k+1;
            end
        end
    end
end
indtable(k:end,:)=[];
ind = indtable(:,1);
ind = [ind;find(merge==0)];
data_part = data(ind,:);
label_part = label(ind);

% run supCPM
ydata_part = supCPM(data,label,no_dims,compel_force,geodist,degree,ratio,k,change,niter,seed,factor);

% recover removed points
n = size(data_part,1);
ydata = zeros(n,2);
ind = indtable(:,1);
ind = [ind;find(merge==0)];
D = pdist2(ydata_part,ydata_part);
rsmall = min(D(D~=0));

ydata(ind,:) = ydata_part;
for i = 1:size(indtable,1) %number of merge clusters
    cl = find(merge==i);
    center = indtable(i,1);
    ydata(cl,:) = repmat(ydata(center,:),[length(cl),1])...
                                  + rand(length(cl),2) * sqrt(rsmall/2);
end
end