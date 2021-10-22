function [grad,TraceSw,TraceSb] = sup_grad(ydata,label,cluster_num,cluster_size)
          [n,no_dims] = size(ydata);
          TraceSw = 0;
          endid = cumsum(cluster_size);
          startid = endid-cluster_size+1;
          for i = 1:cluster_num
              ypart = ydata(startid(i):endid(i),:);
              TraceSw = TraceSw + sum(sum(ypart).^2)/cluster_size(i);
          end
          temp = TraceSw;
          TraceSw = TraceSw - sum(sum(ydata).^2)/n;
          
          TraceSb = trace(ydata'* ydata) - temp;
          y_mean = mean(ydata);
          colsum = zeros(cluster_num,no_dims);
          for i = 1:no_dims
             colsum(:,i) = accumarray(label,ydata(:,i));
          end
          for i = 1:cluster_num
              colsum(i,:) = colsum(i,:)/cluster_size(i);
          end
          rep = repelem(1,no_dims);
          temp = repelem(colsum,cluster_size,rep);
          yb = ydata - temp;
          yw = temp - y_mean;
          grad = (TraceSw * yb- TraceSb * yw)/TraceSw^2;
end