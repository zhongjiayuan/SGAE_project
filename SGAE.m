
clear;
[mydata,gene_names] = xlsread('MHC-to-HCC_data.xlsx');
mydata=log(mydata+1);
data_size=size(mydata);
gene_num=data_size(1);
boxsize = 0.1;
cell_num=[54,70,41,65,70,77,70];
time_point=7;
nn=1;
sce_matrix=zeros(gene_num,sum(cell_num));
all_count=0;

for kk=1:time_point
    data=mydata(:,nn:sum(cell_num(1:kk)));
    nn=1+sum(cell_num(1:kk));
    [m,n] = size(data);
    upper = zeros(m,n);
    lower = zeros(m,n);
    for i = 1 : m
        [s1,s2] = sort(data(i,:));
        n0 = n - sum(sign(s1));
        h = round(boxsize*(n - n0)/2);
        k = 1;
        while k <= n
            s = sum(s1 == s1(k)) - 1;
            if s >= h
                upper(i,s2(k : k + s)) = s1(k);
                lower(i,s2(k : k + s)) = s1(k);
            else
               upper(i,s2(k : k + s)) = s1(min(n, k + s + h));
               lower(i,s2(k : k + s)) = s1(max(n0*(n0 > h) + 1, k - h));
            end
            k = k + s + 1;
        end
    end
    identity=ones(m,1);
    for k = 1 : n
        tic
        all_count=all_count+1;
        b = bsxfun(@and,(bsxfun(@le,data,upper(:,k)) & ...
            bsxfun(@ge,data,lower(:,k))),data(:,k) > 0);
        a = sum(b,2);
        b = sparse(double(b));
        b = b*b';
        mi_entropy= b/n.*log(n*b./(a*a' + eps) + eps);
        mi_entropy(mi_entropy < 0) = 0;
       count=full(sum(mi_entropy~=0));
       h_entropy= -(b/n).*log(b./(a*identity' + eps) + eps);
       h_entropy(h_entropy < 0) = 0;
       h_entropy(find(mi_entropy==0))=0;
       entropy=mi_entropy+ h_entropy;
       local_net_entropy=sum(full(entropy),2)./count';
       local_net_entropy=fillmissing(local_net_entropy,'constant',0);
       sce_matrix(:,all_count)=local_net_entropy;
       disp(['Cell ' num2str(k) ' is completed']);
       toc
    end
end
save sce_matrix.mat;

