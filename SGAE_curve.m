load sce_matrix.mat;

subplot(2,1,1)
 for i=1:all_count
    [sorted_H,idx]=sort(sce_matrix(:,i),'descend');
    SH(i)=sum(sorted_H(1:gene_num*0.05));
 end
result(1)=mean(SH(1:54));
result(2)=mean(SH(:,55:124));
result(3)=mean(SH(:,125:165));
result(4)=mean(SH(:,166:230));
result(5)=mean(SH(:,231:300));
result(6)=mean(SH(:,301:377));
result(7)=mean(SH(:,378:447));
combineData= [SH(1:54),SH(55:124),SH(125:165),SH(166:230),SH(231:300),SH(301:377),SH(378:447)];        % 组合 
group = [2*ones(1,54),3*ones(1,70),4*ones(1,41),5*ones(1,65),6*ones(1,70),7*ones(1,77),8*ones(1,70)];  % 给每一个对应的值，设定标签
boxplot(combineData,group);
hold on;
profile=mydata;
cell_num=[54,70,41,65,70,77,70];
count=0;
profile_size=size(profile);
for i=1:profile_size(2)
    [sorted_H,idx]=sort(profile(:,i),'descend');
    exp_SH(i)=mean(sorted_H(1:profile_size(1)*0.05));
end
exp_result(1)=mean(exp_SH(1:54));
exp_result(2)=mean(exp_SH(:,55:124));
exp_result(3)=mean(exp_SH(:,125:165));
exp_result(4)=mean(exp_SH(:,166:230));
exp_result(5)=mean(exp_SH(:,231:300));
exp_result(6)=mean(exp_SH(:,301:377));
exp_result(7)=mean(exp_SH(:,378:447));
t=1:7;
y1=result;
plot(t,y1);
set(gca,'ylim',[24,36]);
set(gca,'yTick',[24:2:36]) ;
ylabel('Entropy');
box off
subplot(2,1,2)
profile=mydata;
profile_size=size(profile);
for i=1:profile_size(2)
    [sorted_H,idx]=sort(profile(:,i),'descend');
    exp_SH(i)=mean(sorted_H(1:profile_size(1)*0.05));
end
y2=exp_result;
combineData = [exp_SH(1:54),exp_SH(55:124),exp_SH(125:165),exp_SH(166:230),exp_SH(231:300),exp_SH(301:377),exp_SH(378:447)];        % 组合 
group = [2*ones(1,54),3*ones(1,70),4*ones(1,41),5*ones(1,65),6*ones(1,70),7*ones(1,77),8*ones(1,70)];  % 给每一个对应的值，设定标签
boxplot(combineData,group,'notch','on','sym','mo','Colors',[0.4 0.81 0.66]);
hold on;
plot(t,y2);
ylabel('Gene expression');
set(gca,'ylim',[6,6.8]);
set(gca,'yTick',[6:0.2:6.8]) ;
B={' E10.5' ' E11.5' ' E12.5' ' E13.5' ' E14.5' ' E15.5' ' E17.5'};
set(gca,'XTickLabel',B);
box off








