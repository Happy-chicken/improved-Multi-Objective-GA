function [Chrom,ObjV] = integrate(Chrom,ObjV)
k=3;
MP=length(Chrom);
% fprintf('达到合并周期，正在进行种群合并...\n');
diversity_group=zeros(MP,1);
for i=1:MP
    [~,C,~,~]=kmeans(Chrom{i},k,'distance','cosine','MaxIter',50,'Replicates',5);
    D = pdist(C, 'cosine'); % 计算余弦距离
    diversity = std(D);
    diversity_group(i)=diversity;
end
[~,minI]=min(diversity_group);
[~,maxI]=max(diversity_group);

integrated_pop=[Chrom{minI};Chrom{maxI}];
[integrated_pop, ~, ~] = unique(integrated_pop, 'rows'); % 去除相同行
load scheduleData Jm T JmNumber
[~, objv, ~, ~]=cal(integrated_pop,JmNumber,T,Jm);
if minI < maxI
    Chrom(minI)=[];
    Chrom(maxI-1)=[];
    Chrom{end+1}=integrated_pop;
    ObjV(minI)=[];
    ObjV(maxI-1)=[];
    ObjV{end+1}=objv;
elseif minI > maxI
    Chrom(minI-1)=[];
    Chrom(maxI)=[];
    Chrom{end+1}=integrated_pop;
    ObjV(minI-1)=[];
    ObjV(maxI)=[];
    ObjV{end+1}=objv;
end

end

