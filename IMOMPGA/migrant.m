function [Chrom,ObjV]=migrant(Chrom,ObjV)
%% 移民算子
MP=length(Chrom);
% fprintf('达到移民周期，正在进行种群移民...\n');
for i=1:MP
    [~,minI]=min(ObjV{i});  % 找出第i种群中最优的个体
    next_i=i+1;                % 目标种群（移民操作中）
    if next_i>MP
        next_i=mod(next_i,MP);
    end
    if mean(ObjV{next_i}) < mean(ObjV{i}) 
        [~,maxI]=max(ObjV{next_i});          %  找出目标种群中最劣的个体
        %% 目标种群最劣个体替换为源种群最优个体
        Chrom{next_i}(maxI,:)=Chrom{i}(minI,:);
        ObjV{next_i}(maxI)=ObjV{i}(minI);
    else
        a = randi([1, 50], 1, 2); % 生成 2 个随机数,等待移民交流
%         Chrom{next_i}(a,:)=Chrom{i}(a,:);
%         ObjV{next_i}(a)=ObjV{i}(a);
        Chrom{next_i}=[Chrom{next_i}; Chrom{i}(a,:)];
        ObjV{next_i}=[ObjV{next_i};ObjV{i}(a)];
    end
    
end

