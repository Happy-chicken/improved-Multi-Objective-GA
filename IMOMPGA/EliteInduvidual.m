function [MinObjV,MinChrom]=EliteInduvidual(Chrom,ObjV,MinObjV,MinChrom)
%% 人工选择算子
MP=length(Chrom);  %种群数
for i=1:MP
    [MinO,minI]=min(ObjV{i});   %找出第i种群中最优个体
    if MinO<MinObjV(i)
        MinObjV(i)=MinO;         %记录各种群的精华个体
        MinChrom(i,:)=Chrom{i}(minI,:);  %记录各种群精华个体的编码
    end
end
