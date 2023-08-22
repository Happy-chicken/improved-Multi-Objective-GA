function [PVal, ObjV, P, S]=cal(Chrom,JmNumber,T,Jm)

% 功能说明：       根据基因群,计算出个群中每个个体的调度工序时间，
%                 保存最小时间的调度工序和调度工序时间
% 输入参数：
%       Chrom     为基因种群  
%       T         为各工件各工序使用的时间 
%       Jm        为各工件各工序使用的机器 

% 输出参数:
%       PVal      为最佳调度工序时间 
%       P         为最佳输出的调度工序 
%       ObjV      为群中每个个体的调度工序时间
%       S         为最佳输出的调度基因
%
fuzzy_mat=[40, 55, 60, 70;
           30, 40, 50, 55;
           65, 80, 90, 110;
           45, 50, 65, 75;
           60, 70, 80, 95;
           50, 55, 70, 80];
%初始化
NIND=size(Chrom,1);
ObjV=zeros(NIND,1);
MinVal = 0;
%  工件个数 工序个数 
[PNumber, MNumber]=size(Jm);

for i=1:NIND  
    
    %取一个个体
    S=Chrom(i,:);
    
    %根据基因，计算调度工序
    P= calp(S,PNumber);
    
    %根据调度工序，计算出调度工序时间
    PVal=caltime(S,P,JmNumber,T,Jm); 
        
    %取完成时间
    MT=max(PVal);
    TVal=max(MT);  
    
    for j=1:PNumber* MNumber
        val= P(1,j);
        a=(mod(val,100)); %工序
        b=((val-a)/100); %工件
        finish_time(b)=PVal(2,j);
    end
    content = calcontent(finish_time,fuzzy_mat);
    %保存时间
    if sum(content) == 0
        ObjV(i,1)=TVal + 10000;
    else
        ObjV(i,1)=TVal + 500/sum(content);
    end
%     ObjV(i,1)=TVal;
%     ObjV(i,2)=sum(content);
    
    %初始化
    if i==1
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end
    
    %记录 最小的调度工序时间、最佳调度工序时间 最佳输出的调度工序
    if MinVal> ObjV(i,1)
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end   
end 
 
%最佳调度工序时间 最佳输出的调度工序
 PVal=Val1;
 P=Val2;
 S=STemp;
 
 

 
 
 
 