%% 清空环境
clc;clear
close all
%% 下载数据
load scheduleData Jm T JmNumber
%工序 时间
tic
%% 基本参数
NIND=50;        %个体数目
MAXGEN=100;      %最大遗传代数
MP=10; %种群数目
GGAP=0.7;       %代沟
XOVR=0.6 + (0.9-0.6)*rand(MP,1);       %交叉率0.7-0.9
MUTR=0.6 + (0.6-0.1)*rand(MP,1);       %变异率0.1-0.6
gen=0;          %代计数器
gen0 =0;
%PNumber 工件个数 MNumber  工序个数
[PNumber, MNumber]=size(Jm);
trace=zeros(2, MAXGEN);      %寻优结果的初始值
trace_group={};
WNumber=PNumber*MNumber;     %工序总个数

%% 初始化
Number=zeros(1,PNumber);     % PNumber 工件个数
for i=1:PNumber
    Number(i)=MNumber;         %MNumber工序个数
end
for k=1:MP
    % 代码2层，第一层工序，第二层机器
    Chrom=zeros(NIND,2*WNumber);
    for j=1:NIND
        WPNumberTemp=Number;
        for i=1:WNumber
            
            %随机产成工序
            val=unidrnd(PNumber);
            while WPNumberTemp(val)==0
                val=unidrnd(PNumber);
            end
            
            %第一层代码表示工序
            Chrom(j,i)= val;
            WPNumberTemp(val)=WPNumberTemp(val)-1;
            
            %第2层代码表示机器
            Temp=Jm{val,MNumber-WPNumberTemp(val)};
            SizeTemp=length(Temp);
            %随机产成工序机器
            Chrom(j,i+WNumber)= unidrnd(SizeTemp);
            
        end
    end
    Chrom_group{k}=Chrom;
end
%计算目标函数值
for k=1:MP
    [~, ObjV, ~, ~]=cal(Chrom_group{k},JmNumber,T,Jm);
    ObjV_group{k}=ObjV;
end
MinObjV=ones(MP,1) * 100000;           %记录精华种群
[~,choromlen]=size(Chrom);
MinChrom=zeros(MP,choromlen) ; %记录精华种群的编码
minY=10000; % 当代所有中种群的最佳

% % 初始化帕累托前沿和帕累托集合
% pareto_front = fitnesses(1,:);
% pareto_set = cell(1,1);
% pareto_set{1} = population{1};
%% 循环寻找
keep_gen=0;
while gen<MAXGEN
    %代计数器增加
    gen=gen+1;
    for i=1:length(Chrom_group)
        %分配适应度值
        FitnV{i}=ranking(ObjV_group{i});
        %选择操作
        SelCh{i}=select('rws', Chrom_group{i}, FitnV{i}, GGAP);
        %交叉操作
        SelCh{i}=across(SelCh{i},XOVR(i),Jm,T);
        %变异操作
        SelCh{i}=aberranceJm(SelCh{i},MUTR(i),Jm,T);
        %计算目标适应度值
        [PVal, ObjVSel, P, S]=cal(SelCh{i},JmNumber,T,Jm);
        %重新插入新种群
        [Chrom_group{i}, ObjV_group{i}] =reins(Chrom_group{i}, SelCh{i},1, 1, ObjV_group{i}, ObjVSel);
        %保存最优值
        trace(1, gen)=min(ObjV_group{i});
        trace(2, gen)=mean(ObjV_group{i});
        trace_group{i}=trace;
    end
    [Chrom_group, ObjV_group] = migrant(Chrom_group, ObjV_group); % 移民
    [MinObjV,MinChrom]=EliteInduvidual(Chrom_group,ObjV_group,MinObjV,MinChrom);     % 人工选择精华种群
    
    [minval, index]=min(MinObjV);
    YY(gen)=minval;%最小时间
    [PVal, ObjVSel, P, S]=cal(MinChrom(index,:),JmNumber,T,Jm);
    % 记录保持代数
    if YY(gen) == minY
        keep_gen = keep_gen +1;
    end
     %记录 最小的工序
    if YY(gen)<minY   %判断当前优化值是否与前一次优化值相同
        minY=YY(gen); %更新最优值
        Val1=PVal;
        Val2=P;% 工序
        STemp=S; %基因
    end
    % 连续5代不变，提前停止
%     if keep_gen>5
%         break;
%     end
end

% 当前最佳值
PVal=Val1; %工序时间
P=Val2;  %工序
S=STemp; %调度基因含机器基因
toc
%% 描绘解的变化
figure(1)
for i=1:MP-5
    plot(trace_group{i}(1,:),'lineWidth',1.5);
    hold on;
    plot(trace_group{i}(2,:),'-.','lineWidth',1.5,'marker','^');
    hold on
end
xlabel('进化代数')
ylabel('每个种群的最优解变化')
title(['各个种群的最优解变化'])
xlim([1,gen])


%% 
hold on
plot(1:gen,YY,'lineWidth',1.5)
legend('合并种群1解的变化','合并种群1均值的变化','合并种群2解的变化','合并种群2均值的变化','合并种群3解的变化','合并种群3均值的变化',...
'合并种群4解的变化','合并种群4均值的变化','合并种群5解的变化','合并种群5均值的变化','MPGA变化曲线')
xlabel('MPGA进化代数')
ylabel('MPGA最优解变化')
title('进化过程')
grid;
xlim([1,gen])
%% 显示最优解
figure(3);
MP=S(1,PNumber*MNumber+1:PNumber*MNumber*2);
c_map = [0.57, 0.69, 0.30
         0.89, 0.88, 0.57
         0.76, 0.49, 0.58
         0.47, 0.76, 0.81
         0.21, 0.21, 0.35
         0.28, 0.57, 0.54
         0.07, 0.35, 0.40
         0.41, 0.20, 0.42
         0.60, 0.24, 0.18
         0.76, 0.84, 0.65];
finish_time=zeros(1,6);
for i=1:WNumber  
    val= P(1,i);
    a=(mod(val,100)); %工序
    b=((val-a)/100); %工件
    Temp=Jm{b,a};
    mText=Temp(MP(1,i));
    
    x1=PVal(1,i);
    x2=PVal(2,i);
    
    y1=mText-1 + 0.15;
    y2=mText - 0.15;
    plotRec(x1,x2,mText -0.15);
    
    plotRec(PVal(1,i),PVal(2,i),mText-0.15);
    hold on;
    finish_time(b)=x2;
%     [1-1/b,1/b,b/PNumber]
    fill([x1,x2,x2,x1],[y1,y1,y2,y2],c_map(b));
    text((x1+x2)/2,mText-0.25,num2str(P(i)));
end
xlabel('时间')
ylabel('机器')
title('甘特图')
%% 
fuzzy_mat=[40, 55, 60, 70;
           30, 40, 50, 55;
           65, 80, 90, 110;
           45, 50, 65, 75;
           60, 70, 80, 95;
           50, 55, 70, 80];
content = calcontent(finish_time,fuzzy_mat);
for i=1:length(content)
    fprintf('工件%d的顾客满意度为%f\n',i,content(i));
end
fprintf('最少时间%ds\n',YY(end));
%% GA
