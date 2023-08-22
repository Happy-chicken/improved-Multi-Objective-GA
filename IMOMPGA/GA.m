%% ��ջ���
clc;clear
close all
%% ��������
load scheduleData Jm T JmNumber
%���� ʱ��
tic
%% ��������
NIND=50;        %������Ŀ
MAXGEN=50;      %����Ŵ�����
GGAP=0.9;       %����
XOVR=0.8;       %������
MUTR=0.6;       %������
gen=0;          %��������
%PNumber �������� MNumber  �������
[PNumber, MNumber]=size(Jm);  
trace=zeros(2, MAXGEN);      %Ѱ�Ž���ĳ�ʼֵ
WNumber=PNumber*MNumber;     %�����ܸ���

%% ��ʼ��
Number=zeros(1,PNumber);     % PNumber ��������
for i=1:PNumber
    Number(i)=MNumber;         %MNumber�������
end

% ����2�㣬��һ�㹤�򣬵ڶ������
Chrom=zeros(NIND,2*WNumber);
for j=1:NIND
    WPNumberTemp=Number;
    for i=1:WNumber
        
        %������ɹ���
        val=unidrnd(PNumber);
        while WPNumberTemp(val)==0
            val=unidrnd(PNumber);
        end
        
        %��һ������ʾ����
        Chrom(j,i)= val;
        WPNumberTemp(val)=WPNumberTemp(val)-1;
        
        %��2������ʾ����
        Temp=Jm{val,MNumber-WPNumberTemp(val)};
        SizeTemp=length(Temp);
        %������ɹ������
        Chrom(j,i+WNumber)= unidrnd(SizeTemp);
        
    end
end
 
%����Ŀ�꺯��ֵ
[~, ObjV, ~, ~]=cal(Chrom,JmNumber,T,Jm);  

%% ѭ��Ѱ��

while gen<MAXGEN
    
    %������Ӧ��ֵ
    FitnV=ranking(ObjV);  
    %ѡ�����
    SelCh=select('rws', Chrom, FitnV, GGAP);       
    %�������
    SelCh=across(SelCh,XOVR,Jm,T);          
    %�������
    SelCh=aberranceJm(SelCh,MUTR,Jm,T);            
    
    %����Ŀ����Ӧ��ֵ
    [PVal, ObjVSel, P, S]=cal(SelCh,JmNumber,T,Jm);   
    %���²�������Ⱥ
    [Chrom, ObjV] =reins(Chrom, SelCh,1, 1, ObjV, ObjVSel);       
    %������������
    gen=gen+1;       
    
    %��������ֵ
    trace(1, gen)=min(ObjV);       
    trace(2, gen)=mean(ObjV);  
    
    % ��¼���ֵ
    if gen==1
        Val1=PVal;
        Val2=P;
        MinVal=min(ObjV);%��Сʱ��
        STemp=S;
    end
    %��¼ ��С�Ĺ���
    if MinVal> trace(1,gen)
        Val1=PVal;
        Val2=P;
        MinVal=trace(1,gen);
        STemp=S;
    end
    
end

% ��ǰ���ֵ
PVal=Val1; %����ʱ��
P=Val2;  %���� 
S=STemp; %���Ȼ��򺬻�������
toc
%% ����ı仯
figure(1)
plot(trace(1,:),'lineWidth',1.5,'marker','^');
hold on;
plot(trace(2,:),'-.','lineWidth',1.5,'marker','d');
grid;
legend('Changes in solutions','Changes in population mean');

xlabel('GA generations')
ylabel('Change in optimal solution of GA')
title('Evolutionary Process ')

%% ��ʾ���Ž�
figure(2);
MP=S(1,PNumber*MNumber+1:PNumber*MNumber*2);
c_map = [0.57, 0.69, 0.30
         0.89, 0.88, 0.57
         0.76, 0.49, 0.58
         0.47, 0.76, 0.81
         0.07, 0.35, 0.40
         0.21, 0.21, 0.35
         0.28, 0.57, 0.54
         0.60, 0.24, 0.18
         0.41, 0.20, 0.42
         0.76, 0.84, 0.65];
finish_time=zeros(1,6);
for i=1:WNumber  
    val= P(1,i);
    a=(mod(val,100)); %����
    b=((val-a)/100); %����
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
xlabel('Time')
ylabel('Machine')
title('Gantt Chart')
%% 
% fuzzy_mat=[40, 55, 60, 70;
%            30, 40, 50, 55;
%            65, 80, 90, 110;
%            45, 50, 65, 75;
%            60, 70, 80, 95;
%            50, 55, 70, 80];
% content = calcontent(finish_time,fuzzy_mat);