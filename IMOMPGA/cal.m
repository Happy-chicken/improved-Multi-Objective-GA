function [PVal, ObjV, P, S]=cal(Chrom,JmNumber,T,Jm)

% ����˵����       ���ݻ���Ⱥ,�������Ⱥ��ÿ������ĵ��ȹ���ʱ�䣬
%                 ������Сʱ��ĵ��ȹ���͵��ȹ���ʱ��
% ���������
%       Chrom     Ϊ������Ⱥ  
%       T         Ϊ������������ʹ�õ�ʱ�� 
%       Jm        Ϊ������������ʹ�õĻ��� 

% �������:
%       PVal      Ϊ��ѵ��ȹ���ʱ�� 
%       P         Ϊ�������ĵ��ȹ��� 
%       ObjV      ΪȺ��ÿ������ĵ��ȹ���ʱ��
%       S         Ϊ�������ĵ��Ȼ���
%
fuzzy_mat=[40, 55, 60, 70;
           30, 40, 50, 55;
           65, 80, 90, 110;
           45, 50, 65, 75;
           60, 70, 80, 95;
           50, 55, 70, 80];
%��ʼ��
NIND=size(Chrom,1);
ObjV=zeros(NIND,1);
MinVal = 0;
%  �������� ������� 
[PNumber, MNumber]=size(Jm);

for i=1:NIND  
    
    %ȡһ������
    S=Chrom(i,:);
    
    %���ݻ��򣬼�����ȹ���
    P= calp(S,PNumber);
    
    %���ݵ��ȹ��򣬼�������ȹ���ʱ��
    PVal=caltime(S,P,JmNumber,T,Jm); 
        
    %ȡ���ʱ��
    MT=max(PVal);
    TVal=max(MT);  
    
    for j=1:PNumber* MNumber
        val= P(1,j);
        a=(mod(val,100)); %����
        b=((val-a)/100); %����
        finish_time(b)=PVal(2,j);
    end
    content = calcontent(finish_time,fuzzy_mat);
    %����ʱ��
    if sum(content) == 0
        ObjV(i,1)=TVal + 10000;
    else
        ObjV(i,1)=TVal + 500/sum(content);
    end
%     ObjV(i,1)=TVal;
%     ObjV(i,2)=sum(content);
    
    %��ʼ��
    if i==1
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end
    
    %��¼ ��С�ĵ��ȹ���ʱ�䡢��ѵ��ȹ���ʱ�� �������ĵ��ȹ���
    if MinVal> ObjV(i,1)
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end   
end 
 
%��ѵ��ȹ���ʱ�� �������ĵ��ȹ���
 PVal=Val1;
 P=Val2;
 S=STemp;
 
 

 
 
 
 