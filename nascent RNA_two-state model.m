function  pp3=hebingonoff()
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
%data1 = on_off1();
%data = on_off2();



% k_01=0.435; %OFF到ON的切换率
% k_10=5;%ON到OFF的切换率
% r=5;%转录率，生成�?
% L=0.4436;%基因长度
% v=0.08;%聚合酶移动�?�率
% T=L/v;
k_10=5;%ON��OFF���л���
k_01=10/11;%OFF��ON���л���0.435
r=5;%5������
L=0.4436;%���򳤶�
v=0.08;%�ۺ�ø�ƶ�����
% T=L/v;
T=1;

m1=(k_01*r*T)/(k_01+k_10);%稳�?�的均�??
state=1;
initiation = 0; 
time=0; 
num=100;
count=0;
N=1000;

for T= 0.01:0.1:1 
    count=count+1;
for j=1:N
    store=zeros(num+1,2);
    time=0; 
    initiation=0; 
    state=1;
for i=2:num+1
react=(k_10+r)*state-k_01*(state-1);
if state==1
rand_num = rand(2,1);
tau = -log(rand_num(1))/react;
tauu(i,:) = tau;
  if rand_num(2)<k_10/react
    state=0;
  else
    initiation=initiation+1;
  end
else
    rand_num = rand(2,1);
    tau = -log(rand_num(1))/react;
    state=1;
end
time=time+tau;
store(i,:)=[initiation time];
end %for i的循环结�?
for k = 1:size(store,1)
  diff(k,:) = abs(T- store(k,2));
  [~, I] = min(diff);%位置
end

if store(I,2)<T
    mRNA=store(I,1);
else
    mRNA=store(I-1,1);
end
nascent(j,:)=mRNA;
end
data1(count,:)=[mean(nascent) T];
dataF(count,:)=[var(nascent)/mean(nascent) T];
end



% k_01=0.435; %OFF到ON的切换率
% k_10=5;%ON到OFF的切换率
% r=5;%转录率，生成�?
% L=0.4436;%基因长度
% v=0.08;%聚合酶移动�?�率
% T=L/v;
k_10=5;%ON��OFF���л���
k_01=10/11;%OFF��ON���л���0.435
r=5;%5������
L=0.4436;%���򳤶�
v=0.08;%�ۺ�ø�ƶ�����
% T=L/v;
T=1;
m1=(k_01*r*T)/(k_01+k_10)
state=1;
initiation=0;
time=0;
num=100;
count=0;
N=1000;

for t1=0.01:0.1:9
    count=count+1;
for j=1:N
    store=zeros(num+1,2);
    time=0;initiation=0;
    state=1;
for i=2:num+1
react=(k_10+r)*state-k_01*(state-1);

if state==1
rand_num = rand(2,1);
tau = -log(rand_num(1))/react;
if rand_num(2)<k_10/react
    state=0;
else
    initiation=initiation+1;
end
else
    rand_num = rand(2,1);
    tau = -log(rand_num(1))/react;
    state=1;
end
time=time+tau;
store(i,:)=[initiation time];
end
diff=zeros(size(store,1),1);diff1=zeros(size(store,1),1);
for k = 1:size(store,1)
  diff(k,:) = abs(t1- store(k,2));%On a sliding time window, this defines the initial number of nascent
  diff1(k,:) = abs(T+t1- store(k,2));
  [~, I1] = min(diff);
   [~, I2] = min(diff1);
end

if store(I1,2)<t1
    mRNA=store(I1,1);
else
    mRNA=store(I1-1,1);
end

if store(I2,2)<T+t1
    mRNA1=store(I2,1);
else
    mRNA1=store(I2-1,1);
end
nascent(j,:)=mRNA1-mRNA;
end


data(count,:)=[mean(nascent) T+t1];
dataF2(count,:)=[var(nascent)/mean(nascent) T+t1];
end

t_x=[data1(:,2);
    data(:,2)];
t_y=[data1(:,1);
    data(:,1)];
plot(t_x,t_y);

end

