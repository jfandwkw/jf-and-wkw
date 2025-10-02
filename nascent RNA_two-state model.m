function  pp3=hebingonoff()
%UNTITLED3 姝ゅ鏄剧ず鏈夊叧姝ゅ嚱鏁扮殑鎽樿
%   姝ゅ鏄剧ず璇︾粏璇存槑
%data1 = on_off1();
%data = on_off2();



% k_01=0.435; %OFF鍒癘N鐨勫垏鎹㈢巼
% k_10=5;%ON鍒癘FF鐨勫垏鎹㈢巼
% r=5;%杞綍鐜囷紝鐢熸垚鐜?
% L=0.4436;%鍩哄洜闀垮害
% v=0.08;%鑱氬悎閰剁Щ鍔ㄩ?熺巼
% T=L/v;
k_10=5;%ON到OFF的切换率
k_01=10/11;%OFF到ON的切换率0.435
r=5;%5生成率
L=0.4436;%基因长度
v=0.08;%聚合酶移动速率
% T=L/v;
T=1;

m1=(k_01*r*T)/(k_01+k_10);%绋虫?佺殑鍧囧??
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
end %for i鐨勫惊鐜粨鏉?
for k = 1:size(store,1)
  diff(k,:) = abs(T- store(k,2));
  [~, I] = min(diff);%浣嶇疆
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



% k_01=0.435; %OFF鍒癘N鐨勫垏鎹㈢巼
% k_10=5;%ON鍒癘FF鐨勫垏鎹㈢巼
% r=5;%杞綍鐜囷紝鐢熸垚鐜?
% L=0.4436;%鍩哄洜闀垮害
% v=0.08;%鑱氬悎閰剁Щ鍔ㄩ?熺巼
% T=L/v;
k_10=5;%ON到OFF的切换率
k_01=10/11;%OFF到ON的切换率0.435
r=5;%5生成率
L=0.4436;%基因长度
v=0.08;%聚合酶移动速率
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

