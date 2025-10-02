function  hebingthreestate()

% k_20=5;%ON到OFF0的切换率5
% k_01=20/11;%OFF0到OFF1的切换率0.8 20/11
% k_12=20/11;%OFF1到ON的切换率 （OFF停留时间为2.5min）0.8 20/11
% r=8;%5生成率
% T=5.545;%5.545
k_20=5;%ON到OFF0的切换率
k_01=20/11;%OFF0到OFF1的切换率 20/11 0.8
k_12=20/11;%OFF1到ON的切换率 （OFF停留时间为2.5min） 20/11  0.8
r=5;
L=0.4436;
v=0.08;
% T=L/v;
T=1;
m1=(k_01*k_12*r*T)/(k_20*(k_01+k_12)+k_01*k_12)
doff1=0; doff2=0; don=1;
initiation = 0; %0;3
time=0; %0; 5.5
num=200;
count=0;
N=1000;

for T= 0.01:0.1:1 %5.5:0.2:10   %0.05:0.2:5.5
    count=count+1;
for j=1:N
    store=zeros(num+1,2);
    time=0; %0; 5.5
    initiation=0; %0; 3
    don=1; doff1=0; doff2=0;
for i=2:num+1
    h=[don doff1 doff2 don];
    c=[k_20 k_01 k_12 r];
    a1=h.*c;
    react=sum(a1);
    rand_num = rand(2,1);
    tau = -log(rand_num(1))/react;
    tauu(i,:) = tau;
if  rand_num(2)<=k_20*don/react
    don=0;doff1=1;doff2=0; 
else if k_20*don/react<rand_num(2)<=(k_20*don+k_01*doff1)/react
     don=0;doff1=0;doff2=1; 
    else if (k_20*don+k_01*doff1)/react<rand_num(2)<=(k_20*don+k_01*doff1+k_12*doff2)/react
            don=1; doff1=0; doff2=0;
        else
            initiation=initiation+1;
            don=1;doff1=0;doff2=0;
        end
    end
end
time=time+tau;
store(i,:)=[initiation time];
end
for k = 1:size(store,1)
  diff(k,:) = abs(T- store(k,2));
  [~, I] = min(diff);
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

% plot(data1(:,2),data1(:,1));
%hold on;
%figure;
%plot(dataF(:,2),dataF(:,1));
% hold on;


% k_20=5;%ON到OFF0的切换率5
% k_01=20/11;%OFF0到OFF1的切换率0.8 20/11
% k_12=20/11;%OFF1到ON的切换率 （OFF停留时间为2.5min）0.8 20/11
% r=8;%5生成率
% T=5.545;%5.545
k_20=5;%ON到OFF0的切换率
k_01=20/11;%OFF0到OFF1的切换率 20/11 0.8
k_12=20/11;%OFF1到ON的切换率 （OFF停留时间为2.5min） 20/11  0.8
r=5;
L=0.4436;
v=0.08;
% T=L/v;
T=1;

m1=(k_01*k_12*r*T)/(k_20*(k_01+k_12)+k_01*k_12);
doff1=0; doff2=0; don=1;
initiation = 0; %0;3
time=0; %0; 5.5
num=200;
count=0;
N=1000;


for t1= 0.01:0.1:9 %5.5:0.2:10   %0.05:0.2:5.5
    count=count+1;
for j=1:N
    store=zeros(num+1,2);
    time=0; %0; 5.5
    initiation=0; %0; 3
    don=1; doff1=0; doff2=0;
for i=2:num+1
    h=[don doff1 doff2 don];
    c=[k_20 k_01 k_12 r];
    a1=h.*c;
    react=sum(a1); rand_num = rand(2,1);
    tau = -log(rand_num(1))/react;
    tauu(i,:) = tau;
if  rand_num(2)<=k_20*don/react
    don=0;doff1=1;doff2=0; 
else if k_20*don/react<rand_num(2)<=(k_20*don+k_01*doff1)/react
     don=0;doff1=0;doff2=1; 
    else if (k_20*don+k_01*doff1)/react<rand_num(2)<=(k_20*don+k_01*doff1+k_12*doff2)/react
            don=1; doff1=0; doff2=0;
        else
            initiation=initiation+1;
            don=1;doff1=0;doff2=0;
        end
    end
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

% plot(data(:,2),data(:,1));
%hold on;
%figure;
%plot(dataF2(:,2),dataF2(:,1));
% hold on;

t_x=[data1(:,2);
    data(:,2)];
t_y=[data1(:,1);
    data(:,1)];
plot(t_x,t_y);

end
