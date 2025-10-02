function  hebingcrosstalking()

% k_10=5; %on到OFF的切换率 5
% q_1=0.1; %第一条路径打开的概率0.5 强路径0.2
% q_2=1-q_1;%第二条路径打开的概率
% lamdba_1=0.2;%第一条路径打开时OFF到ON的切换率0.4 强路径0.2
% lamdba_2=8; %0.4 第一条路径打开时OFF到ON的切换率0.4 强路径8
% r=5;%5生成率
% T=5.545;%5.545
k_10=5; %on到OFF的切换率
q_1=0.2; %第一条路径打开的概率0.5  强路径0.2
q_2=1-q_1;%第二条路径打开的概率
lamdba_1=0.2;%第一条路径打开时OFF到ON的切换率0.4 强路径0.2
lamdba_2=8; %0.4 第一条路径打开时OFF到ON的切换率0.4  强路径8
r=5;
L=0.4436;
v=0.08;
% T=L/v;
T=1;

m1=(lamdba_1*lamdba_2*r*T)/((lamdba_1*q_2+lamdba_2*q_1)*k_10+lamdba_1*lamdba_2);
% don=1;doff1=0;doff2=0;

initiation = 0; %0;3
time=0; %0; 5.5
num=200;
count=0;
N=1000;

for T= 0.01:0.1:1 %5.545
    count=count+1;
    for j=1:N
        store=zeros(num+1,2);
        time=0; %0; 5.5
        initiation=0; %0; 3
        don=1;
        if rand<=q_1 %  初始状态
            doff1=1;doff2=0;
        else
            doff1=0;doff2=1;
        end
        % don=0;
        %don=1;doff1=0;doff2=0;
        for i=2:num+1
            h=[don don doff1 doff2 don];
            c=[q_1*k_10 q_2*k_10 lamdba_1 lamdba_2 r];
            a1=h.*c;
            react=sum(a1);
            %if don==1
            rand_num = rand(2,1);
            tau = -log(rand_num(1))/react;
            tauu(i,:) = tau;
            if rand_num(2)<=q_1*k_10*don/react
                don=0;doff1=1;doff2=0;
            else if rand_num(2)<=(q_1*k_10*don+q_2*k_10*don)/react
                    don=0;doff1=0;doff2=1;
            else if rand_num(2)<=(q_1*k_10*don+q_2*k_10*don+lamdba_1*doff1)/react
                    don=1;doff1=0;doff2=0;
            else if rand_num(2)<=(q_1*k_10*don+q_2*k_10*don+lamdba_1*doff1+lamdba_2*doff2)/react
                    don=1;doff1=0;doff2=0;
            else
                initiation=initiation+1;
                don=1;doff1=0;doff2=0;
            end
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
% %hold on;
% %figure;
% %plot(dataF(:,2),dataF(:,1));
% hold on;



% k_10=5; %on到OFF的切换率 5
% q_1=0.1; %第一条路径打开的概率0.5 强路径0.2
% q_2=1-q_1;%第二条路径打开的概率
% lamdba_1=0.2;%第一条路径打开时OFF到ON的切换率0.4 强路径0.2
% lamdba_2=8; %0.4 第一条路径打开时OFF到ON的切换率0.4 强路径8
% r=5;%5生成率
% T=5.545;%5.545
k_10=5; %on到OFF的切换率
q_1=0.2; %第一条路径打开的概率0.5  强路径0.2
q_2=1-q_1;%第二条路径打开的概率
lamdba_1=0.2;%第一条路径打开时OFF到ON的切换率0.4 强路径0.2
lamdba_2=8; %0.4 第一条路径打开时OFF到ON的切换率0.4  强路径8
r=5;
L=0.4436;
v=0.08;
% T=L/v;
T=1;

m1=(lamdba_1*lamdba_2*r*T)/((lamdba_1*q_2+lamdba_2*q_1)*k_10+lamdba_1*lamdba_2);
% don=1;;doff1=0;,doff2=0;
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
        don=1;
        if rand<=q_1 %  初始状态
            doff1=1;doff2=0;
        else
            doff1=0;doff2=1;
        end
        % don=0;
        %don=1; doff1=0; doff2=0;
        for i=2:num+1
            h=[don don doff1 doff2 don];
            c=[q_1*k_10 q_2*k_10 lamdba_1 lamdba_2 r];
            a1=h.*c;
            react=sum(a1);
            rand_num = rand(2,1);
            tau = -log(rand_num(1))/react;
            tauu(i,:) = tau;
            if rand_num(2)<=q_1*k_10*don/react
                don=0;doff1=1;doff2=0;
            else if rand_num(2)<=(q_1*k_10*don+q_2*k_10*don)/react
                    don=0;doff1=0;doff2=1;
            else if rand_num(2)<=(q_1*k_10*don+q_2*k_10*don+lamdba_1*doff1)/react
                    don=1;doff1=0;doff2=0;
            else if rand_num(2)<=(q_1*k_10*don+q_2*k_10*don+lamdba_1*doff1+lamdba_2*doff2)/react
                    don=1;doff1=0;doff2=0;
            else
                initiation=initiation+1;
                don=1;doff1=0;doff2=0;
            end
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
% %hold on;
% %figure;
% %plot(dataF2(:,2),dataF2(:,1));
% hold on;

t_x=[data1(:,2);
    data(:,2)];
t_y=[data1(:,1);
    data(:,1)];
plot(t_x,t_y);

end

