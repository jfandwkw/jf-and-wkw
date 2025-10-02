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
Delta=(lamdba_1+lamdba_2+k_10)^2-4*((lamdba_1*q_2+lamdba_2*q_1)*k_10+lamdba_1*lamdba_2);
alpha=(1/2)*(lamdba_1+lamdba_2+k_10+sqrt(Delta));
beta=(1/2)*(lamdba_1+lamdba_2+k_10-sqrt(Delta));
A1=alpha+beta;
A2=alpha*beta;
N=lamdba_1*lamdba_2/(A2);
P=(lamdba_1-alpha)*(lamdba_2-alpha)/(alpha^2*(beta-alpha));
Q=(lamdba_1-beta)*(lamdba_2-beta)/(beta^2*(alpha-beta));
m=-P-Q;
t=0:0.001:10;
%M=zeros(size(t));
%Var=zeros(size(t));
Fano=zeros(size(t));
for i=1:length(t)
    if t(i)<=T
         M(i)=r*(m+N*t(i)+P*exp(-alpha*t(i))+Q*exp(-beta*t(i)));%mean
         twomoment(i)=2*r^2*(m*N*t(i)+m*P*(exp(-alpha*t(i))-1)+m*Q*(exp(-beta*t(i))-1)+(1/2)*N^2*t(i)^2-N*P*(t(i)*exp(-alpha*t(i))+(1/alpha)*(exp(-alpha*t(i))-1))+N*P*t(i)*(exp(-alpha*t(i))-1)+N*Q*t(i)*(exp(-beta*t(i))-1)-N*Q*(t(i)*exp(-beta*t(i))+(1/beta)*(exp(-beta*t(i))-1))+(1/alpha)*P*N*(1-exp(-alpha*t(i)))-P^2*alpha*t(i)*exp(-alpha*t(i))-(1/(alpha-beta))*P*Q*beta*(exp(-beta*t(i))-exp(-alpha*t(i)))+(1/beta)*N*Q*(1-exp(-beta*t(i)))-Q^2*beta*t(i)*exp(-beta*t(i))-(1/(beta-alpha))*P*Q*alpha*(exp(-alpha*t(i))-exp(-beta*t(i))))+M(i);
        %second moment
         Var(i)=2*r^2*(m*N*t(i)+m*P*(exp(-alpha*t(i))-1)+m*Q*(exp(-beta*t(i))-1)+(1/2)*N^2*t(i)^2-N*P*(t(i)*exp(-alpha*t(i))+(1/alpha)*(exp(-alpha*t(i))-1))+N*P*t(i)*(exp(-alpha*t(i))-1)+N*Q*t(i)*(exp(-beta*t(i))-1)-N*Q*(t(i)*exp(-beta*t(i))+(1/beta)*(exp(-beta*t(i))-1))+(1/alpha)*P*N*(1-exp(-alpha*t(i)))-P^2*alpha*t(i)*exp(-alpha*t(i))-(1/(alpha-beta))*P*Q*beta*(exp(-beta*t(i))-exp(-alpha*t(i)))+(1/beta)*N*Q*(1-exp(-beta*t(i)))-Q^2*beta*t(i)*exp(-beta*t(i))-(1/(beta-alpha))*P*Q*alpha*(exp(-alpha*t(i))-exp(-beta*t(i))))+M(i)-M(i)^2;
        %varance
         Fano(i)=Var(i)/M(i);%Fano factor
    else
          M(i)=r*(N*T+P*exp(-alpha*t(i))*(1-exp(alpha*T))+Q*exp(-beta*t(i))*(1-exp(beta*T)));%mean
          twomoment(i)=2*r^2*(m*N*T+m*P*exp(-alpha*t(i))*(1-exp(alpha*T))+m*Q*exp(-beta*t(i))*(1-exp(beta*T))+(1/2)*N^2*T^2-N*P*exp(-alpha*t(i))*(T*exp(alpha*T)+(1/alpha)*(1-exp(alpha*T)))-N*Q*exp(-beta*t(i))*(T*exp(beta*T)+(1/beta)*(1-exp(beta*T)))+(1/alpha)*P*N*(1-exp(-alpha*T))-P^2*alpha*T*exp(-alpha*t(i))-(1/(alpha-beta))*P*Q*beta*(1-exp(-(alpha-beta)*T))*exp(-beta*t(i))+(1/beta)*N*Q*(1-exp(-beta*T))-Q^2*beta*T*exp(-beta*t(i))-(1/(beta-alpha))*P*Q*alpha*(1-exp(-(beta-alpha)*T))*exp(-alpha*t(i)))+M(i);
         %second moment
          Var(i)=2*r^2*(m*N*T+m*P*exp(-alpha*t(i))*(1-exp(alpha*T))+m*Q*exp(-beta*t(i))*(1-exp(beta*T))+(1/2)*N^2*T^2-N*P*exp(-alpha*t(i))*(T*exp(alpha*T)+(1/alpha)*(1-exp(alpha*T)))-N*Q*exp(-beta*t(i))*(T*exp(beta*T)+(1/beta)*(1-exp(beta*T)))+(1/alpha)*P*N*(1-exp(-alpha*T))-P^2*alpha*T*exp(-alpha*t(i))-(1/(alpha-beta))*P*Q*beta*(1-exp(-(alpha-beta)*T))*exp(-beta*t(i))+(1/beta)*N*Q*(1-exp(-beta*T))-Q^2*beta*T*exp(-beta*t(i))-(1/(beta-alpha))*P*Q*alpha*(1-exp(-(beta-alpha)*T))*exp(-alpha*t(i)))+M(i)-M(i)^2;
         %varance
          Fano(i)=Var(i)/M(i);%Fano factor
    end
end

plot(t,M,'g','linewidth',2)
% plot(t,Fano,'r','linewidth',2)
% plot(t,twomoment,'r','linewidth',2)
% plot(t,Var,'r',[0 1 0],'linewidth',2)
hold on;





k_10=5;%ON到OFF的切换率
k_01=10/11;%OFF到ON的切换率0.435
r=5;%5生成率
L=0.4436;%基因长度
v=0.08;%聚合酶移动速率
% T=L/v;
T=1;
A=r*k_10/(k_01+k_10)^2;
B=r*k_01/(k_01+k_10);
C=k_01+k_10;
t=0:0.001:10;
%M=zeros(size(t));
%Var=zeros(size(t));
Fano=zeros(size(t));
for i=1:length(t)
    if t(i)<=T
        M(i)=-A*exp(-C*t(i))+B*t(i)+A; %mean
        Var(i)=B^2*t(i)^2+2*(A^2)*exp(-C*t(i))*(-1+exp(C*t(i))-C*t(i))+(1/C)*4*A*B*(-1+exp(-C*t(i))+C*t(i))+M(i)-M(i)^2;
        %variance
        M2(i)=Var(i)+M(i)^2;%second moment
        Fano(i)=Var(i)/M(i);%fano
    else
        M(i)=A*exp(-C*t(i))*(-1+exp(T*C))+B*T;
        Var(i)=2*(-A^2*C*T*exp(-C*t(i))-A*B*t(i)*exp(-C*t(i))*(1-exp(C*T))+A*B*(t(i)*exp(-C*t(i))-(t(i)-T)*exp(-C*(t(i)-T))+(1/C)*exp(-C*t(i))*(1-exp(C*T)))-A^2*exp(-C*t(i))*(1-exp(C*T))-(1/C)*A*B*(1-exp(-C*T))+(1/2)*B^2*T^2+A*B*T)++M(i)-M(i)^2;
        M2(i)=Var(i)+M(i)^2;
        Fano(i)=Var(i)/M(i);
    end
end

plot(t,M,'y','linewidth',2)
% plot(t,Fano,'r','linewidth',2)
% plot(t,M2,'r','linewidth',2)
% plot(t,Var,'r',[0 1 0],'linewidth',2)
hold on;






k_20=5;%ON到OFF0的切换率
k_01=20/11;%OFF0到OFF1的切换率 20/11 0.8
k_12=20/11;%OFF1到ON的切换率 （OFF停留时间为2.5min） 20/11  0.8
r=5;
L=0.4436;
v=0.08;
% T=L/v;
T=1;
Delta=(k_01+k_12+k_20)^2-4*((k_01+k_12)*k_20+k_01*k_12);
alpha=(1/2)*(k_01+k_12+k_20+sqrt(Delta));
beta=(1/2)*(k_01+k_12+k_20-sqrt(Delta));
A1=alpha+beta;
A2=alpha*beta;
N=k_01*k_12/(A2);
P=(k_01-alpha)*(k_12-alpha)/(alpha^2*(beta-alpha));
Q=(k_01-beta)*(k_12-beta)/(beta^2*(alpha-beta));
m=-P-Q;
t=0:0.001:10;
%M=zeros(size(t));
%Var=zeros(size(t));
Fano=zeros(size(t));
for i=1:length(t)
    if t(i)<=T
        M(i)=r*(m+N*t(i)+P*exp(-alpha*t(i))+Q*exp(-beta*t(i)));%mean
        Var(i)=2*r^2*(m*N*t(i)+m*P*(exp(-alpha*t(i))-1)+m*Q*(exp(-beta*t(i))-1)+(1/2)*N^2*t(i)^2-N*P*(t(i)*exp(-alpha*t(i))+(1/alpha)*(exp(-alpha*t(i))-1))+N*P*t(i)*(exp(-alpha*t(i))-1)+N*Q*t(i)*(exp(-beta*t(i))-1)-N*Q*(t(i)*exp(-beta*t(i))+(1/beta)*(exp(-beta*t(i))-1))+(1/alpha)*P*N*(1-exp(-alpha*t(i)))-P^2*alpha*t(i)*exp(-alpha*t(i))-(1/(alpha-beta))*P*Q*beta*(exp(-beta*t(i))-exp(-alpha*t(i)))+(1/beta)*N*Q*(1-exp(-beta*t(i)))-Q^2*beta*t(i)*exp(-beta*t(i))-(1/(beta-alpha))*P*Q*alpha*(exp(-alpha*t(i))-exp(-beta*t(i))))+M(i)-M(i)^2;
        %variance
        M2(i)=Var(i)+M(i)^2;%second moment
        Fano(i)=Var(i)/M(i);
    else
        M(i)=r*(N*T+P*exp(-alpha*t(i))*(1-exp(alpha*T))+Q*exp(-beta*t(i))*(1-exp(beta*T)));
         Var(i)=2*r^2*(m*N*T+m*P*exp(-alpha*t(i))*(1-exp(alpha*T))+m*Q*exp(-beta*t(i))*(1-exp(beta*T))+(1/2)*N^2*T^2-N*P*exp(-alpha*t(i))*(T*exp(alpha*T)+(1/alpha)*(1-exp(alpha*T)))-N*Q*exp(-beta*t(i))*(T*exp(beta*T)+(1/beta)*(1-exp(beta*T)))+(1/alpha)*P*N*(1-exp(-alpha*T))-P^2*alpha*T*exp(-alpha*t(i))-(1/(alpha-beta))*P*Q*beta*(1-exp(-(alpha-beta)*T))*exp(-beta*t(i))+(1/beta)*N*Q*(1-exp(-beta*T))-Q^2*beta*T*exp(-beta*t(i))-(1/(beta-alpha))*P*Q*alpha*(1-exp(-(beta-alpha)*T))*exp(-alpha*t(i)))+M(i)-M(i)^2;
         M2(i)=Var(i)+M(i)^2;%second moment
         Fano(i)=Var(i)/M(i);
    end
end
plot(t,M,'b','linewidth',2)
% plot(t,Fano,'r','linewidth',2)
% plot(t,M2,'r','linewidth',2)
% plot(t,Var,'r',[0 1 0],'linewidth',2)
hold on;



