%% the error objective function G(sigma,alpha)
clear;
clc;
global xa xb J h kesi
kesi=1;
xa=-1;
xb=1;
J=120;
h=(xb-xa)/2/J;

Nd=0.1:0.1:2;
Nalpha=0.1:0.1:2;
n1=length(Nd);
n2=length(Nalpha);
G=zeros(n1,n2);
Uob=MET(0.5,0.6,0); % A data set of observations on mean exit time from the original systems
for i=1:n1
    tic
    sigma=Nd(i);
    for j=1:n2
        alpha=Nalpha(j);
        Lu=MET(sigma,alpha,1); % compute the MET for all combinations of (sigma,alpha,drift)
        G(i,j)=norm(Lu-Uob,2).^2/norm(Uob,2).^2;
    end
    toc
    [n1 i]
end
[posd,posalpha]=find(G==min(min(G)));
Lsigma=Nd(posd)
Lalpha=Nalpha(posalpha) % the learned sigma and alpha

figure;
mesh(Nd,Nalpha,G')
xlabel('$\sigma$','Interpreter','latex');
ylabel('$\alpha$','Interpreter','latex');
zlabel('$G(\sigma,\alpha)$','Interpreter','latex')
%% Comparation between the learned drift and the true drift
Lsigma=0.5;
Lalpha=0.6;
x=xa:h:xb;
n=length(x);
Lf1=zeros(1,n);
Lf0=zeros(1,n);
for i=1:n
   Lf0(i)=Lff(x(i),0.5,0.6,0); % the true drift
   Lf1(i)=Lff(x(i),Lsigma,Lalpha,2); % The learned drift
end

figure;
plot(x,Lf0,'r',x,Lf1,'b')
xlabel('$x$','Interpreter','latex');
ylabel('Drift $f(x)$','Interpreter','latex');
legend('True','Learn') 
%% Comparation between the mean exit time from the learned SDE and the observations
Uob=MET(0.5,0.6,0);
Lu=MET(Lsigma,Lalpha,2);
figure;
plot(x,Uob,'o',x,Lu,'b')
xlabel('$x$','Interpreter','latex');
ylabel('Mean Exit Time $U(x)$','Interpreter','latex');
legend('Observations','Learn') 
