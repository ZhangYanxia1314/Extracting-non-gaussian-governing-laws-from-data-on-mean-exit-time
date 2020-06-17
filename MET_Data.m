%% A data set of observations on mean exit time from the original systems
clear;
clc;
syms x1
global xa xb J h

kesi=1;
xa=-1;
xb=1;
sigma=0.5;
alpha=0.6;
f1=x1-5*x1.^3;
f=eval(['@(x1)',vectorize(f1)]);
J=120;
h=(xb-xa)/2/J;
x=xa:h:xb;
n=length(x);
Ca=alpha*gamma(1/2+alpha/2)/2.^(1-alpha)/sqrt(pi)/gamma(1-alpha/2);
Ch=sigma/2-kesi*Ca*zeta(alpha-1)*h^(2-alpha);
U=zeros(n,1);
A1=zeros(2*J-1,1);
A2=zeros(2*J-1,1);
A3=zeros(2*J-1,1);
A=zeros(2*J-1,2*J-1);
B=zeros(2*J-1,2*J-1);
for j=-J+1:J-1
    xj=j*h;
    A1(j+J)=Ch/h.^2-f(xj)/2/h;% Uj-1
    A2(j+J)=-2*Ch/h.^2-kesi*Ca/alpha*(1/(-xa+xj)^alpha+1/(xb-xj)^alpha);%Uj
    A3(j+J)=Ch/h^2+f(xj)/2/h;% Uj+1
    s=0;
    for k=-J-j:J-j
        xk=k*h;
        if k==0
            s=s+0;
        elseif k==-J-j || k==J-j
            s=s+1/2/abs(xk)^(1+alpha);
        else
            s=s+1/abs(xk)^(1+alpha);
        end
    end
    A2(j+J)=A2(j+J)-kesi*Ca*h*s;%Uj
    
    for k=-J-j+1:J-j-1
        xk=k*h;
        if k==0
            B(j+J,j+J)=0;
        else
            B(j+J,j+J+k)=1/abs(xk)^(1+alpha);
        end
    end
    
end
A=diag(A2)+diag(A3(1:end-1), 1)+diag(A1(2:end),-1);
B=kesi*Ca*h*B;
Y=-ones(n-2,1);
U(2:n-1)=(A+B)\Y;
%U(2:n-1)=cgstab(A+B,Y); % choose this way in the case of singularity
figure;
plot(x,U,'r')
xlabel('$x$','Interpreter','latex');
ylabel('$MET$','Interpreter','latex');
title({['alpha=',num2str(alpha),' xa=',num2str(xa),' xb=',num2str(xb),' kesi=',num2str(kesi),' d=',num2str(sigma)]})
