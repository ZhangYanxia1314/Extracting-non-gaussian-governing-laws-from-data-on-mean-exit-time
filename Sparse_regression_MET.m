% Sparse regression with observations on MET
% load('Data-MET.mat') % Uob
global xa xb h
x=xa:h:xb;
n=length(x);
Ncoef=16;% the number of basis function coefficient
A=zeros(n,Ncoef+1);
A(:,1)=1*(x'-xa).*(xb-x');
for i=1:Ncoef
    A(:,i+1)=(x'-xa).*(xb-x').*x'.^i;% construct the basis function
end

B=U;
pos=0:1:Ncoef;% 1 x x^2...
for k=1:Ncoef
    c=(A'*A)\(A'*B);
    I=abs(c)<0.01; % Sparse Learning: delete the coefficient less than 0.01
    A(:,I)=[];
    pos(I)=[];
    if isempty(c(I))
        break;
    end
end

n2=length(pos);
UL=zeros(n,1);
for i=1:n2
    s=c(i)*(x'-xa).*(xb-x').*x'.^pos(i); % composite the learnt function
    UL=UL+s;
    strcat('c',num2str(pos(i)+1),'=',num2str(c(i)))
end
c'
figure;
plot(x,U,'r',x,UL,'b')
xlabel('$x$','Interpreter','latex');
ylabel('$MET$','Interpreter','latex');

%%% UL as a function of x
syms x1
UL1=0;
for i=1:n2
    s1=c(i)*(x1-xa).*(xb-x1).*x1.^pos(i);
    UL1=UL1+s1;
end
digits(4)
UL2=vpa(collect(UL1,x1))
r=coeffs(UL2,x1) % the coeffs of U
