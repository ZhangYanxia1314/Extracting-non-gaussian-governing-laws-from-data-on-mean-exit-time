%% the approximated drift based on polynomial basis function
function f=Lfn(sigma,alpha,TF)
global xa xb h
x=xa:h:xb;
n=length(x);
f=zeros(n,1);
for i=1:n
    f(i)=Lff(x(i),sigma,alpha,TF); % Calculate the drift 
end
if TF==1
    Ncoef=4;% the number of basis function coefficient
    A=zeros(n,Ncoef+1);
    A(:,1)=1;
    for i=1:Ncoef
        A(:,i+1)=x'.^i; % construct the basis function
    end
    A1=A(2:end-1,:);
    B=f(2:end-1);
    pos=0:1:Ncoef;% 1 x x^2...
    for k=1:Ncoef
        c=(A1'*A1)\(A1'*B);
        I=abs(c)<0.1; % Sparse Learning: delete the coefficient less than 0.1
        A1(:,I)=[];
        pos(I)=[];
        if isempty(c(I))
            break;
        end
    end
    n2=length(pos);
    Lf2=zeros(n,1);
    for i=1:n2
        s=c(i)*x'.^pos(i); % composite the learnt function
        Lf2=Lf2+s;
    end
    f=Lf2;
    % where c with pos is the sparse coefficient of the learned drift with leared sigma and learned alpha
end
end
