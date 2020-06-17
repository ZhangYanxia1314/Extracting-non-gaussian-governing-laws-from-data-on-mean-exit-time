% calculating the MET with any sigma and alpha
function U=MET(sigma,alpha,TF)

global xa xb J h kesi
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
f=Lfn(sigma,alpha,TF); % Calculate the approximated drift 
for j=-J+1:J-1 
    xj=j*h;
    A1(j+J)=Ch/h^2-f(j+J+1)/2/h;% Uj-1
    A2(j+J)=-2*Ch/h^2-kesi*Ca/alpha*(1/(-xa+xj)^alpha+1/(xb-xj)^alpha);%Uj
    A3(j+J)=Ch/h^2+f(j+J+1)/2/h;%Uj+1
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
A=diag(A2)+diag(A3(1:end-1),1)+diag(A1(2:end),-1);
B=kesi*Ca*h*B;
Y=-ones(n-2,1);
%U(2:n-1)=(A+B)\Y;
U(2:n-1)=cgstab(A+B,Y);
end