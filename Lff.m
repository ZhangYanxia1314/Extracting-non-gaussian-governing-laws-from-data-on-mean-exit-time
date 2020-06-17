% The drift by solving the inverse problem Au=-1
function f=Lff(x,sigma,alpha,TF)
if TF==0
    f=x-5*x.^3;
elseif TF==2
    f=1.0254*x-5.0400*x.^3;
else
    N=18;
    r=zeros(N,1);
    r0=0.9382;
    r=[0,-0.3199,0,0.2981,0,-2.071,0,10.62,0,-36.0,0,72.19,0,-85.53,0,55.07,0,-15.2]; % the coeffs of U
    M1=0; M2=0; dM1=0; dM2=0; U0=r0; dU=0; ddU=0; d3U=0;
    for j=1:N
        for k=1:j
            if j>=2 && mod(k,2)==0
                M1=M1+2*nchoosek(j,k)*r(j)*x^(j-k)*(x+1).^(k-alpha)/(k-alpha);
                if j-k>=1
                    dM1=dM1+2*nchoosek(j,k)*r(j)/(k-alpha)*...
                        ((j-k)*x^(j-k-1)*(x+1)^(k-alpha)+(k-alpha)*x^(j-k)*(x+1)^(k-alpha-1));
                end
            end
            if k-alpha==0
                M2=M2+nchoosek(j,k)*r(j)*x^(j-k)*log((1-x)/(x+1));
            else
                M2=M2+nchoosek(j,k)*r(j)*x^(j-k)*((1-x).^(k-alpha)-(x+1).^(k-alpha))/(k-alpha);
            end
            if j-k>=1 && k-alpha==0
                dM2=dM2+nchoosek(j,k)*r(j)*((j-k)*x^(j-k-1)*log((1-x)/(x+1))+2*x^(j-k)/(x^2-1));
            elseif j-k>=1
                dM2=dM2+nchoosek(j,k)*r(j)/(k-alpha)*...
                    ((j-k)*x^(j-k-1)*((1-x).^(k-alpha)-(x+1).^(k-alpha))+...
                    x^(j-k)*(k-alpha)*(-(1-x).^(k-alpha-1)-(x+1).^(k-alpha-1)));
            end
        end
        dU=dU+j*r(j)*x^(j-1);
        if j>=2;
            ddU=ddU+j*(j-1)*r(j)*x^(j-2);
        end
        if j>=3
            d3U=d3U+j*(j-1)*(j-2)*r(j)*x^(j-3);
        end
        U0=U0+r(j)*x^j;
    end
    
    Ca=alpha*gamma(1/2+alpha/2)/2.^(1-alpha)/sqrt(pi)/gamma(1-alpha/2);
    W=-sigma/2*ddU+Ca/alpha*(1/(x+1)^alpha+1/(1-x)^alpha)*U0-Ca*(M1+M2)-1;
    dW=-sigma/2*d3U+Ca*(-(x+1)^(-1-alpha)+(1-x)^(-1-alpha))*U0+...
        Ca/alpha*(1/(x+1)^alpha+1/(1-x)^alpha)*dU-Ca*(dM1+dM2);
    if abs(dU)<0.05
        f=dW/ddU;
    else
        f=W/dU;
    end
end
end

