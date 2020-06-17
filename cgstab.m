function varargout=cgstab(varargin)
%CGSTAB  BiConjugate Gradients Stabilized(ell) for multiple right hand-sides.
%  X = CGSTAB(A,B) attempts to solve the system of linear equations A*X = B 
%  for X. The coefficient matrix A must be square and the right hand side 
%  B must by N-by-P, where A is N-by-N.
% 
%  CGSTAB will start iterating from an initial guess which by default is 
%  an array of size N-by-P of all zeros. Iterates are produced until the 
%  method either converges, fails, or has computed the maximum number of 
%  iterations. Convergence is achieved when an iterate X has relative 
%  residual NORM(B(:,I)-A*X(:,I))/NORM(B(:,I)) less than or equal to the 
%  tolerance of the method for all I=1:P. The default tolerance is 1e-8. 
%  The default maximum number of iterations is 300. No preconditioning is used.
%  CGSTAB uses BiCGstab(ell) with ell = 4 as default.
%
%  [X,HIST] = CGSTAB(A,B) returns also the convergence history
%  (the norms of the subsequential residuals).
%
%  ... = CGSTAB('Afun',B)
%  The first input argument is either a square matrix (which can be full or 
%  sparse, symmetric or nonsymmetric, real or complex), or a string 
%  containing the name of an M-file which applies a linear operator to a 
%  given array of size SIZE(B). In the latter case, N = SIZE(B,1). 
%
%  The remaining input arguments are optional and can be given in
%  practically any order:
%  ... = CGSTAB(...,X0,TR0,OPTIONS,M) 
%  where
%
%      X0        An N-by-P array of doubles: the initial guess.
%      TR0       An N-by-1 array of doubles: the shadow residual.
%      OPTIONS   A structure containing additional parameters.
%      M         String(s) or arrays(s): the preconditioner.
%
%  If X0 is not specified then X0 = ZEROS(N,P).
%  If TR0 is not specified then TR0 = (B-A*X0)*RAND(P,1).
%  The OPTIONS structure specifies certain parameters in the algorithm.
%
%   Field name          Parameter                              Default
% 
%   OPTIONS.Tol         Convergence tolerance:                 1e-8 
%   OPTIONS.MaxIt       Maximum number of iterations.          300
%   OPTIONS.Disp        Shows size of intermediate residuals.  1
%   OPTIONS.ell         ell for BiCGstab(ell).                 4 
%   OPTIONS.Scaled_acc                                         0
%            specifies in what sense the systems 
%            of equations are to be solved         
%            0: each equation with an accuracy of Tol:
%                NORM(A*X(:,I)-B(:,I))<Tol              all I
%            1: each equation with a relative accuracy Tol:
%                NORM(A*X(:,I)-B(:,I))<Tol*norm(B(:,I)) all I
%            2: each equation in the space with a 
%            relative accuracy Tol:
%                NORM(A*X*MU-B*MU)<Tol*NORM(B*MU) 
%                                   for all MU of size P-by-1
%   OPTIONS.Simultane                                          1
%   OPTIONS.TypePrecond                                        'left'
%            Specifies in what way the preconditioner is to
%            be applied. The precodnitioning will always  
%            be explicit, but it can be left, right, central 
%            (=two-sided if two matrices have been specified), 
%            or via the Eisenstat trick (if A is matrix and 
%            the preconditioner is a diagonal matrix).
%   OPTIONS.Omega                                              0.97
%            If A is a matrix (not a string), 
%            and no preconditioner has been specified, but
%            a type of preconditioning has been specied
%            then d_RILU(Omega) will be used, that is
%            a preconditoner M of the form 
%            M = (TRIL(A,-1)+D)*INV(D)*(TRIU(A,1)+D)
%            with D diagonal with a specific property:
%            if Omega==0, DIAG(A-M)=ZEROS(N,1)        (d_ILU)
%            if Omega==1, (A-M)*ONES(N,1)=ZEROS(N,1)  (d_MILU)
%            else mixure of these.
% 
%  If M is not specified then M = I (no preconditioning).
%  A preconditioner can be specified in the argument list:
%   ... = CGSTAB(...,M,...) 
%   ... = CGSTAB(...,L,U,...) 
%   ... = CGSTAB(...,L,U,P,...) 
%   ... = CGSTAB(...,'M',...) 
%   ... = CGSTAB(...,'L','U'...) 
%  as an N-by-N matrix M (then M is the preconditioner), 
%  or as N-by-N matrices L and U (then  L*U is the preconditioner),
%  or as N-by-N matrices L, U, and P (then P\(L*U) is the preconditioner),
%  or as one or two strings containing the name of  M-files ('M', or 
%  'L' and 'U') which apply a linear operator to a given N-by-P array.
%
%  CGSTAB (without input arguments) lists the options and the defaults.


%   Gerard Sleijpen.
%   Copyright (c) april 99


  global  n_operator Operator_params

  if nargin == 0, ShowOptions, return, end

  %%% set default parameters for the iteration method
  defaults=struct('x0',zeros(2,1),'tr0',zeros(2,1),...
    'tol',1.e-8,'disp',0,'maxit',300,'ell',4,...
    'scaled_acc',1,'simultane',1,'angle',0.7);

  %%% set default parameters for the preconditioner
  defaultsPC=struct('typeprecond','ll','omega',-0.5,...
    'Precond',sparse(2,2),'L',sparse(2,2),'U',sparse(2,2),'P',sparse(2,2));

  %%% Find dimension problem, check definition A and b
  [n,p0,r]=FindDimension(varargin{1:2});

  %%% Find parameters for iterative method
  [x,tr,tol,show,kmax,L,scale,simul,kappa,Pre_params]=...
     setparameters(0,n,p0,defaults,varargin{3:nargin});
  show=Str2Num(show); scale=Str2Num(scale,1); simul=Str2Num(simul,1);

  %%% Find parameters for preconditioning
  [tp,omega,Mpc,Lpc,Upc,Ppc,Operator_params]=...
     setparameters(0,n,p0,defaultsPC,Pre_params{:});

  %%% Set preconditioning
  n=SetPrecond(n,Mpc,Lpc,Upc,Ppc);
  n=SetTypePrecond(n,tp,omega);

  if show>0 & n>0
    ShowChoices(x,r,tr,n,p0,kmax,tol,L,scale,kappa,show,simul)
  end

  %%% handle trivial cases
  if n<1, [varargout{1:2}]=output([],[]); return, end
  if n<L, A=MV(eye(n)); [varargout{1:2}]=output(A\r,[]); return, end
  if p0==0, [varargout{1:2}]=output(zeros(n,0),[]); return, end

  tol=tol/sqrt(p0);

  str='';%'%2d accurate solution%s at step %3d\n';
  str1='\rmax(|r_%i|/|r_0|) = %0.1e;  ';
  str2='\r    |r_%i|/|r_0|  = %0.1e;  ';

  t0=clock; if show, b=r;  end

  %%% non trivial initial guess?
  ntx = (norm(x,1) ~= 0); 
  
  %%% check dimensions initial guess 
  x = CheckDimV(x,n,p0);  
  
  %%% Scale residuals: make the residuals as independent as possible
  [r,x,d0,p,nrm0] = ScaleResid(r,x,p0,tol,scale,ntx);

  %%% Shift the problem and precondition residual
  [r,x,x0,ntx] = PreProcess(r,x,n,p0,ntx);

  J=1:p; J0=1:p0;
  if isempty(tr)
    rt=rand(p,1); rt=rt/norm(rt); tr=r(:,J)*rt;
  else
    tr=tr/norm(tr);
  end

  warning off

%%% Iteration loops
if simul %%%%  banded BiCGstab(ell)

  hist=[mynorm(r),0,0];

  r=[r,zeros(n,p0*L)]; u=zeros(n,L+1); S=reshape(1:p0*(L+1),p0,L+1);
  sigma=1; omega=1; Q=eye(p0);

  j=1; k=0;  K=0; 

  while k<kmax 
    sigma=-omega*sigma;

    for l=1:L     
      [w,rho,tau]=HHolder(r(:,S(J,l))'*tr); v=w/tau; rho=rho';
      for i=1:l
        r(:,S(J,i))=r(:,S(J,i))-(r(:,S(J,i))*w)*v'; 
      end
      Q(:,J)=Q(:,J)-(Q(:,J)*w)*v'; 
      x(:,J)=x(:,J)-(x(:,J)*w)*v'; 

      beta=rho/sigma; 
      I=1:l; jI=S(j,I);
      u(:,I)=r(:,jI)-beta*u(:,I);   
      u(:,l+1)=PMV(u(:,l));
 
      sigma=tr'*u(:,l+1); alpha=rho/sigma; 
      x(:,j) = x(:,j) +alpha*u(:,1);
      r(:,jI)= r(:,jI)-alpha*u(:,I+1); 
            
      r(:,S(J,l+1))=PMV(r(:,S(J,l)));
    end

    [mx,q]=max(mynorm(r(:,J))); jI=S(J(q),:); 
    [z,omega]=Gamma(r(:,jI),kappa,tol);
    u(:,1)=u(:,1)-u(:,I+1)*z;
    for q=1:p, i=J(q);
      x(:,i)=x(:,i)+r(:,S(i,I))  *z;
      r(:,i)=r(:,i)-r(:,S(i,I+1))*z;
    end

    k=k+1;

    nrm=mynorm(r(:,J0));                  
                             hist=[hist;nrm,n_operator,0];
                             if show, fprintf(str1,k,max(nrm)), end
    if min(nrm(J))<tol
      J=J(find(nrm(J)>=tol));
      p1=p-length(J); K=K+p1*k; st=' '; if p1>1, st='s'; end
      fprintf(str,p1,st,k),                 hist(end,end)=p1;
      p=p-p1; if p==0, break, end
     
      if 0
        L1=2; [gamma,I,J1]=Res(r,S(J,1:L1+1),tol);
        if ~isempty(I), I=J(I); 
          fprintf(', plus %d extra',length(I))
          x(:,I)=x(:,I)-r(:,reshape(S(J,1:L1),1,L1*p))*gamma;
          r(:,I)=r(:,I)-r(:,reshape(S(J,2:L1+1),1,L1*p))*gamma;
        end, J=J(J1);     
        p1=length(J); K=K+(p-p1)*k; p=p1;
        if p==0, break, end
      end

     j=J(1);
     end

     if p>1 & mod(k,1)==0
       [rr,s,v]=svd(r(:,J),0); r(:,J)=rr*s;
       x(:,J)=x(:,J)*v; Q(:,J)=Q(:,J)*v;
     end

  end
  x=x*Q'; if show r=r(:,J0)*Q'; end

else %%% standard BiCGstab(ell)

  hist=[]; K=0;

  for j=1:p, 
    sigma=1; omega=1; k=0; nrm=1; hist0=[nrm,n_operator,0];
    rr=[r(:,j),zeros(n,L)]; u=zeros(n,L+1); 

    while k<kmax & nrm>=tol
      sigma=-omega*sigma;

      for l=1:L,  
        I=1:l;   
        rho=tr'*rr(:,l); beta=rho/sigma; 
        u(:,I)=rr(:,I)-beta*u(:,I);   
        u(:,l+1)=PMV(u(:,l));
 
        sigma=tr'*u(:,l+1); alpha=rho/sigma; 
        x(:,j) = x(:,j) +alpha*u(:,1);
        rr(:,I)= rr(:,I)-alpha*u(:,I+1); 
        rr(:,l+1)=PMV(rr(:,l));
      end

      [gamma,omega,nrm]=Gamma(rr,kappa,tol);
      x(:,j)=x(:,j)+rr(:,I)*gamma;
      k=k+1;
                           hist0=[hist0;nrm,n_operator,0];
                           if show, fprintf(str2,k,nrm), end
      if nrm<tol
        fprintf(str,1,'',k)
        if j<p, tol=tol*sqrt(p+1-j)/sqrt(p-j); end 
        if show, rr(:,1)=rr(:,1)-rr(:,I+1)*gamma; end
        %%% sufficient accuracy, no need to update u and r
        break
      end 

      u(:,1)=u(:,1)-u(:,I+1)*gamma;
      rr(:,1)=rr(:,1)-rr(:,I+1)*gamma;

    end
    hist=vec2mat(hist,hist0(:,1));
    r(:,j)=rr(:,1); K=K+k; 
  end
  if p<p0, hist=vec2mat(zeros(1,p0-p)+tol,hist); end
  k=K;
  k1=size(hist,1)-1; K1=(0:k1)'/k1; hist(:,p+1)=K1*n_operator;

end

warning on
%%% Undo preconditioning and shift
x = PostProcess(x,x0,ntx); 

%%% Undo scaling
x=x*d0;   

t0=etime(clock,t0);
hist(1,p0+2)=t0;

if show
  DisplayStatistics(x,r*d0,hist,b,nrm0,k,K,t0)
  PlotHistory(hist,simul,L)
end

[varargout{1:2}]=output(x,hist);

return
%====================================================================
% AUXILLARY SUBROUTINES
%====================================================================
function typ=plt(p,gamma)

  kleur=['b';'r';'g';'k';'m'];
  plotsign=['-o';'-+';'-*';'-s';'-p';'-h';'-+']; 
  plk=size(kleur,1); plk=mod(p-1,plk)+1;
  if nargin==1 | gamma==1
    pls=size(plotsign,1); pls=mod(p-1,pls)+1;
    typ=[kleur(plk,:),plotsign(pls,:)];
  end
  if gamma==2
    typ=[kleur(plk,:),'.:'];
  end

return
%====================================================================
function tol=tolerance(tol,nrm,J)

  p=length(nrm); I=setdiff(1:p,J); nrm=nrm(I);
  tol=sqrt((tol^2-sum(nrm.*nrm))/length(J));

return
%====================================================================
function n=mynorm(r)

  n=sqrt(sum(conj(r).*r));

return
%====================================================================
function x=vec2mat(x,y)

  [n,p]=size(x); [m,q]=size(y);
  if n==0, x=y; return, end
  if m==0, return, end
  if n<m
    x=[[x;ones(m-n,1)*x(n,:)],y]; return
  end
  if n>m
    x=[x,[y;ones(n-m,1)*y(m,:)]]; return
  end
  x=[x,y];

return
%====================================================================
function [gamma,I1,I2]=Res(r,S,tol)

  [p,m]=size(S); I=1:p; m=m*p; J=p+1:m; tol=tol*tol;
  S=reshape(S,1,m);
  G=r(:,S)'*r(:,S);
  gamma=G(J,J)\G(J,I);
  nrm=diag(G(I,I)-G(I,J)*gamma);
  I1=find(nrm<tol); I2=find(nrm>=tol);
  gamma=gamma(:,I1); 
  

return
%====================================================================
function [w,rho,tau]=HHolder(w)
% [v,rho,tau]=HHolder(w)
%   w is a k vector (k=size(w,1))
%   Computes Housholder transformation H, H=I-(1/tau)*v*v',
%   for which H*w=rho*I(:,1) with I=eye(k)

  rho=w'*w; tau=w(1,:);
  rho1=rho-tau'*tau; rho=-sqrt(rho)*sign(tau);
  tau=tau-rho; w(1,:)=tau;
  tau=(rho1+tau'*tau)/2; 

return
%====================================================================
function [z,omega,rnrm]=Gamma(R,kappa,tol)
% kappa=0 is equivalent to the original BiCGstab approach
% Default kappa=0.7

  l=size(R,2)-1; I=2:l; J=[1,l+1];

  for i=1:l+1 
    G(i,1:i)=R(:,i)'*R(:,1:i); G(1:i,i)=G(i,1:i)'; 
  end
  Gamma=G(I,I)\G(I,[1,l+1]);
  Ng=G(J,J)-G(I,J)'*Gamma; 
  %%% Ng should be pos.def. If Ng(1,1) is too small
  %%% follow a slightly more expensive, but more stable computation
  if Ng(1,1)<10*eps*G(1,1)
    [Q,U]=qr([R(:,I),R(:,J)],0); W=U([l,l+1],[l,l+1]);
    Ng=W'*W; 
    Gamma=U(I-1,I-1)\U(I-1,[l,l+1]);
  end
  %%% Ng(1,1)=norm(r-R(:,I)*Gamma(:,1))^2.
  %%% Stop if Ng(1,1) is small enough
  %if abs(Ng(1,1))<0.5*tol^2
  %  z=[Gamma(:,1);0]; omega=0; rnrm=sqrt(abs(Ng(1,1))); return
  %end

  cosine=abs(Ng(1,2))/sqrt(abs(Ng(1,1)*Ng(2,2))); omega=Ng(2,1)/abs(Ng(2,2));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Increase |omega| for enhancing the stability of the BiCG part     %%%
  %%%   Sleijpen & Van der Vorst, Numerical Algor., 10 (1995) 203--223  %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if cosine<kappa, omega=(kappa/cosine)*omega; end

  Gamma=[-1,0;Gamma;0,-1];
  z=Gamma(:,1)-omega*Gamma(:,2); 
  %%% compute norm residual at low dim costs
  rnrm=sqrt(abs(z'*G*z)); 
  z=z(2:l+1);

return
%====================================================================
%  OUTPUT
%====================================================================
function varargout=output(x,hist)

  varargout{1}=x;
  varargout{2}=hist;

return
%====================================================================
function  DisplayStatistics(x,r,hist,b,nrm0,k,K,t0)

  global n_operator type_preconditioning

  p=size(b,2);
  str=''; 
  if max(type_preconditioning==[1,3,4]); str='preconditioned '; end
  fprintf('\n\n r is the recursively updated %sresidual:',str)
  if isempty(find(nrm0~=1))
    accuracy=[mynorm(r)',mynorm(b-MV(x))'];
    fprintf('\n %18s %18s','norm(r)  ','norm(b-A*x)');
  else
    accuracy=[(mynorm(r)./nrm0)',(mynorm(b-MV(x))./nrm0)'];
    fprintf('\n %18s %18s','norm(r)./norm(b)  ','norm(b-A*x)./norm(b)')
  end
  fprintf('\n   %11.4e           %11.4e',accuracy'),fprintf('\n')

  if p == 1, 
    str0='\nSolved in %4i steps, %4i MVs, %0.3g sec.'; 
  else, 
    str0='\nAll solutions in %4i steps, %4i MVs, %0.3g sec.'; 
  end
  fprintf('\n'), 

  fprintf(str0,k,n_operator,t0);
  if p>1, 
    str0='\nAverage/sol      %4i      , %4i    , %0.3g     '; 
    fprintf(str0,ceil(K/p),ceil(n_operator/p),t0/p)
  end
  fprintf('\n\n')
return
%====================================================================
function PlotHistory(hist,simul,L)

  p=size(hist,2)-2;
  K=hist(:,p+1)/p; plot(K,log10(hist(:,1)),plt(p,1))
  if p>1, hold on, plot(K,log10(hist(:,2:p)),plt(p,2)), end
  if simul, bb='banded'; else, bb=''; end
  title(sprintf('%s BiCGstab(%i)',bb,L))
  xlabel(sprintf('# MVs/%i',p))
  ylabel('log_{10} ||r||_2')
  drawnow, zoom on, hold off

return
%====================================================================
function ShowOptions
 fprintf('\nOptions for cgstab. Defaults in {}.\n');

 fprintf('\nPROBLEM\n')
 fprintf('            A: [ square matrix | string ]\n');
 fprintf('            b: [ n=size(A,1) by p array of scalars ]\n');

 fprintf('\nOPTIONAL ARGUMENTS\n');
 fprintf('           x0: [ n by p array of scalars {zeros(n,p)} ]\n');
 fprintf('          tr0: [ n by 1 array of scalars {(b-A*x0)*rand(p,1)} ]\n');
 fprintf('Parameterlist: [ struct {see below}]\n');

 fprintf('\nPARAMETERLIST\n');
 fprintf('          Tol: [ positive scalar {1e-8} ]\n');
 fprintf('         Disp: [ ''yes'' | {''no''} ] or [ 1 | {0}]\n');
 fprintf('        MaxIt: [ positive integer {300} ]\n');
 fprintf('          ell: [ positive integer {4} ]\n');
 fprintf('   Scaled_acc: [ {1} | 0 ]\n');
 fprintf('    Simultane: [ {1} | 0 ]\n');
 fprintf('        Angle: [ double in [0,1] {0.7} ]\n');

 fprintf('\nPRECONDITIONING\n');
 fprintf('  TypePrecond: [ no | {left} | right | central | Eisenstat trick ]\n');
 fprintf('        Omega: [ double {[]} ]\n');
 fprintf('      Precond: [ n by n or n by 2n matrix | string {identity} ]\n');
 fprintf('            L: [ n by  n matrix | string {identity} ]\n');
 fprintf('            U: [ n by  n matrix | string {identity} ]\n');
 fprintf('            P: [ n by  n matrix | string {identity} ]\n');

 fprintf('\n')

return
%====================================================================
function ShowChoices(x,r,tr,n,p,kmax,tol,L,scale,kappa,show,simul)

  global A_operator ...
         type_preconditioning type_preconditioner ...
         L_precond U_precond P_precond RILU_omega

  fprintf('\n\nBiCGstab(%i)\n==================\n',L)

  fprintf('\nPROBLEM\n')
  fprintf('%14s: %s\n','A',TypeA(A_operator,n))
  fprintf('%14s: [%ix%i double]\n','b',n,p);

  fprintf('\nOPTIONAL ARGUMENTS\n');
  fprintf('%14s: ','x0')
  if norm(x)==0
    fprintf('zeros(%i,%i)\n',n,p);
  else
    fprintf('[%ix%i double]\n',n,p);
  end
  fprintf('%14s: ','tr0');
  if isempty(tr)
    fprintf('(b-A*x0)');
    if p>1, fprintf('*rand(%i,1)',p); end, fprintf('\n');
  else
    fprintf('[%ix1 double]\n',n);
  end

  fprintf('\nPARAMETERS\n'); str='%14s: %g\n';
  fprintf(str,'Tol',tol);
  fprintf(str,'Disp',show);
  fprintf(str,'MaxIt',kmax);
  fprintf(str,'Scaled_acc',scale);
  if p>1, fprintf(str,'Simultane',simul); end
  fprintf(str,'Angle',kappa);

  fprintf('\nPRECONDITIONING\n');
  fprintf('%14s:','TypePrecond');
  switch type_preconditioning
    case 0, fprintf(' no preconditioning\n');
    case 1, fprintf(' left preconditioning\n');
    case 2, fprintf(' right preconditioning\n');
    case 3, fprintf(' central preconditioning\n');
    case 4, fprintf(' Eisenstat trick\n');
  end
  if type_preconditioning
    if ~isempty(RILU_omega)
      if RILU_omega<inf
        fprintf('%14s: d_RILU(%g)\n','D',RILU_omega);
      else
        fprintf('%14s: Diag(A)\n','D');
      end
      if type_preconditioning < 4
        if ~isempty(U_precond) 
          fprintf('%14s: (tril(A,-1)+D)/D\n','L_precond');
          fprintf('%14s: triu(A,+1)+D\n','U_precond');
        else
          fprintf('%14s: D\n','L_precond');
        end
      end
    else
      switch type_preconditioner
        case {0,1,6},               
          fprintf('%14s: %s\n','preconditioner',TypeA(L_precond,n))
        case {2,3,4,5,8},           
          fprintf('%14s: %s\n','L_precond',TypeA(L_precond,n))
          fprintf('%14s: %s\n','U_precond',TypeA(U_precond,n))
        case 9
          fprintf('%14s: %s\n','L_precond',TypeA(L_precond,n))
          fprintf('%14s: %s\n','U_precond',TypeA(U_precond,n))
          fprintf('%14s: %s\n','P_precond',TypeA(P_precond,n))
        case 10
          fprintf('%14s: [%ix%i diagonal]\n','L_precond',n,n);
      end
    end
  end
 
  fprintf('\n')

return
%====================================================================
function str=TypeA(A,n)

  if ischar(A)
    str=sprintf('''%s.m''',A);
  elseif issparse(A)
    str=sprintf('[%ix%i sparse]',n,n);
  else
    str=sprintf('[%ix%i double]',n,n);
  end

return
%====================================================================
%  POSTPROCESS final approximation
%====================================================================
function x=PostProcess(x,x0,ntx)

  global type_preconditioning 

  switch type_preconditioning   %% unprecondition approximation in case of
    case 2,     x=Solve(x);     %% right
    case {3,4}, x=Solve(x,'U'); %% two-sided, Eisenstat
  end

  if ntx, x=x0+x; end  %% undo shift if problem was shifted

return
%====================================================================
%  PREPROCESS initial guess, residual
%====================================================================
function [r,x,x0,ntx0]=PreProcess(r,x,n,p,ntx)

  global type_preconditioning 


  x0=zeros(n,p); ntx0=0;
  if ntx
    switch type_preconditioning   %% shift problem in case of
    case {2,3,4}                  %% right, two-sided, Eisenstat 
      r=r-MV(x); x0=x; ntx0=1; x=zeros(n,p);
    end
  end

  switch type_preconditioning     %% precondition residual in case of
    case 1,     r=Solve(r);       %% left 
    case {3,4}, r=Solve(r,'L');   %% two-sided, Eisenstat
  end

return
%====================================================================
function [r,x,d,p,nrm0]=ScaleResid(r,x,p,tol,scale,ntx)
% A*(x/d)=b/d; (b/d)*v=r*s; A*((x/d)*v)=r*s; di=v'*d:  
% If A*y=r*s then x=y*di.
% Suppose  norm(A*y(:,j)-r(:,j)*s(j,j),2) < eps(j). 
% Then     norm(A*x(:,j)-b(:,j),2) < tol*d(j,j)   with   tol^2=sum(eps.^2)
% If, in addition, d=eye(p), then norm(A*x-b,F)<tol.
%
% A*(x/d)=b/d; (b/d)*v=r*s; A*((x/d)*v)=r*s; di=v'*d:  
% A*((x/d)*v)/s=r; di=s*v'*d:   A*y=r then x=y*di


  str='';%'\n%2d accurate solution%s at step %3d';

  nrm=mynorm(r); nrm0=nrm;  d=eye(p); 
  if scale==1, 
    d=diag(nrm); r=r/d; if ntx, x=x/d; end, nrm=ones(1,p);
  else
    nrm0=nrm*0+1;
  end
  if p>1
    [r,s,v]=svd(r,0); d=v'*d; if ntx, x=x*v; end
    switch scale
    case {0,1}, r=r*s; nrm=diag(s)'; 
    otherwise, d=s*d; x=x/s; nrm=ones(1,p);
    end
  end
          
  if min(nrm)<tol
    %%% The scalars in nrm decrease. Therefore to be solved 1:p
    p0=p; p=sum(nrm>=tol); p0=p0-p; st=' '; if p0>1, st='s'; end
    fprintf(str,p0,st,0)
  end 

return
%====================================================================
function x=CheckDimV(x,n,p)

    [nx,px]=size(x);  
    if (nx==n & px==p), return, end 
    if nx~=n
      if nx<n, x=[x;zeros(n-nx,px)]; else, x=x(1:n,:); end
    end
    if px ~= p
      if px<p, x=[x,zeros(n,p-px)]; else, x=x(:,1:p); end
    end

return
%====================================================================
%  ACTIONS OPERATORS
%====================================================================
function y=PMV(y)
% preconditioned matrix-vector multiplication

  global n_operator type_preconditioning

  n_operator=n_operator+size(y,2);

  switch type_preconditioning
    case 0, y=MV(y);                 % no precondioner
    case 1, y=MV(y);                 % left
            y=Solve(y);   
    case 2, y=Solve(y);              % right
            y=MV(y);   
    case 3, y=Solve(y,'U');          % two-sided 
            y=MV(y); 
            y=Solve(y,'L');
    case 4, y=EisenstatOperation(y); % Eisenstat trick
  end

return
%====================================================================
function y=Solve(y,flag)
% Action preconditioner

  global type_preconditioner L_precond U_precond P_precond Operator_params

  if nargin==1
    switch type_preconditioner
      case 0,     y=feval(L_precond,y,Operator_params{:});
      case 1,     y=feval(L_precond,y,'preconditioner',Operator_params{:});
      case 2,     y=feval(L_precond,y,Operator_params{:}); 
                  y=feval(U_precond,y,Operator_params{:});
      case 3,     y=feval(L_precond,y,'preconditioner',Operator_params{:}); 
                  y=feval(U_precond,y,Operator_params{:});
      case {4,5}, y=feval(L_precond,y,'L',Operator_params{:}); 
                  y=feval(L_precond,y,'U',Operator_params{:});
      case 6,     y=L_precond\y;
      case 8,     y=U_precond\(L_precond\y);
      case 9,     y=U_precond\(L_precond\(P_precond*y));
    end
  elseif strmatch(flag,'U')
    switch type_preconditioner
      case {2,3},  y=feval(U_precond,y,Operator_params{:});
      case {4,5},  y=feval(L_precond,y,'U',Operator_params{:});
      case {8,10}, y=U_precond\y;
    end
  elseif strmatch(flag,'L')
    switch type_preconditioner
      case 2,     y=feval(L_precond,y,Operator_params{:});
      case 3,     y=feval(L_precond,y,'preconditioner',Operator_params{:});
      case {4,5}, y=feval(L_precond,y,'L',Operator_params{:});
      case 8,     y=L_precond\y;
      case 10,    y=L_precond\y;
                  SetMatricesEisenstat
                  y=L_precond\y; 
    end
  end
    
return
%====================================================================
function y=MV(y)
% Action matrix/operator A

  global A_operator Operator_params

  if ischar(A_operator)
    y=feval(A_operator,y,Operator_params{:});
  else 
    y=A_operator*y;
  end

return
%====================================================================
function y=EisenstatOperation(y)

  global D_precond L_precond U_precond

  y1=U_precond\y;
  y=y+D_precond*y1; 
  y=y1+L_precond\y;

return
%====================================================================
%   CONSTRUCTION PRECONDITIONER
%====================================================================
function SetMatricesEisenstat
% set matrices for Eisenstat

  global A_operator D_precond L_precond U_precond

  n=size(A_operator,1);

  D_precond=L_precond; 
  L_precond=D_precond\tril(A_operator,-1) + speye(n); 
  U_precond=D_precond\triu(A_operator, 1) + speye(n);
  D_precond=D_precond\spdiags(spdiags(A_operator,0),0,n,n);
  D_precond=D_precond-2*speye(n); 

return
%====================================================================
function D=d_rilu(omega,gamma)
% D for the preconditioner M=(L_A+D)*inv(D)*(D+U_A), where A=L_A+D_A+U_A
% with D such that 
%      diag(D-A)=0      if gamma==1
%      diag(M-A)=0      if omega==0:  ILU
%     (M-A)*ones(n,1)=0 if omega==1: MILU

  global A_operator L_precond U_precond

  t0=clock;

  n=size(A_operator,1); D=diag(diag(A_operator)); 
  if gamma, return, end

  L_precond=tril(A_operator,-1);  U_precond=triu(A_operator,1);

  if omega<inf
    fprintf('Constructing RILU(%g)',omega)
    x=U_precond*(omega*ones(n,1)); omega=1-omega;
    %[I0,J]=find(L_precond); [I,I1]=sort(I0); J=J(I1); nn=length(I);
    [I,J]=find(L_precond); nn=length(I);
    for k=1:nn
      i=I(k); j=J(k);
      D(i,i)=D(i,i)-L_precond(i,j)*(D(j,j)\(omega*U_precond(j,i)+x(j,1)));
    end
    fprintf('\r                       \r')
  end

  L_precond=L_precond/D+speye(n);
  U_precond=U_precond+D;

  % fprintf('\rConstruction RILU(%g) took %g sec.',1-omega,etime(clock,t0))

return
%====================================================================
%        INPUT
%====================================================================
function [n,p,b]=FindDimension(A,b);

  global A_operator n_operator ...

  n_operator = 0;
  A_operator = A;

  if ischar(A_operator)
    n=-1;
    if exist(A_operator) ~=2
      msg=sprintf('  Can not find the M-file ''%s.m''  ',A_operator);
      errordlg(msg,'MATRIX'),n=-2; return
    end
  else
    [m,n] = size(A_operator);
  end

  %%% find the dimension n of the problem
  [nb,p]=size(b);
  if n > 0 
    if nb ~= n | p == 0
      msg=sprintf(' The right hand side vector b');
      msg=sprintf('%s\n is %i by %i, b must be %i by p, p>0.',msg,nb,p,n);
      errordlg(msg,'INPUT VECTOR'),n=-3;  return
    end
  else
    n=nb; 
  end

return
%======================================================================
%================== FIND PRECONDITIONER ===============================
%======================================================================
function n=SetPrecond(n,M,L,U,P)
% finds out how the preconditioners are defined (type_preconditioner)
% and checks consistency of the definition
%
% If M is the preconditioner then P*M=L*U. Defaults: L=U=P=I.
%
% type_preconditioner
%       0:   L M-file, no U,     L ~= A
%       1:   L M-file, no U,     L == A
%       2:   L M-file, U M-file, L ~= A, U ~= A, L ~=U
%       3:   L M-file, U M-file, L == A, U ~= A, 
%       4:   L M-file, U M-file, L ~= A,         L ==U
%       5:   L M-file, U M-file, L == A,         L ==U
%       6:   L matrix, no U
%       8:   L matrix, U matrix  no P
%       9:   L matrix, U matrix, P matrix
%      10:   if Eisenstat & L diag. matrix; set in SetTypePrecond.m

  global A_operator Operator_params ...
         type_preconditioning type_preconditioner ...
         L_precond U_precond P_precond

  if iscell(Operator_params) & ~isempty(Operator_params) & ...
     length(Operator_params) == 1
    Operator_params=Operator_params{1};
  end

  if ~isempty(M), 
     P_precond=U; U_precond=L; L_precond=M;
  else
     P_precond=P; U_precond=U; L_precond=L;
  end

  if isempty(L_precond), return, end

  if ~isempty(U_precond) & ischar(L_precond)~=ischar(U_precond)
     msg=sprintf('  L and U should both be strings or matrices');
     errordlg(msg,'PRECONDITIONER'), n=-1; return
  end
  if ~isempty(P_precond) & (ischar(P_precond) | ischar(L_precond))
      msg=sprintf('  P can be specified only if P, L and U are matrices'); 
      errordlg(msg,'PRECONDITIONER'), n=-1; return
  end  
  tp=6*~ischar(L_precond)+2*~isempty(U_precond)+~isempty(P_precond);
  if n>0 & tp<6
    tp=tp+strcmp(L_precond,A_operator);
    if tp>1, tp=tp+2*strcmp(L_precond,U_precond); end
  end
  type_preconditioner=tp;

  % Check consistency definitions
  if tp<6
    ok=1; 
    if exist(L_precond) ~=2
      msg=sprintf('  Can not find the M-file ''%s.m''  ',L_precond); 
      errordlg(msg,'PRECONDITIONER'), n=-1; return
    end
    if mod(tp,2)==1
      eval('v=feval(A_operator,zeros(n,1),''preconditioner'',Operator_params{:});','ok=0;')
      if ~ok
         msg='Preconditioner and matrix use the same M-file';
         msg1=sprintf(' %s.   \n',L_precond);
         msg2='Therefore the preconditioner is called';
         msg3=sprintf(' as\n\n\tw=%s(v,''preconditioner'')\n\n',L_precond); 
         msg4='Put this "switch" in the M-file.';
         msg=[msg,msg1,msg2,msg3,msg4];
      end
    end
    if tp>3 | ~ok
      ok1=1;
      eval('v=feval(L_precond,zeros(n,1),''L'',Operator_params{:});','ok1=0;')
      eval('v=feval(L_precond,zeros(n,1),''U'',Operator_params{:});','ok1=0;')
      if ~ok1 
        if ok
          msg='L and U use the same M-file';
          msg1=sprintf(' %s.m   \n',L_precond);
          msg2='Therefore L and U are called';
          msg3=sprintf(' as\n\n\tw=%s(v,''L'')',L_precond); 
          msg4=sprintf(' \n\tw=%s(v,''U'')\n\n',L_precond); 
          msg5=sprintf('Check the dimensions and/or\n');
          msg6=sprintf('put this "switch" in %s.m.',L_precond);
          msg=[msg,msg1,msg2,msg3,msg4,msg5,msg6]; 
        end
        errordlg(msg,'PRECONDITIONER'), n=-1; return
      elseif ~ok
        tp=5; type_preconditioner=5; U_precond=L_precond; ok=1; msg=[];
      end
    end
    if tp==0 | tp==2
      eval('v=feval(L_precond,zeros(n,1),Operator_params{:});','ok=0')
      if ~ok
         msg=sprintf('''%s'' should produce %i-vectors',L_precond,n); 
         errordlg(msg,'PRECONDITIONER'), n=-1; return 
      end
    end
    if tp == 2 | tp==3
      if exist(U_precond) ~=2
        msg=sprintf('  Can not find the M-file ''%s.m''  ',U_precond);
        errordlg(msg,'PRECONDITIONER'), n=-1; return
      else
        eval('v=feval(U_precond,zeros(n,1),Operator_params{:});','ok=0')
        if ~ok
          msg=sprintf('''%s'' should produce %i-vectors',U_precond,n);
          errordlg(msg,'PRECONDITIONER'), n=-1; return 
        end
      end
    end
  end

  msg=sprintf('L, U, and P should all be %iX%i.',n,n);
  if tp>7 & ~min([n,n]==size(L_precond) &  [n,n]==size(U_precond))
    if tp==8
      msg=sprintf('Both L and U should be %iX%i.',n,n);
    end
    errordlg(msg,'PRECONDITIONER'), n=-1; return 
  end

  if tp==9 & ~min([n,n]==size(P_precond))
    errordlg(msg,'PRECONDITIONER'), n=-1; return 
  end

  if tp==6
    if min([n,2*n]==size(L_precond)) 
      U_precond=L_precond(:,n+1:2*n); L_precond=L_precond(:,1:n); 
      type_preconditioner=8;
    elseif min([n,3*n]==size(L_precond)) 
      U_precond=L_precond(:,n+1:2*n); P_precond=L_precond(:,2*n+1:3*n);
      L_precond=L_precond(:,1:n); type_preconditioner=9;
    elseif ~min([n,n]==size(L_precond)) 
      msg=sprintf('The preconditioning matrix\n');
      msg2=sprintf('should be %iX%i or %ix%i ([L,U]).\n',n,n,n,2*n); 
      errordlg([msg,msg2],'PRECONDITIONER'), n=-1; return
    end
  end

return
%======================================================================
function n=SetTypePrecond(n,tp,omega)
% finds out how the preconditioning is to be implemented (type_preconditioning)
% and checks consistency of the definition
%
% type_preconditioning
%             0:   no preconditioning
%             1:   explicit left 
%             2:   explicit right 
%             3:   explicit central (= twosided)
%             4:   explicit twosided with Eisenstat trick


  global A_operator Operator_params ...
         type_preconditioning type_preconditioner ...
         L_precond U_precond RILU_omega

  type_preconditioning=...
     strmatch(lower(tp(1:2)),['ll';'no';'le';'ri';'ce';'ei'])-2;

  RILU_omega=[]; if omega==-0.5, omega=[]; end
  % Check consistency definitions
  
  if type_preconditioning==-1
    type_preconditioning=~isempty(L_precond);
  end 
  if type_preconditioning==0, return, end

  if isempty(L_precond) & type_preconditioning>0; 
    msg=sprintf('You request preconditioning, but\n');
    msg=sprintf('%sthere is no preconditioner.',msg);
    msg=sprintf('%s\n\nDo you want to continue with',msg);
    msg=sprintf('%s\n   1) no preconditioner,\n   2) ',msg);
    msg2=sprintf('stop?'); button='ri';
    if ischar(A_operator)
      button=questdlg([msg,msg2],'PRECONDITIONER',...
            'no precond.','stop','stop');
    elseif isempty(omega) % construct preconditioner from A
      msg1=sprintf('a preconditioner constructed from A, or  \n   3) ');
      button=questdlg([msg,msg1,msg2],'PRECONDITIONER',...
             'no precond.','precond.','stop','precond.');
    end
    button=lower(button(1:2));
    if strcmp(button,'st'),  n=-1; return, end
    if strcmp(button,'no'), type_preconditioning=0; return, end

    if isempty(omega)
      msg=sprintf('\nPrecondition with');
      button=questdlg(msg,'PRECONDITIONER',...
          'Diag(A)','SOR','RILU(omega)','RILU(omega)');
      button=lower(button(1:2));
      omega=inf;
      if strcmp(button,'ri')
        helpdlg('Give omega (in Matlab window)','RILU(omega)')
        while omega==inf
          o=input('RILU(0)=ILU; RILU(1)=MILU.\nGive omega (default 0.97):');
          if isempty(o), o=0.97; end
          if isnumeric(o) & all(size(o)==1), omega=o; end 
        end
      end
    end
    type_preconditioner=8; RILU_omega=omega;
    D=d_rilu(omega,strcmp(button,'di')); 
    if type_preconditioning==4, L_precond=D; end
    if strcmp(button,'di'), L_precond=D; type_preconditioner=6; end  
  end

  if max(type_preconditioner==[0;1;6]) & type_preconditioning==3;
    if max(type_preconditioner==[0;1]) 
    %% Check whether L_precond accepts 'L' and 'U'
       ok=1;
       eval('v=feval(L_precond,zeros(n,1),''L'',Operator_params{:});','ok=0;')
       eval('v=feval(L_precond,zeros(n,1),''U'',Operator_params{:});','ok=0;')
       if ok
         type_preconditioner=type_preconditioner+4;
       else
         msg=sprintf('There is no L_preconditioner and U_preconditioner');
         msg=sprintf('%s\nDo you want to continue with left preconditioning?',msg);
         button=questdlg(msg,'PRECONDITIONER','Yes','No','Yes');
         if strcmp(button,'No'),  n=-1; return, end
         type_preconditioning=1;
       end
    end
  end

  %% Check input for Eisenstat trick
  if type_preconditioning==4 
    if ischar(A_operator)        % is A a matrix?
        msg=sprintf('Eisenstat trick requires %s to be a matrix.',A_operator);
        n=-1; errordlg(msg,'PRECONDITIONER'), return
    end 
    if type_preconditioner < 6   % is L a matrix?
        msg=sprintf('Eisenstat trick requires a diagonal matrix %s',L_precond);
        n=-1; errordlg(msg,'PRECONDITIONER'), return
    end 
    if  norm(triu(L_precond,1),1)+norm(tril(L_precond,-1),1)~=0  
                                 % is L diagonal?
        msg=sprintf('Eisenstat trick requires a diagonal matrix');
        n=-1; errordlg(msg,'PRECONDITIONER'), return
    end
    type_preconditioner = 10;
  end

return
%======================================================================
%============== FIND PARAMETERS =======================================
%======================================================================
function varargout=setparameters(N,n,p,dopt,varargin);
%[PAR1,PAR2,...,OTHERPARS]=...
%      setparameters(J,N,P,DEFAULT_OPTIONS,INPAR1,INPAR2,INPAR3,...);
%
% First DEFAULT_OPTIONS values are replaced by new values as extracted
% from the input parameter list INPAR1,INPAR2,INPAR3,.... Then the values 
% DEFAULT_OPTIONS are listed in PAR1,PAR2, ... in order of occurence in 
% DEFAULT_OPTIONS.
%
% If INPARi is a structure then:
% INPARi is checked for exceptable OPTIONS and the appropriate 
% DEFAULT_OPTIONS values are replaced by the new values from OPTIONS.
% Field names are case insensative and only the first 2 characters 
% are checked.
%
% If INPARi is not a structure then:
% Input parameters are classified as operators of dimension N, vectors of
% size N by P, strings, booleans, pos. integers, doubles (in this order), 
% matrices, and within a class, in order of occurence, they replace
% the appropriate values in DEFAULT_OPTIONS if they have not been replaced
% before.
%
% The parameters that are excessive with respect to DEFAULT_OPTIONS
% are collected in OTHERPARS.
%

% Matlab does not distinguish between 1.0 and 1. As a consequence
% the double 1.0 will be interpreted as boolean. 
% The same remark aplies to integers. 


  fdopt = fieldnames(dopt);
  ndo=length(fdopt); J=1:ndo;
  na=length(varargin); JJ=1:na;

  % replace defaults by parameters that are defined by a struct
  for jj=1:na
    ar = varargin{jj};
    if isstruct(ar)
      [dopt,ar,J]=mergeopt(dopt,fdopt,ar,J);
      if isempty(ar), JJ=JJ(find(JJ~=jj)); else, varargin{jj}=ar; end
    end
  end 

  % replace defaults by the first N parameters from list
  for jj=1:min(N,na)
    ar = varargin{jj};
    if (~isempty(find(J==jj)) & (isnumeric(ar) | ischar(ar)))
      dopt=setfield(dopt,char(fdopt(jj,:)),ar);
      J=J(find(J~=jj)); JJ=JJ(find(JJ~=jj));
    end
  end

  if ~isempty(JJ)
  jd = {[],[],[],[],[],[],[],[],[],[]};
  for jj=J
    ar=getfield(dopt,char(fdopt(jj,:)));
    if isempty(ar),           jd{1}=[jd{1},jj]; % empty
    elseif isoperator(ar,n),  jd{2}=[jd{2},jj]; % operator
    elseif ischar(ar),        jd{3}=[jd{3},jj]; % strings
    elseif isboolean(ar),     jd{4}=[jd{4},jj]; % booleans
    elseif isinteger(ar),     jd{5}=[jd{5},jj]; % positive integers
    elseif iscnumber(ar),     jd{6}=[jd{6},jj]; % complex number
    elseif isvector(ar,n,p),  jd{7}=[jd{7},jj]; % vectors
    else,                     jd{8}=[jd{8},jj]; % matrices
    end
  end

  if ~isempty(jd{1}), 
    str=char(fdopt(jd{1}(1),:));
    str=sprintf('The option ''%s'' is empty.',str);
    str=[str,...
  sprintf('\nEmpty strings and empty non-strings are indistinguishable!')];
    error(str); 
  end

  j = {[0],[0],[0],[0],[0],[0],[0],[0],[0],[0]};
  for jj=JJ
    ar=varargin{jj};
    if isempty(ar)
    elseif isoperator(ar,n)                          % operator
      [j{2},dopt,JJ]=CR(j{2},jd{2},jj,ar,dopt,fdopt,JJ);
    elseif ischar(ar)                                % strings
      [j{3},dopt,JJ]=CR(j{3},jd{3},jj,ar,dopt,fdopt,JJ);
    elseif ~isnumeric(ar)                 % rule out cells, structs
    elseif isboolean(ar) %& j{4}<length(jd{4})        % boolean
      [j{4},dopt,JJ]=CR(j{4},jd{4},jj,ar,dopt,fdopt,JJ);
    elseif isinteger(ar) %& j{5}<length(jd{5})        % positive integers
      [j{5},dopt,JJ]=CR(j{5},jd{5},jj,ar,dopt,fdopt,JJ);
    elseif iscnumber(ar)                             % complex number
      [j{6},dopt,JJ]=CR(j{6},jd{6},jj,ar,dopt,fdopt,JJ);
    elseif isvector(ar,n,p)                          % vectors
      [j{7},dopt,JJ]=CR(j{7},jd{7},jj,ar,dopt,fdopt,JJ);
    else                                             % matrices
      [j{8},dopt,JJ]=CR(j{8},jd{8},jj,ar,dopt,fdopt,JJ);
    end
  end
  end

  for jj=1:ndo
   [ar,ok]=CheckT(getfield(dopt,char(fdopt(jj,:))));
   if ok, dopt=setfield(dopt,char(fdopt(jj,:)),ar); end      
   varargout{jj}=ar;
  end

  %%% collect the additional parameters in one cell
  varargout{ndo+1}=varargin(JJ);

return
%======================================================================
function  [j,dopt,J]=CR(j,jd,jj,ar,dopt,fdopt,J)
  if j<length(jd)     
    j=j+1; jo=jd(j); J=J(find(J~=jj));
    dopt=setfield(dopt,char(fdopt(jo,:)),ar);
  end
return
%======================================================================
function [dopt,options,J]=mergeopt(dopt,fdopt,options,J)
% N is the number of first characters to be compared in field names; 

  N=2;

  foptions=fieldnames(options);
  for j=1:length(foptions)
    name=char(foptions(j,:)); l=min(N,length(name));
    jj=min(find(strncmpi(lower(name(1:l)),fdopt,l)));
    if ~isempty(jj) 
      dname=char(fdopt(jj,:));
      ar=getfield(options,name); 
      options=rmfield(options,name);
      dopt=setfield(dopt,dname,ar); J=J(find(J~=jj)); 
    end
  end

return
%======================================================================
function [M,ok]=CheckT(M)
  ok = (istrivialM(M) | istrivialV(M)); if ok, M=[]; end
return
%======================================================================
function ok = istrivialM(M)
  ok = (issparse(M) & min(size(M)==[2,2]) & isempty(find(M~=0)));
return
%======================================================================
function ok = istrivialV(x)
  ok = (isnumeric(x) & min(size(x)==[2,1]) & isempty(find(x~=0)));
return
%======================================================================
function ok=isboolean(ar)
  ok = (length(ar)==1 & (ar == 0 | ar == 1));
return
%======================================================================
function ok=isinteger(ar)
  ok = (length(ar)==1 & round(ar) == ar & ar > 0);
return
%======================================================================
function ok=iscnumber(ar)
  ok = (length(ar)==1);
return
%======================================================================
function ok=isvector(ar,n,p)
  ok = ((size(ar,1)==n & size(ar,2)<=p) | istrivialV(ar));
return
%======================================================================
function ok=isoperator(ar,n)
  if ischar(ar), ok = (exist(ar)==2);
  elseif isnumeric(ar)
    ok = (min(size(ar))==n | istrivialM(ar)) ; 
  else, ok = 0; end
return
%======================================================================
%======================================================================
function x = Str2Num(x,gamma,string)
%Y = Str2Num(X,GAMMA,STRING)
%  GAMMA(1) is the default. If GAMMA is not specified, GAMMA = 0.
%  STRING is a matrix of accepted strings. 
%  If STRING is not specified STRING = ['no ';'yes']
%  STRING(I,:) and GAMMA(I) are accepted expressions for X 
%  If X=GAMMA(I) then Y=X. If X=STRING(I,:), then Y=GAMMA(I+1).
%  For other values of X, Y=GAMMA(1);

  l=1;
  if nargin < 2, gamma=0; end
  if nargin < 3, 
    string=strvcat('no','yes'); 
    l=length(gamma);
    gamma=[gamma,0,1]; 
  end

  if ischar(x)
    i=strmatch(lower(x),string,'exact'); 
    if isempty(i),i=1; else, i=i+l; end, x=gamma(i);
  elseif max((gamma-x)==0)
  elseif gamma(end) == inf
  else, x=gamma(1);
  end
  
return
%======================================================================
%======================================================================
