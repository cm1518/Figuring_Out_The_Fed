function [probspce, lcpce, ldpce, p2pce, x2pce, x1pce, E_dpce, E_cpce, M,C,Md,Cd,Fd]=filter_final(param,pipce,y,i_vec,T, z_initial, g_initial,p2_initial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function that solves the model for given parameter values and solves for unobservables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SETTING UP PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa=param(1)/1.99;
sigma=param(2);
rho_z=param(3);
rho_g=param(4);
pi_bar=param(5);
r_bar=param(6);
lambda=param(7);
lambda_i=param(8);
sigma_measurement=param(9);
lambda_d=param(10);
lambda_id=param(11);
pi_bar_d=param(12);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calibrated parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_c=.5;
p_d=.5;
y_bar=0;
beta=.99;
gamma_f=.99/1.99;
gamma_b=1/1.99;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SETTING MODEL MATRICES (SEE PAPER BY SOEDERLIND THAT I CITE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pi_target=pi_bar;
pi_target_d=pi_bar_d;

i_bar=r_bar-1+pi_bar; %net interest rate, but r_bar is gross real rate

A_11=zeros(6,6);
A_11(1,1)=rho_z;
A_11(2,2)=rho_g;
A_11(3,3)=1;


A_12=zeros(6,2);
A_12(4,1)=1;
A_12(5,2)=1;
B_1=zeros(6,1);
B_1(6,1)=1;

H=[gamma_f 0; 1/sigma 1];

A_21=[1 0 (gamma_f-1+gamma_b)*pi_bar -gamma_b 0 0; 0 -1 1/sigma*(pi_bar-i_bar) 0 0 0];

A_22=[1 -kappa; 0 1];

B_2=[0;(1/sigma)];

A=[A_11 A_12;H\A_21 H\A_22];

B=[B_1;H\B_2];
Q=zeros(8,8);

Q(3,3)=pi_target^2+lambda*y_bar^2+lambda_i*i_bar^2;
Q(3,7)=-pi_target;
Q(3,8)=-y_bar*lambda;
Q(7,3)=-pi_target;
Q(7,7)=1;
Q(8,3)=-y_bar*lambda;
Q(8,8)=lambda;

R=lambda_i;

U=zeros(8,1);
U(3,1)=-lambda_i*i_bar;
i_bar_d=r_bar-1+pi_bar_d;
A_21_d=[1 0 (gamma_f-1+gamma_b)*pi_bar_d -gamma_b 0 0;0 -1 1/sigma*(pi_bar_d-i_bar_d) 0 0 0];
Q_d=zeros(8,8);


Q_d(3,3)=pi_target_d^2+lambda_d*y_bar^2+lambda_id*i_bar_d^2;
Q_d(3,7)=-pi_target_d;
Q_d(3,8)=-y_bar*lambda_d;
Q_d(7,3)=-pi_target_d;
Q_d(7,7)=1;
Q_d(8,3)=-y_bar*lambda_d;
Q_d(8,8)=lambda_d;
R_d=lambda_id;
U_d=zeros(8,1);
U_d(3,1)=-lambda_id*i_bar_d;

A_d=[A_11 A_12;H\A_21_d H\A_22];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SETTING UP STORAGE MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x1pce=zeros(6,T+1);
x1pce(1,1)=z_initial;
x1pce(2,1)=g_initial;
x1pce(3,:)=1;
x1pce(4,1)=pipce(1);
x1pce(4,2)=pipce(1);
x1pce(4,3:end)=pipce(1:end-1);
x1pce(5,1)=y(1);
x1pce(5,2)=y(1);
x1pce(5,3:end)=y(1:end-1)';
x1pce(6,1)=i_vec(1);
x1pce(6,2:end)=i_vec(1:end-1);

x2pce=zeros(2,T);
x2pce(1,:)=pipce';
x2pce(2,:)=y';
p2pce=zeros(2,T+1);
p2pce(:,1)=p2_initial;
probspce=zeros(2,T+1);
lcpce=zeros(1,T);
ldpce=zeros(1,T);


E_cpce=zeros(2,T);
E_dpce=zeros(2,T);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SETTING PRIOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

probspce(1,1)=p_c;
probspce(2,1)=p_d;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOLVING FOR OPTIMAL POLICIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [M,C] = ComItAlg( A,B,Q,R,U,beta,6,2,1.01 );

[Md,Cd,Vd,Fd] = DiscAlg( A_d,B,Q_d,R_d,U_d,beta,6,2,ones(6,6),ones(2,6),...
    (1e-6)*ones(2,1),1,1,2,0.5,0,10000 );
for ii=1:T
     %calculating likelihoods for both models of monetary policymaking

    ldpce(ii)=-1/(2*sigma_measurement^2)*(i_vec(ii)+Fd*x1pce(:,ii))^2-log(sigma_measurement)-.5*log(2*pi);
    lcpce(ii)=-1/(2*sigma_measurement^2)*(i_vec(ii)-C(3,:)*[x1pce(:,ii);p2pce(:,ii)])^2-log(sigma_measurement)-.5*log(2*pi);

    if ii<40





denom=real((exp((lcpce(ii)+log(probspce(1,ii))))+exp((ldpce(ii)+log(probspce(2,ii))))));
if denom==0;
    probspce(1,ii+1)=probspce(1,ii);
    probspce(2,ii+1)=probspce(2,ii);
else
probspce(1,ii+1)=real(exp(lcpce(ii)+log(probspce(1,ii)))/denom);
probspce(2,ii+1)=real(exp(ldpce(ii)+log(probspce(2,ii)))/denom);
end

    else
       ldtemp=sum(ldpce(ii-39:ii)); 
       lctemp=sum(lcpce(ii-39:ii)); 
       
       
       denom=real((exp((lctemp+log(p_c)))+exp((ldtemp+log(p_d)))));
       
       
       if denom==0;
    probspce(1,ii+1)=probspce(1,ii);
    probspce(2,ii+1)=probspce(2,ii);
    
       else
probspce(1,ii+1)=real(exp(lctemp+log(p_c))/denom);
probspce(2,ii+1)=real(exp(ldtemp+log(p_d))/denom);

      end

       
       
    end
%calculate conditional expectations
E_dpce(:,ii)=Cd*Md^2*x1pce(:,ii);
E_cpce(:,ii)=C(1:2,:)*M^2*[x1pce(:,ii);p2pce(:,ii)];




%calculate new shocks
temp=probspce(1,ii+1)*(A_21(1:2,1:2)\(H*(E_cpce(:,ii))-B_2*i_vec(ii+1)-A_22*x2pce(:,ii)-A_21(1:2,3:end)*x1pce(3:end,ii+1)))+probspce(2,ii+1)*(A_21_d(1:2,1:2)\(H*(E_dpce(:,ii))-B_2*i_vec(ii+1)-A_22*x2pce(:,ii)-A_21_d(1:2,3:end)*x1pce(3:end,ii+1)));
x1pce(1:2,ii+1)=temp;
%calculate new end. variables
p2pce(:,ii+1)=M(7:end,:)*[x1pce(:,ii+1);p2pce(:,ii)];


end




