function [M,C,V,F] = DiscAlg( A,B,Q,R,U,bet,n1,n2,Vt1,Ct1,...
                              ConvCrit,Vweight,Fweight,CritLags,step,PrintIt,MaxIter );
%DiscAlg    Solves the LQ problem under discretion, iterating backwards in time.
%
%
%
%  Usage:     [M,C,V,F] = DiscAlg( A,B,Q,R,U,bet,n1,n2,Vt1,Ct1,...
%                                  ConvCrit,Vweight,Fweight,CritLags,step,PrintIt,MaxIter );
%
%  Input:     A,B,Q,R,U,bet,n1,n2:   - see ComItAlg
%             Vt1        n1xn1 matrix: initial guess of value matrix
%             Ct1        n2xn1 matrix, initial guess of C in x2(t)=C*x1(t)
%             ConvCrit   2x1 convergence criteria for abs(V-Vt1)|abs(F-Ft1)
%             Vweight    scalar or n1xn1 matrix with weights for
%                        difference criterion for V
%             Fweight    scalar or (n1+n2)x1 vector with weights for
%                        difference criterion for F
%             CritLags   #lags of CritVar compared with ConvCrit
%             step       scalar in (0,1): factor of updating of V and F as
%                        in Vt1 = step*V + (1-step)*Vt1
%             PrintIt    1: printing iteration number
%                        2: printing  iteration number and convergence criteria
%             MaxIter    scalar, maximum number of iterations (eg. 10000). NEW (March 2003)
%
%  Output:    M        n1x1 matrix, x1(t+1) = M*x1(t) + e(t+1)
%             C        n2xn1 matrix, x2(t)  = C*x1(t)
%             V        n1xn1 matrix, value function is x1(t)'*V*x1(t)
%             F        kxn1 matrix, decision rule is u(t) = -F*x1(t), where
%                      k is number of elements in u(t)
%
%
%  Remark:    (For Octave users) If PrintIt is 1 or 2, then set
%             page_output_immediately = 1;
%             This seems to make Octave work better.
%
%  Calls on:  DiscAlg2
%
%
%  Paul Söderlind, Paul.Soderlind@unisg.ch, Aug 2000, Mar 2003
%-----------------------------------------------------------------------

Q = (Q + Q')/2;                %to make symmetric
R = (R + R')/2;

n = n1 + n2;
Ft1 = 1000;




A11 = A(1:n1,1:n1);
A12 = A(1:n1,(n1+1):n);
A21 = A((n1+1):n,1:n1);
A22 = A((n1+1):n,(n1+1):n);
Q11 = Q(1:n1,1:n1);
Q12 = Q(1:n1,(n1+1):n);
Q21 = Q((n1+1):n,1:n1);
Q22 = Q((n1+1):n,(n1+1):n);

B1 = B(1:n1,:);
B2 = B(n1+1:n,:);
U1 = U(1:n1,:);
U2 = U(n1+1:n,:);



Cdiff = 1000*ones(1,2);
iter = 1;
while any( max(Cdiff) > (ConvCrit')) & (iter < MaxIter);   %iterations

  [M,Ct1,F,V] = DiscAlg2( A11,A12,A21,A22,B1,B2,Q11,Q12,Q21,Q22,R,U1,U2,bet,Ct1,Vt1 ); %solve period t

  Vdiff = max(  max( Vweight.*abs(V-Vt1) )  );        %changes t+1 -> t
  Fdiff = max(  max( Fweight.*abs(F-Ft1) )  );
  Cdiff = [ Vdiff, Fdiff ];
                  


  Vt1 = step*V + (1-step)*Vt1;                        %"downdating"
  
  Ft1 = step*F + (1-step)*Ft1;


  if PrintIt == 1;
    disp(iter);
  elseif PrintIt == 2;
    disp([iter, (max(Cdiff))]);
    elseif PrintIt == 0;
  end;

  iter = iter + 1;

end;                                 %end iterations

if iter >= MaxIter;
  warning('Maximum number of iterations reached');
end;
C=Ct1;
%-----------------------------------------------------------------------

