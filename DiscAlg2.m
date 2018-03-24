function [M,C,Fs,Vs] = DiscAlg2( A11,A12,A21,A22,B1,B2,Q11,Q12,Q21,Q22,R,U1,U2,bet,Ct1,Vst1 );
%DiscAlg2   Solves the LQ problem under discretion in t, given C(t+1) and V(t+1).
%
%
%
%
%  Usage:     [M,C,Fs,Vs] = DiscAlg2( A,B,Q,R,U,bet,n1,n2,Ct1,Vst1 );
%
%  Input:     A,B,Q,R,U,bet,n1,n2: see ComItAlg
%             Ct1            n2xn1  matrix, C(t+1) in x2(t+1) = C(t+1)x1(t+1)
%             Vst1           n1xn1  matrix, V(t+1) value function matrix in t+1
%
%  Output:    M              n1xn1 matrix,
%             C              n2xn1 matrix,
%             Fs             kxn1 matrix,
%             Vs             n1xn1 matrix,
%
%
%
%
%
%  Paul Söderlind, Paul.Soderlind@unisg.ch, Aug 2000
%-----------------------------------------------------------------------





d1mat = eye(size(A22))/(A22-Ct1*A12);
D = d1mat * (Ct1*A11-A21);
G = d1mat * (Ct1*B1-B2);

As = A11 + A12*D;
Bs = B1 + A12*G;

Qs = Q11 + Q12*D + D'*Q21 + D'*Q22*D;
Us = Q12*G + D'*Q22*G + U1 + D'*U2;
Rs = R + G'*Q22*G + G'*U2 + U2'*G;

Fs = eye(size(Rs))/(Rs + bet*Bs'*Vst1*Bs) * (Us' + bet*Bs'*Vst1*As); %u(t) = -F*x1(t)
Vs = Qs - Us*Fs - Fs'*Us' + Fs'*Rs*Fs + bet*(As-Bs*Fs)'*Vst1*(As-Bs*Fs);

C = D - G*Fs;                   %x2(t)=C*x1(t)
M = A11 + A12*C - B1*Fs;        %x1(t+1) = M*x1(t) + e(t+1)
%-----------------------------------------------------------------------

