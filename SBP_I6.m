function [S,MM,BD,QQ,H,x,h] = SBP_I6(m,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 6:te Implicit    SBP Finita differens   %%%
%%% operatorer framtagna av Ken Mattsson    %%%
%%%                                         %%%
%%% Datum: 2018 06 18                       %%% 
%%%                                         %%%
%%% 4 randpunkter, bandad   norm            %%%
%%%                                         %%%
%%%                                         %%%
%%% H           (Normen)                    %%%
%%% D1          (approx första derivatan)   %%%
%%% D2          (approx andra derivatan)    %%%
%%% S           Artificiell Dissipation     %%%
%%%                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% D1 har noggrannhet 3-6-3
% D2 har noggrannhet 2-6-2 (samma inre som vanliga SBP D2)

h=L/(m-1);
x=linspace(0,L,m);

H_U=[0.3707e4 / 0.17280e5 0.3689e4 / 0.17280e5 -0.467e3 / 0.3456e4 0.5e1 / 0.128e3; 0.3689e4 / 0.17280e5 0.6209e4 / 0.5760e4 0.83e2 / 0.1920e4 -0.671e3 / 0.17280e5; -0.467e3 / 0.3456e4 0.83e2 / 0.1920e4 0.4609e4 / 0.5760e4 0.533e3 / 0.3456e4; 0.5e1 / 0.128e3 -0.671e3 / 0.17280e5 0.533e3 / 0.3456e4 0.13627e5 / 0.17280e5;];

h2=-1/30;h1=2/15;h0=4/5;

H=h2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+h1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+h0*diag(ones(m,1),0);

H(1:4,1:4)=H_U;
H(m-3:m,m-3:m)=rot90(H_U,2);
%HI=1/h*inv(H);
H=H*h;


% First derivative SBP operator, 3rd order accurate at first 4 boundary points

q2=-1/12;q1=2/3;
Q=q2*(diag(ones(m-2,1),2) - diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));

Q_U = [0 0.263e3 / 0.360e3 -0.211e3 / 0.720e3 0.1e1 / 0.16e2; -0.263e3 / 0.360e3 0 0.43e2 / 0.48e2 -0.119e3 / 0.720e3; 0.211e3 / 0.720e3 -0.43e2 / 0.48e2 0 0.247e3 / 0.360e3; -0.1e1 / 0.16e2 0.119e3 / 0.720e3 -0.247e3 / 0.360e3 0;];
Q(1:4,1:4)=Q_U;
Q(m-3:m,m-3:m)=rot90(-Q_U,2);

e_l=zeros(m,1);e_l(1)=1;
e_r=zeros(m,1);e_r(m)=1;

QQ=Q-1/2*e_l*e_l'+1/2*e_r*e_r';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Second derivative, 2nd order accurate at first 6 boundary points
m3=1/45;m2=-1/20;m1=-1;m0=37/18;

M=m3*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3))+m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);


% Below with 5 boundary points There was one free parameters left that
% was optimised for leading truncation error
M_U=[0.37627e5 / 0.30240e5 -0.9551e4 / 0.6048e4 0.887e3 / 0.2016e4 -0.3613e4 / 0.30240e5 0.109e3 / 0.7560e4; -0.9551e4 / 0.6048e4 0.105151e6 / 0.30240e5 -0.22703e5 / 0.10080e5 0.2357e4 / 0.6048e4 -0.67e2 / 0.1890e4; 0.887e3 / 0.2016e4 -0.22703e5 / 0.10080e5 0.32917e5 / 0.10080e5 -0.15241e5 / 0.10080e5 0.23e2 / 0.630e3; -0.3613e4 / 0.30240e5 0.2357e4 / 0.6048e4 -0.15241e5 / 0.10080e5 0.14075e5 / 0.6048e4 -0.1999e4 / 0.1890e4; 0.109e3 / 0.7560e4 -0.67e2 / 0.1890e4 0.23e2 / 0.630e3 -0.1999e4 / 0.1890e4 0.15649e5 / 0.7560e4;];
M(1:5,1:5)=M_U;
M(m-4:m,m-4:m)=rot(M_U,2);

M=M/h;
    
d1_U=[-0.25e2 / 0.12e2 4 -3 0.4e1 / 0.3e1 -0.1e1 / 0.4e1;]/h;
d1_l=zeros(1,m);
d1_l(1:5)=d1_U;
d1_r=zeros(1,m);
d1_r(m-4:m)=fliplr(-d1_U);

BD=-e_l*d1_l+e_r*d1_r;

MM=-M+BD;

% Artificial dissipation = -HI*S (in RHS)
% and will lead to order 3-7-3

d4=0*[1 -4 6 -4 1];
DD_4=(diag(ones(m-2,1),2)-4*diag(ones(m-1,1),1)+6*diag(ones(m,1),0)-4*diag(ones(m-1,1),-1)+diag(ones(m-2,1),-2));
DD_4(1:2,1:5)=[d4;d4];DD_4(m-1:m,m-4:m)=[d4;d4];

DD=DD_4'*DD_4;
aa=max(max(DD));
S=-1/aa*DD;

