%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                      %%%
%%% Wave equation 1-D, Boundary cond.    %%%
%%%                                      %%%
%%% To test Implicit SBP                 %%%
%%%                                      %%%
%%% Date: Jan 10 2022                    %%%
%%% Author: Ken Mattsson                 %%%
%%%                                      %%%
%%%                                      %%%
%%%                                      %%%
%%%  u_tt     =  c^u_xx,  x_l<= x <=x_r  %%%
%%%                                      %%%
%%%  Test Imlicit SBP, Neumann BC        %%%
%%%                                      %%%
%%% Initial data                         %%%
%%% u(x,0)=f_1   u_t(x,0)=0              %%%
%%%                                      %%%
%%% As an analytic solution we use a     %%%
%%% Gaussian profile                     %%%
%%%                                      %%%
%%% Solves eq (19) in JCP paper          %%% 
%%%                                      %%%
%%% Se chapter 4.3 for details           %%%
%%%                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;


rr=0.1;               % Width of Gaussian
x_l=-1;x_r=1;         % The boundaries of the domain
c=1;                  % Wave speed

L=x_r-x_l;
theAxes=[x_l x_r -0.7 1.1]; % Regarding the figure

m_start=51;                 % Number of grid-points, first grid
antal=3;                    % Number of grid-refinements
t_1=2.2;                    % End time

m_all=zeros(antal,1);
h_all=zeros(antal,1);
m_all(1)=m_start;
h_all(1)=L/(m_all(1)-1);

% Convergens result
conv=zeros(antal-1);
differens=zeros(antal-1);


for j=2:antal
    if (j==2)
        m_all(j)=50+m_all(j-1);
    elseif (j==3)
        m_all(j)=100+m_all(j-1);
    elseif (j==4)
        m_all(j)=200+m_all(j-1);
    else
        m_all(j)=400+m_all(j-1);
    end
    h_all(j)=L/(m_all(j)-1);
end

for i=1:antal

    m=m_all(i);
    h=h_all(i);
    dt=0.1/c*h;

    [ST,MM,BD,QQ,H,xx,h] = SBP_SL6(m, L); % Construct SL6 SBP
    ordningstyp=' Sixth order Spectral';


    [LL, UU] = lu(sparse(H));

    HL=sparse(MM-BD); % Neumann BC

    max_itter=floor(t_1/dt);        % Antal itterationermax_itter

    n=2*m;

    V=zeros(n,1);                       % Numerical solution

    %x=linspace(x_l,x_r,m)';	        % discrete x values
    x=xx+x_l;

    %%% Initialize

    t=0.0;

    uc1=exp(-((x-c*t)/rr).^2);

    V(1:m)=uc1;



    %%%%%%%%%%%%% RK4 time integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for nr_itter=1:max_itter

        w1=RHS(LL,UU,HL,V,m);
        w2=RHS(LL,UU,HL,V+dt/2*w1,m);
        w3=RHS(LL,UU,HL,V+dt/2*w2,m);
        w4=RHS(LL,UU,HL,V+dt*w3,m);

        V=V+dt/6*(w1+2*w2+2*w3+w4);

        t=t+dt;

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%% Compute the error
    tt=rem(t,2)-2; % This analytic holds for t>1, t<2;

    uc1=exp(-((x-c*t+L)/rr).^2)';
    uc2=-exp(-((x+c*t-L)/rr).^2)';

    exact=uc1/2-uc2/2;

    differens(i)=sqrt(h)*norm(V(1:m)-exact');


end


for j=1:antal-1
    conv(j,:)=log10(differens(j)./differens(j+1))/log10(h_all(j)/h_all(j+1));
end


disp(' ');
disp(' ');
disp('____________________________________________________________________');
disp(' ');
disp([sprintf('m =%3d   ', m_all),'   ',ordningstyp]);
disp('____________________________________________________________________');
disp(' ');
disp(' ');

disp(' ');
disp('------------------');
disp('  m       L2-error     Convergence ');
for j=1:antal
    if j==1
        disp(['$',sprintf('%3d',m_all(j)),'$ ',...
            sprintf('& %1.2f ',log10(differens(j))),sprintf('& %1.2f ',0),'\\']);
    else
        disp(['$',sprintf('%3d',m_all(j)),'$ ',...
            sprintf('& %1.2f ',log10(differens(j))),sprintf('& %1.2f ',conv(j-1)),'\\']);
    end

end
disp(' ____________________________________________________________________');
disp(' ');
disp(' ');


figure(2);
plot(x,exact,'r',x,V(1:m),'b--','LineWidth',1);
%plot(x,V(1:m,1),'r','LineWidth',1);
title(['Numerical solution at t = ',num2str(t_1)]);
axis(theAxes);
grid;xlabel('x');
ax = gca;          % current axes
ax.FontSize = 16;

function ut = RHS(L,U,M,V,m)
y=L\(M*V(1:m,1));
ut2=U\y;
ut=[V(m+1:2*m,1);ut2];
end



