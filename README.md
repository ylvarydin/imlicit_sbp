# imlicit_sbp
The SBP operators are constructed by the Matlab-functions



SBP_I6	6th order Implicit

SBP_I8	8th order Implicit

SBP_I10	10th order Implicit

SBP_I12	12th order Implicit

SBP_SL4	SL order 4

SBP_SL6	SL order 6



You call with the number of grid-points (m) and width of domain (L)  

Return parameters are 



ST 	Artificial dissipation (S=S’<=0)

MM	= (-M + BD) in second derivative SBP (M=M’>=0)

BD	Boundary derivative in second derivative SBP

QQ	Q+1/2*B first derivative SBP (Q+Q’=0)

H	The norm

xx	Grid-points (between 0 and L)

h	Grid-step

  

This is an example of a function call to construct the SL 6 operators



[ST,MM,BD,QQ,H,xx,h] = SBP_SL6(m, L); % Construct SL6 SBP



In the Matlab script “Wave_Implicit_JCP.m” we solve the second order wave equation with Neumann boundary conditions. Here using SBP_SL6

And RK4 for time-integration.



D2 operator = inv(H)*MM

D1 operator = inv(H)*QQ
