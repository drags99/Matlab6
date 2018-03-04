%solves 2d poisson equation 

%defining paramaters
Lx=pi;
Ly=pi;
u0=0;
uL=0;
v0=0;
vL=0;
M=10;%equispaced points in x
N=10;%equispaced points in y
%space between
delx=Lx/(M+1);
dely=Ly/(N+1);
x=0:delx:Lx;
y=0:dely:Ly;

