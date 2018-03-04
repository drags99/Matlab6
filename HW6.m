%solves 2d poisson equation 
clear;
clc;

%defining paramaters
iterate=20000;
Lx=pi;
Ly=pi;
u0=0;
uL=0;
v0=0;
vL=0;
M=1;
N=100;
%deltas
delx=Lx/(N+1);
dely=Ly/(N+1);
x=0:delx:Lx;
y=0:dely:Ly;

f=zeros(N+2,N+2);%preallocate f(x,y)

for i=1:N+2
    for j=1:N+2
    f(i,j)=-2*M*sin(M*x(i))*cosh(M*y(j));
    end
end

% initial values of matrix u
u=zeros(N+2,N+2);
u(1,:)=u0;
u(N+2,:)=uL;
u(:,1)=v0;
u(:,N+2)=vL;


%Gauss Siedel
for z=1:iterate
    for i=2:N+1
        for j=2:N+1
            u(i,j)=(1/4)*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1))-((delx^2)/4)*f(i,j);
        end
    end
end


%exact solution for comparison
exactu=zeros(N+2,N+2);
for i=1:N+2
    for j=1:N+2
        exactu(i,j)=(Lx-y(j))*sin(M*x(i))*sinh(M*y(j));
    end
end

%initial L1 value
L1=0;
numpoints=0;

%sum errors
for i=2:N+1
    for j=2:N+1
        L1=L1+abs(u(i,j)-exactu(i,j));
        numpoints=numpoints+1;
    end
end
%dividing by number of points
L1=L1/(numpoints)

%graphs
for z=1:((N+2)*(N+2))
    for i=2:N+1
        for j=2:N+1
            vectoraproxu(z)=u(i,j);
            vectorexactu(z)=exactu(i,j);
        end
    end
end


