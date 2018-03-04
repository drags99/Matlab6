%solves 2d poisson equation 

%defining paramaters
iterate=10;
Lx=pi;
Ly=pi;
u0=0;
uL=0;
v0=0;
vL=0;
M=1;
N=10;%equispaced points in x and y
%deltas
delx=Lx/(N+1);
dely=Ly/(N+1);
x=delx:delx:Lx-delx;
y=dely:dely:Ly-dely;

f=zeros(N,N);%preallocate f(x,y)

for i=1:N
    for j=1:N
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

