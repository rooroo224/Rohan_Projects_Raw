N = 1001;
x = linspace(-3,2,N);
dx = 5/(N-1);

T_end = 0.2;
dt = dx/2;
t = 0:dt:T_end; %%as per CFL

if (mod(length(t),2) == 0)
    t_sample = length(t)/2;
else t_sample = (length(t)+1)/2;
end

u = sin(pi*x);

% fun = @(f) 2*f; %%linear advection
fun = @(f) 0.5*f*f; %%burgers equation

% flux = @(l,r) 0.5*(fun(l)+fun(r)); %%naive
% flux = @(l,r) (0.5*(fun(l)+fun(r))) - (dx*(r-l)/2/dt); %%Lax-Friedrichs
% flux = @(l,r) 0.5*(fun(l)+fun(r)) - dt*4*(r-l)/2/dx; %%Lax-Wendroff advec
flux = @(l,r) 0.5*(fun(l)+fun(r)) - dt*((l+r)^2)*(r-l)/8/dx; %%Lax-Wendroff for burgers

F(length(t))=struct('cdata',[],'colormap',[]);
for i=1:1:length(t)
    u(i+1,1)=u(i,1)+ (dt*(flux(0,u(i,1))-flux(u(i,1),u(i,2)))/dx);
    for j=2:1:length(x)-1
        u(i+1,j)=u(i,j) + (dt*(flux(u(i,j-1),u(i,j))-flux(u(i,j),u(i,j+1)))/dx);
    end
    u(i+1,N)=u(i,N)+ (dt*(flux(u(i,N-1),u(i,N))-flux(u(i,N),0))/dx);
    plot(x,u(i,:)),axis([-1 1 -1 2])    
    F(i)=getframe;
end

% F1(length(t))=struct('cdata',[],'colormap',[]);
% u1_ex=zeros(length(t),length(x));  %%exact advection solution
% for i=1:1:length(t)
%     for j=1:1:length(x)
%     u1_ex(i,j)=sin(pi*(x(j)-2*t(i)));
%     end
%     plot(x,u1_ex(i,:)),axis([-1 1 -1 2])
%     F1(i)=getframe;
% end

u1_ex = zeros(length(t),length(x));
[xmin,u0min]=fminbnd(@(x) sin(pi*x),-1,1);
[xmax,u0max]=fminbnd(@(x) sin(pi*x),-1,1);
u0max=-u0max;
uEx=@(x,t) fminbnd(@(u) (u-sin(pi*(x-u*t))).^2,u0min,u0max);

F2(length(t))=struct('cdata',[],'colormap',[]);
for i=1:1:length(t)
    for j=1:1:length(x)
        u1_ex(i,j)=uEx(x(j),t(i));
    end
    plot(x,u1_ex(i,:)),axis([-1 1 -1 2])
    F2(i)=getframe;
end

L1_error = 0;
    for j =1:1:length(x)
        L1_error = L1_error + abs(u(t_sample,j)-u1_ex(t_sample,j));
    end

L1_error = L1_error*dx;

err2 = 0;

    for j =1:1:length(x)-1
        err2 = err2 + abs(u(t_sample,j)-((u1_ex(t_sample,j)+u1_ex(t_sample,j+1))/2));
    end

err2 = err2*dx;

