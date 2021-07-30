dx = 0.01;
x = 0:dx:1;
dt = dx/5;
t = 0:dt:0.25;

%%initialisation
u = zeros(3,length(x),length(t));

true_vars = zeros(3,length(x),length(t));
for i=1:1:length(x)
    true_vars(:,i,1)=init1(x(i));
end
E = zeros(length(x),length(t));
for i=1:1:length(x)
    E(1,i) = (true_vars(3,i,1)/0.4) + (true_vars(1,i,1)*(true_vars(2,i,1)^2)/2);
end

for i=1:1:length(x)
    u(:,i,1)=[true_vars(1,i,1);true_vars(2,i,1)*true_vars(1,i,1);E(1,i)];
end
%%end of initialisation

%%time-stepping
for i=1:1:length(t)-1
    u(:,1,i+1)=u(:,1,i)-(dt/dx)*(hll_solv(u(:,1,i),u(:,2,i))-hll_solv(init1(0),u(:,1,i)));
    for j=2:1:length(x)-1
        u(:,j,i+1)=u(:,j,i)-(dt/dx)*(hll_solv(u(:,j,i),u(:,j+1,i))-hll_solv(u(:,j-1,i),u(:,j,i)));
    end
    u(:,length(x),i+1)=u(:,length(x),i)-(dt/dx)*(hll_solv(u(:,length(x),i),init1(1))-hll_solv(u(:,length(x)-1,i),u(:,length(x),i)));
end
%%end of time-stepping

true_vars(1,:,:) = u(1,:,:); %%density(rho)
true_vars(2,:,:) = u(2,:,:)./u(1,:,:);  %%velocity(v)
true_vars(3,:,:) = 0.4*(u(3,:,:)-((u(2,:,:).^2)/2./u(1,:,:)));  %%pressure(p)
E = u(3,:,:);   %%energy(E)

%%plotting
for i=1:1:length(t)
    plot(x,true_vars(1,:,i)),axis([0 1 -2 2])
    F1(i)=getframe;
end

for i=1:1:length(t)
    plot(x,true_vars(2,:,i)),axis([0 1 -2 2])
    F2(i)=getframe;
end

for i=1:1:length(t)
    plot(x,true_vars(3,:,i)),axis([0 1 -2 2])
    F3(i)=getframe;
end

for i=1:1:length(t)
    plot(x,E(1,:,i)),axis([0 1 -3 3])
    F4(i)=getframe;
end
%%end of plotting

%%initial conditions
function u0 = init1(z)
if (z<=0.5)
    u0 = [1.0;0;1];
else
    u0 = [0.125;0;0.1];
end
end

function u0 = init2(z)
if (z<=0.5)
    u0 = [1.0;-1.0;0.4];
else
    u0 = [1;1;0.4];
end
end
