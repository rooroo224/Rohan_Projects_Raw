function flux = hll_solv(L,R)
rho = [L(1) R(1)];
vel = [L(2)/L(1) R(2)/R(1)];
pre = [0.4*(L(3)-((L(2)^2)/2/L(1))) 0.4*(R(3)-((R(2)^2)/2/R(1)))];
c = sqrt(1.4*pre./rho); %%speed of sound

lambda1 = vel-c;    %%lowest eigenvalue
lambda2 = vel+c;    %%highest eigenvalue

a_l = min(lambda1);
a_r = max(lambda2);

act_flux_l = [L(2);rho(1)*vel(1)^2+pre(1);vel(1)*(L(3)+pre(1))];    %%flux function as per equation
act_flux_r = [R(2);rho(2)*vel(2)^2+pre(2);vel(2)*(R(3)+pre(2))];

%%calculating numerical flux
if (a_l>=0)
    flux = act_flux_l;
elseif(a_l<0 && a_r>0)
    flux = (a_r*act_flux_l - a_l*act_flux_r + a_l*a_r*(R-L))/(a_r-a_l);
elseif (a_r<=0)
    flux = act_flux_r;
end