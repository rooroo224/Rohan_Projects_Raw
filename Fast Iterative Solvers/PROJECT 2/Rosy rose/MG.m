function [uout,rmaxl,count] = MG(l, uin, fl, gamma, nu1, nu2, count, rmaxl)

count = count + 1;
[len,~] = size(uin);
N = len - 1;
Nc = N / 2;
h = 1 / N;
rh = zeros(size(uin)); 

ubar = GS(uin, fl, nu1);
rmax = 0;
for i = 2 : N  % Calculate residual rh for fine grid
    for j = 2 : N
        rh(i, j) = fl(i,j) + ((ubar(i-1,j)-2*ubar(i,j)+ubar(i+1,j))/(h^2)) + ((ubar(i,j-1)-2*ubar(i,j)+ubar(i,j+1))/(h^2)) ;
        if rmax < abs(rh(i, j))
            rmax = abs(rh(i, j));
        end
    end
end
rmaxl(count) = rmax;
rH = RESTR(rh, Nc);  % Residual rH for coarse grid. Restriction: fine to coarse.

if l == 1
    etilde0 = zeros(Nc + 1);
    etildeH = GS(etilde0, -rH, 1);
else
    etildeH = zeros(Nc + 1);  % Error etilda_H for coarse grid
    for i = 1 : gamma
        [etildeH, rmaxl,count] = MG(l-1, etildeH, -rH, gamma, nu1, nu2, count, rmaxl);
    end
end

etildeh = PROLONG(etildeH, Nc); % Error etilda_h for fine grid. Prolongation: coarse to fine.
uout = GS(ubar - etildeh, fl, nu2);

end

