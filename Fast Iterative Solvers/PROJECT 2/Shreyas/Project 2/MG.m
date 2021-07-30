function [u] = MG(u, f, g, v1, v2)

u = GS(u,f,v1);

r_f = Res(u,f);

r_c = RESTR(r_f);

Nc = size(r_c,1);
e_c = zeros(Nc,Nc);

if Nc==3        % for coarsest level, only 1 interior node to solve exactly
    e_c = GS(e_c, -r_c, 1);
else
    for i=1:1:g
        e_c = MG(e_c, -r_c, g, v1, v2);
    end
end

e_f = PROLONG(e_c);

u = GS(u-e_f,f,v2);

end
                                        