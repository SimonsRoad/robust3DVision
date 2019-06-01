function [v] = vmap(hat_v)
%VMAP convert [t]x to t
v=[-hat_v(2,3);
    hat_v(1,3);
    -hat_v(1,2)];
end

