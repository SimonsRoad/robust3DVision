function [result] = kronDelta(x,y)
%KRONDELTA Summary of this function goes here
%   kronecker delta function for doubles
if x==y
    result = 1;
else
    result = 0;
end

end

