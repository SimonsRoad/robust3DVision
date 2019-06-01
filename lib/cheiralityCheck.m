function pass = cheiralityCheck(R,t,bearing1,bearing2)
%CHEIRALITYCHECK Given beairng vectors 1 and 2, check cheirality constraint

if size(bearing1,2) ~= size(bearing2,2)
    error('Bearing vectors must have same # of columns!')
end
NUM_POINTS = size(bearing1,2);
pass=true;
for i=1:NUM_POINTS
    f1=bearing1(:,i);
    f2=bearing2(:,i);
    A = [f1, -R*f2];
    d = A\t;
%     disp(d)
    if d(1)*d(2)<0
        pass=false;
    end
end
if R(1,1) < 0
    pass=false;
end

end

