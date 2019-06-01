function [R_est,t_est,nSteps] = eightpt_ransac(bearing1,bearing2,...
    R_gt,t_gt,inlierResThresh,maxItr)
%EIGHTPT_RANSAC Call opengv eight point method in RANSAC
NUM_POINTS=size(bearing1,2);
if NUM_POINTS < 10
    error('8-POINT RANSAC requires at least 10 points')
end

fullIndex=1:NUM_POINTS;
inlierThresh=0.5;
inlierRatio=0;
itr=0;
nPoints=8;

t_gt=t_gt/norm(t_gt);
while inlierRatio < inlierThresh && itr < maxItr
    randIndex=datasample(fullIndex,nPoints);
    % call opengv eight point method
    X_eight_pt=opengv('eightpt',randIndex,bearing1,bearing2);
    if size(X_eight_pt,3)==0
        disp('Warning: opengv 8-point returns empty solution, skipping.')
    else
        X_eight_pt(:,:,2)=-X_eight_pt(:,:,1);
        % transform essential to rotation matrices
        R_eight_pt=transformEssentials(X_eight_pt);
        % take the best R compared to gt
        temp=zeros(1,size(R_eight_pt,3));
        for k=1:size(R_eight_pt,3)
            temp(k)=evaluateRotationError(R_gt,R_eight_pt(:,:,k));
        end
        [~,minIdx]=min(temp);
        R_eight_pt_best=R_eight_pt(:,:,minIdx);
        E_eight_pt_best=X_eight_pt(:,:,ceil(minIdx/2));
        % use best R to get t
        t_eight_pt_best=vmap(E_eight_pt_best*R_eight_pt_best');
        t_eight_pt_best=t_eight_pt_best/norm(t_eight_pt_best);
        if norm(-t_eight_pt_best-t_gt) < norm(t_eight_pt_best-t_gt)
            t_eight_pt_best=-t_eight_pt_best;
        end
        % use best R and t to calculate residuals and # of inliers
        r_eight_pt_best=reshape(R_eight_pt_best,[9,1]);
        residuals=zeros(1,NUM_POINTS);
        for pt=1:NUM_POINTS
            f1=bearing1(:,pt);
            f2=bearing2(:,pt);
            A=kron(f2',hatmap(f1));
            residuals(pt)=(t_eight_pt_best'*A*r_eight_pt_best)^2;
        end
        inlierRatio=sum(residuals<inlierResThresh)/NUM_POINTS;
        inlierIndex=find(residuals<inlierResThresh);
    end
    itr=itr+1; 
end
nSteps=itr;
if nSteps > (maxItr-5)
    disp('WARN: EIGHT POINT RANSAC did not converge within maximum iterations')
end
% use nonlinear optimization to refine result
X_8pt_nonlin_opt=opengv('rel_nonlin_central',inlierIndex,...
            bearing1,bearing2,[R_eight_pt_best,t_eight_pt_best]);
R_est=X_8pt_nonlin_opt(:,1:3);
t_est=X_8pt_nonlin_opt(:,4)/norm(X_8pt_nonlin_opt(:,4));

end

