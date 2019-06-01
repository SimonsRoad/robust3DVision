function v_outlier = addOutlier(v_clean,focal_length,outlier_level,theta)
%addOutlier ADD OUTLIER GIVEN BEARING VECTOR
%find good vector in normal plane based on good conditioning
    %this inPlaneVector1 is perpendicular to v_clean
    inPlaneVector1 = zeros(3,1);
    
    if v_clean(1,1) > v_clean(2,1) && v_clean(1,1) > v_clean(3,1)
        inPlaneVector1(2,1) = 1.0;
        inPlaneVector1(3,1) = 0.0;
        inPlaneVector1(1,1) = 1.0/v_clean(1,1) * (-v_clean(2,1));
    else
        if v_clean(2,1) > v_clean(3,1) && v_clean(2,1) > v_clean(1,1)
            inPlaneVector1(3,1) = 1.0;
            inPlaneVector1(1,1) = 0.0;
            inPlaneVector1(2,1) = 1.0/v_clean(2,1) * (-v_clean(3,1));
        else
            %v_clean(3,1) > v_clean(1,1) && v_clean(3,1) > v_clean(2,1)
            inPlaneVector1(1,1) = 1.0;
            inPlaneVector1(2,1) = 0.0;
            inPlaneVector1(3,1) = 1.0/v_clean(3,1) * (-v_clean(1,1));
        end
    end
    
    %normalize the in-plane vector
    inPlaneVector1 = inPlaneVector1 / norm(inPlaneVector1);
    inPlaneVector2 = cross(v_clean,inPlaneVector1);
    
    noiseX = outlier_level * cos(theta);% / sqrt(2);
    noiseY = outlier_level * sin(theta);% / sqrt(2);
    
    v_outlier = focal_length * v_clean + noiseX * inPlaneVector1 + noiseY * inPlaneVector2;
    
    v_outlier_norm = norm(v_outlier);
    v_outlier = v_outlier ./ v_outlier_norm;
end

