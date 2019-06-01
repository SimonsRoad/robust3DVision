function [v1, v2, t, R, outliers_idx] = create2ViewCentralExperiment(num_points, varargin)
%   Author: Heng Yang hankyang@mit.edu
%   Last update: Feb 02, 2019
%   Adapted from opengv matlab function create2D2DExperimen.m
%   Sample only 3D points that are within FOV of the first camera

params = inputParser;
params.CaseSensitive = false;
params.addParameter('Noise', 0, @(x) isscalar(x) && isfinite(x)); % in px
params.addParameter('FOV', 100, @(x) isscalar(x) && isfinite(x)); %in deg
params.addParameter('MaxParallax', 2, @(x) isscalar(x) && isfinite(x)); %in meters
params.addParameter('FocalLength', 800, @(x) isscalar(x) && isfinite(x)); % in pixel
params.addParameter('OutliersRatio', 0, @(x) isscalar(x) && ...
    isfinite(x) && x>=0 && x<1); 
params.addParameter('OutlierThreshold', 0.01, @(x) isscalar(x) && ...
    isfinite(x) && x>=0); 
params.addParameter('PlotCameras', false, @(x) isscalar(x) && islogical(x));
params.addParameter('PlotResiduals', false, @(x) isscalar(x) && islogical(x));
params.parse(varargin{:});

% Aliases
noise = params.Results.Noise;
fov = params.Results.FOV;
max_parallax = params.Results.MaxParallax;
outlier_fraction = params.Results.OutliersRatio;
focal_length = params.Results.FocalLength;
outlier_threshold = params.Results.OutlierThreshold;

%% generate random view-points
max_rotation = 0.5;

% position and orientation of the first camera (w.r.t. global)
position1 = zeros(3,1);
rotation1 = eye(3);
% position and orientation of the second camera (w.r.t. global)
% t_c2^c1; R_c2^c1
position2 = max_parallax * 2.0 * (rand(3,1) - repmat(0.5,3,1));
rotation2 = generateBoundedR(max_rotation);

%% Generate random point-cloud within FOV
minDepth = 1.0;
maxDepth = 5.0;
deg2rad = pi/180.0;

nSamples = 0;
points_c1 = zeros(3, num_points);
points_c2 = zeros(3, num_points);
pose_c2_c1 = [rotation2, position2; zeros(1,3), 1]; % convert point from c2 to c1: T_c2^c1
while nSamples < num_points
    % sample inside a bounded cone within FOV of the first camera
    alpha = -pi + 2*pi*rand(1); % rotation about z
    theta = (-fov/2 + fov*rand(1))*deg2rad; % deviation from z
    direction = [sin(theta)*cos(alpha); sin(theta)*sin(alpha); cos(theta)];
    depth = minDepth + (maxDepth - minDepth) * rand(1);
    point_c1 = direction * (depth / direction(end)); % make point's z = depth
    % convert point from first camera to second camera
    point_c2_h = pose_c2_c1 \ [point_c1; 1];
    point_c2 = point_c2_h(1:3);
    % check if point in the cone of second camera 
    theta_2 = atan2(sqrt(point_c2(1)^2 + point_c2(2)^2), point_c2(3));
    if (-fov/2*deg2rad < theta_2 && theta_2 < fov/2*deg2rad)
        nSamples = nSamples + 1;
        points_c1(:, nSamples) = point_c1;
        points_c2(:, nSamples) = point_c2;
    end
end

if params.Results.PlotCameras
%     disp('points')
%     disp(points)
%     disp('points in c2')
%     disp(points_c2)
    figure
    scatter3(points_c1(1,:), points_c1(2,:), points_c1(3,:), 'red', 'filled')
    hold on
    scatter3(position1(1), position1(2), position1(3), 'black', '*')
    quiver3(position1(1), position1(2), position1(3), ...
    0.5, 0, 0, 0, 'green', 'LineWidth', 3)
    quiver3(position1(1), position1(2), position1(3), ...
    0, 0.5, 0, 0, 'blue', 'LineWidth', 3)
    quiver3(position1(1), position1(2), position1(3), ...
    0, 0, 0.5, 0, 'red', 'LineWidth', 3)
    text(position1(1), position1(2), position1(3), 'Camera 1')
    
    c2_coords = [0.5,0,0;0,0.5,0;0,0,0.5];
    c2_coords_in_c1 = pose_c2_c1 * [c2_coords;1,1,1];
%     disp(c2_coords)
%     disp(c2_coords_in_c1)
%     disp(position2)
    scatter3(position2(1), position2(2), position2(3), 'black', '*')
    quiver3(position2(1), position2(2), position2(3), ...
    c2_coords_in_c1(1,1)-position2(1), ...
    c2_coords_in_c1(2,1)-position2(2), ...
    c2_coords_in_c1(3,1)-position2(3), ...
    0, 'green', 'LineWidth', 3)
    quiver3(position2(1), position2(2), position2(3), ...
    c2_coords_in_c1(1,2)-position2(1), ...
    c2_coords_in_c1(2,2)-position2(2), ...
    c2_coords_in_c1(3,2)-position2(3), ...
    0, 'blue', 'LineWidth', 3)
    quiver3(position2(1), position2(2), position2(3), ...
    c2_coords_in_c1(1,3)-position2(1), ...
    c2_coords_in_c1(2,3)-position2(2), ...
    c2_coords_in_c1(3,3)-position2(3),...
    0, 'red', 'lineWidth', 3)
    text(position2(1), position2(2), position2(3), 'Camera 2')
    axis equal
end

%% Now create the correspondences by looping through the cameras

v1 = zeros(3,num_points);
v2 = zeros(3,num_points);

for i=1:num_points
    bearingVector1 = points_c1(:,i)/norm(points_c1(:,i));
    bearingVector2 = points_c2(:,i)/norm(points_c2(:,i));
    
    % add noise to the bearing vectors here
    bearingVector1_noisy = addNoise(bearingVector1,focal_length,noise);
    bearingVector2_noisy = addNoise(bearingVector2,focal_length,noise);
    
    % store the normalized bearing vectors along with the cameras they are
    % being seen (we create correspondences that always originate from the
    % same camera, you can change this if you want)
    bearingVector1_noisy_norm = norm(bearingVector1_noisy);
    bearingVector2_noisy_norm = norm(bearingVector2_noisy);
    
    v1(:,i) = bearingVector1_noisy/bearingVector1_noisy_norm;
    v2(:,i) = bearingVector2_noisy/bearingVector2_noisy_norm;
end

if params.Results.PlotCameras
    scale = 5;
    for i=1:num_points
        v2_rot = rotation2 * v2(:,i);
        quiver3(position1(1), position1(2), position1(3),...
            scale*v1(1,i), scale*v1(2,i), scale*v1(3,i), 0, 'yellow')
        quiver3(position2(1), position2(2), position2(3),...
            scale*v2_rot(1), scale*v2_rot(2), scale*v2_rot(3), 0, 'yellow')
    end
end

%% compute relative translation and rotation
R = rotation2;
t = (position2 - position1);
t_gt=t/norm(t);
r_gt=reshape(R,[9,1]);

%% Add outliers
MAX_TENTATIVE = 10000;
tentative = 1;
number_outliers = floor(outlier_fraction*num_points);
outliers_idx = shuffle(1:num_points);
outliers_idx = outliers_idx(1:number_outliers);
failed_outliers = [];
if number_outliers > 0
    for i = outliers_idx
        outlier_residual = 0;
        while outlier_residual < params.Results.OutlierThreshold
            direction = randn(3,1); 
            direction = direction / norm(direction);
            f1=v1(:,i);
            A=kron(direction',hatmap(f1));
            outlier_residual=(t_gt'*A*r_gt)^2; 
            tentative = tentative + 1;
            if(tentative > MAX_TENTATIVE)
                warning('Failed to create an outlier')
                failed_outliers = [failed_outliers i];
            break
            end
        end        
        tentative = 0;
        v2(:,i) = [direction];
    end
end
outliers_idx = setdiff(outliers_idx, failed_outliers);

%% calculate residuals
residuals = zeros(1, num_points);
for i=1:num_points
    f2=v2(:,i); f1=v1(:,i);
    A=kron(f2',hatmap(f1));
    residuals(i)=(t_gt'*A*r_gt)^2;
end

if params.Results.PlotResiduals
    figure
    plot(residuals, 'linewidth', 2);
    title('Residuals from pairs of bearing vectors.')
end

end
