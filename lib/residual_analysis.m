
NOISE_LEVEL = 10; % in pixel
FIELD_OF_VIEW = 100; % in degress
NUM_POINTS = 6;
MAX_PARALLAX = 2; % in meters
OUTLIER_RATIO = 0.2;
DEBUG_SYNTHETIC = false; % true to turn plot on

NUM_ITERATIONS = 100;
residual=zeros(NUM_ITERATIONS,NUM_POINTS);
for i=1:NUM_ITERATIONS
    [bearing1, bearing2, t_gt, R_gt] = create2ViewCentralExperiment(...
        NOISE_LEVEL, NUM_POINTS, FIELD_OF_VIEW, MAX_PARALLAX,...
        OUTLIER_RATIO, DEBUG_SYNTHETIC);
    t_gt=t_gt/norm(t_gt);
    r_gt=reshape(R_gt,[9,1]);
    
    for k=1:NUM_POINTS
        f1=bearing1(:,k);
        f2=bearing2(:,k);
        A=kron(f2',hatmap(f1));
        residual(i,k)=(t_gt'*A*r_gt)^2;
    end
end
