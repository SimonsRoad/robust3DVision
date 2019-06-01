function plotGNC(options)
% plot graduated non-convexity
% supported functions:
% Geman-McClure function
% ...

if ~isfield(options,'costFunction')
	options.costFunction='Geman-McClure';
end
if ~isfield(options,'initialMu')
	options.initialMu=100;
end
if ~isfield(options,'residualRange')
	options.residualRange=-1:0.1:1;
end
if ~isfield(options,'divFactor')
	options.divFactor=1.1;
end
if ~isfield(options,'numSteps')
	options.numSteps=10;
end

figure;
hold on
x=options.residualRange;
mu=options.initialMu;
divFactor=options.divFactor;
for i=1:options.numSteps
    mu=mu/(divFactor^(i-1));
    switch options.costFunction
    	case 'Geman-McClure'
    		y=mu*x.^2 ./ (mu + x.^2);
    		plot(x,y,'linewidth',2);
    	otherwise
    		error('Cost function not supported.');
    end    
end
grid on
title(sprintf('Graduated non-convexity: %s',options.costFunction))
xlabel('Residual')
ylabel('Cost')
hold off

end