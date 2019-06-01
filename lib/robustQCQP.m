function [RStar,tStar,EStar,nrItrs,weights,mu,options] = ...
    robustQCQP(bearing1,bearing2,options)

 % Robust version of QCQP
 % using Geman-McClure robust cost function
 % solve using graduated non-convexity, with QCQP relaxation as base solver
 %
 % Author: Heng Yang
 % Massachusetts Institute of Technology
 % Feb. 01, 2019

if ~isfield(options,'method')
	options.method='Geman-McClure';
end
if ~isfield(options,'QCQPDisp')
	options.QCQPDisp=false;
end
if ~isfield(options,'CVXDisp')
	options.CVXDisp=false;
end
if ~isfield(options,'robustDisp')
	options.robustDisp=false;
end
if ~isfield(options,'stopThresh')
	options.stopThresh=1e-6;
end
if ~isfield(options,'initialMu')
	options.initialMu=50;
end
if ~isfield(options,'maxSteps')
	options.maxSteps=100;
end
if ~isfield(options,'divFactor')
	options.divFactor=1.1;
end

if options.robustDisp
	fprintf('Robust QCQP, initial mu=%g, divFactor=%g, maxSteps=%g.\n',...
		options.initialMu,options.divFactor,options.maxSteps);
end

nrCorrs = size(bearing1,2);
% initial weights: take all measurements into optimization
weights=ones(1,nrCorrs);
residuals=zeros(1,nrCorrs);

cost=inf;
itr=1;
mu=options.initialMu;

QCQPOptions.CVXDisp=options.CVXDisp;
QCQPOptions.QCQPDisp=options.QCQPDisp;

while cost>options.stopThresh && itr<options.maxSteps
	% use weights to solve for E (R and t)
	weightedBearing1=weights.^(1/4) .* bearing1;
	weightedBearing2=weights.^(1/4) .* bearing2;
	[RHat,tHat,RStar,tStar,ZStar,recoverStatus,cvxStatus,pStarRelax,tightness,rankZ] =...
	qcqprelaxation(weightedBearing1,weightedBearing2,QCQPOptions);

	% check if returned four pairs of R and t result in same essential matrix
	errEThresh=1e-3;
	essential=hatmap(tHat{1}) * RHat{1};
	for i=2:length(RHat)
		tempE=hatmap(tHat{i}) * RHat{i}; % QCQP returned either E or -E
		errorE1=norm(tempE-essential, 'fro');
		errorE2=norm(tempE+essential, 'fro');
		if min(errorE1,errorE2)>errEThresh
			warning('QCQP returned R and t pairs are not consistent (lead to different essential)');
			warning(sprintf('Error between essential %g and essential 1 are %g and %g.\n',i,errorE1,errorE2));
			warning()
		end
	end

	% use RHat and tHat to solve for weights
	for i=1:nrCorrs
		f1=bearing1(:,i);
		f2=bearing2(:,i);
		residuals(i)=(f1' * essential * f2)^2;
	end
	weights=( mu ./ (mu + residuals) ).^2;

	% update cost
	cost=weights*residuals';

	% divide mu and move to next step
	mu = mu * options.divFactor;

	if options.robustDisp
		fprintf('Robust approach at step %g: cost=%g, mu=%g, rankZ=%g.\n',itr,cost,mu,rankZ);
	end

    % log data
    rankZStep(itr)=rankZ;
    tightnessStep(itr)=tightness;
    pStarRelaxStep(itr)=pStarRelax;
    costStep(itr)=cost;
    weightsStep(itr,:)=weights;
	% update iteration
	itr=itr+1;
end

nrItrs=itr;
EStar=essential;

if options.robustDisp
	figure;
	plot(costStep,'linewidth',2);
	set(gca, 'YScale', 'log');
	xlabel('Step')
	ylabel('Cost')
	title('Cost function w.r.t step')
	hold off
end

options.rankZStep=rankZStep;
options.tightnessStep=tightnessStep;
options.pStarRelaxStep=pStarRelaxStep;
options.costStep=costStep;
options.weightsStep=weightsStep;






