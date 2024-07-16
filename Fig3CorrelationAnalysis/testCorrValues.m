tic
%% read in data of KC activities
data = readtable('all KC values.csv');
allKCvalues = data(3:end,2:end); % remove header rows, and first column
allKCvalues = table2array(allKCvalues);
linearized = allKCvalues(:); % linearize the array
linearized(isnan(linearized)) = []; % remove NaNs
numOdors = 7;
numSamples = length(linearized)/numOdors;
if floor(numSamples)~=numSamples
    error('number of measurements not evenly divisible by numOdors');
end
rectData = reshape(linearized,numOdors,numSamples); % make a 7 x numSamples array - each column is one brain

%% read in true correlation values
realCorrValues = table2array(readtable('allKCAPLcorrValues.csv'));

% % Gradient descent to find the correct noise value
% mu_std = -2; %-1;
% sigma_std = 0; %0.2742; %0.5;
% 
% numTrials = 1;
% converged = false;
% targetMedian = 0.7;
% targetStdCorr = 0.3245;
% epsilonMedian = 0.01;
% epsilonStdCorr = 0.01;
% alpha = 0.001;
% medians = zeros(10000,1);
% numIters = 0;
% while ~converged
%     numIters = numIters+1;
%     corrValues = zeros(numSamples,numTrials);
%     for i=1:numTrials
%         for j=1:numSamples
%             x = rectData(:,j);
%             stdev = exp(mu_std+randn(1,1)*sigma_std); % sample a random number from a log-normal distribution
%             y = x+randn(numOdors,1).*stdev;
%             corrValues(j,i) = corr(x,y);
%         end
%     end
%     medians(numIters) = median(corrValues,"all");
%     medianDiff = median(corrValues,"all") - targetMedian;
%     stdDiff = std(corrValues,0,"all") - targetStdCorr;
%     medianConverged = abs(medianDiff) < epsilonMedian;
%     stdConverged = abs(stdDiff) < epsilonStdCorr;
%     if medianConverged && stdConverged
%         converged = true;
%     else
%         if ~medianConverged
%             % if the correlation is too low (medianDiff<0), need to
%             % decrease the error (decrease mu_std)
%             mu_std = mu_std + medianDiff*alpha;
%         end
%         if ~stdConverged
%             sigma_std = sigma_std - stdDiff*alpha;
%         end
%     end
%     [mu_std, sigma_std, medianDiff, stdDiff]
% 
% end

%% run the simulation
numTrials = 1000;
stdev = 0.1353;
corrValues = zeros(numSamples,numTrials);
for i=1:numTrials
    for j=1:numSamples
        x = rectData(:,j);
        y = x + randn(numOdors,1)*stdev;
        corrValues(j,i) = corr(x,y);
    end
end

%% plot ecdfs of true correlation values vs results of simulation
figure
hold on
bins = -0.995:0.01:0.995;
realECDF = interpECDF(realCorrValues,bins);
stairs(bins-0.005,realECDF);
simHistValues = zeros(numTrials,length(bins));
for t=1:numTrials
    try
        [simHistValues(t,:),~] = interpECDF(corrValues(:,t),bins);
    catch
        error('here we are')
    end
end
errorRange = prctile(simHistValues,[2.5 97.5],1);
meanSimHistValues = mean(simHistValues,1,"omitnan");
stairs(bins-0.005,meanSimHistValues,'k');
binEdges = [-1 repelem((-.99:0.01:0.99),2) 1];
fillx = [binEdges, fliplr(binEdges)];
filly = repelem([errorRange(1,:) fliplr(errorRange(2,:))],2);

fillx(isnan(fillx)) = 0;
filly(isnan(filly)) = 0;
fill(fillx, filly,'k','FaceAlpha',0.2,'Edgecolor','none')
ylim([0 1])

xlabel('Correlation','FontSize',6)
ylabel('Cumulative distribution','FontSize',6)
axis square
set(gca,'FontSize',6)
set(gca,'TickLength',[.03 .03],'TickDir','out')
set(gca,'units','centimeters','position',[1,1,2.4,2.4])
xticks(-1:0.5:1);
yticks(0:0.5:1);
xtickangle(0)

%% run K-S test
[pvalue, withinSample, KSdiv, meanWithin] = MonteCarloKS2test(realCorrValues, corrValues);
pvalue
toc