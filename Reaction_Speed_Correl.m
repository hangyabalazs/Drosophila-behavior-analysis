function Reaction_Speed_Correl(Red_mutant_groups_13tun)
%REACTION_SPEED_CORREL   Speed-modulation of responses.
%   REACTION_SPEED_CORREL(DATA) calculates the reaction type proportions
%   and probabilities at a each speed and performs ANOVA tests for
%   comparison between groups.
%
%   Required input arguments:
%       RED_MUTANT_GROUPS_13TUN: preprocessed results struct calculated by
%       the DROSOPHILA_SPEED funtion.
%
%   See also DROSOPHILA_SPEED.

% Inicialize variables
ar = 0.002 * 1000;
br = 0.02 * 1000;
stophist = [];
stbefstopsmut = [];
stbefstopsval = [];
slbefslowsmut = [];
slbefslowsval = [];
subefspeedupsmut = [];
subefspeedupsval = [];
nrbefnoreactssmut = [];
nrbefnoreactsval = [];

% Extract speed-reaction histgorams for each group
for mut = 1:size(Red_mutant_groups_13tun,2)
        stoppReactionTime{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedStoppReactionTime;
        stoppReactionTime{mut}(stoppReactionTime{mut} < 0) = 0;
        slowingFirst{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSlowingFirst==1;
        speedupFirst{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSpeedupFirst==1;
        meanSpeedBefore200{mut} = Red_mutant_groups_13tun(mut).age_groups.Shadow_sorted_mean_speed_before_200*1000;
        speedupFirst{mut}(isnan(speedupFirst{mut})) = 0;
        slowingFirst{mut}(isnan(slowingFirst{mut})) = 0;
        stoppingFirst{mut} = ((~isnan(stoppReactionTime{mut}) - slowingFirst{mut}) - speedupFirst{mut}) == 1;
        noReaction{mut} = (((ones(size(stoppReactionTime{mut})) - stoppingFirst{mut}) - slowingFirst{mut}) - speedupFirst{mut}) == 1; 
        onesHist = ones(size(stoppingFirst{mut}));
        [speedCounts{mut},speedCenter{mut}] = hist(meanSpeedBefore200{mut}(onesHist == 1),0:ar:br);
        [stopCounts{mut},stopCenter{mut}] = hist(meanSpeedBefore200{mut}(stoppingFirst{mut}),0:ar:br);
        [speedupCounts{mut},speedupCenter{mut}] = hist(meanSpeedBefore200{mut}(speedupFirst{mut}),0:ar:br);
        [slowCounts{mut},slowCenter{mut}] = hist(meanSpeedBefore200{mut}(slowingFirst{mut}),0:ar:br);
        [noReactCounts{mut},noReactCenter{mut}] = hist(meanSpeedBefore200{mut}(noReaction{mut}),0:ar:br);

        stopCounts{mut}(speedCounts{mut}(:) < 5) = [];
        speedupCounts{mut}(speedCounts{mut}(:) < 5) = [];
        slowCounts{mut}(speedCounts{mut}(:) < 5) = [];
        noReactCounts{mut}(speedCounts{mut}(:) < 5) = [];
        speedCenter{mut}(speedCounts{mut}(:) < 5) = [];
        stopCenter{mut}(speedCounts{mut}(:) < 5) = [];
        speedupCenter{mut}(speedCounts{mut}(:) < 5) = [];
        slowCenter{mut}(speedCounts{mut}(:) < 5) = [];
        noReactCenter{mut}(speedCounts{mut}(:) < 5) = []; 
        speedCounts{mut}(speedCounts{mut}(:) < 5) = [];
        
        medmeanSpeedBefStops(mut) = median(meanSpeedBefore200{mut}(stoppingFirst{mut}));
        meanSpeedsBefStops = meanSpeedBefore200{mut}(stoppingFirst{mut});
        stbefstopsval = vertcat(stbefstopsval,meanSpeedsBefStops);
        stmutsvec = repmat(mut, length(meanSpeedsBefStops),1);
        stbefstopsmut = vertcat(stbefstopsmut,stmutsvec);
        
        medmeanSpeedBefSlows(mut) = median(meanSpeedBefore200{mut}(slowingFirst{mut}));
        meanSpeedsBefslow = meanSpeedBefore200{mut}(slowingFirst{mut});
        slbefslowsval = vertcat(slbefslowsval,meanSpeedsBefslow);
        slmutsvec = repmat(mut, length(meanSpeedsBefslow),1);
        slbefslowsmut = vertcat(slbefslowsmut,slmutsvec);
        
        medmeanSpeedBefSpeedups(mut) = median(meanSpeedBefore200{mut}(speedupFirst{mut}));
        meanSpeedsBefsppedups = meanSpeedBefore200{mut}(speedupFirst{mut});
        subefspeedupsval = vertcat(subefspeedupsval,meanSpeedsBefsppedups);
        sumutsvec = repmat(mut, length(meanSpeedsBefsppedups),1);
        subefspeedupsmut = vertcat(subefspeedupsmut,sumutsvec);
        
        medmeanSpeedBefNoreacts(mut) = median(meanSpeedBefore200{mut}(noReaction{mut}));
        meanSpeedsBefnoreact = meanSpeedBefore200{mut}(noReaction{mut});
        nrbefnoreactsval = vertcat(nrbefnoreactsval,meanSpeedsBefnoreact);
        nrmutsvec = repmat(mut, length(meanSpeedsBefnoreact),1);
        nrbefnoreactssmut = vertcat(nrbefnoreactssmut,nrmutsvec);
        
        [SpedCountsedBefStops{mut}, SpedCentersedBefStops{mut}] = hist(meanSpeedBefore200{mut}(stoppingFirst{mut}),0:ar:br);
        [SpedCountsedBefslowing{mut}, SpedCentersedBefslowing{mut}] = hist(meanSpeedBefore200{mut}(slowingFirst{mut}),0:ar:br);
        [SpedCountsedBefspeedup{mut}, SpedCentersedBefspeedup{mut}] = hist(meanSpeedBefore200{mut}(speedupFirst{mut}),0:ar:br);
        [SpedCountsedBefnoReact{mut}, SpedCentersedBefnoReact{mut}] = hist(meanSpeedBefore200{mut}(noReaction{mut}),0:ar:br);
        
        SpedCountsedBefStops{mut}(SpedCountsedBefStops{mut}(:) < 5) = 0;
        SpedCountsedBefslowing{mut}(SpedCountsedBefslowing{mut}(:) < 5) = 0;
        SpedCountsedBefspeedup{mut}(SpedCountsedBefspeedup{mut}(:) < 5) = 0;
        SpedCountsedBefnoReact{mut}(SpedCountsedBefnoReact{mut}(:) < 5) = 0;
        
        NorSpedCountsedBefStops{mut} = SpedCountsedBefStops{mut}./sum(SpedCountsedBefStops{mut});
        NorSpedCountsedBefslowing{mut} = SpedCountsedBefslowing{mut}./sum(SpedCountsedBefslowing{mut});
        NorSpedCountsedBefspeedup{mut} = SpedCountsedBefspeedup{mut}./sum(SpedCountsedBefspeedup{mut});
        NorSpedCountsedBefnoReact{mut} = SpedCountsedBefnoReact{mut}./sum(SpedCountsedBefnoReact{mut});
end
col = ['r';'k';'k';'k';'g';'y';'m';'b';]; 

close all
stopcents = [];
slowcents = [];
speedupcents = [];
noreactscents = [];
for counts = 17 : -2 : 3
    stopcents(stbefstopsval < counts) = counts/2;
end
for counts = 19 : -2 : 1
    slowcents(slbefslowsval < counts) = counts/2;
end
for counts = 15 : -2 : 1
    speedupcents(subefspeedupsval < counts) = counts/2;
end
for counts = 9 : -2 : 1
    noreactscents(nrbefnoreactsval < counts) = counts/2;
end
stopcents(stopcents == 0) = nan;
slowcents(slowcents == 0) = nan;
speedupcents(speedupcents == 0) = nan;
noreactscents(noreactscents == 0) = nan;

% ANOVA tests
[~,~,statsst] = anovan(stbefstopsval' ,{stbefstopsmut,stopcents'},'model','interaction','varnames',{'stbefstopsmut','stopcents'});
figure;[resultsst,~,~,gnamesst] = multcompare(statsst,'Dimension',[1]);title('stops');

[~,~,statssl] = anovan(slbefslowsval' ,{slbefslowsmut,slowcents'},'model','interaction','varnames',{'slbefslowsmut','slowcents'});
figure;[resultssl,~,~,gnamessl] = multcompare(statssl,'Dimension',[1]);title('slows');

[~,~,statssp] = anovan(subefspeedupsval' ,{subefspeedupsmut,speedupcents'},'model','interaction','varnames',{'subefspeedupsmut','speedupcents'});
figure;[resultssu,~,~,gnamessu] = multcompare(statssp,'Dimension',[1]);title('speedup');

[~,~,statsnr] = anovan(nrbefnoreactsval' ,{nrbefnoreactssmut,noreactscents'},'model','interaction','varnames',{'nrbefnoreactssmut','noreactscents'});
figure;[resultsnr,~,~,gnamesnr] = multcompare(statsnr,'Dimension',[1]);title('noreactions');

% Ratio plots
figure;  % stops
for sbs = 1 : length(NorSpedCountsedBefStops)
    plot(SpedCentersedBefStops{sbs},NorSpedCountsedBefStops{sbs},[col(sbs,1) '-']);hold on;
end
title('Mean speed counts before stops');
xlabel('speed in mm/s');
ylabel('counts');
legend(Red_mutant_groups_13tun.names);
setmyplot_balazs

figure; % slows
for sbs = 1 : length(NorSpedCountsedBefslowing)
    plot(SpedCentersedBefslowing{sbs},NorSpedCountsedBefslowing{sbs},[col(sbs,1) '-']);hold on;
end
title('Mean speed counts before slows');
xlabel('speed in mm/s');
ylabel('ratio');
legend(Red_mutant_groups_13tun.names);
setmyplot_balazs

figure; % speedups
for sbs = 1 : length(NorSpedCountsedBefspeedup)
    plot(SpedCentersedBefspeedup{sbs},NorSpedCountsedBefspeedup{sbs},[col(sbs,1) '-']);hold on;
end
title('Mean speed counts before speedups');
xlabel('speed in mm/s');
ylabel('ratio');
legend(Red_mutant_groups_13tun.names);
setmyplot_balazs

figure; % noreacts
for sbs = 1 : length(NorSpedCountsedBefnoReact)
    plot(SpedCentersedBefnoReact{sbs},NorSpedCountsedBefnoReact{sbs},[col(sbs,1) '-']);hold on;
end
title('Mean speed counts before no reactions');
xlabel('speed in mm/s');
ylabel('ratio');
legend(Red_mutant_groups_13tun.names);
setmyplot_balazs


% Probability plots
figure;
for sbs = 1 : length(stopCounts)
    plot(stopCenter{sbs},(stopCounts{sbs}./speedCounts{sbs}),[col(sbs,1) '-']);hold on;
end
title('stops % distribution');
xlabel('speed in mm/s');
ylabel('probability');
legend(Red_mutant_groups_13tun.names);
setmyplot_balazs

figure;
for sbs = 1 : length(slowCounts)
    plot(slowCenter{sbs},(slowCounts{sbs}./speedCounts{sbs}),[col(sbs,1) '-']);hold on;
end
title('slows % distribution');
xlabel('speed in mm/s');
ylabel('probability');
legend(Red_mutant_groups_13tun.names);
setmyplot_balazs

figure;
for sbs = 1 : length(speedupCounts)
    plot(speedupCenter{sbs},(speedupCounts{sbs}./speedCounts{sbs}),[col(sbs,1) '-']);hold on;
end
title('speedup % distribution');
xlabel('speed in mm/s');
ylabel('probability');
legend(Red_mutant_groups_13tun.names);
setmyplot_balazs

figure;
for sbs = 1 : length(noReactCounts)
    plot(noReactCenter{sbs},(noReactCounts{sbs}./speedCounts{sbs}),[col(sbs,1) '-']);hold on;
end
title('noreact % distribution');
xlabel('speed in mm/s');
ylabel('probability');
legend(Red_mutant_groups_13tun.names);
setmyplot_balazs