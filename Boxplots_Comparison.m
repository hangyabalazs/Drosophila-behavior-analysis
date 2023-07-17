function Boxplots_Comparison(Red_mutant_groups_13tun)
%BOXPLOTS_COMPARISON   Statistics.
%   BOXPLOTS_COMPARISON creates boxplots to compare the behavioural
%   differences between fly genotypes.
%
%   Required input arguments:
%       RED_MUTANT_GROUPS_13TUN: preprocessed results struct calculated by
%       the DROSOPHILA_SPEED funtion.
%
%   See also DROSOPHILA_SPEED.

% Calculate reaction times, stop durations and reaction probabilities for
% each group
for mut = 1 : size(Red_mutant_groups_13tun,2)
        groupName = [ Red_mutant_groups_13tun(mut).names]; 
        plotId{mut} = groupName;
        stoppReactionTime{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedStoppReactionTime;
        stoppReactionTime{mut}(stoppReactionTime{mut} < 0) = 0;
        stoppingTime{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedStoppingTime;
        stoppingTime{mut}(stoppingTime{mut} < 0) = 0;
        stopDuration{mut} = Red_mutant_groups_13tun(mut).age_groups.Shadow_sorted_stopdur/1000;
        slowingReactionTime{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSlowingReactionTime;
        slowingReactionTime{mut}(slowingReactionTime{mut} < 0) = 0;
        slowingTime{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSlowingTime;
        slowingTime{mut}(slowingTime{mut} < 0) = 0;
        slowingFirst{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSlowingFirst==1;
        speedingReactionTime{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSpeedingReactionTime;
        speedingReactionTime{mut}(speedingReactionTime{mut} < 0) = 0;
        speedingTime{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSpeedingTime;
        speedingTime{mut}(speedingTime{mut} < 0) = 0;
        speedupFirst{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSpeedupFirst==1;
        meanSpeedBefore200{mut} = Red_mutant_groups_13tun(mut).age_groups.Shadow_sorted_mean_speed_before_200*1000;
        stantBefore{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedStandBefore;
        runBefore{mut} = ones(size(stantBefore{mut})) - stantBefore{mut};
        speedupFirst{mut}(isnan(speedupFirst{mut})) = 0;
        slowingFirst{mut}(isnan(slowingFirst{mut})) = 0;
        stoppingFirst{mut} = ((~isnan(stoppReactionTime{mut}) - slowingFirst{mut}) - speedupFirst{mut}) == 1;
        FirststoppReactionTime{mut} = nan(size(stoppReactionTime{mut}));
        FirststoppReactionTime{mut}(stoppingFirst{mut}) = stoppReactionTime{mut}(stoppingFirst{mut});
        noReaction{mut} = (((ones(size(stoppReactionTime{mut})) - stoppingFirst{mut}) - slowingFirst{mut}) - speedupFirst{mut}) == 1;
        allReactionTime{mut} = nan(size(stoppReactionTime{mut}));
        allReactionTime{mut}(stoppingFirst{mut}) = stoppReactionTime{mut}(stoppingFirst{mut});
        allReactionTime{mut}(slowingFirst{mut}) = slowingReactionTime{mut}(slowingFirst{mut});
        allReactionTime{mut}(speedupFirst{mut}) = speedingReactionTime{mut}(speedupFirst{mut});
        allActioTime{mut} = nan(size(stoppReactionTime{mut}));
        allActioTime{mut}(slowingFirst{mut}) = slowingTime{mut}(slowingFirst{mut});
        allActioTime{mut}(speedupFirst{mut}) = speedingTime{mut}(speedupFirst{mut});
        
        stopDuration{mut}(stoppingFirst{mut}==0) = NaN;
        stoppingTime{mut}(stoppingFirst{mut}==0) = NaN;
        stopPercDrosi{mut} = sum(stoppingFirst{mut},2)./sum(runBefore{mut},2);
        slowPercDrosi{mut} = sum(slowingFirst{mut},2)./sum(runBefore{mut},2);
        speedupPercDrosi{mut} = sum(speedupFirst{mut},2)./40;
        noreactPercDrosi{mut} = sum(noReaction{mut},2)./40;
        
        meanstopPercDrosi{mut} = mean(stopPercDrosi{mut});
        meanslowPercDrosi{mut} = mean(slowPercDrosi{mut});
        meanspeedupPercDrosi{mut} = mean(speedupPercDrosi{mut});
        meannoreactPercDrosi{mut} = mean(noreactPercDrosi{mut});
end

for round = 1 : 7
    cont = [3 4 3 4 9 9 8];
    muta = [1 1 2 2 5 6 7];
    
    % Perform Wilcoxon ranksum test for statistical comparisons
    siglev = 0.05;
    [manWPAllReact(round),manWHAllReact(round)] = ranksum(nanmean(allReactionTime{cont(round)},2),nanmean(allReactionTime{muta(round)},2),'alpha',siglev);
    [manWPAllAct(round),manWHAllAct(round)] = ranksum(nanmean(allActioTime{cont(round)},2),nanmean(allActioTime{muta(round)},2),'alpha',siglev);
    [manWPStopPec(round),manWHStopPec(round)] = ranksum(stopPercDrosi{cont(round)},stopPercDrosi{muta(round)},'alpha',siglev);
    [manWPSlowPerc(round),manWHSlowPerc(round)] = ranksum(slowPercDrosi{cont(round)},slowPercDrosi{muta(round)},'alpha',siglev);
    [manWPSpeedupPerc(round),manWHSpeedupPerc(round)] = ranksum(speedupPercDrosi{cont(round)},speedupPercDrosi{muta(round)},'alpha',siglev);
    [manWPNoReact(round),manWHNoReact(round)] = ranksum(noreactPercDrosi{cont(round)},noreactPercDrosi{muta(round)},'alpha',siglev);
    [manWPstopReact(round),manWHstopReact(round)] = ranksum(nanmean(FirststoppReactionTime{cont(round)},2),nanmean(FirststoppReactionTime{muta(round)},2),'alpha',siglev);
    [manWPstopAct(round),manWHstopAct(round)] = ranksum(nanmean(stoppingTime{cont(round)},2),nanmean(stoppingTime{muta(round)},2),'alpha',siglev);
    [manWPSlowReact(round),manWHSlowReact(round)] = ranksum(nanmean(slowingReactionTime{cont(round)},2),nanmean(slowingReactionTime{muta(round)},2),'alpha',siglev);
    [manWPSlowAct(round),manWHSlowAct(round)] = ranksum(nanmean(slowingTime{cont(round)},2),nanmean(slowingTime{muta(round)},2),'alpha',siglev);
    [manWPSpeedupReact(round),manWHSpeedupReact(round)] = ranksum(nanmean(speedingReactionTime{cont(round)},2),nanmean(speedingReactionTime{muta(round)},2),'alpha',siglev);
    [manWPSpeedupAct(round),manWHSpeedupAct(round)] = ranksum(nanmean(speedingTime{cont(round)},2),nanmean(speedingTime{muta(round)},2),'alpha',siglev);
    [manWPStopDur(round),manWHStopDur(round)] = ranksum(nanmean(stopDuration{cont(round)},2),nanmean(stopDuration{muta(round)},2),'alpha',siglev);
    [manWPMeanSpeedBef(round),manWHMeanSpeedBef(round)] = ranksum(nanmean(meanSpeedBefore200{cont(round)},2),nanmean(meanSpeedBefore200{muta(round)},2),'alpha',siglev);
      
end

% Plotting results
order = [1 2 3 4 5 6 9 7 8];
color = [1 0 0;0 0 1;0 0 1;1 1 0;1 0 1;0 1 1;0.5 1 1;0 0.5 1; 1 0 0.5];

tagVector =vertcat(repmat(1,size(allActioTime{order(1)},1),1),repmat(2,size(allActioTime{order(2)},1),1),repmat(3,size(allActioTime{order(3)},1),1),...
repmat(4,size(allActioTime{order(4)},1),1),repmat(5,size(allActioTime{order(5)},1),1),repmat(6,size(allActioTime{order(6)},1),1),...
repmat(7,size(allActioTime{order(7)},1),1),repmat(8,size(allActioTime{order(8)},1),1),repmat(9,size(allActioTime{order(9)},1),1));
Lab = [{plotId{order(1)}} {plotId{order(2)}} {plotId{order(3)}} {plotId{order(4)}}];
Col = [color(order(1),:); color(order(2),:); color(order(3),:); color(order(4),:);];

% Freezing
figure;
boxDataStopPerc = vertcat(stopPercDrosi{order(1)},stopPercDrosi{order(2)},stopPercDrosi{order(3)},stopPercDrosi{order(4)},stopPercDrosi{order(5)},...
    stopPercDrosi{order(6)},stopPercDrosi{order(7)},stopPercDrosi{order(8)},stopPercDrosi{order(9)});
boxplot(boxDataStopPerc,tagVector,'labels',Lab,'Color',Col);
setmyplot_balazs
title('Stop ratio');
ylabel('Ratio of reaction');
h=findobj(gca,'tag','Outliers');
delete(h);
ylim auto

% Stop duration
figure;
boxDatastopDuration = vertcat(nanmean(stopDuration{order(1)},2),nanmean(stopDuration{order(2)},2),nanmean(stopDuration{order(3)},2),nanmean(stopDuration{order(4)},2),...
    nanmean(stopDuration{order(5)},2),nanmean(stopDuration{order(6)},2),nanmean(stopDuration{order(7)},2),nanmean(stopDuration{order(8)},2),nanmean(stopDuration{order(9)},2));
boxplot(boxDatastopDuration,tagVector,'labels',Lab,'Color',Col);
setmyplot_balazs
title('stopDuration');
ylabel('Time in s');
h=findobj(gca,'tag','Outliers');
delete(h);
ylim auto

% Mean speed before
figure;
boxDatameanSpeedBefore200 = vertcat(nanmean(meanSpeedBefore200{order(1)},2),nanmean(meanSpeedBefore200{order(2)},2),nanmean(meanSpeedBefore200{order(3)},2),nanmean(meanSpeedBefore200{order(4)},2),...
    nanmean(meanSpeedBefore200{order(5)},2),nanmean(meanSpeedBefore200{order(6)},2),nanmean(meanSpeedBefore200{order(7)},2),nanmean(meanSpeedBefore200{order(8)},2),nanmean(meanSpeedBefore200{order(9)},2));
boxplot(boxDatameanSpeedBefore200,tagVector,'labels',Lab,'Color',Col);
setmyplot_balazs
title('meanSpeedBefore200');
ylabel('Speed in mm/s');
h=findobj(gca,'tag','Outliers');
delete(h);
ylim auto

% Slow down
figure;
boxDataSlowPerc = vertcat(slowPercDrosi{order(1)},slowPercDrosi{order(2)},slowPercDrosi{order(3)},slowPercDrosi{order(4)},slowPercDrosi{order(5)},...
    slowPercDrosi{order(6)},slowPercDrosi{order(7)},slowPercDrosi{order(8)},slowPercDrosi{order(9)});
boxplot(boxDataSlowPerc,tagVector,'labels',Lab,'Color',Col);
setmyplot_balazs
title('Slow ratio');
ylabel('Ratio of reaction');
h=findobj(gca,'tag','Outliers');
delete(h);
ylim auto

% Speed up
figure;
boxDataSpeedupPerc = vertcat(speedupPercDrosi{order(1)},speedupPercDrosi{order(2)},speedupPercDrosi{order(3)},speedupPercDrosi{order(4)},speedupPercDrosi{order(5)},...
    speedupPercDrosi{order(6)},speedupPercDrosi{order(7)},speedupPercDrosi{order(8)},speedupPercDrosi{order(9)});
boxplot(boxDataSpeedupPerc,tagVector,'labels',Lab,'Color',Col);
setmyplot_balazs
title('Speedup ratio');
ylabel('Ratio of reaction');
h=findobj(gca,'tag','Outliers');
delete(h);
ylim auto

% No reaction
figure;
boxDataNoReactPerc = vertcat(noreactPercDrosi{order(1)},noreactPercDrosi{order(2)},noreactPercDrosi{order(3)},noreactPercDrosi{order(4)},noreactPercDrosi{order(5)},...
    noreactPercDrosi{order(6)},noreactPercDrosi{order(7)},noreactPercDrosi{order(8)},noreactPercDrosi{order(9)});
boxplot(boxDataNoReactPerc,tagVector,'labels',Lab,'Color',Col);
setmyplot_balazs
title('No Reactions ratio');
ylabel('Ratio of reaction');
h=findobj(gca,'tag','Outliers');
delete(h);
ylim auto