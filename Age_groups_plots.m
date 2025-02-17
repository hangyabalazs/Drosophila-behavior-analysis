function Age_groups_plots
%AGE_GRROUPS_PLOTS creates the plots for age and genotype group comparisons
%with ANOVA test.

load('F:\Marci\Cikkhez\raw data\age groups raw\Red_mutant_groups_13tun.mat');

GroupNameTag = [];
ageNameTag = [];
AnStopDrosiMeans = [];
AnSlowDrosiMeans = [];
AnSpeedupDrosiMeans = [];
AnNoreactDrosiMeans = [];
AnstopDurPerDrosi = [];
AnMeanSpeedPerDrosi = [];
for mut = 1 : size(Red_mutant_groups_13tun,2)
    groupName = {Red_mutant_groups_13tun(mut).names};
    for age = 1:5
        AgeName = {Red_mutant_groups_13tun(mut).age_groups(age).name};

        stoppReactionTime{mut,age} = Red_mutant_groups_13tun(mut).age_groups(age).ShadowSortedStoppReactionTime;
        stoppReactionTime{mut,age}(stoppReactionTime{mut,age} < 0) = 0;
        % stopDur
        stopDuration{mut,age} = Red_mutant_groups_13tun(mut).age_groups(age).Shadow_sorted_stopdur/1000;
       
        %Slows
        slowingFirst{mut,age} = Red_mutant_groups_13tun(mut).age_groups(age).ShadowSortedSlowingFirst==1;
        speedupFirst{mut,age} = Red_mutant_groups_13tun(mut).age_groups(age).ShadowSortedSpeedupFirst==1;
        meanSpeedBefore200{mut,age} = Red_mutant_groups_13tun(mut).age_groups(age).Shadow_sorted_mean_speed_before_200*1000;
        speedupFirst{mut,age}(isnan(speedupFirst{mut,age})) = 0;
        slowingFirst{mut,age}(isnan(slowingFirst{mut,age})) = 0;
        stoppingFirst{mut,age} = ((~isnan(stoppReactionTime{mut,age}) - slowingFirst{mut,age}) - speedupFirst{mut,age}) == 1;
        noReaction{mut,age} = (((ones(size(stoppReactionTime{mut,age})) - stoppingFirst{mut,age}) - slowingFirst{mut,age}) - speedupFirst{mut,age}) == 1;
        stopDuration{mut,age}(stoppingFirst{mut,age}==0) = NaN;
        
         animNum = size(stoppReactionTime{mut,age},1);
        for fas = 1 : animNum
            GroupNameTag = [GroupNameTag;groupName];
            ageNameTag = [ageNameTag; AgeName];   
        end
        
        stopPercDrosi{mut,age} = sum(stoppingFirst{mut,age},2)./sum(runBefore{mut,age},2);
        slowPercDrosi{mut,age} = sum(slowingFirst{mut,age},2)./sum(runBefore{mut,age},2);
        speedupPercDrosi{mut,age} = sum(speedupFirst{mut,age},2)./40;
        noreactPercDrosi{mut,age} = sum(noReaction{mut,age},2)./40;
        stopDurPerDrosi{mut,age} = nanmean(stopDuration{mut,age},2);
        MeanSpeedPerDrosi{mut,age} = nanmean(meanSpeedBefore200{mut,age},2);
        
        AnStopDrosiMeans = [AnStopDrosiMeans;stopPercDrosi{mut,age}];
        AnSlowDrosiMeans = [AnSlowDrosiMeans;slowPercDrosi{mut,age}];
        AnSpeedupDrosiMeans = [AnSpeedupDrosiMeans;speedupPercDrosi{mut,age}];
        AnNoreactDrosiMeans = [AnNoreactDrosiMeans;noreactPercDrosi{mut,age}];
        AnstopDurPerDrosi = [AnstopDurPerDrosi;stopDurPerDrosi{mut,age}];
        AnMeanSpeedPerDrosi = [AnMeanSpeedPerDrosi;MeanSpeedPerDrosi{mut,age}];
        
        stopsem(mut,age) = std(stopPercDrosi{mut,age})/sqrt(size(stopPercDrosi{mut,age},1));
        slowsem(mut,age) = std(slowPercDrosi{mut,age})/sqrt(size(slowPercDrosi{mut,age},1));
        speedupsem(mut,age) = std(speedupPercDrosi{mut,age})/sqrt(size(speedupPercDrosi{mut,age},1));
        noreactsem(mut,age) = std(noreactPercDrosi{mut,age})/sqrt(size(noreactPercDrosi{mut,age},1));
         
        meanstopPercDrosi(mut,age) = mean(stopPercDrosi{mut,age});       
        meanslowPercDrosi(mut,age) = mean(slowPercDrosi{mut,age});
        meanspeedupPercDrosi(mut,age) = mean(speedupPercDrosi{mut,age});
        meannoreactPercDrosi(mut,age) = mean(noreactPercDrosi{mut,age});
    
        meanMeanSpeedBefore200{mut,age} = mean(meanSpeedBefore200{mut,age},2);
        RealmeanMeanSpeedBefore200(mut,age) = mean(meanMeanSpeedBefore200{mut,age});
        meanspeed200sem(mut,age) = std(meanMeanSpeedBefore200{mut,age})/sqrt(length(meanMeanSpeedBefore200{mut,age}));
        
        stopDurMean{mut,age} = nanmean(stopDuration{mut,age}');
        RealmeanstopDurMean(mut,age) = nanmean(stopDurMean{mut,age},2);
        stopDurstd{mut,age} = nanstd(stopDurMean{mut,age});
        stopDursem(mut,age) = stopDurstd{mut,age}/sqrt(length(stopDurMean{mut,age}));
    end
end
 % ANOVA
[~,~,statsdur] = anovan(AnstopDurPerDrosi,{ageNameTag,GroupNameTag},'model','interaction','varnames',{'ageNameTag','GroupNameTag'});%savefig('stopage');
% figure;[resultsdur,~,~,gnames] = multcompare(statsdur,'Dimension',[2]);title('stopdur');%savefig('noreactage');
[~,~,statmeans] = anovan(AnMeanSpeedPerDrosi,{ageNameTag,GroupNameTag},'model','interaction','varnames',{'ageNameTag','GroupNameTag'});%savefig('stopage');
% figure;[resultmeans,~,~,gnames] = multcompare(statmeans,'Dimension',[2]);title('meanspeed200');%savefig('noreactage');

%% Mann Whitney betw age

MWageNameTag = deal(cell(1,5));
MWstopPercDrosi = deal(cell(1,5));
MWslowPercDrosi = deal(cell(1,5));
MWspeedupPercDrosi = deal(cell(1,5));
MWnoreactPercDrosi = deal(cell(1,5));
MWstopDurPerDrosi = deal(cell(1,5));
MWMeanSpeedPerDrosi = deal(cell(1,5));
for age = 1:5
    for mut = 1:5
        animNum = size(stoppReactionTime{mut,age},1);
        for KW = 1:animNum
            AgeName = {AgeNameI{age}};
            MWageNameTag{age} = [MWageNameTag{age}; AgeName];   
        end
        
        MWstopPercDrosi{age} = [MWstopPercDrosi{age};stopPercDrosi{mut,age}];
        MWslowPercDrosi{age} = [MWslowPercDrosi{age};slowPercDrosi{mut,age}];
        MWspeedupPercDrosi{age} = [MWspeedupPercDrosi{age};speedupPercDrosi{mut,age}];
        MWnoreactPercDrosi{age} = [MWnoreactPercDrosi{age};noreactPercDrosi{mut,age}];
    end
end


for agec = 1:4
   [manWPStopPec(agec),manWHStopPec(agec)] = ranksum(MWstopPercDrosi{agec},MWstopPercDrosi{agec+1},'alpha',siglev);
   [manWPSlowPec(agec),manWHSlowPec(agec)] = ranksum(MWslowPercDrosi{agec},MWslowPercDrosi{agec+1},'alpha',siglev);
   [manWPSpeedPec(agec),manWHSpeedPec(agec)] = ranksum(MWspeedupPercDrosi{agec},MWspeedupPercDrosi{agec+1},'alpha',siglev);
   [manWPnoreact(agec),manWHnoreact(agec)] = ranksum(MWnoreactPercDrosi{agec},MWnoreactPercDrosi{agec+1},'alpha',siglev);
end

% Plots
figure;
colorLine = ['g';'b';'m';'k';'r'];
for i = 1:5
  errorshade(1:5,meanstopPercDrosi(i,:)',stopsem(i,:)','LineColor',colorLine(i),'ShadeColor',colorLine(i));hold on;
end
setmyplot_balazs
legend('Dop1,2 EcR', 'Dop1R', 'Dop2R', 'w1118', 'y1w67c23');
title('Stop ratio');

figure;
for i = 1:5
    errorshade(1:5,meanslowPercDrosi(i,:)',slowsem(i,:)','LineColor',colorLine(i),'ShadeColor',colorLine(i));hold on;
end
setmyplot_balazs
legend('Dop1,2 EcR', 'Dop1R', 'Dop2R', 'w1118', 'y1w67c23');
title('Slow ratio');

figure;
for i = 1:5
    errorshade(1:5,meanspeedupPercDrosi(i,:)',speedupsem(i,:)','LineColor',colorLine(i),'ShadeColor',colorLine(i));hold on;
end
setmyplot_balazs
legend('Dop1,2 EcR', 'Dop1R', 'Dop2R', 'w1118', 'y1w67c23');
title('Speedup ratio');

figure;
for i = 1:5
    errorshade(1:5,meannoreactPercDrosi(i,:)',noreactsem(i,:)','LineColor',colorLine(i),'ShadeColor',colorLine(i));hold on;
end
setmyplot_balazs
legend('Dop1,2 EcR', 'Dop1R', 'Dop2R', 'w1118', 'y1w67c23');
title('Noreact ratio');


figure;
for i = 1:5
    errorshade(1:5,RealmeanMeanSpeedBefore200(i,:)',meanspeed200sem(i,:)','LineColor',colorLine(i),'ShadeColor',colorLine(i));hold on;
end
setmyplot_balazs
legend('Dop1,2 EcR', 'Dop1R', 'Dop2R', 'w1118', 'y1w67c23');
title('Mean speed bef');

figure;
for i = 1:5
    errorshade(1:5,RealmeanstopDurMean(i,:)',stopDursem(i,:)','LineColor',colorLine(i),'ShadeColor',colorLine(i));hold on;
end
setmyplot_balazs
legend('Dop1,2 EcR', 'Dop1R', 'Dop2R', 'w1118', 'y1w67c23');
title('Stop dur ratio');

close 
