function Normality_test
%NORMALITY_TEST tests the normality of the response probability distribution

%   Márton Kajtor
%   Institute of Experimental Medicine, Budapest, Hungary
%   Lendulet Laboratory of Systems Neuroscience 
%   13-Jul-2023

% Loading data
load('F:\Marci\Cikkhez\raw data\2 weeks old raw\Red_mutant_groups_13tun.mat');
speedupFirst{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSpeedupFirst==1;
speedupFirst{mut}(isnan(speedupFirst{mut})) = 0;
slowingFirst{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSlowingFirst==1;
slowingFirst{mut}(isnan(slowingFirst{mut})) = 0;
stoppReactionTime{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedStoppReactionTime;
stoppReactionTime{mut}(stoppReactionTime{mut} < 0) = 0;
stoppingFirst{mut} = ((~isnan(stoppReactionTime{mut}) - slowingFirst{mut}) - speedupFirst{mut}) == 1;
stopPercDrosi{mut} = sum(stoppingFirst{mut},2)./sum(runBefore{mut},2);
noReaction{mut} = (((ones(size(stoppReactionTime{mut})) - stoppingFirst{mut}) - slowingFirst{mut}) - speedupFirst{mut}) == 1;
meanSpeedBefore200{mut} = Red_mutant_groups_13tun(mut).age_groups.Shadow_sorted_mean_speed_before_200*1000;
stopDuration{mut} = Red_mutant_groups_13tun(mut).age_groups.Shadow_sorted_stopdur/1000;
    
% Calculate the number of actions for each fly
stopPercDrosi{mut} = sum(stoppingFirst{mut},2)./sum(runBefore{mut},2);
slowPercDrosi{mut} = sum(slowingFirst{mut},2)./sum(runBefore{mut},2);
speedupPercDrosi{mut} = sum(speedupFirst{mut},2)./40;
noreactPercDrosi{mut} = sum(noReaction{mut},2)./40;

% Normality test (Kolmogorov-Smirnov test)
for normgrp = 1:3
    stopsnorm(normgrp) = kstest(stopPercDrosi{normgrp});
    slownorm(normgrp) = kstest(slowPercDrosi{normgrp});
    speedupnorm(normgrp) = kstest(speedupPercDrosi{normgrp});
    noreactnorm(normgrp) = kstest(noreactPercDrosi{normgrp});
    meanspeedbefnorm(normgrp) = kstest(meanSpeedBefore200{normgrp});
    stopdurnorm(normgrp) = kstest(stopDuration{normgrp});
end
end

