function Scatter_plots
% SCATTER_PLOTS schows the speed of the flies befor and after the looming
% stimulus. 

%   Márton Kajtor
%   Institute of Experimental Medicine, Budapest, Hungary
%   Lendulet Laboratory of Systems Neuroscience 
%   13-Jul-2023

load('C:\Users\User\Documents\muslica\muslica cikk\DUDALTER_mutant_groups_13tun.mat');
for mut = 1 : size(DUDALTER_mutant_groups_13tun,2)
    groupNameI{mut} = ['DUDA ' DUDALTER_mutant_groups_13tun(mut).names];
    groupName = {DUDALTER_mutant_groups_13tun(mut).names};
    plotId{mut} = groupName;

    for age = 1
        speedStPerdrosi2sec{mut,age} = DUDALTER_mutant_groups_13tun(mut).age_groups(age).speed_st_perdrosi_2sec;
        raw_pre_shadow{mut} = speedStPerdrosi2sec{mut,age}(:,1:201);
        raw_post_shadow{mut} = speedStPerdrosi2sec{mut,age}(:,201:401);
        
        row_mean_pre_shadow{mut} =  mean(raw_pre_shadow{mut},2);
        row_mean_post_shadow{mut} = mean(raw_post_shadow{mut},2);
        
%          slowing first index

        slow_first_matrix{mut} = DUDALTER_mutant_groups_13tun(mut).age_groups(age).ShadowSortedSlowingFirst;
        SlowFirstIndex{mut} = find(slow_first_matrix{mut}' == 1);
        raw_pre_shadow_Slow{mut} = speedStPerdrosi2sec{mut,age}(SlowFirstIndex{mut},150:201);
        raw_post_shadow_Slow{mut} = speedStPerdrosi2sec{mut,age}(SlowFirstIndex{mut},201:301); 
        
%         speedup first
        
        speedup_first_matrix{mut} = DUDALTER_mutant_groups_13tun(mut).age_groups(age).ShadowSortedSpeedupFirst;
        SpeedFirstIndex{mut} = find(speedup_first_matrix{mut}' == 1);
        raw_pre_shadow_Speedup{mut} = speedStPerdrosi2sec{mut,age}(SpeedFirstIndex{mut},180:201);
        raw_post_shadow_Speedup{mut} = speedStPerdrosi2sec{mut,age}(SpeedFirstIndex{mut},201:301); 
        
%       stopping first index
        
        stoppReactionTime{mut} = DUDALTER_mutant_groups_13tun(mut).age_groups(age).ShadowSortedStoppReactionTime;
        Stop_first_matrix{mut} = ((~isnan(stoppReactionTime{mut})-~isnan(slow_first_matrix{mut}))-~isnan(speedup_first_matrix{mut})) == 1;
        StopFirstIndex{mut} = find(Stop_first_matrix{mut}' == 1);
        raw_pre_shadow_Stops{mut} = speedStPerdrosi2sec{mut,age}(StopFirstIndex{mut},180:201);
        raw_post_shadow_Stops{mut} = speedStPerdrosi2sec{mut,age}(StopFirstIndex{mut},201:401);  
        
        noreact_first_matrix{mut} = DUDALTER_mutant_groups_13tun(mut).age_groups(age).ShadowSortedSlowingFirst;
        noreactFirstIndex{mut} = ((ones(size(stoppReactionTime{mut})) - (~isnan(stoppReactionTime{mut})-~isnan(slow_first_matrix{mut}))-~isnan(speedup_first_matrix{mut}))) == 1;
        NoReactIndex{mut} = find(noreactFirstIndex{mut}' == 1);
        raw_pre_shadow_noreact{mut} = speedStPerdrosi2sec{mut,age}(NoReactIndex{mut},180:201);
        raw_post_shadow_noreact{mut} = speedStPerdrosi2sec{mut,age}(NoReactIndex{mut},201:301);          

%       Calc mean values
        Duda_PreSh_Stop(mut) = mean(mean(raw_pre_shadow_Stops{mut},2));
        Duda_PostSh_Stop(mut) = mean(mean(raw_post_shadow_Stops{mut},2));
        Duda_PreSh_Slow(mut) = mean(mean(raw_pre_shadow_Slow{mut},2));
        Duda_PostSh_Slow(mut) = mean(mean(raw_post_shadow_Slow{mut},2));
        Duda_PreSh_Speedup(mut) = mean(mean(raw_pre_shadow_Speedup{mut},2));
        Duda_PostSh_Speedup(mut) = mean(mean(raw_post_shadow_Speedup{mut},2));
        Duda_PreSh_Noreact(mut) = mean(mean(raw_pre_shadow_noreact{mut},2));
        Duda_PostSh_Noreact(mut) = mean(mean(raw_post_shadow_noreact{mut},2));
    end 
end


stop_pre_trial_mean =  mean(raw_pre_shadow_Stops{5},2);
stop_post_trial_mean = min(raw_post_shadow_Stops{5}');
figure;
scatter(stop_pre_trial_mean,stop_post_trial_mean);
xlim([0 0.014]);
ylim([0 0.014]);
setmyplot_balazs
xlabel('mean speed before shadow in mm/ms');
ylabel('min speed aftere shadow in mm/ms');
title('stop');

for miert = 1 : size(raw_post_shadow_Slow{5},1)
    smooth_post_shadow_Slow{5}(miert,:) = smooth(raw_post_shadow_Slow{5}(miert,:),'linear',35);
    Duda_PostSh_Slow(miert) = min(smooth_post_shadow_Slow{5}(smooth_post_shadow_Slow{5}(miert,:) > 0)');
end

minzeroindex = find(min(smooth_post_shadow_Slow{5}') == 0);
Duda_PreSh_Slow = mean(raw_pre_shadow_Slow{5},2);
Duda_PostSh_Slow = min(smooth_post_shadow_Slow{5}');
for sm = 1 : size(minzeroindex,2)
    eqzero(sm) = find(smooth_post_shadow_Slow{5}(minzeroindex(sm),:) == 0,1,'first');
    if eqzero(sm) ~= 1
        Duda_PostSh_Slow(1,minzeroindex(sm)) = smooth_post_shadow_Slow{5}(minzeroindex(sm),eqzero(sm)-4);
    end
end

seceqzero = find(Duda_PostSh_Slow == 0); Duda_PostSh_Slow(seceqzero) = []; Duda_PreSh_Slow(seceqzero) = [];

figure;
scatter(Duda_PreSh_Slow,Duda_PostSh_Slow);
xlim([0 0.018]);
ylim([0 0.018]);
setmyplot_balazs
xlabel('mean speed before shadow in mm/ms');
ylabel('min speed aftere shadow in mm/ms');
title('slow');

Duda_PreSh_Speedup = mean(raw_pre_shadow_Speedup{5},2);
Duda_PostSh_Speedup = mean(raw_post_shadow_Speedup{5},2);
figure;
scatter(Duda_PreSh_Speedup,Duda_PostSh_Speedup);
xlim([0 0.014]);
ylim([0 0.014]);
setmyplot_balazs
xlabel('mean speed before shadow in mm/ms');
ylabel('mean speed aftere shadow in mm/ms');
title('speedup');

Duda_PreSh_Noreact = mean(raw_pre_shadow_noreact{5},2);
Duda_PostSh_Noreact = mean(raw_post_shadow_noreact{5},2);
figure;
scatter(Duda_PreSh_Noreact,Duda_PostSh_Noreact);
xlim([0 0.018]);
ylim([0 0.018]);
setmyplot_balazs
xlabel('mean speed before shadow in mm/ms');
ylabel('mean speed aftere shadow in mm/ms');
title('noreact');

end 