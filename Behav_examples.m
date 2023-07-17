function Behav_examples(Red_mutant_groups_13tun)
%BEHAV_EXAMLPES   Example plots of responses.
%   BEHAV_EXAMLPES(DATA) creates the plots for figure 2C The function is
%   using the raw data before and after the shadow presentation to showcase
%   the 4 main behaviour types in single trials, and the average speed
%   characteristic of the sorted behaviour.

%   Required input arguments:
%       RED_MUTANT_GROUPS_13TUN: preprocessed results struct calculated
%       with the drosophila_speed funtion.

% See also DROSOPHILA_SPEED.

% Extract data for each group
for mut = 1:size(Red_mutant_groups_13tun,2)
    
        stoppReactionTime{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedStoppReactionTime;
        stoppReactionTime{mut}(stoppReactionTime{mut} < 0) = 0;
        slowingFirst{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSlowingFirst==1;
        speedupFirst{mut} = Red_mutant_groups_13tun(mut).age_groups.ShadowSortedSpeedupFirst==1;
        speedStPerdrosi2sec{mut} = Red_mutant_groups_13tun(mut).age_groups.speed_st_perdrosi_2sec*1000;
        meanSpeedBefore200{mut} = Red_mutant_groups_13tun(mut).age_groups.Shadow_sorted_mean_speed_before_200;
        stoppAfter{mut} = Red_mutant_groups_13tun(mut).age_groups.collected_freezing;
        speedupFirst{mut}(isnan(speedupFirst{mut})) = 0;
        slowingFirst{mut}(isnan(slowingFirst{mut})) = 0;
        stoppingFirst{mut} = ((~isnan(stoppReactionTime{mut}) - slowingFirst{mut}) - speedupFirst{mut}) == 1;
        noReaction{mut} = (((ones(size(stoppReactionTime{mut})) - stoppingFirst{mut}) - slowingFirst{mut}) - speedupFirst{mut}) == 1;
        
        % Creating the example plots for freezing
        kstopp = speedStPerdrosi2sec{mut}(find(stoppingFirst{mut}' == 1),:);
        figure; plot(mean(kstopp));
        for i = 1 : 1 %size(kstopp,1)
            figure; plot(-2:0.01:2,smooth(kstopp(i,:),'linear',7));
            title('stopping first');
            rectangle('Position',[-0.20,0,0.20,20],'FaceColor',[0.9290 0.6940 0.1250 0.33],'EdgeColor','none','LineWidth',3);
            rectangle('Position',[-0.05,0,0.05,20],'FaceColor',[0.8500 0.3250 0.0980 0.50],'EdgeColor','none','LineWidth',3);
            rectangle('Position',[0,0,1,20],'FaceColor',[0.4660 0.6740 0.1880 0.33],'EdgeColor','none','LineWidth',3);
            xlim([-1 1.5]);
            ylabel('Speed (mm/s)');
            xlabel('Time (s)');
            setmyplot_balazs
        end
        
        % Creating the mean plots for freezing
        figure;
        errorshade(-2:0.01:2,mean(kstopp),std(kstopp)/sqrt(size(kstopp,1)));
        title('stopping first');
        rectangle('Position',[-0.20,0,0.20,20],'FaceColor',[0.9290 0.6940 0.1250 0.33],'EdgeColor','none','LineWidth',3);
        rectangle('Position',[-0.05,0,0.05,20],'FaceColor',[0.8500 0.3250 0.0980 0.50],'EdgeColor','none','LineWidth',3);
        rectangle('Position',[0,0,1,20],'FaceColor',[0.4660 0.6740 0.1880 0.33],'EdgeColor','none','LineWidth',3);
        xlim([-1 1.5]);
        ylabel('Speed (mm/s)');
        xlabel('Time (s)');
        ylim([0 20]);
        setmyplot_balazs
 
        % Creating the example plots for slow down
        kslow = speedStPerdrosi2sec{mut}(find(slowingFirst{mut}' == 1),:);
        for i = 1 : 1 %size(kslow,1)
            figure; plot(-2:0.01:2,smooth(kslow(i,:),'linear',7));
            title('slowing first');
            rectangle('Position',[-0.20,0,0.20,20],'FaceColor',[0.9290 0.6940 0.1250 0.33],'EdgeColor','none','LineWidth',3);
            rectangle('Position',[-0.05,0,0.05,20],'FaceColor',[0.8500 0.3250 0.0980 0.50],'EdgeColor','none','LineWidth',3);
            rectangle('Position',[0,0,1,20],'FaceColor',[0.4660 0.6740 0.1880 0.33],'EdgeColor','none','LineWidth',3);
            xlim([-1 1.5]);
            ylabel('Speed (mm/s)');
            xlabel('Time (s)');
            setmyplot_balazs
        end
        
        % Creating the mean plots for slow down
        figure;
        errorshade(-2:0.01:2,mean(kslow),std(k)/sqrt(size(kslow,1)));
        title('slowing first');
        rectangle('Position',[-0.20,0,0.20,20],'FaceColor',[0.9290 0.6940 0.1250 0.33],'EdgeColor','none','LineWidth',3);
        rectangle('Position',[-0.05,0,0.05,20],'FaceColor',[0.8500 0.3250 0.0980 0.50],'EdgeColor','none','LineWidth',3);
        rectangle('Position',[0,0,1,20],'FaceColor',[0.4660 0.6740 0.1880 0.33],'EdgeColor','none','LineWidth',3);
        xlim([-1 1.5]);
        ylabel('Speed (mm/s)');
        xlabel('Time (s)');
        ylim([0 20]);
        setmyplot_balazs

        % Creating the example plots for speedup
        kspeedup = speedStPerdrosi2sec{mut}(find(speedupFirst{mut}' == 1),:); 

        for i = 1 : 1 %size(kspeedup,1)
            figure; plot(-2:0.01:2,smooth(kspeedup(i,:),'linear',7));
            title('speedup');
            rectangle('Position',[-0.20,0,0.20,20],'FaceColor',[0.9290 0.6940 0.1250 0.33],'EdgeColor','none','LineWidth',3);
            rectangle('Position',[-0.05,0,0.05,20],'FaceColor',[0.8500 0.3250 0.0980 0.50],'EdgeColor','none','LineWidth',3);
            rectangle('Position',[0,0,1,20],'FaceColor',[0.4660 0.6740 0.1880 0.33],'EdgeColor','none','LineWidth',3);
            xlim([-1 1.5]);
            ylabel('Speed (mm/s)');
            xlabel('Time (s)');
            setmyplot_balazs
        end
        
        % Creating the mean plots for speedup
        figure;
        errorshade(-2:0.01:2,mean(kspeedup),std(kspeedup)/sqrt(size(kspeedup,1)));
        title('speedup first');
        rectangle('Position',[-0.20,0,0.20,20],'FaceColor',[0.9290 0.6940 0.1250 0.33],'EdgeColor','none','LineWidth',3);
        rectangle('Position',[-0.05,0,0.05,20],'FaceColor',[0.8500 0.3250 0.0980 0.50],'EdgeColor','none','LineWidth',3);
        rectangle('Position',[0,0,1,20],'FaceColor',[0.4660 0.6740 0.1880 0.33],'EdgeColor','none','LineWidth',3);
        xlim([-1 1.5]);
        ylabel('Speed (mm/s)');
        xlabel('Time (s)');
        ylim([0 20]);
        setmyplot_balazs
        
        % Creating the example plots for no reaction
        kno = speedStPerdrosi2sec{mut}(find(noReaction{mut}' == 1),:);
        for i = 1 : 1 %size(kno,1)
            figure; plot(-2:0.01:2,smooth(kno(i,:),'linear',7));
            title('noreact');
            rectangle('Position',[-0.20,0,0.20,20],'FaceColor',[0.9290 0.6940 0.1250 0.33],'EdgeColor','none','LineWidth',3);
            rectangle('Position',[-0.05,0,0.05,20],'FaceColor',[0.8500 0.3250 0.0980 0.50],'EdgeColor','none','LineWidth',3);
            rectangle('Position',[0,0,1,20],'FaceColor',[0.4660 0.6740 0.1880 0.33],'EdgeColor','none','LineWidth',3);
            xlim([-1 1.5]);
            ylabel('Speed (mm/s)');
            xlabel('Time (s)');
            setmyplot_balazs
        end
        
        % Creating the mean plots for no reaction
        figure;
        errorshade(-2:0.01:2,mean(kno),std(kno)/sqrt(size(kno,1)));
        title('Not reacted');
        rectangle('Position',[-0.20,0,0.20,20],'FaceColor',[0.9290 0.6940 0.1250 0.33],'EdgeColor','none','LineWidth',3);
        rectangle('Position',[-0.05,0,0.05,20],'FaceColor',[0.8500 0.3250 0.0980 0.50],'EdgeColor','none','LineWidth',3);
        rectangle('Position',[0,0,1,20],'FaceColor',[0.4660 0.6740 0.1880 0.33],'EdgeColor','none','LineWidth',3);
        xlim([-1 1.5]);
        ylabel('Speed (mm/s)');
        xlabel('Time (s)');
        ylim([0 20]);
        setmyplot_balazs   

end