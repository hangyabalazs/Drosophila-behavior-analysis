function [Red_mutant_groups_13tun] = Drosophila_speed(testname,mainpath,sfilename)
%DROSIPHILA_SPEED   Quantifies drosophila fear behaviour.
%   DROSOPHILA_SPEED(TESTNAME,MAINPATH,SFILENAME) analyzes Drosophila
%   spatial coordinate in the function of the timastamps, in relation to
%   shadow onset and offset times loaded% from spearate files. The code
%   detects and quantifies the beavioural properties of each individual fly
%   and saves the results.
%
%   Required input arguments:
%       TESTNAME: Possible inputs:
%                 '*2 weeks old raw' -- run the analysis on the 2 weeks old
%                                       flies for comparion between
%                                       genotypes (data loaded from the
%                                       folder named same as the TESTNAME
%                                       input)
%                 '*age groups raw' -- run the analysis on the 5 age groups
%                                      dataset for comparioson between age
%                                      groups and genotypes
%                                      (data loaded from the folder named
%                                      same as the TESTNAME input)
%       MAINPATH: Main data path SFILENAME: Name of the results file

%   Márton Kajtor
%   Institute of Experimental Medicine, Budapest, Hungary
%   Lendulet Laboratory of Systems Neuroscience
%   13-Jul-2023

% Sampling rate
sr = 100;   % sampling rate

% Speed test parameters
coordinate_fail_tolerance = 10;   % detect discontinuity for position correction
shadow_onset_correction_factor = 0;    % 0 ms
nTunnels = 13;   % number of parallel drosi runways
instability = 130;  % duration of the movement by instabil stops

% Tunnel length measurements
xtun = 53 * ones(1,13);  % the real lenght of the tunnel
ytun = 5 * ones(1,13); % the real width of the tunnel
age_groups = []; % age group struct preallocation
Red_mutant_groups = []; % storage structure preallocation
dir_ch_dens_trh = 0.21; % constans for detecting standing animal detection instability

% Data directory
com_dir = fullfile(mainpath);  %  directory
com_fold = dir(strcat(com_dir, testname)); % select the two existing folder
n_c_Folders = numel(com_fold);  % the number of mutant folders

for cA = 1 : n_c_Folders
    c_foldername = com_fold(cA).name;
    mut_dir = fullfile([com_dir c_foldername '\']);  % mutant directory
    mut_fold = dir(strcat(mut_dir, '*,')); % select the folders with specific symbole
    n_m_Folders = numel(mut_fold);  % the number of mutant folders
    for iM = 1:n_m_Folders  % mutant groups
        m_foldername = mut_fold(iM).name;                % current folder's name
        my_m_Folder = fullfile(mut_dir,m_foldername);    % combines the path with the current folder name
        cd(my_m_Folder); % open mutant folder
        age_dir = fullfile(([mut_dir '\' m_foldername '\']));
        age_fold = dir(strcat(age_dir, '*old'));   % find the folders of age groups
        n_a_Folders = numel(age_fold); % the number of the age group folders
        
        for iA = 1:n_a_Folders              % age groups
            
            Shadow_sorted_stop_latency = [];
            Shadow_sorted_stopdur = [];
            shadow_sorted_freezing = [];
            shadow_sorted_slow = [];
            Shadow_sorted_speedup = [];
            Shadow_sorted_turn = [];
            Shadow_sorted_mean_speed_before_200 = [];
            vercat_speed_st_perdrosi_2sec = [];
            vercat_speed_st_perdrosi_4sec = [];
            ShadowSortedSlowingReactionTime = [];
            ShadowSortedSlowingTime = [];
            ShadowSortedStoppingTime = [];
            ShadowSortedStoppReactionTime = [];
            ShadowSortedSlowingFirst = [];
            ShadowSortedSpeedingReactionTime = [];
            ShadowSortedSpeedingTime = [];
            ShadowSortedSpeedupFirst = [];
            ShadowSortedSlowStop = [];
            ShadowSortedFastStop = [];
            ShadowSortedisstop_before = [];
            
            a_foldername = age_fold(iA).name; % age folders folder name
            my_a_Folder = fullfile(age_dir,a_foldername); % age group folder directory
            cd(my_a_Folder); % open age groups folder
            
            directory = fullfile([mut_dir '\' m_foldername '\' a_foldername '\']);
            folders = dir(strcat(directory,'*_1'));  % first recording of the day; if only relevant folders present, use '*\'
            nFolders = numel(folders);      % number of folders in the directory
            
            for iF = 1:nFolders
                foldername = folders(iF).name;                % current folder's name
                myFolder = fullfile(directory,foldername);    % combines the path with the current folder name
                cd(myFolder); % open the folder
                
                % Preallocation
                [sum_freezing_sess,sum_run_after_sess,sum_turn_after_sess,...
                    sum_souldstop_before_sess,avg_stop_latency_sess,...
                    avg_start_latency_sess,avg_stop_duration_sess, avg_stop_ratio_sess,...
                    avg_walk_ratio_sess, sum_start_sess, sum_slow_after_sess,...
                    avg_speed_wout_stop_permut] = deal(zeros(1,nFolders));
                
                [sum_freezing_person,sum_run_after_person,sum_turn_after_person,...
                    sum_shouldstop_before_person, avg_stop_latency_person,...
                    avg_start_latency_person,avg_stop_duration_person,...
                    NumShadows, sum_speed_wout_stop_persess, collected_data, sum_isstart_after,...
                    stop_ratio, walk_ratio, sum_slow_after] = deal(zeros(1,nTunnels));
                
                [remd_nan_start_lateny,remd_nan_stop_lateny,remd_nan_stop_duration,...
                    freezing,remd_nan_rstart_lateny,remd_nan_rstop_lateny,all_avg_speed_perdrosi,...
                    remd_nan_rstop_duration, speed_st_perdrosi_2sec, speed_st_perdrosi_4sec,...
                    Manual_isstop_before,Manual_Instabil_isstop_after,Manual_isstop_after,...
                    Manual_isslow_after, Manual_stop_duration,manual_wn,Manual_check,...
                    ShadowOffset_forcalc, ShadowOnset_forcalc,xmax_tun, ymax_tun, xmin_tun,...
                    ymin_tun, max_matrix, min_matrix,size_matrix,calc_dist_per_coordinate,...
                    error_list, all_start_lat_persess,all_stop_dur_persess, timediff,...
                    all_stop_lat_persess,all_rstart_lat_persess,all_rstop_dur_persess,...
                    all_rstop_lat_persess, rmean_speed_before_200, rmean_speed_after_200,...
                    risstop_after,rstop_latency, rstart_latency, rstop_duration,...
                    rmean_speed_after_600, risrun_after,avg_speed_person,risturn_after,...
                    rmean_direction_before_600, rmean_direction_after_600,  instabil_stop,...
                    mean_speed_before_200, mean_speed_before_600, mean_speed_after_200, isstop_before,...
                    stop_latency, start_latency, stop_duration,real_size, Position, PosDiff, Speed...
                    mean_speed_after_600, isrun_after, isslow_after,Shadow, ShadowOnset, ShadowOffset,...
                    ShadowIndicator, rShadowOnset, isstop_after,risstop_before,Position_b, Position_s...
                    mean_direction_before_600, mean_direction_after_600, isturn_after,isstart_after,...
                    slowingReactionTime,slowingTime,stoppingTime,stoppReactionTime,slowingFirst,speedingReactionTime,...
                    speedingTime,speedupFirst,slowStop,fastStop] = deal(cell(1,nTunnels));
                
                % Calculate accurate shadow arrival with the correction of tunnel width
                % The shadow onset is based on calculation (with the speed of shadow) not just on detection
                % Shadow arrival detection and calculation
                try
                    for iT = 1:2
                        det_filename = ['det_',num2str(iT)];  % the detection file name
                        Shadow{iT} = load([det_filename '.csv']);  % shadow detection data
                        
                        mstime = Shadow{iT}(:,1) * 3600000 + Shadow{iT}(:,2) * 60000 + ...
                            Shadow{iT}(:,3) * 1000 + Shadow{iT}(:,4);         % time to ms
                        Shadow{iT}(end-1:end,2:3) = NaN(2,2);
                        Shadow{iT} = [mstime Shadow{iT}(:,5) Shadow{iT}(:,6)];   % Shadow: timstamps, x pos., y pos. (shadow != NaN)
                        ShadowIndicator{iT} = ~isnan(Shadow{iT}(:,2));   % shadow detection
                        invers_ShadowIndicator = (ShadowIndicator{iT} * -1) + 1;
                        invers_corrected_ShInd = bwareaopen(invers_ShadowIndicator,25);
                        bright_corrected_ShadowIndicator = (invers_corrected_ShInd * -1) + 1;
                        Shadow_corrected_ShadowIndicator = bwareaopen(bright_corrected_ShadowIndicator,sr*4);
                        
                        ShadowOnset_forcalc{iT} = find(diff(Shadow_corrected_ShadowIndicator)>0.5) + 1;   % shadow onset
                        ShadowOffset_forcalc{iT} = find(diff(Shadow_corrected_ShadowIndicator)<-0.5);   % shadow offset
                        
                        if length(ShadowOnset_forcalc{iT}) ~= 40 % if there is less than 40 shadow detected give error
                            error('Error Shadow detection');
                        end
                        
                        if length(ShadowOffset_forcalc{iT}) ~= length(ShadowOnset_forcalc{iT})   % if the session finished during shadow detection (wrong detection)
                            ShadowOnset_forcalc{iT}(end) = [];
                        end
                    end
                    shadow_onset_diff_cell = (ShadowOnset_forcalc{2}(:,1) - ShadowOnset_forcalc{1}(:,1 ));
                    cell_per_mm = median(shadow_onset_diff_cell)/120; % 120 is the distance between the two extrenal tunnel
                    ShadowOnset{1}(:,1) = round(ShadowOnset_forcalc{1}(:,1) - (cell_per_mm * 12)); % calculate the shadow detection error (4 mm from the tunnel edge to the det point)
                    
                    for iT = 1:nTunnels
                        
                        ShadowOnset{iT}(:,1) = round(ShadowOnset{1}(:,1) + (iT-1) * cell_per_mm * 9); % calculating the shadow arrival from the first shadow detection and shadow speed
                        error_list{iT} = {};
                        real_size{iT} = [xtun(iT); ytun(iT)];   % actual tunnel dimensions
                        try
                            if nTunnels == 0
                                error('There is no file in the folder');
                            end
                            try
                                
                                % Load tracking data
                                brighness_filename = ['tun_',num2str(iT) '_B'];   % csv file names
                                shadow_filename = ['tun_',num2str(iT) '_S'];
                                Position_b{iT} = load([brighness_filename '.csv']);  % position data for bright condition
                                Position_s{iT} = load([shadow_filename '.csv']);  % position data for dark condition
                                
                                % Convert timestamps to ms
                                mstime_b = time2ms(Position_b{iT});         % time in ms
                                mstime_s = time2ms(Position_s{iT});         % time in ms
                                Position_b{iT} = [mstime_b Position_b{iT}(:,5) Position_b{iT}(:,6)];   % Position (bright): timstamps, x pos., y pos.
                                Position_s{iT} = [mstime_s Position_s{iT}(:,5) Position_s{iT}(:,6)];   % Position (shadow): timstamps, x pos., y pos.
                                srp = 1 / mean(diff(Position_b{iT}(:,1)/1000));  % sampling rate for position
                                if abs(sr-srp) > 3      % sampling rate check
                                    error('drosophila_speed: Sampling rate error.')
                                end
                                
                                % Merge bright and dark positions
                                snansx = isnan(Position_s{iT}(:,2));
                                snanx = find(snansx);
                                for Col = 2:3
                                    for iNs = 1:length(snanx)
                                        scinx = snanx(iNs);
                                        Position_s{iT}(scinx,Col) = Position_b{iT}(scinx,Col);
                                    end
                                end
                                Position{iT} = Position_s{iT};
                            catch
                                filename = ['tun_',num2str(iT)];
                                Position{iT} = load([filename '.csv']);
                                
                                mstime = time2ms(Position{iT});
                                Position{iT} = [mstime Position{iT}(:,5) Position{iT}(:,6)];
                                srp = 1 / mean(diff(Position{iT}(:,1)/1000));  % sampling rate for position
                                if abs(sr-srp) > 3      % sampling rate check
                                    error('drosophila_speed: Sampling rate error.')
                                end
                            end
                            
                            % Interpolate missing position values
                            pnans = isnan(Position{iT}(:,2));
                            if sum(pnans)*2 > length(Position{iT}(:,2))
                                error('drosophila_speed: Too much missing point')
                            end
                            while isnan(Position{iT}(end,2))
                                Position{iT}(end,:) = [];
                            end
                            nanx = find(pnans);
                            if isempty(nanx) == 0
                                first_nonan = 1;
                                if nanx(1,1) == 1
                                    nonanh = find(~isnan(Position{iT}(:,2)));
                                    first_nonan = nonanh(1,1);
                                    for nantozero = 1 : first_nonan - 1
                                        Position{iT}(nantozero,2) = 0;
                                        Position{iT}(nantozero,3) = 0;
                                    end
                                end
                                for iN = first_nonan:length(nanx)
                                    cinx = nanx(iN);
                                    nanpre = cinx - 1;   % sequential updates ensure that prev. one is non-NaN
                                    nanpost = cinx + find(~pnans(cinx+1:end),1,'first');
                                    Position{iT}(cinx,2) = interp1([nanpre, nanpost],...
                                        [Position{iT}(nanpre,2) Position{iT}(nanpost,2)],cinx);
                                    Position{iT}(cinx,3) = interp1([nanpre, nanpost],...
                                        [Position{iT}(nanpre,3) Position{iT}(nanpost,3)],cinx);
                                end
                            end
                            
                            % Position vector correction
                            if Position{iT}(1,2) == 0                               % if first position is nan or 0 then change it to the first nonzero coocrdinate
                                firstnozero = find(0<Position{iT}(:,2),1,'first');
                                Position{iT}(1,2)= Position{iT}(firstnozero,2);
                                Position{iT}(1,3)= Position{iT}(firstnozero,3);
                            end
                            
                            dpt = abs(diff(Position{iT}(:,2)));
                            try
                                while any(dpt>coordinate_fail_tolerance(1))
                                    firsterror = find(dpt>coordinate_fail_tolerance(1),1,'first') + 1;   % first missing point
                                    inx = 1;
                                    errorinx = firsterror;
                                    while abs(Position{iT}(firsterror+inx,2)-Position{iT}(firsterror-1,2)) > coordinate_fail_tolerance(1)
                                        errorinx = [errorinx firsterror+inx]; %#ok <AGROW>   % find all consequtive missing points
                                        inx = inx + 1;
                                    end
                                    replacevector = interp1([firsterror-1 errorinx(end)+1],Position{iT}([firsterror-1 errorinx(end)+1],2),errorinx);  % interpolate missing values
                                    Position{iT}(errorinx,2) = replacevector;   % replacement, x
                                    
                                    replacevector = interp1([firsterror-1 errorinx(end)+1],Position{iT}([firsterror-1 errorinx(end)+1],3),errorinx);  % interpolate missing values
                                    Position{iT}(errorinx,3) = replacevector;   % replacement, y
                                    
                                    dpt = abs(diff(Position{iT}(:,2)));
                                end
                            catch
                                Position{iT}(errorinx,:) = [];
                            end
                            
                            % Calibrate positions to physical distances
                            xmax_tun{iT} = max(Position{iT}(:,2));   % NOTE: assumes the flies covered the entire tunnel
                            ymax_tun{iT} = max(Position{iT}(:,3));
                            xmin_tun{iT} = min(Position{iT}(:,2));
                            ymin_tun{iT} = min(Position{iT}(:,3));
                            max_matrix{iT} = [xmax_tun{iT}; ymax_tun{iT}];
                            min_matrix{iT} = [xmin_tun{iT}; ymin_tun{iT}];
                            size_matrix{iT} = max_matrix{iT} - min_matrix{iT};  % the length and width of the tunnels
                            calc_dist_per_coordinate{iT} = real_size{iT} ./ size_matrix{iT};  % distance-tracking data conversion
                            
                            % Calculate walking distance in mm
                            PosDiff{iT}(:,1) = diff(smooth(Position{iT}(:,1),'linear',21));   % time difference in ms (Position(6))
                            PosDiff{iT}(:,2) = abs(diff(smooth(Position{iT}(:,2),'linear',21))) * calc_dist_per_coordinate{iT}(1,1);  % absolute translation in x
                            PosDiff{iT}(:,3) = abs(diff(smooth(Position{iT}(:,3),'linear',21))) * calc_dist_per_coordinate{iT}(2,1); % absolute translation in y
                            PosDiff{iT}(:,4) = sqrt(PosDiff{iT}(:,2).^2+PosDiff{iT}(:,3).^2);   % 2D translation distance
                            PosDiff{iT}(:,5) = diff(smooth(Position{iT}(:,2),'linear',21)) * calc_dist_per_coordinate{iT}(1,1); % directional translation in x
                            
                            % Calculate min_speed
                            value_PosDiff = [];
                            value_PosDiff(:,1) = Position{iT}(:,1);
                            value_PosDiff(end) = [];
                            value_PosDiff(:,1) = value_PosDiff(:,1)-value_PosDiff(1);
                            value_PosDiff(:,2) = diff(smooth(Position{iT}(:,1),'linear',15));
                            value_PosDiff(:,3) = abs(diff(smooth(Position{iT}(:,2),'linear',15))) * calc_dist_per_coordinate{iT}(1,1); % absolute translation in x
                            value_PosDiff(:,4) = abs(diff(smooth(Position{iT}(:,3),'linear',15))) * calc_dist_per_coordinate{iT}(2,1);% absolute translation in y
                            value_PosDiff(:,5) = sqrt(value_PosDiff(:,3).^2+value_PosDiff(:,3).^2);   % 2D translation distance
                            value_PosDiff(:,6) = diff(smooth(Position{iT}(:,2),'linear',15)) * calc_dist_per_coordinate{iT}(1,1);   % directional translation in x
                            Speed{iT} = PosDiff{iT}(:,4) ./ PosDiff{iT}(:,1);  % Speed in mm/ms
                            summmmm = value_PosDiff(:,3) + value_PosDiff(:,4);
                            plus_minus(value_PosDiff(:,6) < 0) = -0.01;
                            plus_minus(value_PosDiff(:,6) > 0) = 0.01;
                            plus_minus(value_PosDiff(:,6) == 0) = 0;
                            [~,cross] = valuecrossing(value_PosDiff(:,1)',value_PosDiff(:,6)',0,'both');
                            sm_cross = smooth(cross,'linear',5);
                            location = [];
                            locationDiff = [];
                            distInTimeWind = [];
                            location(:,1) = (Position{iT}(:,2)-min(Position{iT}(:,2)));
                            location(:,2) = (Position{iT}(:,3)-min(Position{iT}(:,3)));
                            location(location(:,1) == 0,1) = 1;
                            location(location(:,2) == 0,2) = 1;
                            locationDiff(:,1) = diff(Position{iT}(:,2)); % x coordinate diff
                            locationDiff(:,2) = diff(Position{iT}(:,3)); % y coordinate diff
                            for timeshift = 1 : length(locationDiff)-1
                                distInTimeWind(timeshift,1)= abs(locationDiff(timeshift,1) + locationDiff(timeshift + 1,1));
                                distInTimeWind(timeshift,2)= abs(locationDiff(timeshift,2) + locationDiff(timeshift + 1,2));
                            end
                            
                            % Calculate walking speed
                            standCellX = find(distInTimeWind(:,1) <= 0.5);
                            standCellY = find(distInTimeWind(:,2) <= 0.3);
                            standconXY = intersect(standCellX(:,1),standCellY(:,1));
                            Speed1 = Speed{iT};
                            Speed1(standconXY(:,1)) = 0;
                            sm_cross_speed = Speed{iT}(find(sm_cross > dir_ch_dens_trh));
                            sort_sm_cross_speed =  sort(sm_cross_speed);
                            min_speed_limit(iT) = sort_sm_cross_speed(round(length(sort_sm_cross_speed) * 0.85));
                            Speed{iT} = smooth(Speed1,'linear',15);
                            if min_speed_limit(iT) > 0.004
                                error('Too high resolution by standting (wrong bonsai tracking)')
                            end
                            NumShadows(iT) = length(ShadowOnset{iT});
                            nzr = sum(abs(timediff{iT})>2) / length(timediff{iT});   % ratio of shifted timestamps w 2 ms tolerance
                            try
                                if nzr > 0.15 || any(abs(timediff{iT})>50)   % 1.5% tolerance for >15 ms timestamp differences
                                    error('drosophila_speed:timestamp mismatch.')
                                end
                            catch
                            end
                            
                            % Shadow overlay plot
                            sum_speed_wout_stop_persess(iT) = sum(Speed{iT}(find(Speed{iT}> min_speed_limit(iT)),1))...
                                /length(find(find(Speed{iT}> min_speed_limit(iT))));
                            collected_data(iT,10) = sum_speed_wout_stop_persess(iT);
                            
                            non_shadow_index = find(Shadow_corrected_ShadowIndicator == 0);   % walk stand ratio without shadowtime
                            non_shadow_index(end) = [];
                            try
                                stop_ratio(iT) = length(find(Speed{iT}(non_shadow_index) < min_speed_limit(iT)))...
                                    /length(non_shadow_index);
                                walk_ratio(iT) = length(find(Speed{iT}(non_shadow_index) > min_speed_limit(iT)))...
                                    /length(non_shadow_index);
                            catch
                                non_shadow_index_err = non_shadow_index;
                                try
                                    non_shadow_index_err(find(non_shadow_index(lenght(Speed{iT}))):end) = [];
                                catch
                                    non_shadow_index_err(end-1:end) = [];
                                    
                                end
                            end
                            
                            % Stop detect
                            [mean_speed_before_200{iT}, mean_speed_after_200{iT}, isstop_before{iT}, isstop_after{iT},...
                                mean_speed_before_600{iT},stop_latency{iT}, start_latency{iT}, stop_duration{iT}, speed_st_2sec,...
                                mean_speed_after_600{iT}, isrun_after{iT},isslow_after{iT}, avg_speed_sta_2sec,...
                                avg_speed_sta_4sec, mean_direction_before_600{iT}, mean_direction_after_600{iT},...
                                speed_st_4sec, isturn_after{iT},isstart_after{iT}, instabil_stop{iT},...
                                Manual_isstop_before{iT},Manual_Instabil_isstop_after{iT},Manual_isstop_after{iT},...
                                Manual_isslow_after{iT},Manual_stop_duration{iT},manual_wn{iT},Manual_check{iT},...
                                slowingReactionTime{iT},slowingTime{iT},stoppingTime{iT},stoppReactionTime{iT},slowingFirst{iT},speedingReactionTime{iT},...
                                speedingTime{iT},speedupFirst{iT},slowStop{iT},fastStop{iT}] = stop_detection(ShadowOnset{iT},...
                                Speed{iT},Position{iT},PosDiff{iT}(:,5),sr,min_speed_limit(iT), shadow_onset_correction_factor,...
                                instability,sm_cross,dir_ch_dens_trh);
                            
                            speed_st_perdrosi_2sec{iT} = speed_st_2sec;
                            speed_st_perdrosi_2sec{iT}(speed_st_perdrosi_2sec{iT} < min_speed_limit(iT)) = 0;
                            
                            avg_heatmap_speed_persess_2sec(iT,:) = avg_speed_sta_2sec;
                            speed_st_perdrosi_4sec{iT} = speed_st_4sec;
                            speed_st_perdrosi_4sec{iT}(speed_st_perdrosi_4sec{iT} < min_speed_limit(iT)) = 0;
                            avg_heatmap_speed_persess_4sec(iT,:) = avg_speed_sta_4sec;
                            
                            % Remove nans from matrix
                            remd_nan_start_lateny{iT} = rem_nan(start_latency{iT});
                            remd_nan_stop_lateny{iT} = rem_nan(stop_latency{iT});
                            remd_nan_stop_duration{iT} = rem_nan(stop_duration{iT});
                            
                            remd_nan_rstart_lateny{iT} = rem_nan(rstart_latency{iT});
                            remd_nan_rstop_lateny{iT} = rem_nan(rstop_latency{iT});
                            remd_nan_rstop_duration{iT} = rem_nan(rstop_duration{iT});
                            
                            % Sum freezing per fly
                            for sms = 1 : size(isstop_after{iT},2)
                                if isstop_after{iT}(1,sms) == 1
                                    freezing{iT}(1,sms) =  isstop_after{iT}(1,sms) - isstop_before{iT}(1,sms);
                                else
                                    freezing{iT}(1,sms) = 0;
                                end
                            end
                            if size(freezing{iT},2) == 40
                                shadow_sorted_freezing = vertcat(shadow_sorted_freezing,freezing{iT});
                                shadow_sorted_slow = vertcat(shadow_sorted_slow,isslow_after{iT});
                            end
                            avg_speed_person{iT} = sum(Speed{iT}(:,1))/sum(length(Position{iT}(:,1)));
                            
                            Shadow_sorted_stop_latency = vertcat(Shadow_sorted_stop_latency,vertcat(stop_latency{iT}));
                            Shadow_sorted_stopdur = vertcat(Shadow_sorted_stopdur,vertcat(stop_duration{iT}));
                            Shadow_sorted_speedup = vertcat(Shadow_sorted_speedup,vertcat(isrun_after{iT}));    % 1 or 0
                            Shadow_sorted_turn = vertcat(Shadow_sorted_turn,vertcat(isturn_after{iT}));   % 1 or 0
                            Shadow_sorted_mean_speed_before_200 = vertcat(Shadow_sorted_mean_speed_before_200,vertcat(mean_speed_before_200{iT}));
                            vercat_speed_st_perdrosi_2sec = vertcat(vercat_speed_st_perdrosi_2sec,vertcat(speed_st_perdrosi_2sec{iT}));
                            vercat_speed_st_perdrosi_4sec = vertcat(vercat_speed_st_perdrosi_4sec,vertcat(speed_st_perdrosi_4sec{iT}));
                            ShadowSortedSlowingReactionTime = vertcat(ShadowSortedSlowingReactionTime,vertcat(slowingReactionTime{iT}));
                            ShadowSortedSlowingTime = vertcat(ShadowSortedSlowingTime,vertcat(slowingTime{iT}));
                            ShadowSortedStoppingTime = vertcat(ShadowSortedStoppingTime,vertcat(stoppingTime{iT}));
                            ShadowSortedStoppReactionTime = vertcat(ShadowSortedStoppReactionTime,vertcat(stoppReactionTime{iT}));
                            ShadowSortedSlowingFirst = vertcat(ShadowSortedSlowingFirst,vertcat(slowingFirst{iT}));
                            ShadowSortedSpeedingReactionTime = vertcat(ShadowSortedSpeedingReactionTime,vertcat(speedingReactionTime{iT}));
                            ShadowSortedSpeedingTime = vertcat(ShadowSortedSpeedingTime,vertcat(speedingTime{iT}));
                            ShadowSortedSpeedupFirst = vertcat(ShadowSortedSpeedupFirst,vertcat(speedupFirst{iT}));
                            ShadowSortedSlowStop = vertcat(ShadowSortedSlowStop,vertcat(slowStop{iT}));
                            ShadowSortedFastStop = vertcat(ShadowSortedFastStop,vertcat(fastStop{iT}));
                            ShadowSortedisstop_before = vertcat(ShadowSortedisstop_before,vertcat(isstop_before{iT}));
                        catch
                            close all
                            continue
                        end
                        close all
                    end
                catch
                    close all
                    continue
                end
                cd(myFolder);
            end
            
            if nFolders ~= 0 || n_a_Folders ~= 0
                
                % Save data in struct
                age_groups(iA).name = a_foldername;
                age_groups(iA).speed_st_perdrosi_2sec = vercat_speed_st_perdrosi_2sec;
                age_groups(iA).speed_st_perdrosi_4sec = vercat_speed_st_perdrosi_4sec;
                age_groups(iA).collected_freezing = shadow_sorted_freezing;
                age_groups(iA).Shadow_sorted_stop_latency = Shadow_sorted_stop_latency;
                age_groups(iA).Shadow_sorted_stopdur = Shadow_sorted_stopdur;
                age_groups(iA).Shadow_sorted_mean_speed_before_200 = Shadow_sorted_mean_speed_before_200;
                age_groups(iA).Shadow_sorted_speedup = Shadow_sorted_speedup;
                age_groups(iA).Shadow_sorted_slow = shadow_sorted_slow;
                age_groups(iA).Shadow_sorted_turn = Shadow_sorted_turn;
                age_groups(iA).ShadowSortedSlowingReactionTime = ShadowSortedSlowingReactionTime;
                age_groups(iA).ShadowSortedSlowingTime = ShadowSortedSlowingTime;
                age_groups(iA).ShadowSortedStoppingTime = ShadowSortedStoppingTime;
                age_groups(iA).ShadowSortedStoppReactionTime = ShadowSortedStoppReactionTime;
                age_groups(iA).ShadowSortedSlowingFirst = ShadowSortedSlowingFirst;
                age_groups(iA).ShadowSortedSpeedingReactionTime = ShadowSortedSpeedingReactionTime;
                age_groups(iA).ShadowSortedSpeedingTime = ShadowSortedSpeedingTime;
                age_groups(iA).ShadowSortedSpeedupFirst = ShadowSortedSpeedupFirst;
                age_groups(iA).ShadowSortedSlowStop = ShadowSortedSlowStop;
                age_groups(iA).ShadowSortedFastStop = ShadowSortedFastStop;
                age_groups(iA).ShadowSortedStandBefore = ShadowSortedisstop_before;
                
                close all
            else
                age_groups(iA).name = a_foldername;
                age_groups(iA).speed_st_perdrosi_2sec = [];
                age_groups(iA).speed_st_perdrosi_4sec = [];
                age_groups(iA).collected_freezing = [];
                age_groups(iA).Shadow_sorted_stop_latency = [];
                age_groups(iA).Shadow_sorted_stopdur = [];
                age_groups(iA).Shadow_sorted_mean_speed_before_200 = [];
                age_groups(iA).Shadow_sorted_speedup = [];
                age_groups(iA).Shadow_sorted_slow = [];
                age_groups(iA).Shadow_sorted_turn = [];
                age_groups(iA).ShadowSortedSlowingReactionTime = [];
                age_groups(iA).ShadowSortedSlowingTime = [];
                age_groups(iA).ShadowSortedStoppingTime = [];
                age_groups(iA).ShadowSortedStoppReactionTime = [];
                age_groups(iA).ShadowSortedSlowingFirst = [];
                age_groups(iA).ShadowSortedSpeedingReactionTime = [];
                age_groups(iA).ShadowSortedSpeedingTime = [];
                age_groups(iA).ShadowSortedSpeedupFirst = [];
                age_groups(iA).ShadowSortedSlowStop = [];
                age_groups(iA).ShadowSortedFastStop = [];
                age_groups(iA).ShadowSortedStandBefore = [];
                
            end
        end
        Red_mutant_groups_13tun(iM).names = m_foldername;
        Red_mutant_groups_13tun(iM).age_groups = age_groups;
        
    end
    cd(com_fold);
    save(sfilename, 'Red_mutant_groups_13tun');
end

% -------------------------------------------------------------------------
function [mean_speed_before_200, mean_speed_after_200, isstop_before, isstop_after,...
    mean_speed_before_600, stop_latency, start_latency, stop_duration,speed_st_2sec,...
    mean_speed_after_600, isrun_after, isslow_after,avg_speed_sta_2sec,...
    avg_speed_sta_4sec, mean_direction_before_600, mean_direction_after_600,  speed_st_4sec, isturn_after,...
    isstart_after, instabil_stop, Manual_isstop_before,Manual_Instabil_isstop_after,...
    Manual_isstop_after, Manual_isslow_after,Manual_stop_duration,manual_wn,Manual_check,...
    slowingReactionTime,slowingTime,stoppingTime,stoppReactionTime,slowingFirst,speedingReactionTime,...
    speedingTime,speedupFirst,slowStop,fastStop] = stop_detection(ShadowOnset,Speed,Position,...
    DirectionalX,sr,min_speed_limit_iT,shadow_onset_correction_factor,...
    instability,sm_cross,dir_ch_dens_trh)
% STOP_DETECTION qulifies the behaviour after the shadow onset.

% Shadow-triggered speed average
[st, sta, ~, ~] = stacall(ShadowOnset,Speed,sr,sr*4);  % 2sec before and after shadow onset
speed_st_2sec = st;
avg_speed_sta_2sec = sta;

[st, sta, ~, ~] = stacall(ShadowOnset,Speed,sr,sr*8); % 4 sec before and after shadow onset
speed_st_4sec = st;
avg_speed_sta_4sec = sta;

close all

% Stop detection
NumShadows = length(ShadowOnset);

% Preallocation
[mean_speed_before_200, mean_speed_after_200, isstop_before, isstop_after,...
    stop_latency, start_latency, stop_duration,...
    mean_speed_before_600, mean_speed_after_600, isrun_after,...
    mean_direction_before_600, mean_direction_after_600, isturn_after, Manual_isstop_before,...
    Manual_Instabil_isstop_after, Manual_isstop_after, Manual_isslow_after,...
    Manual_stop_duration, isrun_before, instabil_stop, isstart_after, isslow_after,...
    manual_wn,Manual_check,slowingFirst,slowingReactionTime,slowingTime,stoppingTime,...
    stoppReactionTime,speedingReactionTime,speedingTime,speedUpFirst,slowStop,fastStop,speedupFirst] = deal(nan(1,NumShadows));

for iS = 1:NumShadows
    so = ShadowOnset(iS)-shadow_onset_correction_factor*sr;
    wn = 0.2 * sr;    % 200 ms window
    
    % Calculating mean speeds before shadow
    mean_speed_before_200(iS) = mean(Speed(so-wn:so));  % average speed in before-window  0.2 sec
    mean_speed_after_200(iS) = mean(Speed(so:so+wn));   % average speed in after-window  0.2 sec
    mean_speed_before_100(iS) = mean(Speed(so-wn/2:so)); % for slow and runs
    mean_cross_before_200(iS) = mean(sm_cross(so-wn:so));
    mean_speed_before_600(iS) = mean(Speed(so-3*wn:so));
    mean_speed_after_600(iS) = mean(Speed(so:so+3*wn));
    isstop_before(iS) = mean_speed_before_200(iS) <= min_speed_limit_iT;
    isrun_before(iS) = ~isstop_before(iS);
    
    % Detecting stop or start after the shadow arrival
    dvc = valuecrossing(Position(so:so+sr*10,1)', Speed(so:so+sr*10)',min_speed_limit_iT ,'down') - Position(so);   % 10 sec max latency
    uvc = valuecrossing(Position(so:end-1,1)', Speed(so:end)', min_speed_limit_iT ,'up') - Position(so);  % interpolate exact starting time
    % if there is a false spike caused by shadow then delet this (too sort starts 5 ms)
    walk_duration = [];
    if ~isempty(dvc) && length(dvc)>1 && any(uvc<dvc(2))  && ~isstop_before(iS)
        for belsodvc = 2:length(dvc)
            walk_duration(belsodvc) = dvc(belsodvc)-uvc(find(uvc<dvc(belsodvc),1,'last'));
        end
    end
    if ~isempty(dvc) && any(uvc<dvc(1))  &&  any(isstop_before)
        for belsodvc = 1:length(dvc)
            walk_duration(belsodvc) = dvc(belsodvc)-uvc(find(uvc<dvc(belsodvc),1,'last'));
        end
    end
    
    % Instabil stop detection
    instab_ind = find(walk_duration < instability & walk_duration > 15);
    if any(walk_duration < instability) && any(walk_duration>15) && any(instab_ind) && ~isempty(dvc)...
            && any(uvc<dvc(instab_ind(1)))  && isstop_before(iS) == 0
        instabil_stop(iS) = 1;
    else
        instabil_stop(iS) = 0;
    end
    if any(walk_duration < instability & walk_duration ~=0) && ~isstop_before(iS)
        bad_dvc = find(walk_duration < instability & walk_duration ~= 0);
        bad_uvc = find(walk_duration < instability & walk_duration ~= 0)-1;
        bad_uvc(find(bad_uvc == 0)) = [];
        uvc(bad_uvc) = [];
        dvc(bad_dvc) = [];
    end
    
    % Stop latency calculation
    if    isstop_before(iS) == 0 && ~isempty(uvc) && ~isempty(dvc)
        try
            stop_latency(iS) = dvc(1);
        catch
        end
        isstop_after(iS) = 1;
    else
        stop_latency(iS) = NaN;
        isstop_after(iS) = 0;
    end
    if ~isempty(dvc) && ~isempty(uvc)
        try
            start_latency(iS) = uvc(find(uvc>dvc(1),1,'first'));
        catch
        end
        isstart_after(iS) = 1;
    else
        start_latency(iS) = NaN;
        isstart_after(iS) = 0;
    end
    if isstop_before(iS) == 0
        stop_duration(iS) = start_latency(iS) - stop_latency(iS);   % pause duration
    end
    
    % Turn detection
    mean_direction_before_600(iS) = mean(DirectionalX(so-3*wn:so));  % average direction in before-window
    mean_direction_after_600(iS) = mean(DirectionalX(so:so+3*wn));  % average direction in after-window
    isturn_after(iS) = ~isstop_after(iS) && ...
        mean_direction_before_600(iS) * mean_direction_after_600(iS) < 0;   % change direction after shadow
    
    % Speeding up or slowing down detection
    % average speed in before-window 3*wn:so
    % average speed to 600 after-window (there supposed to be a pause before the jump/run)
    reacted = 0;
    for middle = so : so + 100           % If speed equals 0, this cycle will run
        if reacted ~= 1
            slowStart = [];
            speedStart = [];
            slowEnd = [];
            speedEnd = [];
            if mean(diff(Speed(middle : 2 : middle + 10))) < -0.0002  % slowing detection
                slowStart = valuecrossing(Position(middle-5:middle+10,1)',smooth(diff(Speed(middle-6:middle+10)),'linear',3)',-0.0001,'down') - Position(so);
                if isempty(slowStart)
                    slowStart = 0;
                end
                line([(middle-so)+200,(middle-so)+200],[0,0.01],'Color','magenta')
                slowEnd = valuecrossing(Position(middle+10:middle+100,1)',smooth(diff(Speed(middle+10:middle+101)),'linear',3)',-0.000000000000000001,'up') - Position(so);
                slowingReactionTime(iS) = slowStart(1);
                slowingTime(iS) = slowEnd(1) - slowStart(1);
                if isstop_after(iS) == 1 && slowEnd(1) * 1.2 >= stop_latency(iS)
                    stoppingTime(iS) = dvc(1) - slowingReactionTime(iS);
                    stoppReactionTime(iS) = slowStart(1);
                    stopStart = slowStart(1);
                    reacted = 1;
                elseif isstop_after(iS) == 1 && slowEnd(1) * 1.2 < start_latency(iS)
                    slowingFirst(iS) = 1;
                    reacted = 1;
                end
            end
            try
                if mean(diff(Speed(middle : 2 : middle + 10))) > 0.0002  % speed up detection
                    speedStart = valuecrossing(Position(middle-5:middle+10,1)',smooth(diff(Speed(middle-6:middle+10)),'linear',3)',0.0001,'up') - Position(so);
                    speedEnd = valuecrossing(Position(middle+10:middle+100,1)',smooth(diff(Speed(middle+10:middle+101)),'linear',3)',0.000000000000000001,'down') - Position(so);
                    if isempty(speedStart)
                        speedStart = 0;
                    end
                    speedingReactionTime(iS) = speedStart(1);
                    speedingTime(iS) = speedEnd(1) - speedStart(1);
                    reacted = 1;
                    if isstop_after(iS) == 1 && speedEnd(1) * 1.2 < start_latency(iS) || isstop_after(iS) ~= 1
                        speedupFirst(iS) = 1;
                    end
                end
            catch
            end
        end
    end
    
    % Counting fast and slow stops
    try
        if isstop_after(iS) == 1 && (isempty(slowStart) || slowingFirst(iS) == 1)
            stopStart = valuecrossing(Position((so:so + round(dvc(1)/10)),1)',diff(smooth(Speed((so:(so + 1 + round(dvc(1)/10)))),'linear',3)'),0,'down') - Position(so) + 10;
            stoppingTime(iS) = dvc(1) - stopStart(end);  % if the deceleration speed is not 1.3* greather or slower then the current speed
            stoppReactionTime(iS) = stopStart(end);
        end
        if isstop_after(iS) == 1 && isempty(slowStart)
            stopSlope = abs(mean(diff(Speed((so + round(stopStart(end)/10)): 2 :so + round(dvc(1)/10)))));
            if stopSlope(1) < 0.00025
                slowStop(iS) = 1;
            elseif stopSlope(1) >= 0.00025
                fastStop(iS) = 1;
            end
        end
    catch
        
    end
end

% -------------------------------------------------------------------------
function [X2, S] = smooth(X,str,wn)
%SMOOTH   Smooth with a moving average.
%   [O S] = SMOOTH(X,M,WN) performs smoothing on X by averaging in a
%   sliding window of WN size. M describes the data type, i.e. 'linear' or
%   'circular'. Output is returned in O with standard error of the
%   smoothing in S.
%
%   See also CONV.

n = wn;
nn = (n - 1) / 2;
m = length(X);
X2 = zeros(size(X));
S = zeros(size(X));
switch str
    case 'linear'
        for k = 1:m
            ind1 = max(1,k-nn);
            ind2 = min(k+nn,m);
            X2(k) = mean(X(ind1:ind2));
            S(k) = std(X(ind1:ind2));
        end
    case 'circular'
        for k = 1:m
            if k - nn < 1
                X2(k) = mean([X(mod2(k-nn,m):m); X(1:k+nn)]);
                S(k) = std([X(mod2(k-nn,m):m); X(1:k+nn)]) / sqrt(n);
            elseif k + nn > m
                X2(k) = mean([X(k-nn:m); X(1:mod2(k+nn,m))]);
                S(k) = std([X(k-nn:m); X(1:mod2(k+nn,m))]) / sqrt(n);
            else
                X2(k) = mean(X(k-nn:k+nn));
                S(k) = std(X(k-nn:k+nn)) / sqrt(n);
            end
        end
end

% -------------------------------------------------------------------------
function [t,kross] = valuecrossing(x,y,v,opt)
%VALUECROSSING   Interpolate value-crossings of a function.
%   T = VALUECROSSING(X,Y,V) returns the V-crossings of the function (X,Y).
%   VALUECROSSING uses linear interpolation.
%   T = VALUECROSSING(X,Y,V,OPT) determines upwards, downwards or all
%   value-crossings using 'up', 'down' or 'both' options, respectively.

%   Balazs Hangya
%   Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   New York 11724, USA
%   balazs.cshl@gmail.com
%
%   Institute of Experimental Medicine
%   Szigony street 43, Budapest 1083, Hungary
%   hangyab@koki.hu

% Input argument check
narginchk(3,4)
if nargin == 3
    opt = 'both';
end

% Interpolate crossings
lup = y <= v & [y(2:end) v] > v;
ldown = y >= v & [y(2:end) v] < v;
switch opt
    case 'up'
        fnd = find(lup);
    case 'down'
        fnd = find(ldown);
    case 'both'
        fnd = find(lup|ldown);
end
xf = x(fnd);
xff = x(fnd+1);
yf = y(fnd);
yff = y(fnd+1);
mpl = (v - yf) ./ (yff - yf);
kross(lup == 1) = 1;
kross(lup == 0) = 0;
kross(ldown == 1) = 1;
t = xf + (xff - xf) .* mpl;

% -------------------------------------------------------------------------
function tms = time2ms(t)
%TIME2MS converts the date to ms.
numDate = size(t,2);
switch numDate
    case 6
        tms = t(:,1) * 3600000 + t(:,2) * 60000 + t(:,3) * 1000 + t(:,4);         % time in ms
    case 8
        if any(t(:,3)-mean(t(:,3)))
            error('Marci recorded at midnight.')
        end
        if any(t(:,2)-mean(t(:,2)))
            error('Marci recorded at midnight on the last day of the month.')
        end
        if any(t(:,1)-mean(t(:,1)))
            error('Marci recorded at New Year''s Eve')
        end
        tms = t(:,4) * 3600000 + t(:,5) * 60000 + t(:,6) * 1000 + t(:,7);         % time in ms
    otherwise
        error('Unrecognized timestamp format.')
end

% -------------------------------------------------------------------------
function remd_nans = rem_nan(a)
%REM_NAN Removing NaNs from array and givs it back in order without NaNs.

noNans_in_a = ~isnan(a);           % indicate noNans
loc_noNans_in_a = find(noNans_in_a); % location of noNans
[remd_nans] = deal(zeros(1,sum(noNans_in_a))); % preallocation
lo = 0;
for rnan = 1 : length(loc_noNans_in_a')
    itt = loc_noNans_in_a(1,rnan);
    lo = lo + 1;
    remd_nans(1,lo) = a(1,itt);
end