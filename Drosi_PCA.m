function Drosi_PCA(Red_mutant_groups_13tun)
%DROSI_PCA   Principal Component Analysis.
%   DROSI_PCA(DATA) performs Principal Component Analyisis on the stimulus
%   triggered speed data followed by a cluster analysis approach based on
%   the first 3 principal components explaining the most variance, to
%   identify the main characteristic reaction types to the threatening
%   stimuli.
%
%   Required input arguments:
%       RED_MUTANT_GROUPS_13TUN: preprocessed results struct calculated
%       by the DROSOPHILA_SPEED funtion.
%
%   See also DROSOPHILA_SPEED

%   Márton Kajtor
%   Institute of Experimental Medicine, Budapest, Hungary
%   Lendulet Laboratory of Systems Neuroscience
%   13-Jul-2023

for mut = 1 : size(Red_mutant_groups_13tun,2)
    speedStPerdrosi2sec{mut} = Red_mutant_groups_13tun(mut).age_groups.speed_st_perdrosi_2sec;
    prePcaSt2Sec = speedStPerdrosi2sec{mut}(:,180:300);
    PcaSt2Sec{mut} = prePcaSt2Sec;
    
    % PCA and clustering reactions
    % set paramters
    PCA_num = 3 % number of principle components used for clustering
    clustnum = 8 % number of required clusters
    y =PcaSt2Sec{mut}';
    window = 1:120; % window for pca
    [PC,PCA1]=pca(y(window,:)'); % running principal component analysis
    PCA2=PCA1(:,1:PCA_num); % extracting the first two pricinpal components
    Dend = linkage(PCA2,'complete','euclidean'); % calculating the dendrogram of psths
    Clusters = cluster(Dend,'maxclust',clustnum); % clustering based on the tree
    D = pdist(PCA2); % euclidean distrance between point in the artifical space
    cutoff = Dend(end-clustnum+2,3); % cutting the tree to get 'clustnum' clusters
    set(gca,'Ydir','reverse');
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    for pcaind = 1:clustnum
        figure; % figure for heatmap
        imagesc(speedStPerdrosi2sec{mut}(Clusters == pcaind,:));
        ylabel('Trials');
        xlabel('Time');
        xlim([100 350]);
        setmyplot_balazs
        figure;  % figure for mean speed plot
        plot(mean(speedStPerdrosi2sec{mut}(Clusters == pcaind,:)),'b-');
        title(['PCA cluster num: ' num2str(pcaind)]);
        line([195 195],[0 0.03],'Color',[0.6350 0.0780 0.1840]);
        xlim([100 350]);
        ylabel('mm/ms');
        xlabel('time');
        setmyplot_balazs
    end
    
    close all
end