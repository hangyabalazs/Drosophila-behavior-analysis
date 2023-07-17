function Main_Drosophila_Behaviour
%MAIN_DROSOPHILA_BEHAVIOUR   Main wrapper for the figures presented in the
%Kajtor et al. manuscript.
%   This codes quantifies drosophila fear behaviour and creates the
%   processed data files and figures from the raw .csv files containing the
%   spatial coordinates of individual flies in the function of time. The
%   function creates a new folder in the data folder to save the result
%   files.

%   Márton Kajtor
%   Institute of Experimental Medicine, Budapest, Hungary
%   Lendulet Laboratory of Systems Neuroscience 
%   13-Jul-2023

% Genotype comparison
path = 'F:\Marci\Cikkhez\raw data\2 weeks old raw\'; % define path to the dataset
sfilename = 'Red_mutant_groups_13tun.mat'; % define results file name 
Red_mutant_groups_13tun = Drosophila_speed('*old raw',path,sfilename);

% Figure 2
Drosi_PCA(Red_mutant_groups_13tun,1)
Behav_examples(Red_mutant_groups_13tun)

% Figure 3
Boxplots_Comparison(Red_mutant_groups_13tun)

% Figure 4
Reaction_Speed_Correl(Red_mutant_groups_13tun)

% Age comparison
path = 'F:\Marci\Cikkhez\raw data\age groups raw\'; % define path to the dataset
sfilename = 'Red_mutant_groups_13tun.mat'; % define results file name 
Red_mutant_groups_13tun = Drosophila_speed('*groups raw',path,sfilename);

% Figure 5
Age_groups_plots(Red_mutant_groups_13tun)