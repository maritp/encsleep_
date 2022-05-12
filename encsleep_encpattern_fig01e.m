%% enclseep. encoding pattern. 
% identify an encoding pattern by contrasting oscillatory power (1-20Hz)
% during encoding and a control condition (pvt)
% this script (1) creates figure 1e and 
%             (2) saves the encoding topography for each participant 
%                 (in a subject x channel matrix)

% required toolboxes on your path: fieldtrip

clear
close all

github_repo = '/Users/petzka/Documents/GitHub/encsleep_'; %directory of github repo where the scripts live
addpath(genpath(github_repo)) 

datfolder = '/Volumes/MEMTOSH/encsleep'; %directory of data folder

%% where to get ffts
datdir = fullfile(datfolder, 'wake', 'encpvt_fft');
datOI = 'encpvt_fft'; % name of data file

load(fullfile(datdir, datOI))
dat = fft_;
%% where to save
savedir = fullfile(datfolder, 'wake', 'encpvt_topo'); % encoding topographies are saved here ... 
savename = 'encpvt_topo'; %... with this name

%% chose artifact correction method
settings_.perc = 95; % percentage for percentile | should be empty if outlier is defined
settings_.outlier = 0; %define std for upper limit e.g. 2.5 std + mean

%% topo preps
[layout, neighbours] = getting_layout('Bham-64CH-Lay.mat',...
    'Bham-64CH-Neighbours.mat', 0);

%%

for isub = 1:numel(settings_.nsub)
    
    %% exclude x% of trial distribution based on percentile OR outliers OR nothing & average
    for ipart = 1:numel(settings_.parts)
        
        pow_tmp = dat{1,isub}{1,ipart}.powspctrm; %trl x ch x freq
        
        % exclude x% of trial distribution based on percentile
        if ~isempty(settings_.perc) && settings_.outlier == 0
            
            out_tmp = perc_disc(pow_tmp, settings_.perc); %discard part of trial distribution
            pow_mean = squeeze(nanmean(out_tmp));
            
        elseif isempty(settings_.perc) && settings_.outlier ~= 0
            
            pow_mean_tmp = mean(pow_tmp); % calculate mean and std
            pow_std_tmp = std(pow_tmp);
            
            thres_up = pow_mean_tmp + settings_.outlier*pow_std_tmp;
            
            out_tmp = pow_tmp;
            out_tmp(out_tmp > thres_up) = NaN;
            pow_mean = squeeze(nanmean(out_tmp));
            
        else
            
            pow_mean = squeeze(mean(pow_tmp));
            
        end
        
        dat{1,isub}{1,ipart}.powspctrm = pow_mean;
        dat{1,isub}{1,ipart}.dimord = 'chan_freq';
        
        
    end

    %%
    enc{isub} = dat{1,isub}{1,1};
    pvt{isub} = dat{1,isub}{1,2};
    
end

%% grand average, but keep individuals

cfg = [];
cfg.keepindividual = 'yes';
GA_enc = ft_freqgrandaverage(cfg,enc{:});
GA_pvt = ft_freqgrandaverage(cfg,pvt{:});

%% stats

GA = {};

%--- decide here which two conditions to compare
GA{1} = GA_enc; 
GA{2} = GA_pvt; 

conditions = {'enc' 'pvt'};
condcolors = {'k' 'b'};

% if isempty(GA{2})
%     GA{2} = GA{1};
%     GA{2}.powspctrm = zeros(size(GA{1}.powspctrm));
%     conditions{2} = 'zeros';
% end

%% FT stats

cfg                     = [];
cfg.frequency           = [1,20]; 

cfg.channel             = 'all'; 
cfg.statistic           = 'depsamplesT';
cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'cluster'; % 'no', cluster;
cfg.alpha               = .025;
cfg.clusteralpha        = .05;
cfg.tail                = 0;
cfg.correcttail         = 'no'; % alpha prob no
cfg.neighbours          = neighbours; 
cfg.minnbchan           = 3;
cfg.avgovertime         = 'no'; 
cfg.avgoverchan         = 'no';
cfg.avgoverfreq         = 'no';
cfg.computecritval      = 'yes';

cfg.numrandomization    = 1000;

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;

nsub = size(GA{1}.powspctrm,1);
% set up design matrix
design = zeros(2,2*nsub);
for i = 1:nsub
    design(1,i) = i;
end
for i = 1:nsub
    design(1,nsub+i) = i;
end
design(2,1:nsub)        = 1;
design(2,nsub+1:2*nsub) = 2;

cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

% run stats
[Fieldtripstats] = ft_freqstatistics(cfg, GA{:});
length(find(Fieldtripstats.mask))

%% plot TOPO of significant t stats, summed across freq

sig_clusters_neg = unique(Fieldtripstats.negclusterslabelmat .* Fieldtripstats.mask);

if numel(sig_clusters_neg) > 1
    for icluster = 2:numel(sig_clusters_neg)
        
        b = Fieldtripstats.stat .* (Fieldtripstats.negclusterslabelmat==sig_clusters_neg(icluster));
        tsums = squeeze(nansum(b,2));
        tsums = repmat(tsums,[1 size(Fieldtripstats.stat,2)]);
        
        sigmat  = Fieldtripstats.negclusterslabelmat==sig_clusters_neg(icluster);
        sigfreq = round(Fieldtripstats.freq(any(sigmat)));
        
        tmp         = Fieldtripstats;
        tmp.plotme  = tsums;
        
%         figure; set(gcf, 'position', [0 0 1500 1500]) 
        figure
        topo             = [];
        topo.layout      = layout;
        topo.parameter   = 'plotme';
        topo.gridscale   = 360;
        topo.marker      = 'off';
        topo.comment     = 'no';
        topo.style       = 'straight';
        cfg.shading      = 'interp';
        
        ft_topoplotER(topo, tmp);
        cold = colormap('hot');
        cold = fliplr(cold);
        cold = flipud(cold);
        colormap(cold)
        set(gcf, 'Color', 'w');
        colorbar
        
        title(sprintf('sig. neg. cluster %d/%d, %s vs. %s\nsig. @ %0.3f, %s %s corr.\nfreqs:%s',...
            icluster-1,numel(sig_clusters_neg)-1,conditions{1},conditions{2},...
            cfg.alpha, cfg.method,cfg.correctm,num2str(sigfreq)),'interpreter','none');
    end
end

%% plot t-vals of one channel

chOI = {'CP4'};
idx_chOI = find(ismember(GA_enc.label, chOI));

dat_plot = Fieldtripstats.stat(idx_chOI,:);

figure
bar(dat_plot, 'Facecolor', [0.7 0.7 0.7])
rectangle('Position', [sigfreq(1) -4 sigfreq(end)-sigfreq(1) 5],...
    'FaceColor', [0.7 0.7 0.7 0.4], 'EdgeColor', [0.9 0.9 0.9])
ylim([-4 3])
box off

%% write out effect vals into nsubjects x channels matrix
dat = [];
freqOI.fft = {'6 20'}; % define freqOI based on cluster based perm test

for i = 1:numel(freqOI.fft)
    this_FOI = str2num(freqOI.fft{i});
    this_FOI_idx = nearest(GA{1}.freq,this_FOI(1)):nearest(GA{1}.freq,this_FOI(2));
    
    dat{1}{i} = mean(abs(GA{1}.powspctrm(:,:,this_FOI_idx))-abs(GA{2}.powspctrm(:,:,this_FOI_idx)),3);

end

%% save
if ~exist(savedir)
    mkdir(savedir)
end

save(fullfile(savedir,savename),'dat', 'settings_', 'freqOI')
