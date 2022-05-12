%% enclseep. overlap between encoding and sleep topographies & their prediction for memory consolidation
% this script (1) creates Figure 2b & 3a & 3c

clear
close all

github_repo = '/Users/petzka/Documents/GitHub/encsleep_'; %directory of github repo where the scripts live
addpath(genpath(github_repo)) 

datfolder = '/Volumes/MEMTOSH/encsleep'; %directory of data folder

%% where to get
wakedatdir = fullfile(datfolder, 'wake', 'encpvt_topo');
wakedatOI = 'encpvt_topo';

sleepdatdir = fullfile(datfolder, 'sleep', 'detect_topo');
sleepdatOI = 'detect_topo';

behavdir = fullfile(datfolder, 'behaviour');
behavdatOI = 'behav_var';

corr_type = 'Spearman';

%% load data
dat_wake = load([wakedatdir, filesep, wakedatOI]);
dat_sleep = load([sleepdatdir, filesep, sleepdatOI]);
dat_behav = load([behavdir, filesep, behavdatOI]);

%% subjects         
nsub = dat_sleep.settings_.nsub;

%% selection of sleep events&properties wake data
dat_sleep.evtOI 
dat_wake.freqOI

dat_sel = {'fast_spin' 'meanEvtMaxAmp';...
    'fast_spin' 'meanEvtLen';....
    'fast_spin' 'density';...
    'SOs' 'meanEvtMinAmp';...
    'SOs' 'meanEvtDuration';...
    'SOs' 'density';...
    'fft' '6 20'}; % 1st column = events, 2nd colum = properties

%% behavioural data 
dat_behav = dat_behav.out;
behav_sel = {'seq', 'corr'}; 

%% get sleep & wake data of interest
ntopos = size(dat_sel,1); 
sleep_evtOI = fieldnames(dat_sleep.evtOI);
wake_evtOI = fieldnames(dat_wake.freqOI);

for idat = 1:ntopos
    
    % get indices for data
    idx_sleep01 = find(ismember(sleep_evtOI, dat_sel{idat,1})); % index for event
    idx_wake01 = find(ismember(wake_evtOI, dat_sel{idat,1})); % index for type of fft
    
    % get sleep & wake data for corresponding events & properties
    if ~isempty(idx_sleep01)
        sleep_prop = dat_sleep.evtOI.(sleep_evtOI{idx_sleep01}); %possible sleep porperties of evt
        idx_sleep02 = find(ismember(sleep_prop, dat_sel{idat,2}));%get index of selected prop
        
        dat_tmp = dat_sleep.dat{idx_sleep01}{idx_sleep02};
        
        % take abs of min amplitude of SOs (to get same correlation direction as with spindles)
        if strcmpi('SOs', dat_sel{idat,1}) && ...
                strcmpi('meanEvtMinAmp',dat_sel{idat,2})
            dat_tmp = abs(dat_tmp);
        end
        
        label_tmp = [dat_sel{idat,1} dat_sel{idat,2}];
    elseif ~isempty(idx_wake01)
        wake_prop = dat_wake.freqOI.(wake_evtOI{idx_wake01});
        idx_wake02 = find(ismember(wake_prop, dat_sel{idat,2}));
        
        dat_tmp = dat_wake.dat{idx_wake01}{idx_wake02};
        label_tmp = [dat_sel{idat,1} dat_sel{idat,2}];
    end
    
    datOI{idat} = dat_tmp;
    label_{idat} = label_tmp;
end

%% correlation matrix. overlap between encoding & sleep topo
corr_mat = [];

for i = 1:ntopos 
    for j = 1:ntopos
        for isub = 1:numel(nsub)
            corr_mat(i,j,isub) = corr(datOI{i}(isub,:)',datOI{j}(isub,:)',...
                'type', corr_type);
        end
    end
end

%% Fisher z
corr_mat_z = 0.5*(log(1+corr_mat) - log(1-corr_mat)); % same as atanh(corr_mat)

%% plot figure 2b

matOI = corr_mat_z;
datOI_ = [];
datOI_{1,1} = squeeze(matOI(1,7,:));
datOI_{1,2} = squeeze(matOI(2,7,:));
datOI_{1,3} = squeeze(matOI(3,7,:));
datOI_{2,1} = squeeze(matOI(4,7,:));
datOI_{2,2} = squeeze(matOI(5,7,:));
datOI_{2,3} = squeeze(matOI(6,7,:));

datOI_ = datOI_';
col_plot = [0, 0.4470, 0.7410;
    0.5 0.5 0.5];
figure
ax(1) = gca;
ax(1).XGrid = 'on';
ax(1).GridLineStyle = '-';
raincloud_edit(datOI_, col_plot)

%% behaviour. 
dat_tmp = dat_behav.(behav_sel{1,1}).(behav_sel{1,2});

dat_tmp(dat_tmp == 1) = NaN; % r = 1. atanh(1) = Inf
max_corr = max(dat_tmp(:));
replace_ = 1 - (1-max_corr)/2; % replacement for r = 1.
dat_tmp(isnan(dat_tmp)) = replace_;

dat_tmp = atanh(dat_tmp);

behav_var = dat_tmp(:,2)./(dat_tmp(:,1)./100);

%% select eeg variables to correlate with behaviour
var1_idx = [1:3]; % sleep spindles 
var2_idx = [7]; % wake

count_idx = 1;
eeg_vars = [];

for ivar1 = 1:numel(var1_idx)
    for ivar2 = 1:numel(var2_idx)
        eeg_vars(:,count_idx) = squeeze(matOI(var1_idx(ivar1),var2_idx(ivar2),:));
        eeg_vars_label{count_idx} = [label_{var1_idx(ivar1)} '.' label_{var2_idx(ivar2)}];
        count_idx = count_idx + 1;
    end
end

%% corr behav & eeg
all_vars = [behav_var eeg_vars];
behav_var_label = 'seq';

corr_mat02 = [];
pvals_mat02 =[];

for i = 1:size(all_vars,2)
   for j = 1:size(all_vars,2)
       [corr_mat02(i,j), pvals_mat02(i,j)] = corr(all_vars(:,i), all_vars(:,j),...
           'type', corr_type);   
   end
end

%% correlation matrix
% figure;
% imagesc(corr_mat02)
% set(gca,'xtick',1:size(all_vars,2),'xticklabel',...
%     [behav_var_label eeg_vars_label], 'TickLabelInterpreter','none')
% set(gca,'ytick',1:size(all_vars,2),'yticklabel',...
%     [behav_var_label eeg_vars_label])
% title(sprintf('%s correlations',corr_type))
% sig_mask = pvals_mat02<.05;
% colorbar
% caxis([-1 1])

%% correlation with behaviour
behav4corr = behav_var;

ivar_eeg = 1; % spindle amplitude
eeg4corr = eeg_vars(:,ivar_eeg);

new_eeg = eeg4corr;
new_behav = behav4corr;

[corr_val, p_val] =corr(behav4corr, eeg4corr, 'type', corr_type);
[corr_val, p_val] =corr(new_eeg, new_behav, 'type', corr_type);

%% plot figure 3a

X = new_eeg;
Y = new_behav;
figure
plot(eeg4corr, behav4corr, 'k.',...
    'MarkerSize',19); %lsline
box('off')

%% permutation test

d_sleep = dat_sleep.dat{1}{1};
d_wake = dat_wake.dat{1}{1};
d_behav = behav4corr;

%% corr original
corr_mat_orig = [];
for isub = 1:size(d_sleep,1)
    corr_mat_orig(isub) = corr(d_sleep(isub,:)',d_wake(isub,:)','type',corr_type);
end

[corr_val, p_val] =corr(d_behav, corr_mat_orig', 'type', corr_type);

fprintf('\n--- perm.test.original. r = %2.2f, p = %2.3f\n\n', corr_val, p_val)

%% shuffle encoding & sleep patterns & plot permutation distribution
% to plot permutation distribution for shuffled encoding patterns set in
% line 227. sel3 = sel1
% to plot permutation distribution for shuffled sleep spindle patterns set in
% line 227. sel3 = sel2

shams = 1000;
meanmat = nan(1,shams);
sigmat = nan(1,shams);
behavcorr = nan(1,shams);
p_vals = nan(1,shams);

cnt = 0;

while 1
    sel1 = randperm(numel(nsub)); %sleep
    sel2 = randperm(numel(nsub)); %task
    sel3 = sel2; %behav
    
   if not(any(sel1==sel2))
    
        cnt = cnt+1
        corr_mat_sh = [];
        for isub = 1:numel(nsub)
            corr_mat_sh(isub) = corr(d_sleep(sel1(isub),:)',...
                d_wake(sel2(isub),:)','type',corr_type);
        end

        meanmat(cnt) = mean(corr_mat_sh);
        corrmat_z = 0.5*(log(1+corr_mat_sh) - log(1-corr_mat_sh));
        sigmat(cnt)  = ttest(corrmat_z);
        [behavcorr(cnt), p_vals(cnt)] = corr(d_behav(sel3),...
            corr_mat_sh', 'type', corr_type);
        
   end
    if cnt==shams
        break
    end
end

fprintf('\n-----%3.0f%% iterations significant\n',100*mean(sigmat));

% plot figure 3C

figure
hist(behavcorr,100)
vline(corr_val, 'r')
title(sprintf('p = %0.3f',mean(abs(behavcorr) > abs(corr_val))))
% 
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = 'k';
box off

