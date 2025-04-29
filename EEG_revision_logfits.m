filedir = '/Users/ryszard/Downloads/etad_files';
cd(filedir)

spmdir = '/Users/ryszard/Documents/spm12_mac/';
addpath(spmdir)
spm('Defaults','EEG')

addpath('/Users/ryszard/Downloads/hannah_elife_matlab_scripts/analysis')

%% some settings:

% some trials show strong linear trends (e.g. increasing amplitude over
% time) - this can remove these trends
do_detrend = 0; % detrend single-trial data? 1: yes, 0: no

% here we decide if we do decoding based on single EEG channels, or on
% principal components grouping several channels together
do_pca = 1;     % 1: PCA over channels; 0: original channels

% here we decide if we use all available channels/components or if we select
% only a subset for decoding
do_selchan = 1; % 1: select channels/components based on SNR; 0: analyse all channels/components; 2: based on F stat (how strongly each channel/component differentiates between faces, houses and chairs)

% if do_selchan == 1, we can set an SNR threshold (cut off threshold)
snr_db = 8;    % SNR threshold for channel selection (dB)

% if do_selchan == 2, we can set a number of most sensitive channels
fstat_nchan = 5; % number of channels selected based on F stat

pps = setdiff(1:31,20) ;

for s= 1:length(pps) % select participants of interest
    try
        % go to their folders
        if pps(s)<10
            datadir = strcat('/Users/ryszard/Downloads/etad_files/pp0',num2str(pps(s)));
        else
            datadir = strcat('/Users/ryszard/Downloads/etad_files/pp',num2str(pps(s)));
        end
        cd(datadir)


        pps(s)

        clear accuracy_leading* accuracy_trailing* % clean up variables for saving later

        for erp = 2 % 1: leading, 2: trailing

            temp = dir('eTad*.mat'); % single-trial EEG files to be loaded

            %% load  data

            D = spm_eeg_load(temp(1).name); % load the file
            bdtrls = D.badtrials;

            data = D(:,:,:);
            data(D.badchannels,:,:) = NaN; % replace bad channels
            data(indchantype(D,'Other'),:,:) = NaN; % replace non-EEG channels (ECG, EOG)

            % in the EEG files, stimulus labels are:

            orig_labels = D.conditions';
            unique_labels = sort(D.condlist)';

            stim = orig_labels;

            % due to merged blocks etc. it is possible that two
            % consecutive images are trailing (or leading). let's
            % delete them. this is only a problem for s=25, pps(s)=26

            double_trailing = [];
            for i=2:length(stim)
                if length(find(strfind(stim{i},'trailing')))>0 & length(find(strfind(stim{i-1},'trailing')))>0
                    double_trailing = [double_trailing i];
                end
            end

            double_leading = [];
            for i=2:length(stim)
                if length(find(strfind(stim{i},'leading')))>0 & length(find(strfind(stim{i-1},'leading')))>0
                    double_leading = [double_leading i];
                end
            end

            if length(double_trailing)>0 | length(double_leading)>0 %% s=25, pps(s)=26
                double_trials = [double_leading double_trailing];
                stim(double_trials) = [];
                data(:,:,double_trials) = [];
                bdtrls = setdiff(bdtrls, double_trials);
                bdtrls(find(bdtrls>=min(double_trials))) = bdtrls(find(bdtrls>=min(double_trials))) - length(double_trials);
            end

            % due to battery issues etc. it is possible that the entire
            % recording starts with a trailing image or ends with a leading
            % image

            if length(find(strfind(stim{1},'trailing')))>0 % if for some reason it still ends with leading
                stim(1) = [];
                data(:,:,1) = [];
            end
            if length(find(strfind(stim{end},'leading')))>0 % if for some reason it still ends with leading
                stim(end) = [];
                data(:,:,end) = [];
            end

            % now make sure to exclude not just single bad trials but
            % entire pairs of leading and trailing images (so e.g. if a
            % leading image is bad, it will also mark the consecutive
            % trailing image as bad; and vice versa)

            bdtrls_leading = bdtrls(find(rem(bdtrls,2)==1)); % bad trials that are leading (i.e. odd)
            bdtrls_trailing = bdtrls(find(rem(bdtrls,2)==0)); % bad trials that are trailing (i.e. even)

            bdtrls = [bdtrls_leading bdtrls_leading+1 bdtrls_trailing bdtrls_trailing-1]; % exclude matching leading/trailing trials

            bdtrls(find(bdtrls<1)) = []; % in case the first trial is trailing and bad, remove it from the list
            bdtrls(find(bdtrls>size(data,3))) = []; % in case the last trial is leading and bad, remove it from the list

            bdtrls = unique(bdtrls);

            data(:,:,bdtrls) = NaN; % get rid of bad trials

            if do_detrend==1 % if you detrend, this will simply remove the linear trend
                for j=1:size(D,1)
                    j
                    to_detrend = squeeze(data(j,:,:));
                    data(j,:,:)=detrend(to_detrend',1,'omitnan')';
                end
            end

            data(find(isnan(mean(nanmean(data,3),2))),:,:) = []; % remove bad channels
            stim(find(isnan(mean(nanmean(data,2),1)))) = [];
            data(:,:,find(isnan(mean(nanmean(data,2),1)))) = []; % remove bad trials



            if erp == 1 % leading
                pick_trials = find(cellfun(@numel, strfind(stim,'leading'))>0);
            else % trailing
                pick_trials = find(cellfun(@numel, strfind(stim,'trailing'))>0);
            end

            data = data(:,:,pick_trials);
            leading_labels = stim(find(cellfun(@numel, strfind(stim,'leading'))>0));
            trailing_labels = stim(find(cellfun(@numel, strfind(stim,'trailing'))>0));

            if do_pca == 1 % replace original channels with principal (temporal) components explaining 99% variance
                [U,S,V] = svd(reshape(data,[size(data,1) size(data,2)*size(data,3)]),'econ');
                no_comp = find(cumsum(diag(S).^2/sum(diag(S).^2))<.99, 1, 'last' ); % find those components that, taken together, explain 99% variance
                data = reshape(V(:,1:no_comp)',[no_comp size(data,2) size(data,3)]); % replace original data with principal components
            end

            if do_selchan == 1 % select channels/components with SNR > threshold
                if do_pca == 0
                    data_erp = nanmean(data,3);
                    snr_perchannel = 10*log10((rms(data_erp(:,indsample(D,0.05):indsample(D,0.15)),2)./rms(data_erp(:,1:indsample(D,0)),2)).^2);
                    selcomps = find(snr_perchannel > snr_db);
                    if length(selcomps) == 0
                        selcomps = find(snr_perchannel > 3);
                    end
                    data = data(selcomps,:,:);
                else
                    snr_perchannel = 10*log10((rms(nanmean(data(:,indsample(D,0.05):indsample(D,0.15),:),3),2)./rms(nanmean(data(:,1:indsample(D,0),:),3),2)).^2);
                    selcomps = find(snr_perchannel > snr_db);
                    if length(selcomps) < 2 % workaround for participant 7, trailing analysis
                        selcomps = find(snr_perchannel > snr_db/2);
                    end
                    data = data(selcomps,:,:);
                end
            end

            if do_selchan == 2 % select channels/components based on F statistic (differences between diff types of tones)
                if do_pca == 0
                    data_erp = cat(3,D(:,:,setdiff(indtrial(D,'1'),D.badtrials)),D(:,:,setdiff(indtrial(D,'2'),D.badtrials)),D(:,:,setdiff(indtrial(D,'3'),D.badtrials)));
                    data_erp([D.badchannels indchantype(D,'Other')],:,:) = [];
                    stimlabels = [ones(length(setdiff(indtrial(D,'1'),D.badtrials)),1)*1; ones(length(setdiff(indtrial(D,'2'),D.badtrials)),1)*2; ones(length(setdiff(indtrial(D,'3'),D.badtrials)),1)*3];
                    for c=1:size(data_erp,1)
                        for t=1:size(data_erp,2)
                            % run an ANOVA
                            [p,anovatab]=anova1(squeeze(data_erp(c,t,:)),stimlabels,'off');
                            fs(c,t)=cell2mat(anovatab(2,5));
                        end
                    end
                    fs=mean(fs,2);
                    % collect all F values
                    fs=sortrows([fs,[1:length(fs)]'],'descend');
                    % select top channels
                    selcomps = fs(1:fstat_nchan,2);
                    data = data(selcomps,:,:);
                else
                    error('option not implemented yet')
                end
            end

            clear trndat testdat covdat distance meandistance

            legal_pairs = {{'leading_Barn' 'trailing_church'} ... % valid 75%
                {'leading_Barn' 'trailing_conference_room'} ... % invalid 25%
                {'leading_beach' 'trailing_church'} ... % valid 75%
                {'leading_beach' 'trailing_conference_room'} ... % invalid 25%
                {'leading_library' 'trailing_conference_room'} ... % valid 75%
                {'leading_library' 'trailing_church'} ... % invalid 25%
                {'leading_restaurant' 'trailing_conference_room'} ... % valid 75%
                {'leading_restaurant' 'trailing_church'} ... % invalid 25%
                {'leading_cave' 'trailing_castle'} ... % control 50%
                {'leading_cave' 'trailing_forest'}}; % control 50%


            %% do decoding

            trial_count = [];
            for j = 1:length(legal_pairs)
                trial_count(j) = length(intersect(find(strcmp(leading_labels,legal_pairs{j}{1})),find(strcmp(trailing_labels,legal_pairs{j}{2}))));
            end

            subsample_trials = min(trial_count); 
            subsample_trials = 47; % same for everyone
            no_samples(erp) = subsample_trials;

            for j = 1:length(legal_pairs)
                % temptrials = intersect(find(strcmp(leading_labels,legal_pairs{j}{1})),find(strcmp(trailing_labels,legal_pairs{j}{2})));
                % my_chosen_trials{j} = temptrials(1:subsample_trials);
                my_chosen_trials{j} = sort(randsample(intersect(find(strcmp(leading_labels,legal_pairs{j}{1})),find(strcmp(trailing_labels,legal_pairs{j}{2}))),subsample_trials));
                trndat{j} = data(:,:,my_chosen_trials{j}); % select the remaining trials as "train data" and average per feature and stimulus label across trials
            end

            strls = NaN(size(data,2),size(data,3));

            for k = 1:size(data,2) % per time point (no sliding)

                if erp==2 % trailing image analysis

                    if rem(k,5) == 0
                        display(strcat(['decoding trailing trials, finished ',num2str(round(100*k/size(data,2))), '%']))
                    end

                    % visual category decoding
                    X = [];
                    y = [];
                    for c = 1:8
                        X = [X; squeeze(trndat{c}(:,k,:))']; % trials x features
                        y = [y; ones(subsample_trials,1)*length(find(strfind(legal_pairs{c}{2},'church')))]; % binary labels: church vs. conference room
                    end
                    testIndices = repmat([1:length(y)/8]',[1 8]);
                    c = cvpartition("CustomPartition",testIndices(:));
                    svmModel = fitcsvm(X, y, 'CrossVal', 'on', 'Leaveout', 'on');
                    accuracy_trailing_visual(k) = mean(svmModel.kfoldPredict == y);
                    accuracy_trailing_visual_decoderoutput(k,:) = svmModel.kfoldPredict == y;
                    
                    strl_lab = reshape(y,[subsample_trials 8]);
                    strl_pred = reshape(svmModel.kfoldPredict,[subsample_trials 8]);
                    strl_corr = strl_lab == strl_pred;

                    strl_valid = strl_corr(:,1:2:8)-strl_corr(:,2:2:8);

                    conds = 1:2:8;

                    for c=1:length(conds)
                        strls(k,my_chosen_trials{conds(c)}) = strl_valid(:,c);
                    end

                  
                end

            end


        end

        strls = strls(:,find(~isnan(strls(1,:))));


        save decoding_nodetrend_pca_snr8db_svm.mat strls

        % end
    catch
    end
end

%% pool data
strls_all = nan(length(pps),180,188);

for s= 1:length(pps) % select participants of interest
    
    % go to their folders
    if pps(s)<10
        datadir = strcat('/Users/ryszard/Downloads/etad_files/pp0',num2str(pps(s)));
    else
        datadir = strcat('/Users/ryszard/Downloads/etad_files/pp',num2str(pps(s)));
    end
    cd(datadir)
    load decoding_nodetrend_pca_snr8db_svm.mat
    strls_all(s,:,1:length(strls)) = strls;
end

%% fit trial-by-trial time series

twin1 = [123 180]; % first significant time window (ms)
twin2 = [280 296]; % second time window

% convert ms to samples/indices
tind1 = [min(find(D.time>=twin1(1)/1000)) max(find(D.time<=twin1(2)/1000))]; 
tind2 = [min(find(D.time>=twin2(1)/1000)) max(find(D.time<=twin2(2)/1000))];

% extract mean decoding accuracy (valid minus invalid) within each time window 
tser1 = squeeze(nanmean(strls_all(:,tind1(1):tind1(2),:),2));
tser2 = squeeze(nanmean(strls_all(:,tind2(1):tind2(2),:),2));

% define fit type (e.g. logarithmic)
rng(12345)
ft_model = fittype('A*log(B*x) + C', 'independent', 'x'); % logarithmic
startPoints = [0 1 0]; % starting points for A, B, C
% ft_model = fittype('A*exp(-B*x) + C', 'independent', 'x'); % exponential
% startPoints = [0 .01 0]; % starting points for A, B, C
smoothf = 5; % smooth data over N trials (for fits)
toler = .001; % tolerance of derivative over trials to determine when decoding plateaus

% fit per participant
log_fit1 = tser1*0;
log_fit2 = tser2*0;

for i=1:size(tser1,1)
    % early time window
    tempfit = fit([1:size(tser1,2)]', smooth(tser1(i,:),smoothf), ft_model, 'StartPoint', startPoints);
    log_fit1(i,:) = tempfit(1:size(tser1,2));
    d_exp = differentiate(tempfit, 1:size(tser1,2)); % First derivative of the exponential fit
    try
        plateau_idx_exp1(i) = find(abs(d_exp) < toler, 1);  % Tolerance can be adjusted
    catch
        plateau_idx_exp1(i) = NaN;
    end

    y = smooth(tser1(i,:),smoothf);
    y_fit = tempfit(1:size(tser1,2));
    residuals = y - y_fit;
    SST = sum((y - mean(y)).^2);
    SSE = sum(residuals.^2);
    R_squared1(i) = 1 - (SSE / SST);

    % late time window
    tempfit = fit([1:size(tser2,2)]', smooth(tser2(i,:),smoothf), ft_model, 'StartPoint', startPoints);
    log_fit2(i,:) = tempfit(1:size(tser2,2));
    d_exp = differentiate(tempfit, 1:size(tser2,2)); % First derivative of the exponential fit
    try
        plateau_idx_exp2(i) = find(abs(d_exp) < toler, 1);  % Tolerance can be adjusted
    catch
        plateau_idx_exp2(i) = NaN;
    end

    y = smooth(tser2(i,:),smoothf);
    y_fit = tempfit(1:size(tser2,2));
    residuals = y - y_fit;
    SST = sum((y - mean(y)).^2);
    SSE = sum(residuals.^2);
    R_squared2(i) = 1 - (SSE / SST);
 
end

figure;

subplot(2,4,1)
smooth_tser1 = tser1*0;
for i=1:30
    smooth_tser1(i,:) = smooth(tser1(i,:),smoothf);
end
shadedErrorBar(1:size(tser1,2),mean(smooth_tser1,1),std(smooth_tser1,[],1)/sqrt(30))
xlabel('trials')
ylabel('rel. decoding acc. (val-inv)')
title('data, early time window')
ylim([-.2 .2])
line([1 size(tser1,2)],[0 0],'linestyle','--','color','k')

subplot(2,4,5)
smooth_tser2 = tser2*0;
for i=1:30
    smooth_tser2(i,:) = smooth(tser2(i,:),smoothf);
end
shadedErrorBar(1:size(tser2,2),mean(smooth_tser2,1),std(smooth_tser2,[],1)/sqrt(30))
xlabel('trials')
ylabel('rel. decoding acc. (val-inv)')
title('data, late time window')
ylim([-.2 .2])
line([1 size(tser1,2)],[0 0],'linestyle','--','color','k')

subplot(2,4,2)
hold on
means = squeeze(mean(mean(reshape(tser1,[size(tser1,1) size(tser1,2)/4 4]),2),1));
sems = squeeze(std(mean(reshape(tser1,[size(tser1,1) size(tser1,2)/4 4]),2),[],1))/sqrt(size(tser1,1));
bar(means)
errorbar(means,sems,'color','blue','linestyle','none')
xticks(1:4)
xlabel('bin')
ylabel('rel. decoding acc. (val-inv)')
title('early time window')
ylim([-.06 .06])
line([0 5],[0 0],'linestyle','--','color','k')

subplot(2,4,6)
hold on
means = squeeze(mean(mean(reshape(tser2,[size(tser2,1) size(tser2,2)/4 4]),2),1));
sems = squeeze(std(mean(reshape(tser2,[size(tser2,1) size(tser2,2)/4 4]),2),[],1))/sqrt(size(tser2,1));
bar(means)
errorbar(means,sems,'color','blue','linestyle','none')
xticks(1:4)
xlabel('bin')
ylabel('rel. decoding acc. (val-inv)')
title('late time window')
ylim([-.06 .06])
line([0 5],[0 0],'linestyle','--','color','k')

subplot(2,4,3)
shadedErrorBar(1:size(log_fit1,2),mean(log_fit1,1),std(log_fit1,[],1)/sqrt(30))
xlabel('trials')
ylabel('rel. decoding acc. (val-inv)')
title('fit, early time window')
ylim([-.06 .06])
line([1 size(tser1,2)],[0 0],'linestyle','--','color','k')

subplot(2,4,7)
shadedErrorBar(1:size(log_fit2,2),mean(log_fit2,1),std(log_fit2,[],1)/sqrt(30))
xlabel('trials')
ylabel('rel. decoding acc. (val-inv)')
title('fit, late time window')
ylim([-.06 .06])
line([1 size(tser1,2)],[0 0],'linestyle','--','color','k')

subplot(2,4,4)
hold on
nanm = nanmedian(plateau_idx_exp1/size(tser1,2));
violinplot(plateau_idx_exp1/size(tser1,2))
line([0 2],[nanm nanm],'linestyle','-','color','k')
title('fit, early, plateau')
ylabel('trials')
ylim([0 1])
prctile(plateau_idx_exp1/size(tser1,2),[1 50 99])

subplot(2,4,8)
hold on
nanm = nanmedian(plateau_idx_exp2/size(tser1,2));
violinplot(plateau_idx_exp2/size(tser1,2))
line([0 2],[nanm nanm],'linestyle','-','color','k')
title('fit, late, plateau')
ylabel('trials')
ylim([0 1])
prctile(plateau_idx_exp2/size(tser1,2),[1 50 99])