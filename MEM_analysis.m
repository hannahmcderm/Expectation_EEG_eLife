exgaussdir = 'C:\Users\admin\Documents\MATLAB\bramzandbelt-exgauss-7cf8c4c';
addpath(exgaussdir)

clear all

pps = setdiff(9:31,[12 20]);

for s = 1:length(pps)
    datadir = 'D:\ES_data_Hannah\data';
    cd(datadir)
    if pps(s)<10
        cd(strcat('P0',num2str(pps(s))));
    else
        cd(strcat('P',num2str(pps(s))));
    end

    cd MEM
    load('ID_B.mat');

    valid_pairs = {'barn' 'church'
        'beach' 'church'
        'library' 'conference_room'
        'restaurant' 'conference_room'};

    invalid_pairs = {'barn' 'conference_room'
        'beach' 'conference_room'
        'library' 'church'
        'restaurant' 'church'};

    paircode = [];
    for i=1:length(logfile.leading_image_list)
        valid_invalid = 0;
        for j=1:size(valid_pairs,1)
            if length(find(strfind(logfile.leading_image_list{i},valid_pairs{j,1}))) && ...
                length(find(strfind(logfile.trailing_image_list{i},valid_pairs{j,2}))) 
                valid_invalid = 1;
            end
        end

        for j=1:size(invalid_pairs,1)
            if length(find(strfind(logfile.leading_image_list{i},invalid_pairs{j,1}))) && ...
                length(find(strfind(logfile.trailing_image_list{i},invalid_pairs{j,2}))) 
                valid_invalid = -1;
            end
        end
        paircode = [paircode valid_invalid];
    end

    logfile.responses = logfile.responses - logfile.stimulustimes' - .8;
    logfile.responses = nansum(logfile.responses')'; % getting rid of indoor vs outdoor columns


    valid_RTs =  logfile.responses(find(paircode == 1));
    invalid_RTs =  logfile.responses(find(paircode == -1));
    neutral_RTs =  logfile.responses(find(paircode == 0));

    valid_RTs(find(valid_RTs>median(valid_RTs)+2*std(valid_RTs) | valid_RTs<median(valid_RTs)-2*std(valid_RTs)))=[]
    invalid_RTs(find(invalid_RTs>median(invalid_RTs)+2*std(invalid_RTs) | invalid_RTs<median(invalid_RTs)-2*std(invalid_RTs)))=[]
    neutral_RTs(find(neutral_RTs>median(neutral_RTs)+2*std(neutral_RTs) | neutral_RTs<median(neutral_RTs)-2*std(neutral_RTs)))=[]

    rts(s,1) = nanmean(log(neutral_RTs));
    rts(s,2) = nanmean(log(valid_RTs));
    rts(s,3) = nanmean(log(invalid_RTs));

end

rts(isinf(rts)) = NaN;

% repeated measures ANOVA
subjects = length(pps);
conditions = {'Neutral', 'Valid', 'Invalid'};
tbl = array2table(rts, 'VariableNames', conditions);
tbl.Subject = (1:subjects)';
Meas = table(conditions', 'VariableNames', {'Condition'});
rm = fitrm(tbl, 'Neutral-Invalid ~ 1', 'WithinDesign', Meas);
ranova_results = ranova(rm)
p_anova = ranova_results.pValue(1)

% eta squared for RM ANOVA
SS_condition = ranova_results.SumSq(1); % SS for Condition
SS_error = ranova_results.SumSq(2);     % SS for Error
SS_total = SS_condition + SS_error;     % Total SS
eta_squared = SS_condition / SS_total;  

% paired t-tests 
[h,p,ci,stats] = ttest(rts(:,1),rts(:,2))
p
stats.tstat

diff_scores = rts(:,1) - rts(:,2);
mean_diff = nanmean(diff_scores);
std_diff = nanstd(diff_scores);
effect_size = mean_diff / std_diff

rts = exp(rts);
validity_effect_beh = rts(:,2) - rts(:,3);
meanRTs = nanmean(rts,1);

% plot within-subjects SEM instead of a regular SEM
rts_for_sem = rts;
for i=1:size(rts_for_sem,1) % per person
    rts_for_sem(i,:) = rts_for_sem(i,:) - nanmean(rts_for_sem(i,:)); % remove that person's average across all 3 conditions
    rts_for_sem(i,:) = rts_for_sem(i,:) + nanmean(rts_for_sem,1); % add the grand average across all 3 conditions
end

stdRTs = nanstd(rts_for_sem,[],1);
semRTs = stdRTs / sqrt(size(rts_for_sem,1));

colours = ['y'; 'b'; 'r'];


figure;
violin(rts, 'xlabel', {'Neutral', 'Valid', 'Invalid'}, 'facecolor', colours, 'edgecolor', 'k', 'facealpha', 0.5, 'mc', 'k', 'medc', []);
hold on
errorbar(1:3, meanRTs, semRTs, 'k.', 'LineWidth', 1.5);
ylabel('Reaction Time (s)');
ylim([0.0 1.2]);
title('Reaction Times across Conditions');
hold off
