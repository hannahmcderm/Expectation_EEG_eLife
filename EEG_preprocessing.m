spmdir = 'C:\Users\admin\Documents\MATLAB\spm12';
noisetooldir = 'C:\Users\admin\Documents\MATLAB\NoiseTools';
addpath(spmdir)
spm('Defaults','EEG')
addpath('C:\Users\admin\Documents\MATLAB\analysis')

datadir = 'D:\ES_data_Hannah\data';

for s = setdiff(1:31, 20);

  %  try

        if s<10
            subdir = fullfile(datadir,strcat('P0',num2str(s),'/'));
        else
            subdir = fullfile(datadir,strcat('P',num2str(s),'/'));
        end

        cd(subdir)
        cd EEG

        temp = dir('*.bdf');

        for i=1:length(temp)

            % convert
            clear S
            S.dataset = temp(i).name;
            S.mode='continuous';
            S.outfile=strcat('spmeeg_',num2str(i),'.mat');
            D = spm_eeg_convert(S);

            % filter
            clear S
            S.D=D.fname;
            S.band='high';
            S.freq=[.1];
            D = spm_eeg_filter(S);

            % notch
            clear S
            S.D=D.fname;
            S.band='stop';
            S.freq=[48 52];
            D = spm_eeg_filter(S);

            % downsample
            clear S
            S.D=D.fname;
            S.fsample_new = 300;
            D = spm_eeg_downsample(S);

            %% remove eye blink artefacts

            clear S
            S.D=D.fname;
            S.task='defaulteegsens';
            D=spm_eeg_prep(S);
            save(D);

            clear matlabbatch
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.D={D.fname};
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.val=1;
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.comment='';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.template=1;
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshres=2;
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).fidname='spmlpa';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).specification.select='lpa';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).fidname='spmnas';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).specification.select='nas';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).fidname='spmrpa';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).specification.select='rpa';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.useheadshape=1;
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.eeg='EEG BEM';
            matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.meg='Single Shell';
            spm_jobman('run',matlabbatch,cell(0,1));

            clear S
            S.D=D.fname;
            badchans = D.badchannels;
            which_chan='AF8';
            S.methods.channels=which_chan;
            S.methods.settings.chanind=indchannel(D,which_chan);
            S.methods.settings.threshold=3;
            S.methods.fun='eyeblink';
            S.mode='mark';
            S.append=1;
            S.methods.settings.excwin=400;
            S.badchanthresh=.2;
            S.prefix='a';
            D = spm_eeg_artefact(S);

            clear S
            S.D= D.fname;
            S.bc=0;
            S.timewin=[-200 400];
            S.trialdef(1).conditionlabel='eyeblink';
            S.trialdef(1).eventtype='artefact_eyeblink';
            S.trialdef(1).eventvalue=which_chan;
            S.reviewtrials=0;
            S.save=1;
            D = spm_eeg_epochs(S);

            clear S
            S.D=D.fname;
            S.method='SVD'; % singular value decomposition (same as principal component analysis)
            S.timewin=[-200 400];
            S.ncomp=2;
            D=spm_eeg_spatial_confounds(S);

            clear S
            fn = D.fname;
            S.D= fn(2:end);
            S.method='SPMEEG';
            S.conffile=D.fname;
            D=spm_eeg_spatial_confounds(S);

            clear S
            S.D=D.fname;
            S.correction='Berg';
            S.save=1;
            D=spm_eeg_correct_sensor_data(S);

            %% denoise (includes re-referencing and epoching)
           % try
                D = update_triggers_hannah(D,fullfile(subdir,'Behavioural'));

                addpath(noisetooldir)
                [y stimcode]=spm_denoise(D); % this includes re-referencing and epoching
                rmpath(noisetooldir)

                cd(fullfile(subdir,'EEG'))

                clear D2
                D2=clone(D,strcat('e',D.fname),size(y));
                D2(indchantype(D,'EEG'),:,:)=y;
                D2 = chantype(D2,1:size(D2,1),'EEG');
                D2 = chanlabels(D2,1:size(D2,1),D.chanlabels(indchantype(D,'EEG')));
                D2 = coor2D(D2,1:size(D2,1),D.coor2D(indchantype(D,'EEG')));
                D2 = units(D2,1:size(D2,1),D.units(indchantype(D,'EEG')));
                clear conditions
                D2=conditions(D2,1:length(D2.conditions),num2cell(stimcode(1:length(D2.conditions))));
                D2=timeonset(D2,-.1);
                save(D2)

                %% baseline correction
                clear S
                S.D=D2.fname;
                D=spm_eeg_load(S.D);
                D(:,:,:)=D(:,:,:)-repmat(mean(D(:,indsample(D,-.1):indsample(D,.0),:),2),[1 length(D.time) 1]);
                save(D);

                %% bad trials
                clear S
                std_per_trial=squeeze(mean(std(D(setdiff(indchantype(D,'EEG'),D.badchannels),:,:),[],2),1));
                badtrls=[find(std_per_trial>nanmedian(std_per_trial)+2*nanstd(std_per_trial));...
                    find(std_per_trial<nanmedian(std_per_trial)-2*nanstd(std_per_trial))];
                if length(badtrls)>0
                    D=badtrials(D,badtrls,ones(1,length(badtrls)));
                end
                save(D);

                %% average

                clear S
                S.D=D.fname;
                S.robust.savew=0;
                S.robust.byconditions=1;
                S.robust.ks=3;
                S.prefix='m';
                D = spm_eeg_average(S);

                %% low-pass filter

                clear S
                S.D=D.fname;
                S.band='low';
                S.freq=48;
                D = spm_eeg_filter(S);
           % catch
           % end

        end

        %% clean up files

        % delete spmeeg* p*spmeeg* fp*spmeeg* afp*spmeeg* Tafp*spmeeg* eafp*spmeeg*
        cd(subdir)
        cd EEG
        mkdir preprocessed
        movefile *eT* preprocessed/

    % catch
    % end


end