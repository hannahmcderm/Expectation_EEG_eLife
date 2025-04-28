%%% BEHAVIOURAL TASK TO DETERMINE IF LEGAL PAIRS FROM ES_EEG HAVE BEEN
%%% LEARNED
%%% Hannah McDermott, NOVEMBER 2023 %%%


function [dt, savethis, logfile] = ES_MEM(ID) % ID should be a number

    datapath = 'C:\Users\CCNB-EEG1\Desktop\Hannah\ES_EEG\experiment_logs';
    mkdir(datapath)

    temp = dir(fullfile(datapath,strcat('ID_',num2str(ID),'*')));
    if length(temp)>0
        error('ID already taken')
    end

    clc;
    Screen('Preference', 'Verbosity', 0);
    Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference', 'VisualDebugLevel', 0);

    ISI = .80; % between leading and trailing images
    pic = .10; % presentation time of each image
    ITI = 2.0; % in between trials

    % define colours

    gray = [200 200 200];
    white = [255 255 255];
    black = [0 0 0];
    fontsize = 25;

    % files directories of stimuli
    categoryFolders = {
        'Cat1', 'C:\Users\CCNB-EEG1\Desktop\Hannah\ES_EEG\barn';
        'Cat2', 'C:\Users\CCNB-EEG1\Desktop\Hannah\ES_EEG\beach';
        'Cat3', 'C:\Users\CCNB-EEG1\Desktop\Hannah\ES_EEG\cave';
        'Cat4', 'C:\Users\CCNB-EEG1\Desktop\Hannah\ES_EEG\library';
        'Cat5', 'C:\Users\CCNB-EEG1\Desktop\Hannah\ES_EEG\restaurant';
        'Cat6', 'C:\Users\CCNB-EEG1\Desktop\Hannah\ES_EEG\Church';
        'Cat7', 'C:\Users\CCNB-EEG1\Desktop\Hannah\ES_EEG\conference_room';
        'Cat8', 'C:\Users\CCNB-EEG1\Desktop\Hannah\ES_EEG\castle';
        'Cat9', 'C:\Users\CCNB-EEG1\Desktop\Hannah\ES_EEG\forest'
        };

    for i=1:size(categoryFolders,1)
        addpath(categoryFolders{i,2});
    end

    numCategories = size(categoryFolders,1);
    imagesByCategory = cell(1, numCategories);

    % get a list of all images
    for categ = 1:numCategories
        filelist = dir(fullfile(categoryFolders{categ,2},'*.*'));
        filelist(1:2) = []; % empty paths
        imagesByCategory{categ} = {filelist.name};
    end

    % create the image pairings

    % Map 1 = C1 LEADING >> C6 TRAILING valid, C7 invalid
    % Map 2 = C2 LEADING >> C6 TRAILING valid, C7 invalid
    % Map 3 = C4 LEADING >> C7 TRAILING valid, C6 invalid
    % Map 4 = C5 LEADING >> C7 TRAILING valid, C6 invalid
    % Map 5 = C3 LEADING >> C8 OR C9 TRAILING
    % Map 6 = C3 LEADING >> C9 OR C8 TRAILING

    leading_maps = [1 2 3 4 5 6];
    leading_categories = [1 2 4 5 3 3];
    leading_probs = [1 1 1 1 1 1];

    trailing_maps = [1 1 2 2 3 3 4 4 5 5 6 6];
    trailing_categories = [6 7 6 7 7 6 7 6 8 9 8 9];
    trailing_probs = [.75 .25 .75 .25 .75 .25 .75 .25 .5 .5 .5 .5];

    no_blocks = 1;
    no_trials_per_block = 24;
    total_no_trials = no_blocks * no_trials_per_block;

    leading_image_list = cell(0);
    trailing_image_list = cell(0);
    for b = 1:no_blocks
        for i = 1:length(leading_maps)
            how_many_images = leading_probs(i)*no_trials_per_block / 6; % divided by 4 possible maps

            which_images = randsample(1:length(imagesByCategory{leading_categories(i)}),how_many_images);
            leading_image_list = [leading_image_list imagesByCategory{leading_categories(i)}(which_images)];
            imagesByCategory{leading_categories(i)}(which_images) = [];
        end

        for i = 1:length(trailing_maps)
            how_many_images = trailing_probs(i)*no_trials_per_block / 6; % divided by 4 possible maps

            which_images = randsample(1:length(imagesByCategory{trailing_categories(i)}),how_many_images);
            trailing_image_list = [trailing_image_list imagesByCategory{trailing_categories(i)}(which_images)];
            imagesByCategory{trailing_categories(i)}(which_images) = [];
        end
    end

    all_trials_reshuffled = Shuffle(1:length(leading_image_list));

    % shuffle the order per block
    for b = 1:no_blocks
        reshuffled_order = all_trials_reshuffled((b-1)*no_trials_per_block+1 : b*no_trials_per_block);
        leading_image_list((b-1)*no_trials_per_block+1 : b*no_trials_per_block) = leading_image_list(reshuffled_order);
        trailing_image_list((b-1)*no_trials_per_block+1 : b*no_trials_per_block) = trailing_image_list(reshuffled_order);
    end

    KbName('UnifyKeyNames');
    kKey = KbName('k');
    sKey = KbName('s');

    [mw, rect]=Screen('OpenWindow', 0, gray,[0 0 800 600]);

    centerX=rect(3)/2;
    centerY=rect(4)/2;

    crossLines=[-50, 50, 0 , 0
        0, 0, -50, 50];

    HideCursor();
    Screen('TextSize', mw, fontsize);
    Screen('TextFont', mw, 'Times');
    myText = 'In the following block, you will see the same pairs of images as before. You should indicate via button press whether the second image in each pair is indoors or outdoors. \n\n If the image is indoors, please press "S". \n\n If the image is outdoors please press "K". \n\n When you are ready, press space to begin.';
    DrawFormattedText(mw, myText, 'center', 'center', white);
    Screen('Flip', mw);
    while true
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown && keyCode(KbName('Space'))
            break;
        end
    end

    Screen('Flip', mw);
    WaitSecs(3);

    k = 1:no_blocks;
    savethis=[];
    logfile = [];
    ListenChar(2);

    for i = 1:no_trials_per_block
        leading_image = imread(leading_image_list{(k-1)*no_trials_per_block+i});
        trailing_image = imread(trailing_image_list{(k-1)*no_trials_per_block+i});

        % INTERPOLATING IMAGES TO MAKE THEM BIGGER

        leading_image = imresize(leading_image,2);
        trailing_image = imresize(trailing_image,2);

        % SHOW LEADING IMAGE

        TexturePointer = Screen(mw, 'MakeTexture', leading_image);
        Screen(mw, 'DrawTexture', TexturePointer);

        Screen('DrawLines', mw, crossLines, 5, white, [centerX centerY]);
        Screen('Flip', mw);

        iti_start = GetSecs;
        logfile.stimulustimes(i) = iti_start;
        tnow = GetSecs;

        responseTime_s = NaN;
        responseTime_k = NaN;

        while tnow < iti_start+pic
            tnow = GetSecs;
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(KbName('s')) == 1
                % Save response time for 's'
                responseTime_s = secs;
                KbReleaseWait;
            elseif keyCode(KbName('k')) == 1
                % Save response time for 'k'
                responseTime_k = secs;
                KbReleaseWait;
            end
        end

        Screen('DrawLines', mw, crossLines, 5, white,[centerX centerY]);
        Screen('Flip', mw);

        while tnow < iti_start+ISI
            tnow = GetSecs;
            [keyIsDown, secs, keyCode]=KbCheck;
            if keyCode(KbName('s')) == 1
                % Save response time for 's'
                responseTime_s = secs;
                KbReleaseWait;
            elseif keyCode(KbName('k')) == 1
                % Save response time for 'k'
                responseTime_k = secs;
                KbReleaseWait;
            end
        end

        % SHOW TRAILING IMAGE

        TexturePointer = Screen(mw, 'MakeTexture', trailing_image);

        Screen(mw, 'DrawTexture', TexturePointer);
        Screen('DrawLines', mw, crossLines, 5, white,[centerX centerY]);
        Screen('Flip', mw);

        iti_start = GetSecs;
        tnow = GetSecs;

        while tnow < iti_start+pic
            tnow = GetSecs;
            [keyIsDown, secs, keyCode]=KbCheck;
            if keyCode(KbName('s')) == 1
                % Save response time for 's'
                responseTime_s = secs;
                KbReleaseWait;
            elseif keyCode(KbName('k')) == 1
                % Save response time for 'k'
                responseTime_k = secs;
                KbReleaseWait;
            end
        end

        Screen('DrawLines', mw, crossLines, 5, white,[centerX centerY]);
        Screen('Flip', mw);

        iti_start = GetSecs;
        tnow = GetSecs;

        while tnow < iti_start+ITI
            tnow = GetSecs;
            [keyIsDown, secs, keyCode]=KbCheck;
            if keyCode(KbName('s')) == 1
                % Save response time for 's'
                responseTime_s = secs;
                KbReleaseWait;
            elseif keyCode(KbName('k')) == 1
                % Save response time for 'k'
                responseTime_k = secs;
                KbReleaseWait;
            end
        end

        savethis = [savethis; responseTime_s, responseTime_k];

        dt=datestr(now,'dd.mm.YY_HH:MM:SS');
        name=['ID_' num2str(ID)];
        logfile.name = name;
        logfile.ID = ID;
        logfile.responses = savethis();
        logfile.leading_image_list = leading_image_list((k-1)*no_trials_per_block+1 : k*no_trials_per_block);
        logfile.trailing_image_list = trailing_image_list((k-1)*no_trials_per_block+1 : k*no_trials_per_block);

        save(fullfile(datapath,name),'logfile')

    end

    % end of experiment
    Screen('TextSize', mw, fontsize);
    [normBoundsRect, ~] = Screen('TextBounds', mw, 'The experiment has ended. Thank you.');
    Screen('DrawText', mw, 'The experiment has ended. Thank you.', centerX-normBoundsRect(3)/2, centerY-normBoundsRect(4)/2, white);
    Screen('Flip', mw);
    WaitSecs(5);
    ShowCursor;
    ListenChar(1);
    Screen('CloseAll')

end