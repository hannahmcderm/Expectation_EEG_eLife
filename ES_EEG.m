% Expectation suppression study, Hannah McDermott

function [dt, blocks, savethis, logfile] = ES_EEG(ID) % ID should be a number

    datapath = 'C:\Users\admin\Documents\MATLAB\experiment_Logs';
    mkdir(datapath)

    %  temp = dir(fullfile(datapath,strcat('ID_',num2str(ID),'*')));
    %  if length(temp)>0
    %     error('ID already taken')
    % end

    clc;
    Screen('Preference', 'Verbosity', 0);
    Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference', 'VisualDebugLevel', 0);


    ISI = .80; % between leading and trailing images
    pic = .10; % presentation time of each image

    ITI_range = [1.3 2.2]; % we will sample uniformly from this ITI range (trailing to leading image)

    % define colours

    gray = [211 211 211];
    white = [255 255 255];
    black = [0 0 0];
    fontsize = 40;

    % EEG port

    %  ioObj = io64;
    % status = io64(ioObj)
    % address = hex2dec('3FF8');

    % sequence

    % condition 1: two-to-one mapping of 75% validity
    % condition 2: two-to-one mapping of 75% validity
    % condition 3: one-to-two mapping of 100% validity (50/50 split)

    % files directories of stimuli
    categoryFolders = {
        'Cat1', 'C:\Users\admin\Documents\MATLAB\Stimuli\forest';
        'Cat2', 'C:\Users\admin\Documents\MATLAB\Stimuli\beach';
        'Cat3', 'C:\Users\admin\Documents\MATLAB\Stimuli\conference_room';
        'Cat4', 'C:\Users\admin\Documents\MATLAB\Stimuli\library';
        'Cat5', 'C:\Users\admin\Documents\MATLAB\Stimuli\restaurant';
        'Cat6', 'C:\Users\admin\Documents\MATLAB\Stimuli\Church';
        'Cat7', 'C:\Users\admin\Documents\MATLAB\Stimuli\cave';
        'Cat8', 'C:\Users\admin\Documents\MATLAB\Stimuli\castle';
        'Cat9', 'C:\Users\admin\Documents\MATLAB\Stimuli\barn'
        };

    numCategories = size(categoryFolders,1);
    imagesByCategory = cell(1, numCategories);

    % get a list of all images
    for categ = 1:numCategories
        filelist = dir(fullfile(categoryFolders{categ,2},'*.*'));
        filelist(1:2) = [];
        imagesByCategory{categ} = {filelist.name};
    end

    % create the image pairings

    % Map 1 = C1 or C2 LEADING >> C6 TRAILING valid, C7 invalid
    % Map 2 = C4 or C5 LEADING >> C7 TRAILING valid, C6 invalid
    % Map 3 = C3 LEADING >> C8 or C9 TRAILING

    leading_maps = [1 1 2 2 3];
    leading_categories = [1 2 4 5 3];
    leading_probs = [.5 .5 .5 .5 1];

    trailing_maps = [1 1 2 2 3 3];
    trailing_categories = [6 7 7 6 8 9];
    trailing_probs = [.75 .25 .75 .25 .5 .5];

    no_blocks = 8;
    no_trials_per_block = 216;
    total_no_trials = no_blocks * no_trials_per_block;

    leading_image_list = cell(0);
    trailing_image_list = cell(0);
    for b = 1:no_blocks
        for i = 1:length(leading_maps)
            how_many_images = leading_probs(i)*no_trials_per_block / 3; % divided by 3 possible maps

            which_images = randsample(1:length(imagesByCategory{leading_categories(i)}),how_many_images);
            leading_image_list = [leading_image_list imagesByCategory{leading_categories(i)}(which_images)];
            imagesByCategory{leading_categories(i)}(which_images) = [];
        end

        for i = 1:length(trailing_maps)
            how_many_images = trailing_probs(i)*no_trials_per_block / 3; % divided by 3 possible maps

            which_images = randsample(1:length(imagesByCategory{trailing_categories(i)}),how_many_images);
            trailing_image_list = [trailing_image_list imagesByCategory{trailing_categories(i)}(which_images)];
            imagesByCategory{trailing_categories(i)}(which_images) = [];
        end
    end

    % shuffle the order per block
    for b = 1:no_blocks
        leading_image_list((b-1)*no_trials_per_block+1 : b*no_trials_per_block) = Shuffle(leading_image_list((b-1)*no_trials_per_block+1 : b*no_trials_per_block));
        trailing_image_list((b-1)*no_trials_per_block+1 : b*no_trials_per_block) = Shuffle(trailing_image_list((b-1)*no_trials_per_block+1 : b*no_trials_per_block));
    end



    % load task

    KbName('UnifyKeyNames');
    escapeKey = KbName('Escape');
    enterKey = KbName('Return');
    spacebarKey = KbName('Space');

    [mw, rect]=Screen('OpenWindow', 0, gray,[0 0 800 600]);
    % [mw, rect]=Screen('OpenWindow', 1, gray);

    centerX=rect(3)/2;
    centerY=rect(4)/2;


    crossLines=[-50, 50, 0 , 0
        0, 0, -50, 50];

    HideCursor();
    Screen('TextSize', mw, fontsize);
    Screen('TextFont', mw, 'Times');
    myText = 'The experiment is about to begin. \n\n In the following block, you will see a series of different images of everyday scenes. \n\n Occasionally, the image will appear upside down. \n\n Please press the space bar when the image is upside down. \n\n There will be 8 blocks, each lasting 10 minutes. \n\n When you are ready, press space to begin.';
    DrawFormattedText(mw, myText, 'center', 'center', white);
    Screen('Flip', mw);
    while true
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown && keyCode(KbName('space'))
            break;
        end
    end

    Screen('Flip', mw);
    WaitSecs(1);

    % for k = 1:no_blocks
    for k = 1:2
        savethis=[];
        logfile = [];
        ListenChar(2);

        % 5% catch trials where image is shown upside down.
        targets = [zeros(1,round(no_trials_per_block*.95)) ones(1,round(no_trials_per_block*.05))];
        targets = randsample(targets,length(targets));

        % for i = 1:no_trials_per_block
        for i = 1:10

            leading_image = imread(leading_image_list{(k-1)*no_trials_per_block+i});
            trailing_image = imread(trailing_image_list{(k-1)*no_trials_per_block+i});

            % TRY INTERPOLATING / RESAMPLING THESE IMAGES TO GET e.g.
            % 1024*1024*3

            % SHOW LEADING IMAGE

            TexturePointer = Screen(mw, 'MakeTexture', leading_image);
            Screen(mw, 'DrawTexture', TexturePointer);

            Screen('DrawLines', mw, crossLines, 5, white, [centerX centerY]);
            Screen('Flip', mw);

            iti_start = GetSecs;
            logfile.stimulustimes(i) = iti_start;
            tnow = GetSecs;

            %         io64(ioObj,address,1); % send EEG trigger

            while tnow < iti_start+pic
                tnow = GetSecs;
                [keyIsDown, secs, keyCode] = KbCheck;
                if keyCode(KbName('Space')) == 1
                    savethis=[savethis;
                        keyIsDown, secs, keycode(KbName('Space'))];
                    KbReleaseWait
                end
            end

            Screen('DrawLines', mw, crossLines, 5, white,[centerX centerY]);
            Screen('Flip', mw);

            %    io64(ioObj,address,0); % reset EEG trigger

            while tnow < iti_start+ISI
                tnow = GetSecs;
                [keyIsDown, secs, keyCode]=KbCheck;
                if keyCode(KbName('Space')) == 1
                    savethis=[savethis;
                        keyIsDown, secs, keyCode(KbName('Space'))]
                    KbReleaseWait;
                end
            end


            % SHOW TRAILING IMAGE

            if targets(i)==1 % show trailing image upside down
                TexturePointer = Screen(mw, 'MakeTexture', flipud(trailing_image));
            else % show normally
                TexturePointer = Screen(mw, 'MakeTexture', trailing_image);
            end

            Screen(mw, 'DrawTexture', TexturePointer);
            Screen('DrawLines', mw, crossLines, 5, white,[centerX centerY]);
            Screen('Flip', mw);

            iti_start = GetSecs;
            tnow = GetSecs;

            %         io64(ioObj,address,2); % send EEG trigger

            while tnow < iti_start+pic
                tnow = GetSecs;
                [keyIsDown, secs, keyCode]=KbCheck;
                if keyCode(KbName('Space')) == 1
                    savethis=[savethis;
                        keyIsDown, secs, keyCode(KbName('Space'))]
                    KbReleaseWait;
                end
            end

            Screen('DrawLines', mw, crossLines, 5, white,[centerX centerY]);
            Screen('Flip', mw);

            trial_ITI = randsample(ITI_range(1):.001:ITI_range(2),1);

            iti_start = GetSecs;
            tnow = GetSecs;

            %          io64(ioObj,address,0); % reset EEG trigger

            while tnow < iti_start+trial_ITI
                tnow = GetSecs;
                [keyIsDown, secs, keyCode]=KbCheck;
                if keyCode(KbName('Space')) == 1
                    savethis=[savethis;
                        keyIsDown, secs, keyCode(KbName('Space'))]
                    KbReleaseWait;
                end
            end

        end

        dt=datestr(now,'dd.mm.YY_HH:MM:SS');
        name=['ID_' num2str(ID) '_Block_' num2str(k)];
        logfile.name = name;
        logfile.ID = ID;
        logfile.block = k;
        logfile.responsetimes = savethis(:,2);
        logfile.targets = targets;
        logfile.leading_image_list = leading_image_list((k-1)*no_trials_per_block+1 : k*no_trials_per_block);
        logfile.trailing_image_list = trailing_image_list((k-1)*no_trials_per_block+1 : k*no_trials_per_block);

        save(fullfile(datapath,name),'logfile')

        if k < 8
            myText = 'Press Space to continue.';

            DrawFormattedText(mw, myText, 'center', 'center', white);
            Screen('Flip', mw);
            while true
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyIsDown && keyCode(KbName('Space'))
                    break;
                end
            end

            Screen('Flip', mw);
            WaitSecs(1);

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
        %  clear io64;
    end
