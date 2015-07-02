% XXX s - XXX TRs @ 2.00 sec/TR

function IEM_hexMap_attn

% mapping task - adapted from wmDrop_hexMap.m
% TCS 6/18/2015
%
% instead of grid of mapping positions, now present stimuli along a
% hexagonal grid which is offset a different amount each run
%
% on a subset of trials, mapping stimulus has contrast increment/decrement
% for a few flicker cycles - responsd as to whether it dimmed/brightened
% (only on those trials)
%




%% get user input about sub name etc...


warning('off','MATLAB:dispatcher:InexactMatch');
%Screen('Preference', 'SkipSyncTests', 1)
%get subject info
prompt = {'Subject Name','Ang offset? (degrees, -30 to 30)', 'Block number', 'Random seed', 'fMRI','Eye tracking','Contrast change? (%, 0-1)'};
%grab a random number and seed the generator (include this as a gui field
%in case we want to repeat an exact sequence)
s = round(sum(100*clock));
%put in some default answers
defAns = {'XXX','XXX','1',num2str(s),'1','0','0.2'};
box = inputdlg(prompt,'Enter Subject Information...', 1, defAns);

p.exptName = 'wmLazSpace_hexMap';

if length(box)==length(defAns)
    p.subName=char(box{1});
    p.sessionNum=1;%str2num(box{2});
    %p.cond = str2num(box{2});
    p.angOffset = str2num(box{2});
    p.nBlocks=eval(box{3});
    p.rndSeed=str2num(box{4});
    p.fMRI=str2num(box{5});
    p.eyeTracking = str2num(box{6});
    p.targContrastChange = str2num(box{7}); % WM difficulty is how far, in radians from center of stim, targ & probe are apart (multiples of pi/2)
    rand('state',p.rndSeed);  %actually seed the random number generator
else
    return
end

%p.keys = 66; %keys that correspond to elements in p.dir respectively (B)

if p.eyeTracking
    disp('setting up eyetracking...');
    try
        write_parallel(0); % initial communication w/ eyetracker?
    catch thisError
        disp('EyeTracking error...continuing with p.eyeTracking = 0;');
        disp(thisError);
        p.eyeTracking = 0;
        p.eyeTrackingStopped = 'at setup';
    end
end

%% --------------------begin user define parameters----------------------------%

% Screen (parameter) setup
p.windowed = ~p.fMRI;                     % if 1 then the display will be in a small window, use this for debugging.

% monitor stuff
p.refreshRate = 60;                 % refresh rate is normally 100 but change this to 60 when on laptop!
if p.fMRI
    p.vDistCM = 375; % CM
    p.screenWidthCM = 120; % CM
    % LUT:
    p.LUTfn = 'C:/Documents and Settings/serenceslab/My Documents/MATLAB/wmDrop/correctGammafMRI.mat';
else % behavioral room w/ eyetracker
    p.vDistCM = 62;
    p.screenWidthCM = 51;
    % LUT:
    p.LUTfn = 'Y:/Monitor Calib/correctGamma_eyetracker_2.28.12.mat';
    
end

% ideally, let's do 2 hex steps = wmDrop's target eccentricity (max ecc of
% target center = 3*steps, then stimRad on top oft that

% maybe 35 if we drop center point, which we should...
p.nLoc = 37; % instead of grid, let's just do number of positions, then have x, y coords for this number of positions

%stimulus geometry (in degrees, note that these vars have the phrase 'Deg'
%in them, used by pix2deg below)


%p.stepSize = 1.625;% p.wmEcc/2;
p.stepSize = 1.5;


p.radDeg = 0.9;


p.jitterRadDeg = 0.5; % jitter each position using draw from uniform circle w/ this radius

%p.radDeg = p.usedScreenSizeDeg/(p.nLoc+2);
%p.sfDeg  = 1.5 * .4532;
p.sfDeg = .6785;  % TODO: check this out...


% font stuff
p.fontSize = 24;
p.fontName = 'ARIAL';
p.textColor = [100, 100, 100];

%Response keys
if ismac
    p.escape = 41;
else
    p.escape = 27;
end

% TODO: also allow left/right arrow...
p.keys = [KbName('b'),KbName('y')]; % 1 and 2 keys
p.space = KbName('space');
p.start = KbName('t');

% TODO: use LUT, make sure these are consistent w/ WM

p.backColor     = [128, 128, 128];      % background color

%fixation point properties
p.fixColor      = [180 180 180];%[155, 155, 155];
% p.fixSizeDeg    = .25 * .4523;                  % size of dot
p.fixSizeDeg    = .2;
p.appSizeDeg    = p.fixSizeDeg * 4;     % size of apperture surrounding dot



% target properties (fix/checks)
p.stimContrast = 0.7;

%p.contrastChange = 0.5; % will change by this*100% difference from bg color
%p.targetContrast = p.stimContrast - p.stimContrastChange;
%p.targetFixColor = p.backColor + (abs(p.fixColor - p.backColor)) * (1-p.fixContrastChange);

% trial setup

p.repetitions = 1;
if p.windowed % DEBUG MODE
    p.repetitions = 1;
end
p.nTrainingTrials = p.repetitions * p.nLoc;    % length p.cond, for now, should always be 1

p.nTargTrials = 10; % 5 increment, 5 decrement (21% target trials overall)
p.nTrials = p.nTrainingTrials + p.nTargTrials; % no stim

p.percentTarg = p.nTargTrials/p.nTrials;

% stimulus timing (in seconds - convert to vid frames, if necessary, later)
p.stimExposeDur = 3; % keep this constant
%p.wmTargetDur   = 0.5; % s
%p.wmProbeDur    = 0.75; % s
p.flickerFreq   = 6; % Hz
p.flickerPer  = 1/p.flickerFreq; % s

p.targOnsetWindow = [0.5 2]; % times at whcih target could start


p.ITIs = linspace(2,6,p.nTrials)';
p.ITIs = p.ITIs(randperm(length(p.ITIs)));
%p.ITI           = 2.0 * ones(p.nTrials,1);
p.trialDur      = p.stimExposeDur + p.ITIs;
p.nTRsWait      = 0;
p.startWait     = 2;
p.passiveDur    = 10; % in sec, how long to wait at end of a block (and beginning?)
p.TR            = 2.00;
p.expDur        = sum(p.trialDur) + p.passiveDur + p.startWait + p.nTRsWait*p.TR;

p.targDur = 0.5;


disp(sprintf('SCAN DURATION: %i',p.expDur));
disp(sprintf('# TRs: %i',ceil(p.expDur/p.TR)));

ListenChar(2);

%% --------------------Screen properties----------------------------%

load(p.LUTfn);
p.LUT = 255*correctGamma;





% correct all the colors
p.fixColor = round(p.LUT(p.fixColor))';
%p.targetFixColor = p.LUT(round(p.targetFixColor))';
p.textColor = round(p.LUT(p.textColor))';
p.backColor = round(p.LUT(p.backColor))';

%Start setting up the display
AssertOpenGL; % bail if current version of PTB does not use

% Open a screen
Screen('Preference','VBLTimestampingMode',-1);  % for the moment, must disable high-precision timer on Win apps

% figure out how many screens we have and pick the last one in the list
s=max(Screen('Screens'));
p.black = BlackIndex(s);
p.white = WhiteIndex(s);

if p.windowed
    Screen('Preference', 'SkipSyncTests', 1);
    [w p.sRect]=Screen('OpenWindow', s, p.backColor, [50,50,800,600]);
else
    [w p.sRect]=Screen('OpenWindow', s, p.backColor);
    HideCursor;
end

% Enable alpha blending with proper blend-function. We need it
% for drawing of smoothed points:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% test the refesh properties of the display
p.fps=Screen('FrameRate',w);          % frames per second
p.ifi=Screen('GetFlipInterval', w);   % inter-frame-time
if p.fps==0                           % if fps does not register, then set the fps based on ifi
    p.fps=1/p.ifi;
end
p.flickerFrames = 1/p.flickerFreq*p.fps;

% make sure the refreshrate is ok
if abs(p.fps-p.refreshRate)>5
    Screen('CloseAll');
    disp('CHANGE YOUR REFRESH RATE')
    ListenChar(0);
    %clear all;
    return;
end

%% HEXAGONAL GRID - IMPORTANT! for now, just manually generating

% 37 points 
xgrid = [-1.5 -0.5 0.5 1.5  -2 -1 0 1 2   -2.5 -1.5 -0.5 0.5 1.5 2.5  -3 -2 -1 0 1 2 3  -2.5 -1.5 -0.5 0.5 1.5 2.5  -2 -1 0 1 2  -1.5 -0.5 0.5 1.5];
ygrid = [-3 -3 -3 -3       -2 -2 -2 -2 -2  -1 -1 -1 -1 -1 -1           0 0 0 0 0 0 0    1 1 1 1 1 1                  2 2 2 2 2         3 3 3 3]*(sqrt(3)/2);

% make sure p.nLoc==length(xgrid);


% cartesian coords

[thgrid, rgrid] = cart2pol(xgrid,ygrid); % CHECK THIS FOR up/down
thgrid = thgrid + deg2rad(p.angOffset);
[xgrid_adj, ygrid_adj] = pol2cart(thgrid,rgrid);

% SCREEN COORDINATES!!!!

p.xGridDeg =  xgrid_adj*p.stepSize;
p.yGridDeg = -ygrid_adj*p.stepSize;

clear xgrid ygrid thgrid rgrid xgrid_adj ygrid_adj;

% convert relevant timing stuff to vid frames for stim presentation
p.stimExpose = round((p.stimExposeDur*1000)./(1000/p.refreshRate));

% if running the real experiment (not debugging), then hide cursor and set
% priority high
if ~p.windowed
    HideCursor; % Hide the mouse cursor
    % set the priority up way high to discourage interruptions
    Priority(MaxPriority(w));
end

% convert from degrees to pixel units
p = deg2pix(p);

% now we have p.usedScreenSizePix, p.radPix

% compute and store the center of the screen: p.sRect contains the upper
% left coordinates (x,y) and the lower right coordinates (x,y)
center = [(p.sRect(3) - p.sRect(1))/2, (p.sRect(4) - p.sRect(2))/2];

Screen('TextSize', w, p.fontSize);
Screen('TextStyle', w, 1);
Screen('TextFont',w, p.fontName);
Screen('TextColor', w, p.textColor);

% file/directory operations
p.root = pwd;

if ~isdir([p.root '/data/'])
    mkdir([p.root '/data/']);
end
p.subNameAndDate = [p.subName '_' datestr(now,30)];

%% start a block loop
for b=p.nBlocks
    
    
    fName = sprintf('%s/data/%s_%s_sess%02.f_run%02.f.mat',p.root, p.subName,p.exptName,p.sessionNum,b);
    %fName=[p.root, '/Subject Data/', p.exptName, num2str(p.cond),'_', p.subNameAndDate, '_sess', num2str(p.sessionNum), '_blkNum', num2str(b), '.mat'];
    if exist(fName,'file')
        Screen('CloseAll');
        msgbox('File name already exists, please specify another', 'modal')
        ListenChar(0);
        return;
    end
    
    % when eyetracking updated, use: p.et_fname = sprintf('%s/%s/%s_%s_sess%02.f_run%01.f_%s.idf',p.exptName,p.subName, p.subName,p.exptName,p.sessionNum,p.runNum, datestr(now,30));
    
    % open eyetracker file
    if p.eyeTracking
        try
            write_parallel(64); % open file
            disp('ET started');
        catch thisError
            disp('Error during file opening');
            p.eyeTracking = 0;
            p.eyeTrackingStopped = 'After 64';
            disp(thisError);
        end
    end
    
    % randomly repeat some positions, though they'll be jittered
    % differently
    targ_locs_tmp = randperm(length(p.xGridDeg));
    
    p.stimLocs = (1:length(p.xGridDeg))';
    p.stimLocs = [p.stimLocs; targ_locs_tmp(1:p.nTargTrials).'];
    clear targ_locs_tmp;
    
    p.stimLocsDeg = [p.xGridDeg(p.stimLocs)' p.yGridDeg(p.stimLocs)'];
    
    p.targPresent = zeros(size(p.stimLocs));
    p.targPresent((p.nTrainingTrials+1):end) = 1;
    
    p.contrastChange = nan(p.nTrials,1);
    p.contrastChange(p.nTrainingTrials + (1:(p.nTargTrials/2))) = 1; % decrement
    p.contrastChange(p.nTrainingTrials + ((1+p.nTargTrials/2):p.nTargTrials)) = 3; % increment
    
    
    
    
    
    % let's just do these as x, y rather than r, theta
    % jitter each position by random draw from uniform disc w/ radius
    % p.jitterRad
    %tang = rand(p.nTrials,1); trad = rand(p.nTrials,1)*2*pi;
    jrad = p.jitterRadDeg*rand(p.nTrials,1); jang = rand(p.nTrials,1)*2*pi;
    p.jitterDeg =  (repmat(jrad,1,2) .* [cos( jang ) sin( jang ) ]);   % CARTESIAN COORDINATES UP HERE, SCREEN COORDS JUST BEFORE STIM PRESENTATION!!!!!
    clear tang trad;
    
    p.stimLocsDeg = p.stimLocsDeg + p.jitterDeg; % in cartesian coords
    p.stimLocsPix = p.stimLocsDeg * p.ppd;
    
    
    
    
    p.rndInd = randperm(p.nTrials);
    
    
    p.stimLocs = p.stimLocs(p.rndInd);
    p.stimLocsDeg = p.stimLocsDeg(p.rndInd,:);
    p.stimLocsPix = p.stimLocsPix(p.rndInd,:);
    
    p.targPresent = p.targPresent(p.rndInd);
    p.contrastChange = p.contrastChange(p.rndInd);
    
    % maybe also targ dimming sequence?
    
    
    
    %allocate some arrays for storing the subject response
    p.responseTime =          nan(p.nTrials,1);       % store the rt on each trial
    p.responseTimeRaw = nan(p.nTrials,1);             % raw "timestamp" to be compared to trial events
    
    p.resp =      nan(p.nTrials,1);
    p.correct =   nan(p.nTrials,1);
    p.falseAlarm =nan(p.nTrials,1);
    
    
    p.trialStart =  nan(p.nTrials,1);    % equal to stimStart
    p.trialEnd =  nan(p.nTrials,1);
    p.stimStart =  nan(p.nTrials,1);
    p.stimEnd =  nan(p.nTrials,1);
    p.targStart = nan(p.nTrials,1);
    
    
    ch = make_checkerboard(p.radPix,p.sfPix,p.stimContrast);
    ch{1} = (ch{1}-min(ch{1}(:)))/(max(ch{1}(:))-min(ch{1}(:)));
    ch{1} = ch{1}-0.5; 
    ch{2} = (ch{2}-min(ch{2}(:)))/(max(ch{2}(:))-min(ch{2}(:)));
    ch{2} = ch{2}-0.5;
    % n_contrasts x n_phases
    stim = nan(3,2);
    
    p.checkContrasts = [p.stimContrast - p.targContrastChange p.stimContrast p.stimContrast + p.targContrastChange];
    p.checkContrasts(p.checkContrasts>1)=1;
    
    % from wmLazSpace_stim1.m
    % we want stim to be n_contrasts x 2
    for cc = 1:length(p.checkContrasts)
        
        for ii = 1:length(ch)
            thisch = (ceil(255*(ch{ii}*p.checkContrasts(cc)))+1)+127;
            %sca;
            stim(cc,ii)=Screen('MakeTexture', w, round(p.LUT(thisch)));
            
            clear thisch;
        end
    end
    
    % when we show stimuli, will be
    % stim(p.stimTargSequ(t,frame),p.flickerSequ(f)) when flickerSequ ~= 2
    
    % pick the stimulus sequence for every trial (the exact grating to be shown)
    p.flickerSequ = repmat([ones(1,p.flickerFrames/2) 0*ones(1,p.flickerFrames/2) 2*ones(1,p.flickerFrames/2) 0*ones(1,p.flickerFrames/2)],1,(0.5)*p.stimExpose/p.flickerFrames);
    
    p.stimTargSequ = 2*ones(p.nTrials,length(p.flickerSequ));
    
    
    
    
    targ_fr = find(diff(p.flickerSequ))+1;
    targ_t = targ_fr*p.ifi;
    %targ_fr = targ_fr(6:(end-11)); % targ starts b/w 0.5 s and end-1 s, is 0.5 s long
    targ_fr = targ_fr(targ_t>=p.targOnsetWindow(1) & targ_t<=p.targOnsetWindow(2));
    
    % time and frame at which target starts
    p.targOnsetT = nan(p.nTrials,1);
    p.targOnsetF = nan(p.nTrials,1);
    
    for i=1:p.nTrials
        
        if p.targPresent(i)==1
            % this will be "2" for all non-targ trials, and "1" for ~500 ms
            % for targ trials w/ decrement, 3 with increment
            p.targOnsetF(i) = targ_fr(randsample(length(targ_fr),1));
            p.stimTargSequ(i,p.targOnsetF(i):(p.targOnsetF(i)+p.targDur*p.fps-1)) = p.contrastChange(i);
            p.targOnsetT(i) = p.ifi*p.targOnsetF(i);
        end
        
        
    end
    
    %
    
    str1 = sprintf('Block %i complete',b);
  
    instr_text = 'Did checkerboard dim (left) or brighten (right)?';
    instr_text2 = '(Only some trials)';
    
    tCenter1 = [center(1)-RectWidth(Screen('TextBounds', w, instr_text))/2 center(2)/2];
    tCenter2 = [center(1)-RectWidth(Screen('TextBounds', w, instr_text2))/2 center(2)/2];
    
    Screen('DrawText', w, instr_text, tCenter1(1), tCenter1(2)-100, p.textColor);
    Screen('DrawText', w, instr_text2, tCenter2(1), tCenter2(2), p.textColor);
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %change fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    
    %after all initialization is done, sit and wait for scanner synch (or
    %space bar)
    resp=0;
    %     disp('checking for response');
    %     while 1
    %         [resp, timeStamp] = checkForResp([p.space,p.start],p.escape);
    %         if resp==p.space || resp==p.start || resp==p.escape
    %             break;
    %         end
    %         if resp==-1; ListenChar(0); return; end;
    %     end
    
    FlushEvents;
    GetChar;
    
    disp('starting block');
    cumTime = GetSecs;
    
    % start eyetracking...
    if p.eyeTracking
        try
            write_parallel(192); % start recording
            disp('recording started');
        catch thisError
            disp('Error during recording initiation');
            p.eyeTracking = 0;
            p.eyeTrackingStopped = 'After 192';
            disp(thisError);
        end
    end
    
    %waited = 1; % we've already waited for 1 TR here...
    
    
    FlushEvents;
    
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %change fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    cumTime = GetSecs;
    p.startExp = GetSecs;%cumTime;
    
    while GetSecs < p.startExp + p.startWait
        [resp, timeStamp] = checkForResp([],p.escape);
        if resp == p.escape
            Screen('CloseAll');
            clear mex;
            return;
        end
        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %change fixation point
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        %continue;
    end
    
    cumTime = GetSecs;
    
    
    save(fName,'p');
    
    %% here is the start of the trial loop
    for t=1:p.nTrials
        
        % FROM HERE BELOW, IN SCREEN COORDINATES!!!!!
        xLoc = center(1) + p.stimLocsPix(t,1);
        yLoc = center(2) - p.stimLocsPix(t,2);
        
        %targXLoc = p.targetLoc(t,1);% + center(1);
        %targYLoc = p.targetLoc(t,2);% + center(2);
        
        
        if ~p.fMRI
            disp(sprintf('dim fix = %i, dim stim = %i',p.dimFix(t),p.dimStim(t)));
        end
        
        p.trialStart(t) = GetSecs;
        
        if p.eyeTracking
            try
                write_parallel(192+t); % event ID will be trial number for now
                disp(sprintf('Start of trial %i',t));
            catch thisError
                disp(thisError);
                p.eyeTracking = 0;
                p.eyeTrackingStopped = sprintf('Start of trial %i (%i)',t,192+t);
            end
        end
        
        % now just use CenterRectOnPoint
        stimRect = CenterRectOnPoint([0 0 2*p.radPix 2*p.radPix],xLoc, yLoc);
        
        
        
        
        frmCnt=1;
        p.stimStart(t) = GetSecs;   % start a clock to get the stim onset time
        
        %% STIMULUS
        %       while frmCnt<=p.stimExpose % if we want multiple exposure durations, add that here
        for ff = 1:length(p.flickerSequ)
            
            
            
            % when we show stimuli, will be
            % stim(p.stimTargSequ(t,frame),p.flickerSequ(f)) when flickerSequ ~= 2
            
            
            
            % checkerboard
            if p.flickerSequ(ff)~=0
                Screen('DrawTexture',w,stim(p.stimTargSequ(t,ff),p.flickerSequ(ff)),Screen('Rect',stim(p.flickerSequ(ff))),stimRect);
            end
            
            % apperture around fixation
            Screen('DrawDots', w, [0 0], p.appSizePix, p.backColor, center, 1);
            
            
            
            Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0);
            Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            Screen('Flip', w);
            
            % if first target frame, save target onset (for response time)
            if p.stimTargSequ(t,ff)~=2 && isnan(p.targStart(t))
                p.targStart(t) = GetSecs;
            end
            
            % check response...
            
            [resp, timeStamp] = checkForResp(p.keys, p.escape);
            if resp==-1; ListenChar(0); return; end;
            
            if resp && find(p.keys==resp) && isnan(p.resp(t))
                p.resp(t) = find(p.keys==resp); % (FIRST response is that recorded)
                p.responseTimeRaw(t) = timeStamp;
            end
            
            
            
        end
        
        p.stimEnd(t) = GetSecs;
        
        
        
        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0);
        Screen('DrawingFinished',w);
        Screen('Flip',w);
        
        % in case there's weird timing issues?
        
        while GetSecs <= cumTime + p.stimExposeDur
            [resp, timeStamp] = checkForResp(p.keys, p.escape);
            if resp==-1; ListenChar(0); return; end;
            
            if resp && find(p.keys==resp) && isnan(p.resp(t))
                p.resp(t) = find(p.keys==resp); % (most recent response is that recorded)
                p.responseTimeRaw(t) = timeStamp;
            end
            
        end
        
        % clear out screen
        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %draw fixation point
        Screen('Flip',w);
        
        % write end of eyetracking event (give time for return to fix)
        if p.eyeTracking
            try
                write_parallel(255); % event ID will be trial number for now
                disp(sprintf('ITI of trial %i',t));
                %clear msg out;
            catch thisError
                disp(thisError);
                p.eyeTracking = 0;
                p.eyeTrackingStopped = sprintf('ITI of trial %i (%i)',t,255);
                % attempt to save
            end
        end
        
        while GetSecs <= cumTime + p.trialDur(t)
            [resp, timeStamp] = checkForResp(p.keys, p.escape);
            if resp==-1; ListenChar(0); return; end;
            
            if resp && find(p.keys==resp) && isnan(p.resp(t))
                p.resp(t) = find(p.keys==resp); % (most recent response is that recorded)
                p.responseTimeRaw(t) = timeStamp;
                
            end
            
        end
        
        
        %% compute response acc
        if p.resp(t)==2   % for simplicity, change our "2" resposne to "3" to match contrast dec/inc
            p.resp(t)=3;
        end
        
        if isnan(p.contrastChange(t)) && ~isnan(p.resp(t))
            p.falseAlarm(t) = 1;
            p.correct(t) = 0;
        elseif isnan(p.contrastChange(t)) && isnan(p.resp(t))
            p.falseAlarm(t) = 0;
            p.correct(t) = 1;
        elseif p.resp(t)==p.contrastChange(t)
            p.correct(t) = 1;
        else
            p.correct(t) = 0;
        end
        
        % compute response times (if there was a contrast change, relative
        % to target change, otherwise, relative to stim start)
        if ~isnan(p.contrastChange(t))
            p.responseTime(t) = p.responseTimeRaw(t) - p.targStart(t);
        else
            p.responseTime(t) = p.responseTimeRaw(t) - p.stimStart(t);
        end
        
        
        
        
        cumTime = cumTime + p.trialDur(t);
        
        
        p.trialEnd(t) = GetSecs;
        
        save(fName, 'p');
        
    end
    %% end trial loop
    
    %   10s passive fixation at end of block
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, center, 0); %change fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    while GetSecs <= cumTime + p.passiveDur
        [resp, timeStamp] = checkForResp(p.escape, p.escape);
        if resp==-1; ListenChar(0); return; end;
    end
    
    p.endExp = GetSecs;
    
    p.accuracy = 100*nanmean(p.correct(~isnan(p.contrastChange)));% p.wmCorrect/(p.wmCorrect+p.wmIncorrect)*100;
    p.nFalseAlarms = sum(p.falseAlarm(~isnan(p.falseAlarm)));
    
    %save trial data from this block
    save(fName, 'p');
    
    str1 = sprintf('Block %i complete',b);
    tCenter1 = [center(1)-RectWidth(Screen('TextBounds', w, str1))/2 center(2)/2];
    str2 = sprintf('Accuracy (targets): %.02f%%',p.accuracy);
    tCenter2 = [center(1)-RectWidth(Screen('TextBounds', w, str2))/2 center(2)/2];
    str3 = sprintf('%i False alarms',p.nFalseAlarms);
    tCenter3 = [center(1)-RectWidth(Screen('TextBounds', w, str3))/2 center(2)/2];
    
    Screen('DrawText', w, str1, tCenter1(1), tCenter1(2)-100, p.textColor);
    Screen('DrawText', w, str2, tCenter2(1), tCenter2(2), p.textColor);
    Screen('DrawText', w, str3, tCenter3(1), tCenter3(2)+100, p.textColor);
    
    
    % put up a message to wait for a space bar press.
    Screen('Flip', w);
    while resp~=p.space
        [resp, timeStamp] = checkForResp(p.space, p.escape);
        if resp==-1; ListenChar(0); return; end;
    end
    
    % stop recording and save
    if p.eyeTracking
        try
            write_parallel(0); % close file, end recording
        catch thisError
            disp(thisError);
            p.eyeTracking = 0;
            p.eyeTrackingStopped = sprintf('End of block %i',b);
        end
    end
    
end
%% end block loop

ListenChar(0);
Screen('CloseAll');
return


%--------------------------------------------------------------------------
function p = deg2pix(p)
% converts degrees visual angle to pixel units before rendering
% with PTB. Needs p.screenWidthCM and p.vDistCM
% js - 10.2007

% figure out pixels per degree, p.sRect(1) is x coord for upper left of
% screen, and p.sRect(3) is x coord for lower right corner of screen
p.ppd = pi * (p.sRect(3)-p.sRect(1)) / atan(p.screenWidthCM/p.vDistCM/2) / 360;

% get name of each field in p
s = fieldnames(p);

% convert all fields with the word 'Deg' from degrees visual angle to
% pixels, then store in a renmaed field name
for i=1:length(s)
    ind = strfind(s{i}, 'Deg');
    if ind
        curVal = getfield(p,s{i});
        tmp = char(s{i});
        newfn = [tmp(1:ind-1), 'Pix'];
        p = setfield(p,newfn,curVal*p.ppd);
    end
end
