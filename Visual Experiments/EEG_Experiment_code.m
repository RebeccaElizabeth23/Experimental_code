
function EEG_Experiment_code

% pilot experiment to determine optimal wedge size

close all;
E.subj = 'test';
useVP = 0;             % 0 means windowed screen, 1 means full screen
trackermode = 0;       % 0 means not tracking, 1 means tracking
atttask = 'H';            % H means high and L means low

disp(atttask)

% This should stop writing to the code when running in eyetracker mode
echo off;

% Sets examples stimuli parameters
Hori = [45 135];
Lori = [0 0];

% Sets out features of the wedges and central stimuli
ST.npixelsperdegree = 36;       % at 57cm
ST.wedgesizesdeg = 10;  % outer radius of wedges
ST.wedgeinnerradius = 2;
ST.duration = 72;   % in seconds
ST.centrefreq = 4;
ST.flankfreqs = 120./[18 16];
ST.flankcontrast = 0.15;
ST.circfreq = 8;
ST.radfreq = 8;
ST.envradfreq = 2;
ST.envsharpness = 0.99;
ST.targetsize = 1*ST.npixelsperdegree;
ST.targetSF = [2 4];
ST.imsize = ST.npixelsperdegree*max(ST.wedgesizesdeg)*2;

E.nreps = 2;
E.ncond = 4;
E.condsdone(1:E.ncond) = 0;
E.ITI = 1; % This should be 10 seconds

% Randomised trial order for hemifield presentation
E.ntrials = 4;
E.trialnumber = 8;
E.triallist = [1,2,3,4,1,2,3,4];
E.trialorder = E.triallist(randperm(length(E.triallist)));
E.trialsdone(1:E.ntrials) = 0;

% The orientation order of the target - 12 orientations with 6 trials 
%(72 gratings for the 72 second block)
% 4 prerandomised cent. target - 6 orientations and of those 3 for each sf
% Done to shop the target appearing in the 1st 2 and ...
% to stop the same type of target happening twice in a row


targetlist1 = [[75 2]; [30 2]; [90 2]; [75 4]; [90 2]; [60 2]; [135 4]; ...
[120 4]; [45 2]; [0 4]; [135 2]; [150 2]; [60 2]; [0 2]; [165 4]; [60 4];...
[75 4]; [165 2]; [60 4]; [150 4]; [60 2]; [120 2]; [45 2]; [90 4]; [0 2];...
[60 4]; [15 2]; [120 2]; [15 4]; [0 4]; [120 4]; [90 4]; [15 4]; [150 4];...
[135 4]; [165 2]; [15 2]; [135 4]; [105 2]; [135 2]; [0 4]; [105 4];...
[30 2]; [45 4]; [150 2]; [105 4]; [75 2]; [165 4]; [30 4]; [135 2]; ...
[45 4]; [120 2];[165 4]; [45 2]; [90 4]; [15 2]; [45 4]; [165 2]; [0 2];...
[120 4]; [15 4]; [150 2]; [105 2]; [30 2]; [75 2]; [30 4]; [90 2]; ...
[30 4]; [105 2]; [150 4]; [105 4]; [75 4]];
        
targetlist2 = [[15 2];[165 2];[120 2];[75 2];[60 2];[15 2];[60 2];[150 4];...
[15 4];[105 2];[120 4];[30 2];[90 4];[165 4];[135 2];[60 4];[45 4];...
[150 2];[30 4]; [90 2];[0 2];[150 4];[30 4];[165 2];[0 4];[105 2];[75 4];...
[15 2];[30 4];[165 4];[30 2];[45 2];[120 4];[15 4];[135 4];[75 4];[105 2];...
[90 4];[45 4];[90 2];[120 2];[45 4];[105 4];[90 4];[0 2];[165 4];[45 2];...
[105 4];[60 4];[135 2];[0 4];[15 4];[30 4];[75 2];[135 4];[75 4];[45 2];...
[60 4];[150 2];[30 2];[0 4];[150 2];[165 2];[60 2];[135 4];[0 2];[120 2];...
[135 2];[75 2];[105 4];[90 2];[120 4]];
        
        
targetlist3 = [[75,4];[15,2];[30,4];[90,2];[60,2];[0,2];[135,2];[45,4];...
[120,2];[165,4];[90,4];[135,2];[0,4];[75,2];[0,2];[75,4];[30,2];[15,2];...
[60,4];[15,2];[60,4];[165,4];[90,4];[135,4];[150,2];[15,4];[75,2];[135,4];...
[45,2];[75,4];[45,4];[135,2];[165,2];[150,2];[30,4];[15,4];[120,4];[60,2];...
[120,4];[0,4];[150,4];[30,2];[105,4];[60,2];[165,2];[150,4];[105,4];[45,2];...
[150,4];[90,2];[150,2];[15,4];[120,4];[90,4];[120,2];[165,4];[75,2];...
[30,2];[45,4];[0,2];[30,4];[135,4];[105,4];[165,2];[105,2];[60,4];[120,2];...
[45,2];[105,2];[90,2];[0,4];[105,2]];
        
        
targetlist4 = [[135,2];[150,2];[45,2];[75,4];[60,4];[135,2];[105,4];[0,2];...
[135,2];[105,2];[120,4];[135,4];[15,2];[60,2];[105,4];[120,2];[90,4];...
[120,2];[90,2];[165,4];[75,2];[0,4];[75,2];[120,4];[15,4];[30,4];[90,2];...
[15,2];[120,4];[45,4];[30,4];[150,2];[0,2];[60,4];[165,4];[30,4];[45,2];...
[75,4];[165,2];[120,2];[60,2];[105,4];[150,4];[45,2];[15,2];[0,4];[105,2];...
[0,4];[15,4];[150,4];[165,2];[75,4];[15,4];[30,2];[165,2];[0,2];[90,4];...
[45,4];[30,2];[150,4];[105,2];[135,4];[30,2];[60,2];[90,2];[135,4];[165,4];...
[45,4];[150,2];[90,4];[60,4];[75,2]];
        
% Randomised order for target list - this will then randomly pick one of
% the 4 target list orders to show on every trial
E.ntlist = 4;
E.listpresentations = 8;
E.tlistlist = [1,2,3,4,1,2,3,4];
E.tlistorder = E.tlistlist(randperm(length(E.tlistlist)));
E.tlistssdone(1:E.ntlist) = 0;



% checks if a participants file exists and names one if not
fname = strcat('Results/',E.subj,'data.mat');
if exist(fname, 'file')
    load(fname);
else
    runnumber = 0;
    save(fname,'runnumber');
end
runnumber = runnumber + 1;


WaitSecs(0.01); % important to load in the MEX file before the expt starts
GetSecs;
InitializePsychSound;
tr = PsychPortAudio('Open',[],[],[],[],1);
PsychPortAudio('FillBuffer', tr, MakeBeep(440*sqrt(2),0.05,44100).*0.5);

try  % start the 'try/catch' loop
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    PsychGPUControl('SetDitheringEnabled', 0);
    
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    if useVP == 1;      % using a ViewPixx or ProPixx
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        Datapixx('EnableVideoScanningBacklight');% Required if VIEWPixx.
        Datapixx('RegWr');
        
        PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
        PsychImaging('AddTask', 'General', 'EnableDataPixxM16OutputWithOverlay');
        [w, winRect] = PsychImaging('OpenWindow', screenNumber, 0, [], [], [], 0, 0, kPsychNeedFastBackingStore);
        % THIS IS THE IMPORTANT THING TO DO, NOTE THE LAST ARGUMENT IS 0.
        Screen('LoadNormalizedGammaTable', w, linspace(0,1,256)'*ones(1,3), 0);
        
        HideCursor;
        ST.greylevel = doimagegamma(0.5);
    else
        % Make window on current screen
        Screen('Preference', 'SkipSyncTests', 1);
        rect = [1 1 1024 1024];
        [w, winRect] = Screen('OpenWindow',screenNumber,128,rect);
        ST.greylevel = 128;
    end
    
    [width, height] = Screen('WindowSize', w);
    ifi=Screen('GetFlipInterval', w);
    ifims = ifi*1000;
    
    
    %Creates the grey background and draws
    Screen('FillRect', w, ST.greylevel);
    drawfixation(w, width, height, 0);
    Screen('Flip', w);

    % no  frames in 1 cycle (should be 12 or 24)
    ST.nframestarget = round((1000/ifims)/ST.centrefreq);
    % no  frames in 1 cycle (should be 12 or 24)
    ST.nframesflankers = round((1000/ifims)./ST.flankfreqs); 
    
    % Creates target gratiing
    targetframes = ST.nframestarget*4;
    targetwaveform = sin(2 .* ST.centrefreq .* (ifims:ifims:targetframes*ifims) .* pi./1000 + (3*pi/2));
    targetwaveform = (targetwaveform + 1)./2;
    
    % Creates the central task stimuli at 2 sf
    targetstim_2 = mkgrating(ST.targetsize, ST.targetSF(1), 90, 90, 1) .* make_soft_window(ST.targetsize,ST.targetsize,0.9);
    for n = 1:targetframes
        comp = targetwaveform(n).*targetstim_2;
        comp= (1+comp)/2;
        if useVP == 1
            comp = doimagegamma(comp);
        end
        targetlist_2(n) = Screen('MakeTexture', w, comp, [], [], 2);
    end
    
   % Creates the central task stimuli at 4 sf
    targetstim_4 = mkgrating(ST.targetsize, ST.targetSF(2), 90, 90, 1) .* make_soft_window(ST.targetsize,ST.targetsize,0.9);
    for n = 1:targetframes
        comp = targetwaveform(n).*targetstim_4;
        comp = (1+comp)/2;
        if useVP == 1
            comp = doimagegamma(comp);
        end
        targetlist_4(n) = Screen('MakeTexture', w, comp, [], [], 2);
    end
    
    
    
    
    flankframes = ST.nframesflankers;
    flankwaveformL = sin(2 .* ST.flankfreqs(1) .* (ifims:ifims:flankframes(1)*ifims) .* pi./1000 + (3*pi/2));
    flankwaveformL = (flankwaveformL + 1)./2;
    flankwaveformR = sin(2 .* ST.flankfreqs(2) .* (ifims:ifims:flankframes(2)*ifims) .* pi./1000 + (3*pi/2));
    flankwaveformR = (flankwaveformR + 1)./2;
    
    
    
    [u, v] = meshgrid((1:ST.imsize)-(ST.imsize+2)/2, (1:ST.imsize)-(ST.imsize+2)/2);
    f = sqrt(u.^2 + v.^2); % radial coordinate
    theta = atan2(v,u);    % angular coordinate
    f = f./max(f(ST.imsize/2,:));
    theta = ((theta./pi)+1)./2;
    
    grating1 = cos(ST.circfreq .* f .* 2 .* pi);
    grating1 = round((grating1+1)./2);
    
    grating2 = cos(ST.radfreq .* theta .* 2 .* pi);
    grating2 = round((grating2+1)./2);
    
    grating3 = cos(ST.envradfreq .* theta .* 2 .* pi);
    grating3 = round((grating3+1)./2);
    
    plaid = (2.*(grating1-0.5)) .* (2.*(grating2-0.5));
    
    innerwidth = ST.wedgeinnerradius*ST.npixelsperdegree*2;
    innerwindow = zeros(ST.imsize);
    innerwindow((1+ST.imsize/2-innerwidth/2):(ST.imsize/2+innerwidth/2),(1+ST.imsize/2-innerwidth/2):(ST.imsize/2+innerwidth/2)) = make_soft_window(innerwidth,innerwidth,ST.envsharpness);
    
    % Determines the size of the flankers
    flanksize = ST.npixelsperdegree*ST.wedgesizesdeg*2;
    outerwindowtemp = make_soft_window(flanksize,flanksize,ST.envsharpness);
    outerwindow = zeros(ST.imsize);
    outerwindow((1+ST.imsize/2-flanksize/2):(ST.imsize/2+flanksize/2),(1+ST.imsize/2-flanksize/2):(ST.imsize/2+flanksize/2)) = outerwindowtemp;
    
    fullwindow = outerwindow - innerwindow;
    
    newplaid = ST.flankcontrast.*plaid.*fullwindow.*grating3;
    
    for n = 1:flankframes(1)
        
        comp = flankwaveformL(n).*newplaid(:,1:ST.imsize/2);
        comp = (1+comp)/2;
        if useVP  == 1
            comp = doimagegamma(comp);
        end
        flanklistL(n) = Screen('MakeTexture', w, comp, [], [], 2);
    end
    for n = 1:flankframes(2)
        comp = flankwaveformR(n).*newplaid(:,(ST.imsize/2+1):end);
        comp = (1+comp)/2;
        if useVP  == 1
            comp = doimagegamma(comp);
        end
        flanklistR(n) = Screen('MakeTexture', w, comp, [], [], 2);
    end
    
    
    r = [1 1 ST.targetsize ST.targetsize];
    targetrect = CenterRectOnPoint(r, width*0.5, height*0.5);
    examplelrect = CenterRectOnPoint(r, (width*0.5)-100, height*0.5);
    examplerrect = CenterRectOnPoint(r, (width*0.5)+100, height*0.5);
    r = [1 1 ST.imsize/2 ST.imsize];
    leftrect = CenterRectOnPoint(r, width*0.5-(ST.imsize/4), height*0.5);
    rightrect = CenterRectOnPoint(r, width*0.5+(ST.imsize/4), height*0.5);
    
    
    mouseloop;  % wait for subject to click mouse before experiment starts
    WaitSecs(0.5);
    
    trialcounter = 0;
    buttonisdown = 0;
    exitcode = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if trackermode == 1              % initialize eyelink
        if EyelinkInit()~= 1;
            return
        end
        el=EyelinkInitDefaults(w);      % initialise with default settings
    end
    
    % first, gather data on stable fixation with no stimulus
    
    if trackermode == 1
        % open up the file to save eyetracker data to
        Eyelink('OpenFile', strcat(E.subj,'r',num2str(runnumber),'F.edf')); 
        EyelinkDoTrackerSetup(el);                        % do calibration
        
        Screen('FillRect', w, 128);    % set to mean luminance
        drawfixation(w, width, height, 0);
        Screen('Flip', w);
        mouseloop;
        
        Screen('FillRect', w, 128);    % set to mean luminance
        
        
        drawfixation(w, width, height, 0);
        vbl = Screen('Flip', w);
        
        Eyelink('StartRecording');
        
        Screen('FillRect', w, 128);    % set to mean luminance
        drawfixation_L(w, width, height, 0)            
        Screen('Flip', w, vbl + 10 - ifi/2);        % wait 10 seconds
        
        Eyelink('ReceiveFile','F.edf', strcat('/Results/', E.subj,'/', E.subj,'r',num2str(runnumber),'F.edf'))
        Eyelink('StopRecording');                   % stop eyelink
        Eyelink('CloseFile')
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    checkeye = 4;
    while ~exitcode
        Screen('FillRect',w, ST.greylevel);
        Screen('Flip', w);
        
        clickcounter = 0;
        trialcounter = trialcounter + 1;
        condition = E.trialorder(trialcounter);
        targetcondition = E.tlistorder(trialcounter);
        
        framecounter = 0;
        orientindex = 0;
        sfindex = 0;
            %disp(sflist(sfindex))
            %disp(orlist(orient        
        
            
            %%%%%%%%%%%%%%%%%%%%%%% calibrate eye tracker again
        if trackermode == 1
            if trialcounter == checkeye;
                % open up the file to save eyetracker data to
                Eyelink('OpenFile', strcat(E.subj,'r',num2str(runnumber),'T',num2str(trialcounter),'.edf'));
                EyelinkDoTrackerSetup(el);               % do calibration
            end
        end
        Screen('FillRect',w, ST.greylevel);
        Screen('Flip', w);
        disp(atttask)
        
        if atttask == 'H';
            Screen('DrawTexture', w, targetlist_2(10), [], examplelrect, Hori(1));
            Screen('DrawTexture', w, targetlist_4(10), [], examplerrect, Hori(2));
        end
        if atttask == 'L';
            Screen('DrawTexture', w, targetlist_2(10), [], examplelrect, Lori(1));
            Screen('DrawTexture', w, targetlist_4(10), [], examplerrect, Lori(2));
        end
       
        % Wait for a mouse click to start, then start recording eye data
        drawfixation_L(w, width, height, 0);
        Screen('Flip', w);
        mouseloop;
        if trackermode == 1;
            Eyelink('StartRecording');
        end
        
        % Keeps stimulus on screen for set duration
        while framecounter < round(ST.duration*(1/ifi))
            framecounter = framecounter + 1;
            targetindex = mod(framecounter-1,targetframes)+1;
            leftframeindex = mod(framecounter-1,flankframes(1))+1;
            rightframeindex = mod(framecounter-1,flankframes(2))+1;
            if targetindex == 1
                orientindex = orientindex + 1;
                sfindex = sfindex + 1;
            end
            
            if targetcondition == 1;
                thistargetlist = targetlist1;
            elseif  targetcondition == 2;
                thistargetlist = targetlist2;
            elseif  targetcondition == 3;
                thistargetlist = targetlist3;
            elseif  targetcondition == 4;
                thistargetlist = targetlist4;
            end   
                
                
            trialorderinfo(runnumber,trialcounter,orientindex,1) = thistargetlist(orientindex,1);
            trialorderinfo(runnumber,trialcounter,orientindex,2) = thistargetlist(orientindex,2); 
            
            Screen('FillRect',w, ST.greylevel);
            
            if condition == 1
                Screen('DrawTexture', w, flanklistL(leftframeindex), [], leftrect);
            elseif condition == 2
                Screen('DrawTexture', w, flanklistR(rightframeindex), [], rightrect);
            elseif condition == 3
                Screen('DrawTexture', w, flanklistL(leftframeindex), [], leftrect);
                Screen('DrawTexture', w, flanklistR(rightframeindex), [], rightrect);
            end
                    
            if thistargetlist(sfindex,2) == 2
                Screen('DrawTexture', w, targetlist_2(targetindex), [], targetrect, thistargetlist(orientindex,1));
            else
                Screen('DrawTexture', w, targetlist_4(targetindex), [], targetrect, thistargetlist(orientindex,1));
            end
        
            lastflip = Screen('Flip', w);
            if framecounter==1
                trialstart = lastflip;
            end
            
            if useVP == 1
                if targetindex == 1
                    if framecounter ==1
                        % first trigger also contains condition code
                        Datapixx('SetDoutValues', transformindex(condition*10)); 
                        Datapixx('RegWrRd');
                    else
                        Datapixx('SetDoutValues', transformindex(1));
                        Datapixx('RegWrRd');
                    end
             trialorderinfo(runnumber,trialcounter,orientindex,3) = lastflip-trialstart;
             disp(lastflip)
                end
                if targetindex == 6  % set to zero after 3frame pairs(50ms)         
                    Datapixx('SetDoutValues', 0);
                    Datapixx('RegWrRd');
                end
            end
            
            [x,y,buttons] = GetMouse;
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(KbName('Escape'))
                exitcode = 1;
                framecounter = round(ST.duration*(1/ifi))*2;
            end
            
            if sum(buttons)
                if ~buttonisdown
                    clickcounter = clickcounter + 1;
                    clicklog(trialcounter,clickcounter) = lastflip-trialstart;
                    buttonisdown = 1;
                end
            else
                buttonisdown = 0;
            end
            
        end
        E.trialsdone(condition) = E.trialsdone(condition) + 1;
        E.tlistssdone(targetcondition) = E.tlistssdone(targetcondition) +1;

        
        if sum(E.trialsdone)==E.nreps*E.ntrials
            exitcode = 1;
        end
        
        Screen('FillRect',w, ST.greylevel);
        lastflip = Screen('Flip', w);
        
        %%%%%%% CLOSE EYETRACKER
        if trackermode
            Eyelink('ReceiveFile', strcat('E',E.subj,'r',num2str(runnumber),'T',num2str(trialcounter),'.edf'));
            Eyelink('StopRecording');                   % stop eyelink
            Eyelink('CloseFile')
        end
        
        if ~exitcode
            WaitSecs(E.ITI);
        end
    end
  
    allclicks{runnumber} = clicklog;
    alltrialorders(runnumber,:) = E.trialorder;
    alltlistorders(runnumber,:) = E.tlistorder;
    save(fname,'runnumber','allclicks','alltlistorders','alltrialorders','trialorderinfo','task');

catch
    
    lasterr

end

Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);

%Screen('Flip', w);
Screen('Close',targetlist_2);
Screen('Close',targetlist_4);
Screen('Close',flanklistL);
Screen('Close',flanklistR);

ShowCursor;

if useVP        % using a ViewPixx or ProPixx
    Datapixx('DisableVideoScanningBacklight');
end
Screen('CloseAll');
if useVP
    Datapixx('Close');
end

PsychPortAudio('Close', tr);

echo on;

end
%--------------------------------------------------------------------------
function drawfixation(w, width, height, drawframe)

ulx = width/2;
uly = height/2;

Screen('DrawLine', w, [0], ulx-4, uly, ulx+4, uly, 2);
Screen('DrawLine', w, [0], ulx, uly-4, ulx, uly+4, 2);

if drawframe
    r = [0 0 110 110];
    r = CenterRectOnPoint(r, ulx, uly);
    Screen('FrameRect', w, [0], r, 2);
end

end
%--------------------------------------------------------------------------
function drawfixation_L(w, width, height, drawframe)

ulx = width/2;
uly = height/2;

Screen('DrawLine', w, [0], ulx-10, uly, ulx+10, uly, 2);
Screen('DrawLine', w, [0], ulx, uly-10, ulx, uly+10, 2);

if drawframe
    r = [0 0 110 110];
    r = CenterRectOnPoint(r, ulx, uly);
    Screen('FrameRect', w, [0], r, 2);
end

end
%--------------------------------------------------------------------------
function mouseloop

exitcode = 0;

while exitcode==0
    [x,y,buttons] = GetMouse;
    
    if sum(buttons)>0
        exitcode = 1;
    end
end

end

%--------------------------------------------------------------------------
function output = doimagegamma(i)

% gamma corrects the stimuli before sending to Bits++
% adapted from Mark's code, DHB 29.01.08

% parameters from last gamma correct, 12/8/14 on VPixx in M16 mode
k = 0.5344;
Lmax = 101.8721;
j0 = 0.1292;
gamma = 1.9252;
%%%%%

i0 = 0;
imax = 1;                     % Bits++ always scaled between 0 and 1
imean = (i0+imax)/2;
jmax = 1;

% Eqn 2, with j set to 0, to get Lmin
Lmin = k + (Lmax-k)*(max(-j0,0)/(jmax-j0) ).^gamma; 
Lmin = max(Lmin,0);                                % ensure Lmin not <0
Lmean = (Lmin+Lmax)/2;
% desired luminance values Eqn 4
L = Lmean + (Lmax-Lmean)*(i-imean)/(imax-imean); 
% Gamma-corrected lut values, j: Eqn 3
j = ((L - k)/(Lmax-k)).^(1/gamma)*(jmax - j0) + j0; 
output = max(j,j0);                                 % Eqn 3 conditional
output = double(output);

end
%--------------------------------------------------------------------------------------------------
function output = transformindex(input)

% fixes binary inputs for EEG amp.- pins in different order from  ViewPixx
% desired numbers must be <256
% DHB 18/8/14

truebits = 2.^(2:2:24);
dn = dec2bin(input,length(truebits));
output = 0;
for m = 1:length(truebits)
    output = output + truebits(m)*str2num(dn(end-m+1));
end

end
%--------------------------------------------------------------------------
function imag1 = mkgrating(Regionsize, f, o, p, c)

%TSM; 26.6.03
% modified by DHB to make single component gratings only, scaled  -1 to 1
% f is spatial frequency, scaled as cycles per image
% o is orientation(deg), p is phase (deg relative to centre), c is contrast

p = p*pi/180;
o = o*2*pi/360;		% convert from degrees to radians
f = f/Regionsize;
x0 = ((Regionsize+1)/2);
y0 = x0;

u = f .* cos(o) * 2 * pi;
v = f .* sin(o) * 2 * pi;

imag1 = zeros(Regionsize, Regionsize);
[xx, yy] = meshgrid(1:Regionsize, 1:Regionsize);

imag1(:,:) = (c .* sin(u .*(xx-x0) + v.*(yy-y0) + p));

end
%--------------------------------------------------------------------------
function mask = make_soft_window(W,H,D)

% Mark's code for making a raised cosine window

% SYNTAX: mask = make_soft_window(W,H,[D])
% ends array 'mask' that is 1 inside  circular window, shading to 0 outside
% W, H are  width and height of whole (rectangular or square) array,in pix.
% Diameter of the soft window at half-height defaults to 0.90 units
%    where 1 unit = image width or height (whichever is smaller)
% Smoothing is by convolution with a cosine half-cycle of width 0.1 units
% Optional parameter D specifies  diameter ( relative units, range 0 -> 1)
% MAG, 27.2.04

%soft window parameters
if nargin<3, D = 0.9; end % sets default diameter to 0.9
radius = min(W*D/2,H*D/2);% radius in pixels
blur = 2*(min(W/2,H/2) - radius);  % blur half-cycle
L = blur;
X1 = [-L/2:L/2];

% 1-D blur function (applied twice, in x and y)
WinKernel = cos(X1*pi/L); % half-cycle cosine
%image coordinates - X and Y arrays
X = [1:W] - W/2;
Y = [1:H] - H/2;
xx = repmat(X,H,1);
yy = repmat(Y',1,W);

% make circular soft window
% logical 0 outside the circle,1 inside it
mask = single((xx.*xx + yy.*yy) < radius^2); 
mask = conv2(WinKernel,WinKernel,mask,'same'); 	% smooth the mask
mask = mask/max(max(mask));						% scale the mask 0-1
% figure(2);plot(X,mask(H/2,:),'r-',Y,mask(:,W/2));
mask = double(mask);
end
%-------------------------------------------------------------------------
