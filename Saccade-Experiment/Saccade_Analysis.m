
clear all

%%% This code 
%%%
%participant = ['19';'27';'28';'32';'46';'56';'66';'67';'71';'75';'89';'94'];
participant = ['01']
% Loop though each of the participants
for i = 1:size(participant,1);
    thisparticipant = participant(i,:);
    
    datalocation = '/Volumes/Untitled/Becca_new/Becca/Saccades/' ;
    
    % Format the path to each participants data file
    raweyefile = strcat(datalocation,num2str(thisparticipant),'/',num2str(thisparticipant),'Asac.edf'),datalocation,thisparticipant,thisparticipant ;
    % Using the Edf2Mat package, process the data into pixel location of
    % gaze.
    % https://github.com/uzh/edf-converter
    edf1 = Edf2Mat(raweyefile)
    
    datafile = strcat(datalocation,num2str(thisparticipant),'/',num2str(thisparticipant),'Asaccadeinfo.dat'),datalocation,thisparticipant,thisparticipant ;
    Asaccadeinfo = load(datafile);
    
    starttimes = zeros(1,6);
    
    clear starttimes;
    % Initial loop though each trial to get condition info
    for trials = 1:66;
    trialtime = Asaccadeinfo(trials,3);
    trialtime = trialtime - 3;
    starttimes(trials,1) = trialtime;
    end

% Loop though each of the 66 conditions
    for condition = 1:66;
        thiscondition = Asaccadeinfo(condition,1);
        % If its a left visual field condition
        if thiscondition == 1;
            % We're only interested in the 2 seconds after stimulus onset
            for samples = 1:2000;
                % stimulus onset
                starttimes(condition,2) = thiscondition;
                start = starttimes(condition,1);
                framestart = (1000 * start);
                nextframe = (framestart-1)+samples;
                starttimes(condition,3) = int64(framestart);
                PosX = edf1.Samples.posX(int64(nextframe));
                % If fixation comes within 1 visual degree of target
                if PosX >= 1368;
                    finalPosX = PosX;
                    finalframe = nextframe;
                    framedifms = (finalframe-framestart);
                    starttimes(condition,4) = int64(finalframe);
                    starttimes(condition,5) = int64(finalPosX);
                    starttimes(condition,6) = framedifms;
                    break
                else    
                    continue
                end

            end
        % If it's a right visual field condition
        else
            for samples = 1:2000;
                starttimes(condition,2) = thiscondition;
                start = starttimes(condition,1);
                framestart = (1000 * start);
                starttimes(condition,3) = int64(framestart);
                nextframe = (framestart-1)+samples;
                PosX = edf1.Samples.posX(int64(nextframe));
                if PosX <= 542;
                    finalPosX = PosX;
                    finalframe = nextframe;
                    framedifms = (finalframe-framestart);
                    starttimes(condition,4) = int64(finalframe);
                    starttimes(condition,5) = int64(finalPosX);
                    starttimes(condition,6) = framedifms;
                    break
                else    
                    continue
                end

            end
        end
    end



    lmean = 0;
    lcount = 0;
    lcond = zeros(1,33);
    rmean = 0;
    rcount =0;
    rcond = zeros(1,33);
    % Loop though conditions to calculate means.
    for row = 1:66;
        if starttimes(row,6) == 0;
           starttimes(row,6) = NaN ;
        elseif starttimes(row,2) == 1;
            lcond(1,row) = starttimes(row,6);
            rcond(1,row) = NaN;
            lcount = lcount + 1;
            lmean = lmean + starttimes(row,6);
        else
            rcond(1,row) = starttimes(row,6);
            lcond(1,row) = NaN;
            rcount = rcount + 1;
            rmean = rmean + starttimes(row,6);
        end
    end

% Calculate Standard deviation of 
lstd = nanstd(lcond(1,:));
rstd = nanstd(rcond(1,:));

% Mean reaction time for L and R
starttimes(68,6) = lmean/lcount;
starttimes(69,6) = rmean/rcount;
% STD reaction time for L and R
starttimes(70,6) = lstd;
starttimes(71,6) = rstd;


% Save file out
savefile = strcat('/Users/rek514/Desktop/PhD/EEG Experiment/Saccade_results/',thisparticipant,'Antisaccadeanalysis.dat');
save(savefile);
dlmwrite(savefile, starttimes, 'delimiter', '\t');
end

