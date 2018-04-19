%% FitnessMemoryAnalysis.m 
% Load in data from Guiseppe's Spatial Configuration Task after unzipping
% folder you download from them and renaming zip file 'ExerciseMemoryData'

% Data from Matt Pontifex fitness study, 

% Compare with variable in database file fitness percentile

%path = M:\Data\FitnessMemory\Spatial Configuration Task Data 2018 04 18\
%filename = FR100_Spatial Configuration Task.csv

% Participants with VO2 Confidence 0 
% - FR125
% - FR116

%% load the data 
clear all
close all
clc

%path where the data is:
main_path = 'M:\Data\FitnessMemory\Spatial Configuration Task Data 2018 04 18\';
filenames = ls(main_path); %get list of all the files
filenames(1:2,:) = []; %remove dots in first two spaces

%% subs in each condition
n_subs = size(filenames,1); %how many subjects
subs = zeros(1,n_subs);
for i_sub = 1:n_subs
    subs(i_sub) = str2double(filenames(i_sub,3:5)); %get a list of the actual subject numbers
end

%% some settings
n_trials = 60;
windowsize = 10; %size of window to smooth data >2 wanter to use pick an even number
blocksize = 10;  %how much trials to include in block
n_blocks = n_trials/blocksize; %how many blocks - make these all divisible

%trials to critereon
window_size = 5;
correct_trials = 4; %picked these because the most people reached this level (least were given value of 60)

%% make a place to save all the data for all these subjects 
% (subjects x trials)
all_hitdata = nan(n_subs,n_trials);
all_RTdata = nan(n_subs,n_trials);
all_window_mean = nan(n_subs,n_trials);
all_window_RT = nan(n_subs,n_trials);
all_block_mean = nan(n_subs,n_blocks);
trials2crit = zeros(1,n_subs);
slope_out = zeros(1,n_subs); %linear fit of the window average 
intercept_out = zeros(1,n_subs);

%% loop through each subjects
for i_sub = 1:n_subs %subjects 
    
    this_sub_num = subs(i_sub); %get the number from the list above    
    fprintf(['Loading data for subject ' num2str(this_sub_num) '. \n']);
    [num,txt,raw] = xlsread([main_path 'FR' num2str(this_sub_num) '_Spatial Configuration Task.csv']); %read in the data as numbers, txt, and both
    %saves the numbers in one variable, and the text in another, 
    
    %% find the hits and RT rows and pull out data
    %so Matt's data files are funny, and are not getting parsed correctly,
    %the hits line has the first trial in the first column with the label,
    %and the last column (trial) has a quotation mark
    %Hits:"0 1 1 0 0 ... 0 0 0"
    
    hitdata = NaN(1,n_trials); 
    hitdata(2:n_trials-1) = num(1,:); %get all the hit data that worked
    hitdata(1) = str2double(txt{11,1}(end));  %get the first trial from the label Hits:"x
    %this looks funny {11,1}(end) because it is a string inside a cell, so the cell
    %address is {11,1} 11th row, 1st column, the string address is (end),
    %the last character
    hitdata(n_trials) = str2double(txt{11,n_trials}(1)); %get the last trial from the last part

    RTdata = NaN(1,n_trials);
    RTdata(2:n_trials-1) = num(2,:); %get the RT data
    RTdata(1) = str2double(txt{12,1}(6:end));
    RTdata(n_trials) = str2double(txt{12,n_trials}(1:end-1));
    
    %save into variable with other subjects
    all_hitdata(i_sub,:) = hitdata;
    all_RTdata(i_sub,:) = RTdata;
    
    %% trials to crit
 
    for i_win=window_size:n_trials % go through each window
        correct_sum = sum(hitdata(i_win-(window_size-1):i_win)); %find the trials correct in window
        if correct_sum == correct_trials || i_win == n_trials  %if reach criterion
            trials2crit(i_sub) = i_win; %assign value
            break %break out of loop (don't keep checking)
        end
    end
    
    
    %% plot the data
%     figure;
%     suptitle(['Sub' num2str(this_sub_num)]);
%     subplot(4,1,1);
%         plot(1:n_trials,hitdata); %plot raw hit data
%             title('Hits');
%             xlabel('Trial Number');
%             ylabel('Accuracy');
%             ylim([-.5 1.5]);  %adjust to fit everything in
%             xlim([1 60]);
% 
%     subplot(4,1,2);
%         plot(1:n_trials,RTdata); %plot raw RT data
%             title('RT');
%             xlabel('Trial Number');
%             ylabel('ReactionTime');
%             xlim([1 60]);
    
    
    %% Moving window average
    window_mean = nan(1,60); %reset the values
    for i_trial = 1:n_trials-windowsize  %go through each window and average and assign value %skips last values so no overhang
        window_mean(i_trial+ceil(windowsize/2)) = mean(hitdata(i_trial:i_trial+windowsize)); %skip first five places to save
    end
    
    window_RT = nan(1,60); %reset the values
    for i_trial = 1:n_trials-windowsize  %go through each window and average and assign value %skips last values so no overhang
        window_RT(i_trial+ceil(windowsize/2)) = mean(RTdata(i_trial:i_trial+windowsize)); %skip first five places to save
    end
     
%     subplot(4,1,3);
%         plot(1:n_trials,window_mean); %plot raw hit data
%             title(['Moving window size ' num2str(windowsize)]);
%             xlabel('Trial Number');
%             ylabel('Accuracy');
%             xlim([1 60]);
%             line([0 60],[.2 .2],'color','k')
    
    all_window_mean(i_sub,:) = window_mean; %save each session
    all_window_RT(i_sub,:) = window_RT; %save each session
    
    
    %% Fit linear trend to window_mean
    
    x = 1:n_trials; 
    is_vals = find(~isnan(window_mean)==1); %only the values not NaNs, edges are nans from the moving window
    P = polyfit(x(is_vals),window_mean(is_vals),1);
    yfit = P(1)*x+P(2); %fit a linear function to the values
   
    slope_out(i_sub) = P(1);
    intercept_out(i_sub) = P(2);
    
    % to see the fit
%     figure;
%     plot(x(is_vals),window_mean(is_vals))
%     hold on;
%     plot(x,yfit,'r-.');
    
    %% blocked trials
    
    block_mean = zeros(1,n_blocks);
    for i_block = 1:n_blocks
        block_mean(i_block) = mean(hitdata(  1+((i_block-1)*blocksize):     i_block*blocksize)); %save the block mean into variable
    end
    
%     subplot(4,1,4);
%         plot(1:n_blocks,block_mean); %plot raw hit data
%             title(['Block size ' num2str(blocksize)]);
%             xlabel('Block Number');
%             ylabel('Accuracy');
%             xlim([1 n_blocks]);
%             ylim([0 1]);
%             line([0 n_blocks],[.2 .2],'color','k')
%             set(gca,'XTick',[1:n_blocks])

    all_block_mean(i_sub,:) = block_mean; %save each session
    
    
    
end %loop for subjects







%% compute some stats
% all_window_mean is a subjects x trials

figure;

grand_window_mean = mean(all_window_mean,1);
grand_window_mean_se = std(all_window_mean,[],1)/sqrt(n_subs);

%plot it
subplot(1,3,1);
plot(1:n_trials,grand_window_mean,'k',...
    1:n_trials,grand_window_mean+grand_window_mean_se,'r',...
    1:n_trials,grand_window_mean-grand_window_mean_se,'r');
    title('Smoothed Hit Rate over trials');
    xlabel('Trial Number');
    ylabel('Proportion correct');

%block means
grand_block_mean = mean(all_block_mean,1);
grand_block_mean_se = std(all_block_mean,[],1)/sqrt(n_subs);

subplot(1,3,2);
plot(1:n_blocks,grand_block_mean,'k',...
    1:n_blocks,grand_block_mean+grand_block_mean_se,'r',...
    1:n_blocks,grand_block_mean-grand_block_mean_se,'r');
    title('Smoothed Hit Rate over blocks');
    xlabel('Block Number');
    ylabel('Proportion correct');
    

%get out the smoothed RT and average over subs
grand_window_RT = mean(all_window_RT,1);
grand_window_RT_se = std(all_window_RT,[],1)/sqrt(n_subs);

%plot it
subplot(1,3,3);
plot(1:n_trials,grand_window_RT,'k',...
    1:n_trials,grand_window_RT+grand_window_RT_se,'r',...
    1:n_trials,grand_window_RT-grand_window_RT_se,'r');    
    title('Smoothed RT over trials');
    xlabel('Trial Number');
    ylabel('RT (ms)');
    
    
%% trials 2 crit
%trials2crit is a subjects vector 

figure;
    subplot(1,3,3);
    hist(trials2crit,20);
    title('Trials2Crit');
    
    subplot(1,3,1);
    hist(slope_out,20);
    title('Slope');
    
    subplot(1,3,2);
    hist(intercept_out,20);
    title('Intercept');
   
    
    
%% load in xls file
database_path = 'M:\Data\FitnessMemory\';
[NUM,TXT,RAW]= xlsread([database_path 'FiRe_DemographicsDataEntry_2018_04_10.xls']);

%grab subject numbers (not the same number)
sub_names = TXT(2:end,1);
n_db_subs = length(sub_names);
database_subs = zeros(1,n_db_subs);
for i_item = 1:n_db_subs
    database_subs(i_item) = str2double(sub_names{i_item}(3:end));
end

%grab fitness info
% vo2 max is column 95 of NUM, vo2 max percent is 96
database_fitness = NUM(:,96);

%% loop through task subjects and lookup fitness

task_fitness = zeros(1,n_subs);
for i_sub = 1:n_subs
    i_db = find(database_subs == subs(i_sub));
    if isempty(i_db)
        task_fitness(i_sub) = NaN;
    else
        task_fitness(i_sub) = i_db;
    end
end

%% remove subjects data with no vo2 max confidence
% Participants with VO2 Confidence 0 
% - FR125
% - FR116

remove_subs = [125, 116];
for i_sub = 1:n_subs
    for i_remsub = 1:length(remove_subs)
        if subs(i_sub) == remove_subs(i_remsub)
            task_fitness(i_sub) = NaN;
        end
    end
end



%% correlate fitness with 3 measures
%intercept_out, slope_out, trials2crit
%with
%task_fitness
%task_fitness has some NaN's of missing values

%don't use subs with no database entry
fitness = task_fitness(~isnan(task_fitness));
intercept = intercept_out(~isnan(task_fitness));
slope = slope_out(~isnan(task_fitness));
TTC = trials2crit(~isnan(task_fitness));


figure;
subplot(1,3,1);
    scatter(fitness,slope);
        P = polyfit(fitness,slope,1);
        yfit = P(1)*fitness+P(2); %fit a linear function to the value
   hold on;
    plot(fitness,yfit,'r-.');
    xlabel('VO2Max %');
    ylabel('LearningSlope');
    xlim([0 100]);
    [r, p] = corr(fitness',slope');
    title(['Slope correlation = ' num2str(round(r,4)) '. pvalue = ' num2str(round(p,4))]);

subplot(1,3,2);
    scatter(fitness,intercept);
        P = polyfit(fitness,intercept,1);
        yfit = P(1)*fitness+P(2); %fit a linear function to the value
   hold on;
    plot(fitness,yfit,'r-.');
    title('Intercept');
    xlabel('VO2Max %');
    ylabel('LearningIntercept');
    xlim([0 100]);
   [r, p] = corr(fitness',intercept');
    title(['Intercept correlation = ' num2str(round(r,4)) '. pvalue = ' num2str(round(p,4))]);

subplot(1,3,3);
    scatter(fitness,TTC);
        P = polyfit(fitness,TTC,1);
        yfit = P(1)*fitness+P(2); %fit a linear function to the value
   hold on;
    plot(fitness,yfit,'r-.');
    title('TTC');
    xlabel('VO2Max %');
    ylabel('LearningTTC');
    xlim([0 100]);
   [r, p] = corr(fitness',TTC');
    title(['TTC correlation = ' num2str(round(r,4)) '. pvalue = ' num2str(round(p,4))]);









   



