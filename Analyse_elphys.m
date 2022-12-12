% function [elphys_complete,elphys_cleaned] = Analyse_elphys(make_figures)
% Decent as it is.

%% Parameters

% When reviewing figures
if ~exist("pausetime","var")
    pausetime = 0.5;          % how long should the script pause before moving on when make_figures is turned on?
end

% Sampling rate
sr = 20000;             % indicate samples per second

% For paired-pulse experiments
prepulse_time = 7;      % indicate ms
postpulse_time = 20;    % indicate ms

% For LTP experiments
baseline_sweeps = 90;

posttet_sweeps = 300;

sweeps_to_mean = 30;

tetanisations_in_LTP = 4;

baseline_time = 15;

sweep_interval = 10;    % indicate seconds

minimal_accepted_LTP = 1.1;

allowed_baseline_noise = 0.1;

% Criteria

minimal_accepted_amplitude = 0.2;

% To remove files before a certain date
DateStringStart = '18-06-2022';

% To remove files after a certain date
DateStringEnd   = '10-08-2022';

%% Conditional parameters

if ~exist("make_figures","var")
    make_figures = false;
end

if make_figures == true
    rising_figs_plonx = true;
else
    rising_figs_plonx = false;
end

sweeps_in_longest_experiment = baseline_sweeps+posttet_sweeps;

prepulse_samples = prepulse_time * (sr/1000);
postpulse_samples = postpulse_time * (sr/1000);
LTP_baseline_means = baseline_sweeps / sweeps_to_mean;
LTP_minutes_per_mean = sweeps_to_mean * sweep_interval / 60;

if mod(LTP_baseline_means,1) ~= 0
    disp("The basline period cannot be divided into equal intervals with the chosen sweeps_to_mean")
    return
end

%% Metadata extraction
% Create a list of the files
pathway_ = 'C:/Users/nikol/Documents/SDC_Teaching/Masters_Thesis/Nordvej/Data/Lettensenter/';
cd(pathway_)
elphys_files = dir(strcat('*/*.abf'));
elphys_files = struct2table(elphys_files);

% Remove files outside the chosen time span
formatInStart   = 'dd-mm-yyyy';
startdate       = datenum(DateStringStart,formatInStart);
formatInEnd     = 'dd-mm-yyyy';
enddate         = datenum(DateStringEnd,formatInEnd);

relevant_files_by_time = elphys_files.datenum(:) >= startdate & elphys_files.datenum(:) <= enddate;
elphys_files = elphys_files(relevant_files_by_time,:);

% Start creating the table
files_count = length(elphys_files.datenum);
% To look at a single file
if ~exist("lookat","var")
    files_for_analysis = 1:files_count;
else
    files_for_analysis = lookat;
end

%% Defining output table variables
Experimental_Setup      = strings(files_count,1);
Section_no              = zeros(files_count,1);
Experimenter            = strings(files_count,1);
Electrode_Resistance    = zeros(files_count,1);
Stimulation_Duration    = zeros(files_count,1);
Stimulation_Intensity   = zeros(files_count,1);
Experiment              = strings(files_count,1);
Animal_no               = zeros(files_count,1);
Event_no                = zeros(files_count,1);
Init_slp_idx            = zeros(files_count,2);
Final_slp_idx           = zeros(files_count,2);

no_mean_sweeps = floor(sweeps_in_longest_experiment/sweeps_to_mean);
means_offset = LTP_minutes_per_mean/2;

data_size = [files_count,floor(no_mean_sweeps)];

Slopes                  = zeros(data_size);
LTP_minutes             = repmat([LTP_minutes_per_mean:LTP_minutes_per_mean:LTP_minutes_per_mean*no_mean_sweeps]-means_offset,[files_count 1]);
Norm_slopes             = zeros(data_size);
Amplitudes              = zeros(data_size);
Norm_amplitudes         = zeros(data_size);
PV_amplitudes           = zeros(data_size);
Rel_slopes              = zeros(data_size);
Norm_rel_slopes         = zeros(data_size);
Rel_amplitudes          = zeros(data_size);
Norm_rel_amplitudes     = zeros(data_size);

elphys_files = addvars(elphys_files,Animal_no,Experimental_Setup,Section_no,Experimenter,Electrode_Resistance,...
    Experiment,Stimulation_Duration,Stimulation_Intensity,Event_no,Init_slp_idx,Final_slp_idx,Slopes,LTP_minutes,...
    Norm_slopes,PV_amplitudes,Rel_slopes,Norm_rel_slopes,Rel_amplitudes,Norm_rel_amplitudes);

table_vars = string(elphys_files.Properties.VariableNames);
for i = 1:length(table_vars)
    clear(table_vars(i))
end
clear table_vars


for file_no = files_for_analysis
    fprintf('Processing file number %.0f \n',file_no)

    if elphys_files.name{file_no}(1:4) == 'Conc'
        continue
    end

    folder_string = elphys_files.folder{file_no};

    sep_ind = find(folder_string == '\');

    setup_note = folder_string(sep_ind(end)+1:end);

    elphys_files.Experimental_Setup(file_no) = setup_note;

    file_location = strcat(elphys_files.folder{file_no},'\',elphys_files.name{file_no});

    [d,si,h] = abfload(file_location);
    event_no = size(d,3);
    elphys_files.Event_no(file_no) = event_no;

    if event_no == 90
        elphys_files.Experiment(file_no) = "Baseline";
    elseif event_no == 1 && (elphys_files.Experiment(file_no-1) == "Baseline" || elphys_files.Experiment(file_no-1) == "Tetanisation")
        elphys_files.Experiment(file_no) = "Tetanisation";
    elseif event_no > 6 && elphys_files.Experiment(file_no-1) == "Tetanisation"
        elphys_files.Experiment(file_no) = "Post-tetanisation";
    elseif event_no == 6
        elphys_files.Experiment(file_no) = "Paired pulse";
    elseif event_no > 1 && event_no < 6
        elphys_files.Experiment(file_no) = "Test paired pulse";
    elseif event_no > 6 && event_no < 90
        elphys_files.Experiment(file_no) = "Abandoned baseline";
    end

    if elphys_files.Experimental_Setup(file_no) == "Outer_setup"
        elphys_files.Experimenter(file_no) = "NAK";
    end

    if elphys_files.Experimental_Setup(file_no) == "Outer_setup"
        active_channel = 2;
    else
        active_channel = 1;
    end

    % Defining a time vector
    time_ms = 1:length(d(:,1,1));
    time_ms = time_ms*1000 ./ sr;

    % Calculate mean sweep
    the_sweeps = squeeze(d(:,active_channel,:));
    mean_sweep = squeeze(mean(the_sweeps,2));
    mean_sweep = mean_sweep - mean_sweep(1);

    if make_figures == true
        if exist("summary_fig","var")
            isFigureHandle = ishandle(summary_fig) && strcmp(get(summary_fig,'type'),'figure');
            if ~isFigureHandle
                summary_fig = figure("Name","Summary of experiment");
            end
        else
            summary_fig = figure("Name","Summary of experiment");
        end
    end


    %% Automated analysis for paired-pulse experiments

    if elphys_files.Experiment(file_no) == "Paired pulse"

        % Find the first pulse based on stimulus artefact
        raw_pulse1 = mean_sweep(1:50*20);
        t_pulse1 = min(time_ms(raw_pulse1 == min(raw_pulse1)));
        id_pulse1 = find(time_ms == t_pulse1);

        % Find the second pulse based on stimulus artefact
        raw_pulse2 = mean_sweep(id_pulse1+40*20+1:end);
        raw_pulse2 = raw_pulse2 - raw_pulse2(1);
        t_pulse2 = min(time_ms(raw_pulse2 == min(raw_pulse2))) + t_pulse1 + 40;
        id_pulse2 = find(time_ms == t_pulse2);

        % Extract the pulse neighborhoods

        %    loc_pulse1 = time_ms < t_pulse1 + postpulse_time & time_ms >= t_pulse1 - prepulse_time;
        pulse1 = mean_sweep(id_pulse1-prepulse_samples:id_pulse1+postpulse_samples);
        pulse1 = pulse1 - mean(pulse1((prepulse_time-1.5)*20:(prepulse_time-0.5)*20));

        %    loc_pulse2 = time_ms < t_pulse2 + postpulse_time & time_ms >= t_pulse2 - prepulse_time;
        pulse2 = mean_sweep(id_pulse2-prepulse_samples:id_pulse2+postpulse_samples);
        pulse2 = pulse2 - mean(pulse2((prepulse_time-1.5)*20:(prepulse_time-0.5)*20));

        collected_mean_sweeps = [pulse1,pulse2];

        time_pulse = linspace(-prepulse_time,postpulse_time,length(pulse1));

        inst_slope_pulse1 = [pulse1(2:end) - pulse1(1:end-1) ; 0];
        inst_slope_pulse2 = [pulse2(2:end) - pulse2(1:end-1) ; 0];
        inst_slope_mat = [inst_slope_pulse1,inst_slope_pulse2] .* (sr/1000);

% May be useful for improving find_rising_phase
%         for inst_cols =  1:width(inst_slope_mat)
%             for inst_rows = 5:(length(inst_slope_mat)-4)
%                 inst_slope_mat(inst_rows,inst_cols) = mean(inst_slope_mat(inst_rows-4:inst_rows+4,inst_cols));
%             end
%         end

        [neg_slope_agreement,neg_slope_mat,selected_idx,init_idx] = find_rising_phase(inst_slope_mat,time_pulse,rising_figs_plonx);
        %         neg_slope_pulse1 = inst_slope_pulse1 < 0;
        %         neg_slope_pulse2 = inst_slope_pulse2 < 0;
        %         neg_slope_agreement = neg_slope_pulse1 .* neg_slope_pulse2;

        slope_start = time_pulse(min(selected_idx));
        slope_end = time_pulse(max(selected_idx));
        %
        % Amplitudes
        amp_idx = [selected_idx; linspace(max(selected_idx)+1,max(selected_idx)+4*20,4*20)'];
        mean_amplitudes = min(collected_mean_sweeps(amp_idx,:));
        norm_amplitudes = mean_amplitudes ./ mean_amplitudes(1);

        % Slope analysis
        sloperow = zeros(1,width(collected_mean_sweeps));
        for i = 1:width(collected_mean_sweeps)
            slope = polyfit(time_pulse(selected_idx),collected_mean_sweeps(selected_idx,i),1);
            sloperow(1,i) = slope(1);
        end

        if elphys_files.Experiment(file_no) == "Post-tetanisation"
            norm_sloperow = sloperow ./ mean(sloperow(ceil(LTP_baseline_means*1/3)+1:LTP_baseline_means));
        elseif elphys_files.Experiment(file_no) == "Paired pulse"
            norm_sloperow = sloperow ./ sloperow(1);
        end

        % Prevolley amplitudes
        PV_search_start = find(time_pulse == 0.5);
        PV_search_end = selected_idx(1);
        PV_amplitudes = min(collected_mean_sweeps(PV_search_start:PV_search_end,:));

        % Slopes relative to PV ampltudes
        rel_slopes = sloperow ./ PV_amplitudes;

        % Slopes relative to PV amplitudes, normalised
        norm_rel_slopes =  rel_slopes ./ mean(rel_slopes(1));

        % Amplitudes relative to PV amplitudes
        rel_amplitudes = mean_amplitudes ./ PV_amplitudes;

        % Amplitudes relative to PV amplitudes, normalised
        norm_rel_amplitudes = rel_amplitudes ./ mean(rel_amplitudes(1));



        % Making figures if make_figures is true
        if make_figures == true

            slope_marker = [-2 2];
            mark_slope_start = [slope_start slope_start];
            mark_slope_end = [slope_end slope_end];

            subplot(2,1,1)
            plot(time_pulse,collected_mean_sweeps,LineWidth=1)
            hold on;
            plot(mark_slope_start,slope_marker,Color=[0 0 0],LineStyle="--")
            plot(mark_slope_end,slope_marker,Color=[0 0 0],LineStyle="-")
            yline(mean(mean_amplitudes),LineStyle=":")
            yline(mean(PV_amplitudes),LineStyle="-.")
            legend(["Pulse 1","Pulse 2","Slope start","Slope end","Mean amplitude","Mean PV"])
            hold off;
            ylim([mean_amplitudes(1)*1.5 0.1])
            xlim([-prepulse_time+5 postpulse_time-10])
            title(["Mean sweeps of file ", elphys_files.name(file_no),"Paired pulse"])
            xlabel("Time after stimulus peak, ms")
            ylabel("Potential, mV")

            subplot(2,1,2)
            plot(time_pulse,inst_slope_mat,Color=[0.75 0.75 0.75])
            hold on;
            plot(time_pulse,mean(inst_slope_mat,2),Color=[0 0 0])
            plot(time_pulse,neg_slope_agreement*0.1,LineWidth=1,Color=[1 0 0])
            hold off;
            ylim([-1 1])
            xlim([-prepulse_time+5 postpulse_time-10])
            title("First derivative with slope interval indicated")
            xlabel("Time after stimulus peak, ms")
            ylabel("Instantaneous slope, mV/ms")
        end
    end

    %% Automated analysis for LTP experiments
    if elphys_files.Experiment(file_no) == "Post-tetanisation"
        backtrack = file_no-tetanisations_in_LTP-1;
        if elphys_files.Experiment(backtrack) ~= "Baseline"
            fprintf('\n Post-tetanisation could not be matched to a baseline recording.')
            continue
        end
        loc_baseline = strcat(elphys_files.folder{backtrack},'\',elphys_files.name{backtrack});
        [baseline_d,baseline_si,baseline_h] = abfload(loc_baseline);
        % d, si, and h are already extracted from the current file

        ltp_post_sweeps = squeeze(d(:,active_channel,:));
        ltp_pre_sweeps = squeeze(baseline_d(:,active_channel,:));
        ltp_sweeps = [ltp_pre_sweeps, ltp_post_sweeps];

        ltp_mean_sweeps = zeros(length(ltp_sweeps),floor(width(ltp_sweeps)/sweeps_to_mean));
        for i = 1:width(ltp_mean_sweeps)
            working_mean = mean(ltp_sweeps(:,i*sweeps_to_mean+1-sweeps_to_mean:i*sweeps_to_mean),2);
            ltp_mean_sweeps(:,i) = working_mean;
        end

        peak_idx = find(ltp_mean_sweeps == min(ltp_mean_sweeps),1);

        ltp_sweeps_neighborhood = ltp_mean_sweeps(peak_idx-prepulse_samples:peak_idx+postpulse_samples,:);

        for i = 1:width(ltp_sweeps_neighborhood)
            sweepi = ltp_sweeps_neighborhood(:,i);
            sweepi = sweepi - mean(sweepi((prepulse_time-1.5)*20:(prepulse_time-0.5)*20));
            ltp_sweeps_neighborhood(:,i) = sweepi;
        end

        collected_mean_sweeps = ltp_sweeps_neighborhood;

        inst_slope_mat = [ltp_sweeps_neighborhood(2:end,:) - ltp_sweeps_neighborhood(1:end-1,:) ; zeros(1,width(ltp_sweeps_neighborhood))] .*(sr/1000);

        if make_figures == true
            rising_figs_plonx = true;
        end
        [neg_slope_agreement,neg_slope_mat,selected_idx,init_idx] = find_rising_phase(inst_slope_mat,time_pulse,rising_figs_plonx);

        slope_start = time_pulse(min(selected_idx));
        slope_end = time_pulse(max(selected_idx));

        % Amplitudes
        amp_idx = [selected_idx; linspace(max(selected_idx)+1,max(selected_idx)+4*20,4*20)'];
        mean_amplitudes = min(collected_mean_sweeps(amp_idx,:));

        norm_amplitudes = mean_amplitudes ./ mean(mean_amplitudes(floor(LTP_baseline_means*1/3)+1:LTP_baseline_means));

        % Slope analysis
        sloperow = zeros(1,width(collected_mean_sweeps));
        for i = 1:width(collected_mean_sweeps)
            slope = polyfit(time_pulse(selected_idx),collected_mean_sweeps(selected_idx,i),1);
            sloperow(1,i) = slope(1);
        end

        if elphys_files.Experiment(file_no) == "Post-tetanisation"
            norm_sloperow = sloperow ./ mean(sloperow(floor(LTP_baseline_means*1/3)+1:LTP_baseline_means));
        elseif elphys_files.Experiment(file_no) == "Paired pulse"
            norm_sloperow = sloperow ./ sloperow(1);
        end

        % Prevolley amplitudes
        PV_search_start = find(time_pulse == 0.5);
        PV_search_end = selected_idx(1);
        PV_amplitudes = min(collected_mean_sweeps(PV_search_start:PV_search_end,:));

        % Slopes relative to PV ampltudes
        rel_slopes = sloperow ./ PV_amplitudes;

        % Slopes relative to PV amplitudes, normalised
        norm_rel_slopes =  rel_slopes ./ mean(rel_slopes(floor(LTP_baseline_means*1/3)+1:LTP_baseline_means));

        % Amplitudes relative to PV amplitudes
        rel_amplitudes = mean_amplitudes ./ PV_amplitudes;

        % Amplitudes relative to PV amplitudes, normalised
        norm_rel_amplitudes = rel_amplitudes ./ mean(rel_amplitudes(floor(LTP_baseline_means*1/3)+1:LTP_baseline_means));


        % Figures
        if make_figures == true
            slope_marker = [-2 2];
            mark_slope_start = [slope_start slope_start];
            mark_slope_end = [slope_end slope_end];

            subplot(2,1,1)
            plot(time_pulse,ltp_sweeps_neighborhood)
            ylim([mean_amplitudes(1) 0.1])
            xlim([-prepulse_time postpulse_time])
            hold on;
            plot(mark_slope_start,slope_marker,Color=[0 0 0],LineStyle="--")
            plot(mark_slope_end,slope_marker,Color=[0 0 0],LineStyle="--")
            yline(mean(mean_amplitudes)*1.5,":")
            yline(mean(PV_amplitudes),"-.")
            hold off;
            title(["Mean sweeps of file ", elphys_files.name(file_no),"LTP"])
            xlabel("Time after stimulus peak, ms")
            ylabel("Potential, mV")

            subplot(2,1,2)
            plot(time_pulse,inst_slope_mat,Color=[0.75 0.75 0.75])
            hold on;
            plot(time_pulse,mean(inst_slope_mat,2),Color=[0 0 0])
            plot(time_pulse,neg_slope_agreement*0.1,LineWidth=1,Color=[1 0 0])
            hold off;
            ylim([-1 1])
            xlim([-prepulse_time postpulse_time])
            title("First derivative with slope interval indicated")
            xlabel("Time after stimulus peak, ms")
            ylabel("Instantaneous slope, mV/ms")
        end
    end

    %
    %     %% Get amplitudes
    %     amp_idx = [selected_idx; linspace(max(selected_idx)+1,max(selected_idx)+4*20,4*20)'];
    %     mean_amplitudes = min(collected_mean_sweeps(amp_idx,:));
    %
    %     if elphys_files.Experiment(file_no) == "Post-tetanisation"
    %         norm_amplitudes = mean_amplitudes ./ mean(mean_amplitudes(ceil(LTP_baseline_means*2/3):LTP_baseline_means));
    %     elseif elphys_files.Experiment(file_no) == "Paired pulse"
    %         norm_amplitudes = mean_amplitudes ./ mean_amplitudes(1);
    %     end

    %% Collecting all sweeps
    % use probe_dataset, collect all mean sweeps throughout the analysis

    %% Write statistics to table
    elphys_files.Init_slp_idx(file_no,:)                            = init_idx;
    elphys_files.Final_slp_idx(file_no,:)                           = [min(selected_idx) max(selected_idx)];
    elphys_files.Slopes(file_no,1:width(sloperow))                  = sloperow;
    elphys_files.Norm_slopes(file_no,1:width(norm_sloperow))        = norm_sloperow;
    elphys_files.Amplitudes(file_no,1:width(mean_amplitudes))       = mean_amplitudes;
    elphys_files.Norm_amplitudes(file_no,1:width(norm_amplitudes))  = norm_amplitudes;
    elphys_files.PV_amplitudes(file_no,1:width(PV_amplitudes))      = PV_amplitudes;
    elphys_files.Rel_slopes(file_no,1:width(rel_slopes))            = rel_slopes;
    elphys_files.Norm_rel_slopes(file_no,1:width(norm_rel_slopes))  = norm_rel_slopes;
    elphys_files.Rel_amplitudes(file_no,1:width(rel_amplitudes))    = rel_amplitudes;
    elphys_files.Norm_rel_amplitudes(file_no,1:width(norm_rel_amplitudes))  = norm_rel_amplitudes;

    %% Pause to view progress

    if make_figures == true && (elphys_files.Experiment(file_no) == "Post-tetanisation" || elphys_files.Experiment(file_no) == "Paired pulse")
        pause()
    end
end

%% Adding data from external files

Mousterkey = readtable('C:\Users\nikol\Documents\SDC_Teaching\Masters_Thesis\Mouster_key.xlsx');
Masterkey = readtable('C:\Users\nikol\Documents\SDC_Teaching\Masters_Thesis\Elphys_inventory.xlsx');


for i = 1:height(elphys_files)
    the_filename = string(elphys_files.name{i});
    file_row = find(string(Masterkey.name) == the_filename);

    elphys_files.Experimenter(i) = Masterkey.Experimenter(file_row);
    elphys_files.Section_no(i) = Masterkey.Section_no(file_row);
    elphys_files.Electrode_Resistance(i) = Masterkey.Electrode_resistance_M__(file_row);
    elphys_files.Animal_no(i) = Masterkey.animal_no(file_row);
    elphys_files.Stimulation_Duration(i) = Masterkey.Stimulation_duration__s_(file_row);
    elphys_files.Stimulation_Intensity(i) = Masterkey.Stimulation_intensity_V_(file_row);
    elphys_files.Measure(i) = Masterkey.Measure(file_row);
    elphys_files.Month(i) = Masterkey.Month(file_row);
    elphys_files.Day(i)     = Masterkey.Day_of_month(file_row);
    elphys_files.manual_ok(i)  = Masterkey.manual_ok(file_row);
    elphys_files.popspike_error(i) = Masterkey.pops(file_row);
    elphys_files.manual_low_signal(i) = Masterkey.manual_low_signal(file_row);
end

for i = 1:height(elphys_files)
    animal_ID = elphys_files.Animal_no(i);
    lur = find(Mousterkey.Animal_ == animal_ID);
    Startage(i,:) = Mousterkey.StartAge(lur);
    Genotype(i,:) = categorical(Mousterkey.Genotype(lur));
    Sex(i,:)  =   categorical(Mousterkey.Sex(lur));
    Treatment(i,:) = categorical(Mousterkey.Treatment(lur));
    Cage(i,:) = categorical(Mousterkey.Cage_(lur));
end

elphys_files.Slopes(elphys_files.Slopes == 0) = nan;
elphys_files.Norm_slopes(elphys_files.Norm_slopes == 0) = nan;
elphys_files.Amplitudes(elphys_files.Amplitudes == 0) = nan;
elphys_files.Norm_amplitudes(elphys_files.Norm_amplitudes == 0) = nan;
elphys_files.PV_amplitudes(elphys_files.PV_amplitudes == 0) = nan;
elphys_files.Rel_slopes(elphys_files.Rel_slopes == 0) = nan;
elphys_files.Norm_rel_slopes(elphys_files.Norm_rel_slopes == 0) = nan;
elphys_files.Rel_amplitudes(elphys_files.Rel_amplitudes == 0) = nan; 
elphys_files.Norm_rel_amplitudes(elphys_files.Norm_rel_amplitudes == 0) = nan;

elphys_files.Rel_amplitudes = elphys_files.Amplitudes ./ elphys_files.PV_amplitudes;
elphys_files.Rel_slopes = elphys_files.Slopes ./ elphys_files.PV_amplitudes;

elphys_complete = addvars(elphys_files,Startage,Genotype,Sex,Treatment,Cage);

elphys_complete.Stim_strength = elphys_complete.Stimulation_Duration/1000000 .* elphys_complete.Stimulation_Intensity;

clear Startage
clear Genotype
clear Sex
clear Treatment
clear Cage

%% Completing the table with animal data and manual selections

% What to remove

useless_rows = false(height(elphys_complete),1);
by_file = count(string(elphys_complete.name),"Conc") == 1;
useless_rows = useless_rows | by_file;
by_experiment = (elphys_complete.Experiment ~= "Post-tetanisation" & elphys_complete.Experiment ~= "Paired pulse");
useless_rows = useless_rows | by_experiment;
by_stim = (isnan(elphys_complete.Stimulation_Duration) | isnan(elphys_complete.Stimulation_Intensity));
by_stim = by_stim | elphys_complete.Stimulation_Duration == 0 | elphys_complete.Stimulation_Intensity == 0;
useless_rows = useless_rows | by_stim;

auto_ok = ~useless_rows;

elphys_complete = addvars(elphys_complete,auto_ok);

% Restrictions by signal REVISE: LTP_success
last_baseline_means = ceil(LTP_baseline_means*1/3)+1:LTP_baseline_means;
elphys_complete.Baseline_succes = sum(abs(elphys_complete.Norm_slopes(:,1:LTP_baseline_means)-1) < allowed_baseline_noise,2) == length(1:LTP_baseline_means);
elphys_complete.Posttet_success = max([elphys_complete.Norm_slopes(:,LTP_baseline_means+1:LTP_baseline_means+2)],[],2) > minimal_accepted_LTP;
elphys_complete.Signal_ok = mean(elphys_complete.Amplitudes,2,'omitnan') < -minimal_accepted_amplitude;    % Consider refactoring the minimal accepted signal

elphys_complete.manual_ok = elphys_complete.manual_ok == 1;
elphys_complete.popspike_error = elphys_complete.popspike_error == 1;
elphys_complete.manual_low_signal = elphys_complete.manual_low_signal == 1;
elphys_complete.LTP_success = elphys_complete.Baseline_succes;




good_files = elphys_complete.manual_ok & ...
    elphys_complete.auto_ok & ...
    elphys_complete.Signal_ok & ...
    (elphys_complete.Experiment == "Post-tetanisation" | ...
    elphys_complete.Experiment == "Paired pulse");

elphys_cleaned = elphys_complete(good_files,:);

%% Extra

% last_file_day = '--';
%
% % Identifies section number - revise
% for file_no = 1:(files_count-5)
%     we_have_a_date = elphys_files.date{file_no};
%     day_so_far = str2double(we_have_a_date(1:2));
%
%     if elphys_files.Experiment(file_no) == "Baseline" && elphys_files.Experiment(file_no+5) == "Post-tetanisation"
%         elphys_files.Section_no(file_no:file_no+5) = repmat(current_section,[6,1]);
%     end
%
%     if day_so_far == last_file_day && elphys_files.Experiment(file_no) == "Baseline"
%         current_section = current_section+1;
%     end
%
%     if day_so_far ~= last_file_day
%         current_section = 1;
%     end
%
%
%     last_file_day = day_so_far;
%
%
%
% end

% 1233 is a great example file
%
if make_figures == true
    sweeps_no = size(d,3);
    time = [0:1600]/20;
    for current_sweep = 1:(sweeps_no/2)
        little_mean = mean(d(:,2,current_sweep*2-1:current_sweep*2),3);
        little_mean = little_mean - mean(little_mean(100:150));
        plot(time,little_mean(1:1601),LineWidth=1)
        hold on
    end
    xlim([0 1600/20])
    xlabel("ms")
    ylim([-1 0.1])
    ylabel("mV")
    lgd = legend(["200µs, 20 V","200µs, 30 V","200µs, 40 V"]);
    lgd.FontSize = 10;
    hold off
end
%

%     function Baseline_Success()
%         polyfit()
% 
%     end
% end