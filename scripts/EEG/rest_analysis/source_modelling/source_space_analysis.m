function [] = source_space_analysis(subject_id, paths, p)
% DFA and 1/f Analysis in sensor and source space

    %% Load the eLORETA inverse operator
    load(paths.eloreta);

    %% Allocate space to save results
    dfa_alpha_sensor = zeros(n_sensors, 1);
    power_alpha_sensor = zeros(n_sensors, 1);
    slope_sensor = zeros(n_sensors, 1);
    
    dfa_alpha_source = zeros(n_sources, 1);
    power_alpha_source = zeros(n_sources, 1);
    slope_source = zeros(n_sources, 1);
    
    %% Load data
    fprintf('Processing subject %s\n', subject_id);

    filepath = [paths.preproc subject_id];
    filename = [subject_id 'post_ICA_rejected_data.set'];
    fprintf('EEG filename: %s\n', filename);

    EEG = pop_loadset('filename', filename, 'filepath', filepath);
    assert(EEG.srate == p.dfa.Fs, "Unexpected sampling frequency");
    assert(numel(EEG.chanlocs) == numel(all_chanlocs), "Unexpected number of channels");
    data = double(EEG.data);

    %% Compute sensor-space values
    fprintf('Computing alpha power, alpha DFA and 1/f slope for each channel...\n');
    for i_sensor = 1:n_sensors
        ch_data = data(i_sensor, :);

        % Alpha power and DFA
        ch_data_filt = filter_fir(ch_data, p.fir.hp, p.fir.lp, ...
            p.fir.Fs, p.fir.fir_order);
        ch_data_env = abs(hilbert(ch_data_filt));

        [~, ~, DFA_exp] = Scaling_DFA(ch_data_env, p.dfa.Fs, ...
            p.dfa.DFA_SmallTime, p.dfa.DFA_LargeTime, p.dfa.DFA_SmallTimeFit, ...
            p.dfa.DFA_LargeTimeFit, p.dfa.DFA_Overlap, p.dfa.DFA_Plot);
        
        dfa_alpha_sensor(i_sensor) = DFA_exp;
        power_alpha_sensor(i_sensor) = mean(ch_data_env .^ 2);

        % PSD and 1/f slope
        [spec_ch, freqs] = pwelch(ch_data, p.psd.nfft, p.psd.noverlap, ...
            p.psd.nfft, EEG.srate);
        fooof_ch = fooof(freqs, spec_ch, p.fooof.fit_range, ...
            p.fooof.fooof_settings, p.fooof.return_model);
        slope_sensor(i_sensor) = fooof_ch.aperiodic_params(2);

        fprintf('.');
    end
    fprintf('\n');

    %% Compute source-space values
    fprintf('Computing alpha power, alpha DFA and 1/f slope for each source...\n');
    for i_source = 1:n_sources
        src_data = A_eloreta_fixed(:, i_source)' * data;
        src_data_filt = filter_fir(src_data, p.fir.hp, p.fir.lp, ...
            p.fir.Fs, p.fir.fir_order);
        src_data_env = abs(hilbert(src_data_filt));
        
        % Alpha power and DFA
        [~, ~, DFA_exp] = Scaling_DFA(src_data_env, p.dfa.Fs, ...
            p.dfa.DFA_SmallTime, p.dfa.DFA_LargeTime, p.dfa.DFA_SmallTimeFit, ...
            p.dfa.DFA_LargeTimeFit, p.dfa.DFA_Overlap, p.dfa.DFA_Plot);
        
        dfa_alpha_source(i_source) = DFA_exp;
        power_alpha_source(i_source) = mean(src_data_env .^ 2);

        % PSD and 1/f slope
        [spec_src, freqs] = pwelch(src_data, p.psd.nfft, p.psd.noverlap, ...
            p.psd.nfft, EEG.srate);
        fooof_src = fooof(freqs, spec_src, p.fooof.fit_range, ...
            p.fooof.fooof_settings, p.fooof.return_model);
        slope_source(i_source) = fooof_src.aperiodic_params(2);

        fprintf('.');
        if mod(i_source, 100) == 0
            fprintf(' %d\n', i_source);
        end
    end
    fprintf('\n');
    
    %% Save the results
    savefile = [paths.data subject_id '_DFA_power_slope.mat'];
    fprintf('Saving the results to %s\n', savefile);
    save(savefile, ...
        "dfa_alpha_sensor", "power_alpha_sensor", "slope_sensor", ...
        "dfa_alpha_source", "power_alpha_source", "slope_source");
end


