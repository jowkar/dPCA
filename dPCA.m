classdef dPCA < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure           matlab.ui.Figure                   % dPCA
        Panel              matlab.ui.container.Panel          % Input
        Button             matlab.ui.control.Button           % Load
        CheckBox           matlab.ui.control.CheckBox         % Header
        CheckBox2          matlab.ui.control.CheckBox         % Row labels
        Button4            matlab.ui.control.Button           % Select file
        Button5            matlab.ui.control.Button           % Load
        Button6            matlab.ui.control.Button           % Load
        Button7            matlab.ui.control.Button           % Select file
        CheckBox3          matlab.ui.control.CheckBox         % Header
        CheckBox4          matlab.ui.control.CheckBox         % Row labels
        Button8            matlab.ui.control.Button           % Select file
        CheckBox5          matlab.ui.control.CheckBox         % Header
        CheckBox6          matlab.ui.control.CheckBox         % Row labels
        LabelEditField     matlab.ui.control.Label            % Data matrix...
        EditField          matlab.ui.control.EditField       
        LabelDropDown2     matlab.ui.control.Label            % Delimiter
        DropDown2          matlab.ui.control.DropDown         % Tab, Comma
        LabelEditField2    matlab.ui.control.Label            % Sample info...
        EditField2         matlab.ui.control.EditField       
        LabelEditField3    matlab.ui.control.Label            % Test file (...
        EditField3         matlab.ui.control.EditField       
        LabelDropDown4     matlab.ui.control.Label            % Delimiter
        DropDown4          matlab.ui.control.DropDown         % Tab, Comma
        LabelDropDown5     matlab.ui.control.Label            % Delimiter
        DropDown5          matlab.ui.control.DropDown         % Tab, Comma
        Button12           matlab.ui.control.Button           % Clear all data
        Panel2             matlab.ui.container.Panel          % PCA
        Button2            matlab.ui.control.Button           % Plot PCA
        ButtonGroup        matlab.ui.container.ButtonGroup    % PCA type:
        RadioButton        matlab.ui.control.RadioButton      % 3D
        RadioButton2       matlab.ui.control.RadioButton      % 2D
        DropDown6          matlab.ui.control.DropDown         % None
        LabelDropDown6     matlab.ui.control.Label            % Color by:
        Panel3             matlab.ui.container.Panel          % Data options
        LabelDropDown      matlab.ui.control.Label            % Normalization:
        DropDown           matlab.ui.control.DropDown         % None, log2(...
        LabelDropDown3     matlab.ui.control.Label            % Cutoff vari...
        DropDown3          matlab.ui.control.DropDown         % Mean, CV, None
        LabelListBox       matlab.ui.control.Label            % Genes inclu...
        ListBox            matlab.ui.control.ListBox         
        Button3            matlab.ui.control.Button           % Export genes
        Slider             matlab.ui.control.Slider           % [0 100]
        NumericEditField2  matlab.ui.control.NumericEditField % [-Inf Inf]
        NumericEditField   matlab.ui.control.NumericEditField % [-Inf Inf]
        Label              matlab.ui.control.Label            % Cutoff:
        CheckBox7          matlab.ui.control.CheckBox         % Reverse
        LabelListBox2      matlab.ui.control.Label            % Active filt...
        ListBox2           matlab.ui.control.ListBox          % New
        Button9            matlab.ui.control.Button           % Add filter
        Button10           matlab.ui.control.Button           % Remove filter
        CheckBox_log       matlab.ui.control.CheckBox         % Log10 scale
        Panel4             matlab.ui.container.Panel          % Messages
        TextArea           matlab.ui.control.TextArea        
        Button11           matlab.ui.control.Button           % Clear messages
        CheckBox_transpose matlab.ui.control.CheckBox         % Transpose
    end


    properties (Access = public)
        data % data storage
        fighandle
    end

    methods (Access = private)
        
        function normalize(app)
            if ~isempty(app.data.meas)
                value = app.DropDown.Value;
                app.data.normalization = value;
                switch app.data.normalization
                    case 'None'
                        app.data.meas = app.data.raw;
                    case 'log2(X+1)'
                        app.data.meas = log2(app.data.raw+1);
                        disp([min(min(app.data.meas(~isinf(app.data.meas)))) max(max(app.data.meas(~isinf(app.data.meas))))])
                    case 'RPM'
                        app.data.meas = app.normalize_rpm;
                    case 'log2(RPM+1)'
                        app.data.meas = log2(app.normalize_rpm+1);
                    case 'zscore(log2(RPM+1))'
                        app.data.meas = zscore(log2(app.normalize_rpm+1));
                    case 'zscore(log2(X+1))'
                        app.data.meas = zscore(log2(app.data.raw+1));
                    case 'zscore(X)'
                        app.data.meas = zscore(app.data.raw);
                    case 'zscore(RPM)'
                        app.data.meas = zscore(app.normalize_rpm);
                    otherwise
                        app.add_message('Unsupported normalization option')
                        app.data.meas = app.data.raw;
                end
                if ismember(app.data.normalization,[{'None'} {'log2(X+1)'} {'RPM'} {'log2(RPM+1)'}])
                    if ~ismember({'CV'},app.DropDown3.Items)
                        app.DropDown3.Items = [{'CV'} app.DropDown3.Items];
                    end
                    if ~ismember({'Mean'},app.DropDown3.Items)
                        app.DropDown3.Items = [{'Mean'} app.DropDown3.Items];
                    end
                end
                app.Slider.Limits = [min(min(app.data.meas(~isinf(app.data.meas)))) max(max(app.data.meas(~isinf(app.data.meas))))];
            end
        end
        
        function rpm = normalize_rpm(app)
            d = app.data.raw';
            s = sum(d);
            ss = repmat(s,size(d,1),1);
            rpm = ((d./ss).*(10^6))';
        end

        function add_message(app,message)
            current_message = char(app.TextArea.Value);
            if size(current_message,1)>1
               for i = 1:size(current_message,1)
                   if i == 1
                       str = current_message(i,:);
                   else
                       str = [str '\n' strrep(current_message(i,:),sprintf('\n'),'')];
                   end
               end
               app.TextArea.Value = sprintf([str '\n' message]);
            else
                app.TextArea.Value = sprintf([char(app.TextArea.Value) '\n' message]);
            end
        end

        function update_included_genes(app)
            if all(strcmp(app.ListBox2.Items,'New'))
                app.data.censored = zeros(1,size(app.data.meas,2));
                if isempty(app.data.genes)
                    genes_included = sprintfc('%i',find(~app.data.censored));
                else
                    genes_included = app.data.genes;
                end
                app.ListBox.Items = genes_included;
                app.NumericEditField2.Value = length(genes_included);
                return
            else
                app.data.censored = zeros(1,size(app.data.meas,2));
                for i = 1:length(app.ListBox2.Items)
                    str = char(app.ListBox2.Items(i));
                    if ~strcmp(str,'New')
                        [cutoff_variable,cutoff_value,rev,islog] = app.parse_filter(str);
                     
                        cases = app.DropDown3.Items;
                        option_nr = ismember(cases,cutoff_variable);
                        switch cutoff_variable
                            case 'Mean'
                                v = mean(app.data.meas,1);
                                if islog
                                    v = log10(v);
                                end
                                if ~rev
                                   censored = v < cutoff_value;
                                else
                                   censored = v > cutoff_value; 
                                end
                            case 'CV'
                                v = mean(app.data.meas,1);
                                s = std(app.data.meas,[],1);
                                v = s./v;
                                if islog
                                    v = log10(v);
                                end
                                if ~rev
                                    censored = v < cutoff_value;
                                else
                                    censored = v > cutoff_value;
                                end
                            case 'None'
                                censored = zeros(1,size(app.data.meas,2));
                            otherwise
                                if ~ismember(cases{option_nr},fields(app.data))
                                    app.add_message(sprintf('Error: missing data for variable %s',cases{option_nr}));
                                    return
                                end
                                v = app.data.(cases{option_nr}); %app.data.p, app.data.fold_change, etc
                                if islog
                                    v = log10(v);
                                end
                                if ~rev
                                    censored = v < cutoff_value;
                                else
                                    censored = v > cutoff_value;
                                end
                        end
                        app.data.censored = app.data.censored | censored;
                    end
                end
                
                if isempty(app.data.genes)
                    genes_included = sprintfc('%i',find(~app.data.censored));
                else
                    genes_included = app.data.genes(~app.data.censored);
                end
                app.ListBox.Items = genes_included;
                app.NumericEditField2.Value = length(genes_included);
            end
        end

        function [cutoff_variable,cutoff_value,rev,islog] = parse_filter(app,str)
            if ~strcmp(str,'New')
                presplit = strsplit(str,' (log10)');
                if length(presplit)>1
                    islog = 1;
                else
                    islog = 0;
                end
                str = char(presplit(1));
                [s,m] = strsplit(str,{' < ',' > '});
                cutoff_variable = char(s(1));
                cutoff_value = str2double(s(2));
                if strcmp(m,' < ')
                    rev = 1;
                elseif strcmp(m,' > ')
                    rev = 0;
                else
                    error('Unknown filter criterion. Expected either > or <.')
                end
            else
               cutoff_variable = str;
               cutoff_value = [];
               rev = [];
            end
        end

        function clear_all_data(app)
            app.data = rmfield(app.data,fields(app.data));
            if ishandle(app.fighandle)
                figure(app.fighandle)
            else
                app.fighandle = figure;
                figure(app.fighandle)
            end
            clf
            app.startupFcn;
        end

        function clear_main_data(app)            
            app.data.main_data_path = '';
            app.data.normalization = 'None';
            app.data.censored = [];
            app.data.classes = [];
            app.data.meas = [];
            app.data.raw = [];
            app.data.species = [];
            ish = ishandle(app.fighandle);
            if isempty(ish)
                ish = 0;
            end
            if ~ish
                app.fighandle = figure;
            end
            app.data.cutoff = app.Slider.Value;
            app.data.cutoff_variable = app.DropDown3.Value;
        end
    end


    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.data.main_data_path = [];
            app.data.normalization = 'None';
            ish = ishandle(app.fighandle);
            if isempty(ish)
                ish = 0;
            end
            if ~ish
                app.fighandle = figure;
            end
            app.data.cutoff = app.Slider.Value;
            app.data.cutoff_variable = app.DropDown3.Value;
            app.data.censored = [];
            app.data.genes = [];
            app.data.classes = [];
            app.data.meas = [];
            app.data.raw = [];
            app.data.fc_file.varnames =[];
            app.data.fc_file.size = [];
            app.data.sample_info = [];
            app.data.sample_info_path = '';
            app.data.fc_path = '';
            app.data.fc_file.varnames = [];
            app.data.fc_file.size = [];
        end

        % Button button pushed function
        function load_data(app)
            app.TextArea.Value = 'Loading data...';
            if ~isempty(app.EditField.Value) && isempty(app.data.main_data_path)
                pth = app.EditField.Value;
            else
                pth = app.data.main_data_path;
            end
            app.clear_main_data;
            if isempty(pth)
                d = load('fisheriris');
                app.data.species = d.species;
                app.data.classes = app.data.species;
                app.data.raw = d.meas;
                app.data.meas = app.data.raw;
                app.Slider.Limits = [min(min(app.data.meas)) max(max(app.data.meas))];
                app.data.censored = zeros(1,size(app.data.meas,2));
                app.NumericEditField2.Value = size(app.data.raw,2);
                app.update_included_genes;
                app.add_message('Demo dataset loaded');
            else
                success = 0;
                try
                    if ~app.CheckBox.Value && ~app.CheckBox2.Value && strcmp(app.DropDown2.Value,'Tab')
                        d = load(pth);
                        app.data.species = sprintfc('%i',1:min(size(d)));
                    else
                        del = 'tab';
                        switch app.DropDown2.Value
                            case 'Tab'
                                del = 'tab';
                            case 'Comma'
                                del = ',';
                        end     
                       t = readtable(pth,'ReadVariableNames',app.CheckBox.Value,'ReadRowNames',app.CheckBox2.Value,'Delimiter',del);
                       d = table2array(t);
                       app.data.species = t.VariableNames;
                       app.data.genes = t.RowNames;
                    end
                    if app.CheckBox_transpose.Value
                        app.data.raw = d;
                    else
                        app.data.raw = d';
                    end
                    app.data.meas = app.data.raw;
                    app.add_message('Data loaded');
                    success = 1;
                catch
                    app.add_message('Failed to load data. Please check the path and format of the file.');
                end
                if success
                    app.normalize;
                    app.Slider.Limits = [min(min(app.data.meas)) max(max(app.data.meas))];
                    app.data.censored = zeros(1,size(app.data.meas,2));
                    
                    % Is the number of samples equal to those in the sample info file?
                    if ~isempty(app.data.sample_info);
                        eq = 0;
                        fnames = fields(app.data.sample_info);
                        for i = 1:length(fnames)
                            eq = eq || (length(app.data.sample_info.(fnames{i})) == size(app.data.raw,1));
                            if eq
                                break
                            end
                        end
                        if ~eq
                            app.add_message('Warning: Number of samples in sample info is not equal to the number of samples (columns) in the data matrix');
                        end
                    end
                    
                    % Is the number of genes equal to those in the test file?
                    if ~isempty(app.data.fc_file.size)
                        if app.data.fc_file.size(1) ~= size(app.data.raw,2)
                            app.add_message('Warning: Number of genes in the test file is not equal to the number of genes (rows) in the data matrix');
                        end
                    end
                    app.NumericEditField2.Value = size(app.data.raw,2);
                    app.update_included_genes;
                end
                app.data.main_data_path = '';
            end
        end

        % Button2 button pushed function
        function plot_pca(app)
            try
                app.add_message('Performing PCA and updating plot...')
                if isempty(app.data.censored)
                    d = app.data.meas;
                    [~,score] = pca(d);
                else                    
                    d = app.data.meas(:,~app.data.censored);
                    [~,score] = pca(d);
                end
                if ~isempty(app.data.classes)
                    u = unique(app.data.classes);
                    elem = app.data.classes;
                else
                    u = unique(app.data.species);
                    elem = app.data.species;
                end
                if ishandle(app.fighandle)
                    figure(app.fighandle)
                else
                    app.fighandle = figure;
                    figure(app.fighandle)
                end
                clf
                if ~isempty(app.data.classes)
                    hold all
                    for i = 1:length(u)
                        idx = ismember(elem,u(i));
                        if app.RadioButton2.Value == 0
                            scatter3(score(idx,1),score(idx,2),score(idx,3))
                            view(40,35)
                        else
                            scatter(score(idx,1),score(idx,2))
                        end
                    end
                    legend(u)
                else
                    if app.RadioButton2.Value == 0
                        scatter3(score(:,1),score(:,2),score(:,3))
                        view(40,35)
                    else
                        scatter(score(:,1),score(:,2))
                    end
                end
                grid on
                xlabel('PC1')
                ylabel('PC2')
                if app.RadioButton2.Value == 0
                    zlabel('PC3')
                end
                app.add_message('Done')
            catch ME
                if ((app.RadioButton2.Value == 0) && size(app.data.meas(:,~app.data.censored),2) < 3) || ((app.RadioButton2.Value == 1) && size(app.data.meas(:,~app.data.censored),2) < 2)
                    app.add_message('Too few variables remaining');
                else
                    app.add_message(ME.message);
                end
            end
        end

        % Slider value changed function
        function change_cutoff(app)
            selected_active_filter = app.ListBox2.Value;
            if ~strcmp(selected_active_filter,'New')
                idx = ismember(app.ListBox2.Items,selected_active_filter);
                [cutoff_variable,~,~] = app.parse_filter(selected_active_filter);
                if app.CheckBox7.Value
                    sig = '<';
                else
                    sig = '>';
                end
                str = sprintf('%s %s %d',cutoff_variable,sig,app.Slider.Value);
                if app.CheckBox_log.Value
                   str = sprintf('%s (log10)',char(str)); 
                end
                app.ListBox2.Items{idx} = str;
                app.ListBox2.Value = str;
                app.update_included_genes;
                app.plot_pca;
            end
            app.NumericEditField.Value = app.Slider.Value;
        end

        % EditField value changed function
        function EditFieldValueChanged(app)
            value = app.EditField.Value;
            app.data.main_data_path = value;
        end

        % DropDown value changed function
        function normalization_changed(app, event)
            value = app.DropDown.Value;
            app.data.normalization = value;
            app.normalize;
            app.ListBox2.Items = {'New'};
            if ismember(app.data.normalization,[{'zscore(X)'},{'zscore(log2(X+1))'},{'zscore(RPM)'},{'zscore(log2(RPM+1))'}])
                app.DropDown3.Items = app.DropDown3.Items(~ismember(app.DropDown3.Items,[{'Mean'} {'CV'}]));
                app.DropDown3.Value = 'None';                    
            end
            app.cutoff_variable_changed;
        end

        % DropDown3 value changed function
        function cutoff_variable_changed(app, event)
            try
                app.data.cutoff_variable = app.DropDown3.Value;
                
                [selected_active_filter_cutoff_variable,~,~] = app.parse_filter(app.ListBox2.Value);
                if ~strcmp(app.data.cutoff_variable,selected_active_filter_cutoff_variable)
                    if ~strcmp(app.ListBox2.Value,'New')
                        app.ListBox2.Value = 'New';
                    end
                end
                
                switch app.data.cutoff_variable
                    case 'None'
                        app.Slider.Limits = [0 1];
                        app.Slider.Value = 0;
                        app.Slider.Enable = 'off';
                        app.Button9.Enable = 'off';
                    case 'Mean'
                        app.Slider.Enable = 'on';
                        app.Button9.Enable = 'on';
                        v = mean(app.data.meas,1);
                        app.Slider.Limits = [min(min(v)) max(max(v))];
                        app.Slider.Value = min(min(v));
                        app.log_scale_changed;
                    case 'CV'
                        app.Slider.Enable = 'on';
                        app.Button9.Enable = 'on';
                        v = mean(app.data.meas,1);
                        s = std(app.data.meas,[],1);
                        v = (s./v);
                        app.Slider.Limits = [min(min(v)) max(max(v))];
                        app.Slider.Value = min(min(v));
                        app.log_scale_changed;
                    otherwise
                        app.Slider.Enable = 'on';
                        app.Button9.Enable = 'on';
                        v = app.data.(app.data.cutoff_variable);
                        app.Slider.Limits = [min(min(v)) max(max(v))];
                        app.Slider.Value = min(min(v));
                        app.log_scale_changed;
                end
            catch ME
                if isempty(app.data.meas)
                    app.add_message('Cannot apply cutoff when no data is loaded')
                else
                    app.add_message(ME.message)
                end
            end
        end

        % Button4 button pushed function
        function select_file(app)
            [filename1,filepath1] = uigetfile({'*.*','All Files'},'Select data file');
            if filename1 == 0
                return
            else
                app.data.main_data_path = [filepath1 filesep filename1];
                app.EditField.Value = app.data.main_data_path;
            end
        end

        % Button8 button pushed function
        function select_fold_change_file(app)
            [filename1,filepath1] = uigetfile({'*.*','All Files'},'Select data file');
            if filename1 == 0
                return
            else
                app.data.fc_path = [filepath1 filesep filename1];
                app.EditField3.Value = app.data.fc_path;
            end
        end

        % Button6 button pushed function
        function load_fold_change_file(app)
            if isempty(app.data.fc_path) && isempty(app.EditField3.Value)
                return
            else
                app.add_message('Loading test file...');
                try
                    del = 'tab';
                    switch app.DropDown5.Value
                        case 'Tab'
                            del = 'tab';
                        case 'Comma'
                            del = ',';
                    end     
                    if isempty(app.data.fc_path) && ~isempty(app.EditField3.Value)
                        app.data.fc_path = app.EditField3.Value;
                    end
                    t = readtable(app.data.fc_path,'ReadVariableNames',app.CheckBox5.Value,'ReadRowNames',app.CheckBox6.Value,'Delimiter',del);
                    varnames = t.Properties.VariableNames;
                    
                    if ~isempty(app.data.fc_file) && ~isempty(app.data.fc_file.varnames)
                        fields_to_remove = app.data.fc_file.varnames(~ismember(app.data.fc_file.varnames,varnames));
                        app.data = rmfield(app.data,fields_to_remove);
                    end
                    
                    for i = 1:length(varnames)
                        isn(i) = all(isnumeric(t.(varnames{i})));
                    end
                    if ~any(isn)
                        app.add_message('No column in the fold change table contained entirely numeric data. Please check whether the file contains a header row.')
                        return
                    end
                    if app.CheckBox6.Value
                        app.data.genes = t.Properties.RowNames;
                    end
    
                    if ~isempty(app.data.meas)
                        app.update_included_genes;
                    end
                    for i = 1:length(varnames)
                        app.data.(varnames{i}) = t.(varnames{i})';
                    end
                    app.DropDown3.Items = [{'None'},{'Mean'},{'CV'},varnames];
                    app.data.fc_file.size = size(t);
                    app.data.fc_file.varnames = varnames;
                    if size(t,1) ~= size(app.data.meas,2)
                        app.add_message('Warning: The number of genes (rows) in the test file is not equal to the number of genes (columns, or rows of transposed) in the data matrix.');
                    end
                    app.add_message('Data loaded');
                catch ME
                    disp(ME.message)
                    app.add_message('Failed to load data. Please check the path and format of the file.');
                end
                app.data.fc_path = '';
            end
        end

        % Button7 button pushed function
        function select_sample_info_file(app)
            [filename1,filepath1] = uigetfile({'*.*','All Files'},'Select file');
            if filename1 == 0
                return
            else
                app.data.sample_info_path = [filepath1 filesep filename1];
                app.EditField2.Value = app.data.sample_info_path;
            end
        end

        % Button5 button pushed function
        function load_sample_info_file(app)
            if isempty(app.data.sample_info_path) && isempty(app.EditField2.Value)
                return
            else
                app.add_message('Loading sample info file...');
                try
                    del = 'tab';
                    switch app.DropDown4.Value
                        case 'Tab'
                            del = 'tab';
                        case 'Comma'
                            del = ',';
                    end     
                    if isempty(app.data.sample_info_path) && ~isempty(app.EditField2.Value)
                       app.data.sample_info_path = app.EditField2.Value; 
                    end
                        
                    t = readtable(app.data.sample_info_path,'ReadVariableNames',app.CheckBox3.Value,'ReadRowNames',app.CheckBox4.Value,'Delimiter',del);
                    varnames = t.Properties.VariableNames;
                    for i = 1:length(varnames)
                       app.data.sample_info.(varnames{i}) = t.(varnames{i});
                    end
                    app.DropDown6.Items = [{'None'},varnames];
                    if size(t,1) ~= size(app.data.meas,1)
                        app.add_message('Warning: The number of samples (rows) in the sample info file is not equal to the number of samples (columns, or rows if transposed) in the data matrix');
                    end
                    app.DropDown6.Enable = 'on';
                    app.add_message('Sample info loaded');
                catch
                    app.add_message('Failed to load sample info. Please check the path and format of the file.');
                end
                app.data.sample_info_path = '';
            end
        end

        % DropDown6 value changed function
        function color_by_changed(app, event)
            value = app.DropDown6.Value;
            if strcmp(value,'None')
               app.data.classes = []; 
            else
                if isempty(app.data.sample_info) || ~ismember(value,fields(app.data.sample_info))
                    app.add_message(sprintf('Error: %s is not a column in the sample info file',value))
                else
                    app.data.classes = app.data.sample_info.(value);
                end
            end
            app.plot_pca;
        end

        % Button3 button pushed function
        function export_genes(app)
            genes = app.ListBox.Items;
            t = table(genes');
            t.Properties.VariableNames = {'ID'};
            [FileName,PathName] = uiputfile('*.txt','Save as','exported_genes.txt');
            if FileName == 0
                return
            else
                writetable(t,[PathName filesep FileName]);
            end
        end

        % NumericEditField value changed function
        function cutoff_value_typed(app)
            try
                value = app.NumericEditField.Value;
                app.Slider.Value = value;
            catch ME
                app.add_message('Value of cutoff not within the limits');
            end
        end

        % Button9 button pushed function
        function add_filter(app)
            if app.CheckBox7.Value
                sig = '<';
            else
                sig = '>';
            end
            str = sprintf('%s %s %d',app.DropDown3.Value,sig,app.Slider.Value);
            if app.CheckBox_log.Value
               str = sprintf('%s (log10)',char(str));
            end
            app.ListBox2.Items = [app.ListBox2.Items str];
            if ~isempty(app.data.meas)
                app.update_included_genes;
                app.plot_pca;
            end
        end

        % Button10 button pushed function
        function remove_filter(app)
            if ~strcmp(app.ListBox2.Value,'New')
                fidx = find(ismember(app.ListBox2.Items,app.ListBox2.Value));
                app.ListBox2.Items(fidx(1)) = [];
                %app.ListBox2.Items = app.ListBox2.Items(~ismember(app.ListBox2.Items,app.ListBox2.Value));
                if ~isempty(app.data.meas)
                    app.update_included_genes;
                    app.plot_pca;
                end
            end
        end

        % CheckBox7 value changed function
        function reverse_cutoff_changed(app)
            try
                selected_active_filter = app.ListBox2.Value;
                if ~strcmp(selected_active_filter,'New')
                    idx = ismember(app.ListBox2.Items,selected_active_filter);
                    [cutoff_variable,cutoff_value,~,islog] = app.parse_filter(selected_active_filter);
                    if app.CheckBox7.Value;
                        sig = '<';
                    else
                        sig = '>';
                    end
                    
                    str = sprintf('%s %s %d',cutoff_variable,sig,cutoff_value);
                    if islog
                       str = sprintf('%s (log10)',char(str)); 
                    end
                    app.ListBox2.Items{idx} = str;
                    app.ListBox2.Value = str;
                    app.NumericEditField.Value = cutoff_value;
                    app.update_included_genes;
                    app.plot_pca;                
                end
            catch ME
               app.add_message(sprintf('Error: %s',ME.message)); 
            end
        end

        % ListBox2 value changed function
        function selected_active_filter_changed(app)
            value = app.ListBox2.Value;
            if ~strcmp(value,'New')
                try
                    [cutoff_variable,cutoff_value,rev,islog] = app.parse_filter(value);
                    app.DropDown3.Value = cutoff_variable;
                    app.cutoff_variable_changed;
                    app.CheckBox_log.Value = islog;
                    app.log_scale_changed;
    
                    if cutoff_value < min(app.Slider.Limits)
                        app.Slider.Limits = [cutoff_value max(app.Slider.Limits)];
                    end
                    app.Slider.Value = cutoff_value;
                    app.CheckBox7.Value = rev;
                catch ME
                    app.add_message(sprintf('Error: %s',ME.message))
                end
            end
        end

        % CheckBox_log value changed function
        function log_scale_changed(app)
            try
                ischecked = app.CheckBox_log.Value;
                cutoff_variable = app.DropDown3.Value;
                
                switch cutoff_variable
                    case 'Mean'
                        v = mean(app.data.meas,1);
                    case 'CV'
                        v = mean(app.data.meas,1);
                        s = std(app.data.meas,[],1);
                        v = s./v;
                    case 'None'
                        v = [];
                    otherwise
                        v = app.data.(cutoff_variable); %app.data.p, app.data.fold_change, etc
                end
                
                if ~isempty(v)
                    if ischecked
                        if min(v)>0
                            v = log10(v);
                            mn = min(v);
                            mx = max(v);
                            app.Slider.Limits = [mn mx];
                            app.Slider.Value = min(v);
                        else
                            app.CheckBox_log.Value = 0;
                            app.add_message('Not possible to log transform non-positive values.')
                        end
                    else
                       app.Slider.Limits = [min(v) max(v)];
                       app.Slider.Value = min(v);
                    end
                end
            catch ME
               app.add_message(sprintf('Error: %s',ME.message)); 
            end
        end

        % Button12 button pushed function
        function clear_all_data_button_pushed(app)
            app.EditField.Value = '';
            app.EditField2.Value = '';
            app.EditField3.Value = '';
            app.ListBox2.Items = {};
            app.ListBox.Items = {};
            app.NumericEditField2.Value = 0;
            app.NumericEditField.Value = 0;
            app.TextArea.Value = '';
            app.DropDown6.Items = {'None'};
            app.DropDown6.Value = 'None';
            app.DropDown3.Items = {'None','Mean','CV'};
            app.DropDown3.Value = 'None';
            app.DropDown.Value = 'None';
            app.ListBox2.Items = {'New'};
            app.clear_all_data;
            app.cutoff_variable_changed;
        end

        % Button11 button pushed function
        function clear_messages(app)
            app.TextArea.Value = '';
        end

        % UIFigure close request function
        function UIFigureCloseRequest(app)
            close all
            delete(app)
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 734 653];
            app.UIFigure.Name = 'dPCA';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest);
            setAutoResize(app, app.UIFigure, true)

            % Create Panel
            app.Panel = uipanel(app.UIFigure);
            app.Panel.BorderType = 'line';
            app.Panel.Title = 'Input';
            app.Panel.FontUnits = 'pixels';
            app.Panel.FontSize = 12;
            app.Panel.Units = 'pixels';
            app.Panel.Position = [-1 431 736 222];

            % Create Button
            app.Button = uibutton(app.Panel, 'push');
            app.Button.ButtonPushedFcn = createCallbackFcn(app, @load_data);
            app.Button.Position = [40 142 100 22];
            app.Button.Text = 'Load';

            % Create CheckBox
            app.CheckBox = uicheckbox(app.Panel);
            app.CheckBox.Text = 'Header';
            app.CheckBox.Position = [275 174 58.109375 16];

            % Create CheckBox2
            app.CheckBox2 = uicheckbox(app.Panel);
            app.CheckBox2.Text = 'Row labels';
            app.CheckBox2.Position = [353 174 77.875 16];

            % Create Button4
            app.Button4 = uibutton(app.Panel, 'push');
            app.Button4.ButtonPushedFcn = createCallbackFcn(app, @select_file);
            app.Button4.Position = [150 171 100 22];
            app.Button4.Text = 'Select file';

            % Create Button5
            app.Button5 = uibutton(app.Panel, 'push');
            app.Button5.ButtonPushedFcn = createCallbackFcn(app, @load_sample_info_file);
            app.Button5.Position = [40 81 100 22];
            app.Button5.Text = 'Load';

            % Create Button6
            app.Button6 = uibutton(app.Panel, 'push');
            app.Button6.ButtonPushedFcn = createCallbackFcn(app, @load_fold_change_file);
            app.Button6.Position = [40.5 14 100 22];
            app.Button6.Text = 'Load';

            % Create Button7
            app.Button7 = uibutton(app.Panel, 'push');
            app.Button7.ButtonPushedFcn = createCallbackFcn(app, @select_sample_info_file);
            app.Button7.Position = [150 110 100 22];
            app.Button7.Text = 'Select file';

            % Create CheckBox3
            app.CheckBox3 = uicheckbox(app.Panel);
            app.CheckBox3.Text = 'Header';
            app.CheckBox3.Position = [275 113 58.109375 16];

            % Create CheckBox4
            app.CheckBox4 = uicheckbox(app.Panel);
            app.CheckBox4.Text = 'Row labels';
            app.CheckBox4.Position = [353 113 77.875 16];

            % Create Button8
            app.Button8 = uibutton(app.Panel, 'push');
            app.Button8.ButtonPushedFcn = createCallbackFcn(app, @select_fold_change_file);
            app.Button8.Position = [150 48 100 22];
            app.Button8.Text = 'Select file';

            % Create CheckBox5
            app.CheckBox5 = uicheckbox(app.Panel);
            app.CheckBox5.Text = 'Header';
            app.CheckBox5.Position = [275 51 58.109375 16];

            % Create CheckBox6
            app.CheckBox6 = uicheckbox(app.Panel);
            app.CheckBox6.Text = 'Row labels';
            app.CheckBox6.Position = [353 51 77.875 16];

            % Create LabelEditField
            app.LabelEditField = uilabel(app.Panel);
            app.LabelEditField.HorizontalAlignment = 'right';
            app.LabelEditField.Position = [56.03125 175 84 15];
            app.LabelEditField.Text = 'Data matrix file:';

            % Create EditField
            app.EditField = uieditfield(app.Panel, 'text');
            app.EditField.ValueChangedFcn = createCallbackFcn(app, @EditFieldValueChanged);
            app.EditField.Position = [150.03125 142 464 22];

            % Create LabelDropDown2
            app.LabelDropDown2 = uilabel(app.Panel);
            app.LabelDropDown2.HorizontalAlignment = 'right';
            app.LabelDropDown2.Position = [451.03125 175 48 15];
            app.LabelDropDown2.Text = 'Delimiter';

            % Create DropDown2
            app.DropDown2 = uidropdown(app.Panel);
            app.DropDown2.Items = {'Tab', 'Comma'};
            app.DropDown2.Position = [514.03125 173 100 20];
            app.DropDown2.Value = 'Tab';

            % Create LabelEditField2
            app.LabelEditField2 = uilabel(app.Panel);
            app.LabelEditField2.HorizontalAlignment = 'right';
            app.LabelEditField2.Position = [4.03125 114 139 15];
            app.LabelEditField2.Text = 'Sample info file (optional):';

            % Create EditField2
            app.EditField2 = uieditfield(app.Panel, 'text');
            app.EditField2.Position = [148.03125 81 466 22];

            % Create LabelEditField3
            app.LabelEditField3 = uilabel(app.Panel);
            app.LabelEditField3.HorizontalAlignment = 'right';
            app.LabelEditField3.Position = [42.03125 52 98 15];
            app.LabelEditField3.Text = 'Test file (optional):';

            % Create EditField3
            app.EditField3 = uieditfield(app.Panel, 'text');
            app.EditField3.Position = [149.03125 14 467 22];

            % Create LabelDropDown4
            app.LabelDropDown4 = uilabel(app.Panel);
            app.LabelDropDown4.HorizontalAlignment = 'right';
            app.LabelDropDown4.Position = [455.03125 113 48 15];
            app.LabelDropDown4.Text = 'Delimiter';

            % Create DropDown4
            app.DropDown4 = uidropdown(app.Panel);
            app.DropDown4.Items = {'Tab', 'Comma'};
            app.DropDown4.Position = [512.03125 109 100 22];
            app.DropDown4.Value = 'Tab';

            % Create LabelDropDown5
            app.LabelDropDown5 = uilabel(app.Panel);
            app.LabelDropDown5.HorizontalAlignment = 'right';
            app.LabelDropDown5.Position = [451.03125 52 48 15];
            app.LabelDropDown5.Text = 'Delimiter';

            % Create DropDown5
            app.DropDown5 = uidropdown(app.Panel);
            app.DropDown5.Items = {'Tab', 'Comma'};
            app.DropDown5.Position = [514.03125 48 100 22];
            app.DropDown5.Value = 'Tab';

            % Create Button12
            app.Button12 = uibutton(app.Panel, 'push');
            app.Button12.ButtonPushedFcn = createCallbackFcn(app, @clear_all_data_button_pushed);
            app.Button12.Position = [625 14 100 22];
            app.Button12.Text = 'Clear all data';

            % Create CheckBox_transpose
            app.CheckBox_transpose = uicheckbox(app.Panel);
            app.CheckBox_transpose.Text = 'Transpose';
            app.CheckBox_transpose.Position = [625 175 75.453125 16];

            % Create Panel2
            app.Panel2 = uipanel(app.UIFigure);
            app.Panel2.BorderType = 'line';
            app.Panel2.Title = 'PCA';
            app.Panel2.FontUnits = 'pixels';
            app.Panel2.FontSize = 12;
            app.Panel2.Units = 'pixels';
            app.Panel2.Position = [356 0 141 432];

            % Create Button2
            app.Button2 = uibutton(app.Panel2, 'push');
            app.Button2.ButtonPushedFcn = createCallbackFcn(app, @plot_pca);
            app.Button2.Position = [9 377 100 22];
            app.Button2.Text = 'Plot PCA';

            % Create ButtonGroup
            app.ButtonGroup = uibuttongroup(app.Panel2);
            app.ButtonGroup.BorderType = 'line';
            app.ButtonGroup.Title = 'PCA type:';
            app.ButtonGroup.FontUnits = 'pixels';
            app.ButtonGroup.FontSize = 12;
            app.ButtonGroup.Units = 'pixels';
            app.ButtonGroup.Position = [9 284 123 79];

            % Create RadioButton
            app.RadioButton = uiradiobutton(app.ButtonGroup);
            app.RadioButton.Text = '3D';
            app.RadioButton.Position = [10 32 34.125 16];
            app.RadioButton.Value = true;

            % Create RadioButton2
            app.RadioButton2 = uiradiobutton(app.ButtonGroup);
            app.RadioButton2.Text = '2D';
            app.RadioButton2.Position = [10 10 34.125 16];

            % Create DropDown6
            app.DropDown6 = uidropdown(app.Panel2);
            app.DropDown6.Items = {'None'};
            app.DropDown6.ValueChangedFcn = createCallbackFcn(app, @color_by_changed, true);
            app.DropDown6.Enable = 'off';
            app.DropDown6.Position = [13 227 100 22];
            app.DropDown6.Value = 'None';

            % Create LabelDropDown6
            app.LabelDropDown6 = uilabel(app.Panel2);
            app.LabelDropDown6.HorizontalAlignment = 'right';
            app.LabelDropDown6.Position = [13.03125 256 49 15];
            app.LabelDropDown6.Text = 'Color by:';

            % Create Panel3
            app.Panel3 = uipanel(app.UIFigure);
            app.Panel3.BorderType = 'line';
            app.Panel3.Title = 'Data options';
            app.Panel3.FontUnits = 'pixels';
            app.Panel3.FontSize = 12;
            app.Panel3.Units = 'pixels';
            app.Panel3.Position = [-1 0 358 432];

            % Create LabelDropDown
            app.LabelDropDown = uilabel(app.Panel3);
            app.LabelDropDown.HorizontalAlignment = 'right';
            app.LabelDropDown.Position = [17.03125 379 78 15];
            app.LabelDropDown.Text = 'Normalization:';

            % Create DropDown
            app.DropDown = uidropdown(app.Panel3);
            app.DropDown.Items = {'None', 'log2(RPM+1)', 'zscore(log2(RPM+1))', 'log2(X+1)', 'zscore(log2(X+1))', 'zscore(RPM)', 'zscore(X)', 'RPM'};
            app.DropDown.ValueChangedFcn = createCallbackFcn(app, @normalization_changed, true);
            app.DropDown.Position = [109.03125 377 100 20];
            app.DropDown.Value = 'None';

            % Create LabelDropDown3
            app.LabelDropDown3 = uilabel(app.Panel3);
            app.LabelDropDown3.HorizontalAlignment = 'right';
            app.LabelDropDown3.Position = [12.03125 344 82 15];
            app.LabelDropDown3.Text = 'Cutoff variable:';

            % Create DropDown3
            app.DropDown3 = uidropdown(app.Panel3);
            app.DropDown3.Items = {'Mean', 'CV', 'None'};
            app.DropDown3.ValueChangedFcn = createCallbackFcn(app, @cutoff_variable_changed, true);
            app.DropDown3.Position = [109.03125 340 100 22];
            app.DropDown3.Value = 'Mean';

            % Create LabelListBox
            app.LabelListBox = uilabel(app.Panel3);
            app.LabelListBox.Position = [216.6875 296 88 15];
            app.LabelListBox.Text = 'Genes included:';

            % Create ListBox
            app.ListBox = uilistbox(app.Panel3);
            app.ListBox.Items = {};
            app.ListBox.Position = [216.6875 47 130 219];
            app.ListBox.Value = {};

            % Create Button3
            app.Button3 = uibutton(app.Panel3, 'push');
            app.Button3.ButtonPushedFcn = createCallbackFcn(app, @export_genes);
            app.Button3.Position = [217 13 100 22];
            app.Button3.Text = 'Export genes';

            % Create Slider
            app.Slider = uislider(app.Panel3);
            app.Slider.Orientation = 'vertical';
            app.Slider.ValueChangedFcn = createCallbackFcn(app, @change_cutoff);
            app.Slider.Position = [36.6875 90 3 208];

            % Create NumericEditField2
            app.NumericEditField2 = uieditfield(app.Panel3, 'numeric');
            app.NumericEditField2.Editable = 'off';
            app.NumericEditField2.Position = [217.03125 268 130 22];

            % Create NumericEditField
            app.NumericEditField = uieditfield(app.Panel3, 'numeric');
            app.NumericEditField.ValueChangedFcn = createCallbackFcn(app, @cutoff_value_typed);
            app.NumericEditField.Position = [19.03125 47 82 22];

            % Create Label
            app.Label = uilabel(app.Panel3);
            app.Label.Position = [19 313 37 15];
            app.Label.Text = 'Cutoff:';

            % Create CheckBox7
            app.CheckBox7 = uicheckbox(app.Panel3);
            app.CheckBox7.ValueChangedFcn = createCallbackFcn(app, @reverse_cutoff_changed);
            app.CheckBox7.Text = 'Reverse';
            app.CheckBox7.Position = [19 16 62.5625 16];

            % Create LabelListBox2
            app.LabelListBox2 = uilabel(app.Panel3);
            app.LabelListBox2.HorizontalAlignment = 'right';
            app.LabelListBox2.Position = [108.6875 272 69 15];
            app.LabelListBox2.Text = 'Active filters:';

            % Create ListBox2
            app.ListBox2 = uilistbox(app.Panel3);
            app.ListBox2.Items = {'New'};
            app.ListBox2.ValueChangedFcn = createCallbackFcn(app, @selected_active_filter_changed);
            app.ListBox2.Position = [108.6875 47 100 219];
            app.ListBox2.Value = 'New';

            % Create Button9
            app.Button9 = uibutton(app.Panel3, 'push');
            app.Button9.ButtonPushedFcn = createCallbackFcn(app, @add_filter);
            app.Button9.Position = [109 292 100 22];
            app.Button9.Text = 'Add filter';

            % Create Button10
            app.Button10 = uibutton(app.Panel3, 'push');
            app.Button10.ButtonPushedFcn = createCallbackFcn(app, @remove_filter);
            app.Button10.Position = [109 13 100 22];
            app.Button10.Text = 'Remove filter';

            % Create CheckBox_log
            app.CheckBox_log = uicheckbox(app.Panel3);
            app.CheckBox_log.ValueChangedFcn = createCallbackFcn(app, @log_scale_changed);
            app.CheckBox_log.Text = 'Log10 scale';
            app.CheckBox_log.Position = [217 343 84.125 16];

            % Create Panel4
            app.Panel4 = uipanel(app.UIFigure);
            app.Panel4.BorderType = 'line';
            app.Panel4.Title = 'Messages';
            app.Panel4.FontUnits = 'pixels';
            app.Panel4.FontSize = 12;
            app.Panel4.Units = 'pixels';
            app.Panel4.Position = [496 0 239 432];

            % Create TextArea
            app.TextArea = uitextarea(app.Panel4);
            app.TextArea.Position = [6.03125 48 223 350];

            % Create Button11
            app.Button11 = uibutton(app.Panel4, 'push');
            app.Button11.ButtonPushedFcn = createCallbackFcn(app, @clear_messages);
            app.Button11.Position = [7 13 100 22];
            app.Button11.Text = 'Clear messages';
        end
    end

    methods (Access = public)

        % Construct app
        function app = dPCA()

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
