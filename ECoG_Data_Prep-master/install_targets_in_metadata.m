function metadata = install_targets_in_metadata(metadata, TARGET_DIR)
t_type_fmt = {'category','embedding','similarity'};
    n_tf = numel(t_type_fmt);
    for i_tf = 1:n_tf
        type_fmt = t_type_fmt{i_tf};
        switch type_fmt
            case 'category'
                t_type_categories = list_files(fullfile(TARGET_DIR,type_fmt));
                n_c = numel(t_type_categories);
                for i_c = 1:n_c
                    category_file = fullfile(TARGET_DIR,type_fmt,t_type_categories{i_c});
                    category_label = strip_extension(t_type_categories{i_c});
                    fid = fopen(category_file);
                    tmp = textscan(fid, '%s %u8','Delimiter',',');
                    category_labels = tmp{1};
                    category_targets = tmp{2};
                    fclose(fid);
                    metadata = installCategoryStructure(metadata, category_targets, category_labels, category_label, []);
                end
            case {'embedding','similarity'}
                t_type_sims = list_dirs(fullfile(TARGET_DIR,type_fmt));
                n_ts = numel(t_type_sims);
                for i_ts = 1:n_ts
                    type_sim = t_type_sims{i_ts};
                    t_sim_sources = list_dirs(fullfile(TARGET_DIR,type_fmt,type_sim));
                    n_s = numel(t_sim_sources);
                    for i_s = 1:n_s
                        source = t_sim_sources{i_s};
                        t_source_metrics = list_files(fullfile(TARGET_DIR,type_fmt,type_sim,source));
                        t_source_metrics = t_source_metrics(~strcmp('labels.txt',t_source_metrics));
                        n_m = numel(t_source_metrics);
                        for i_m = 1:n_m
                            metric_file = t_source_metrics{i_m};
                            metric_label = strip_extension(metric_file);
                            structure_file = fullfile(TARGET_DIR,type_fmt,type_sim,source,metric_file);
                            structure_label_file = fullfile(TARGET_DIR,type_fmt,type_sim,source,'labels.txt');
                            structure_matrix = csvread(structure_file);
                            fid = fopen(structure_label_file);
                            tmp = textscan(fid, '%s');
                            structure_labels = tmp{1};
                            fclose(fid);
                            metadata = installSimilarityStructure(metadata, structure_matrix, structure_labels, type_sim, source, metric_label, type_fmt);
                        end
                    end
                end
        end
    end
end