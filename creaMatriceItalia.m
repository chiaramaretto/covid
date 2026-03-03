function [contact_matrix, labels] = creaMatriceItalia()    

    cont_file = '2008_Mossong_POLYMOD_contact_common.csv';    
    part_file = '2008_Mossong_POLYMOD_participant_common.csv';
   
    opts = detectImportOptions(part_file);
    opts.VariableNamingRule = 'preserve';
    participants = readtable(part_file, opts);

    opts_c = detectImportOptions(cont_file);
    opts_c.VariableNamingRule = 'preserve';
    contacts = readtable(cont_file, opts_c);

    italy_parts = participants(participants.part_id >= 4001 & participants.part_id <= 5000, :);
    italy_conts = contacts(ismember(contacts.part_id, italy_parts.part_id), :);

    c_age = italy_conts.cnt_age_exact;

    if iscell(c_age) || isstring(c_age)
        c_age = str2double(string(c_age)); 
    end

    c_min = italy_conts.cnt_age_est_min;
    if iscell(c_min) || isstring(c_min), c_min = str2double(string(c_min)); end
    
    c_max = italy_conts.cnt_age_est_max;
    if iscell(c_max) || isstring(c_max), c_max = str2double(string(c_max)); end

    nan_idx = isnan(c_age);
    if any(nan_idx)
        c_age(nan_idx) = (c_min(nan_idx) + c_max(nan_idx)) / 2;
    end
    italy_conts.age_clean = c_age;

    edges = [0, 10, 20, 30, 40, 50, 60, 70, 80, 120]; 
    labels = {'0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'};

    part_grp = discretize(italy_parts.part_age, edges);
    
    cont_grp = discretize(italy_conts.age_clean, edges);
    grp_summary = groupsummary(table(part_grp), 'part_grp', 'IncludeEmptyGroups', true);
    
    N_i = grp_summary.GroupCount(1:length(labels));
    [~, loc] = ismember(italy_conts.part_id, italy_parts.part_id);
    part_grp_of_cont = part_grp(loc);

    C_ij = zeros(length(labels));
    
    for k = 1:height(italy_conts)
        i = part_grp_of_cont(k); % Riga (Partecipante)
        j = cont_grp(k);         % Colonna (Contatto)
        
        if ~isnan(i) && ~isnan(j)
            C_ij(i,j) = C_ij(i,j) + 1;
        end
    end

    contact_matrix = 4*C_ij ./ N_i; 
    
    contact_matrix(isnan(contact_matrix) | isinf(contact_matrix)) = 0;
end