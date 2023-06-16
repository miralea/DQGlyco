#### functions ####
fasta_to_df <- function(fasta_file) {
  
  # read fasta file
  fa_content <- seqinr::read.fasta(file = fasta_file, seqtype = "AA")
  # parse as a df
  fa_df <- pblapply(fa_content, function(prot_seq) {
    prot_id <- seqinr::getAnnot(prot_seq) %>%
      str_split("\\|", simplify = TRUE) %>%
      nth(., 2)
    prot_seq_df <- data.frame(uniprot_id = prot_id, 
                              AA = as.character(prot_seq)) %>%
      mutate(position = 1:n())
    return(prot_seq_df)
  }) %>%
    bind_rows()
  return(fa_df)
  
}

# get glycosilable asparagines given a fasta df
get_glycosilable_asparagines <- function(fasta_df) {
  
  glycosilable_df <- split(fasta_df, fasta_df$uniprot_id) %>%
    pblapply(. , function(int_df) {
      sequence <- paste0(int_df$AA, collapse = "")
      int_ptm_starts <- str_locate_all(sequence, "N[^P][S|T]")[[1]][, "start"]
      annotated_df <- int_df %>%
        mutate(glycosilable_asparagine = ifelse(position %in% int_ptm_starts, 1, 0))
      return(annotated_df)
    }) %>%
    bind_rows() 
  return(glycosilable_df)
  
}

# retrieves the position of interesting residues in a given sequence
get_mid_position <- function(seq, int_residues) {
  all_positions <- str_locate_all(seq, int_residues)[[1]][,"start"]
  int_position <- all_positions[ceiling(length(all_positions)/2)]
  int_residue <- str_sub(seq, start = int_position, end = int_position)
  out_df <- data.frame(AA = int_residue, position = int_position)
  if(nrow(out_df) == 0) return(data.frame(AA = "NOTFOUND", position = 0))
  return(out_df)
}

# annotates glycan type given from tidy df
annotate_glycan_type <- function(int_df) {
  
  out_df <- int_df %>%
    separate(glycan_content, into = c("hexnac", "hexnac_content", "hex", "hex_content", "end"), sep = "\\(|\\)", extra = "merge", remove = FALSE) %>%
    mutate(hexnac_content = as.numeric(hexnac_content),
           hex_content = as.numeric(hex_content),
           end = ifelse(is.na(end), "", end),
           AA = AA,
           glycan_all = 1) %>% 
    mutate(glycan_type = case_when(
      grepl("Phospho", glycan_content) ~ "phospho",
      grepl("NeuAc", glycan_content) ~ "sialylated",
      grepl("Fuc", glycan_content) ~ "fucosylated",
      hexnac_content == 2 & hex_content > 3  & end == "" ~ "high mannose",
      hexnac_content == 2 & hex_content %in% 1:3 & end == "" ~ "paucimannose",
      hexnac_content <= 2 & hex == "" & end == "" ~ "small",
      TRUE ~ "complex/hybrid"
    ))
  
  return(out_df)
  
}

process_nglyco_excel <- function(ptm_data, ref_proteome_df, custom_breaks = NULL) {
  
  annotated_data <- ptm_data %>%
    dplyr::filter(modification == "glycosylation") %>%
    dplyr::select(uniprot_id = Protein.ID, mod_pep = Peptide,
                  n_position = n_position, obs_mod = Observed.Modifications,
                  glycan_type = glycan_type) %>%
    mutate(glycan_all = 1, AA = "N") %>%
    separate(obs_mod, into = c("glycan_content", "modification_mass"), sep = " % ") 
  
  # parse proteome as a df and subset only to proteins that are detected
  nglyco_proteome <- ref_proteome_df %>%
    subset(uniprot_id %in% unique(annotated_data$uniprot_id))
  
  # find glycosilable asparagines
  glycosilable_df <- get_glycosilable_asparagines(fasta_df = nglyco_proteome) %>%
    dplyr::filter(glycosilable_asparagine == 1)
  
  # remove glycopeptide level, keeping the glycoform data
  glycoform_data <- annotated_data %>%
    dplyr::transmute(uniprot_id, AA, glycan_content, 
                     modification_mass, position = n_position, glycan_all,
                     glycan_type) %>%
    ## MULTIPLE GLYCOPEPTIDES WITH SAME GLYCOFORM COLLAPSED HERE
    distinct()
  
  # create bins of glycoforms
  per_site_n_forms <- glycoform_data %>%
    group_by(uniprot_id, AA, position) %>%
    summarise(n_forms = n()) %>%
    ungroup()
  
  if(is.null(custom_breaks)) {
    qs <- c(0, quantile(per_site_n_forms$n_forms, seq(0,1,length.out = 4)))
  } else {
    qs <- custom_breaks
  }
  per_site_quant_nforms <- per_site_n_forms %>%
    mutate(n_form_cat = cut(per_site_n_forms$n_forms, breaks = qs)) %>%
    dplyr::transmute(uniprot_id,AA, position, nform_range = paste0("glycan_forms_", n_form_cat), val = 1) %>%
    pivot_wider(names_from = nform_range, values_from = val) %>%
    replace(is.na(.), 0)
  
  # collapse to site level
  site_level_df <- glycoform_data %>%
    dplyr::transmute(uniprot_id, AA, position, 
                     glycan_type = paste0("glycan_type_", glycan_type)) %>%
    pivot_longer(-c(uniprot_id, AA, position)) %>%
    dplyr::select(-name) %>%
    distinct() %>%
    mutate(sub = 1) %>%
    pivot_wider(id_cols = c(uniprot_id, AA, position), names_from = value, values_from = sub, values_fill = 0) %>%
    left_join(., per_site_quant_nforms, by = c("uniprot_id", "AA", "position")) %>%
    mutate(glycan_all = 1) %>%
    full_join(., glycosilable_df, by = c("uniprot_id", "AA", "position")) %>%
    replace(is.na(.), 0) %>%
    mutate(is_asparagine = ifelse(AA == "N", 1, 0)) %>%
    mutate(ng_asparagine = ifelse(glycan_all == 0 & is_asparagine == 1, 1, 0),
           ng_glycosilable_asparagine = ifelse(glycan_all == 0 & glycosilable_asparagine == 1 , 1, 0)) 
  
  # create long format dataframe
  to_exclude_ptms <- c("ng_asparagine", "is_asparagine", "glycosilable_asparagine", "glycosilable_st")
  ref_ptm <- "ng_glycosilable_asparagine"
  site_level_long_df <- site_level_df %>%
    pivot_longer(-c(uniprot_id, AA, position), names_to = "ptm", values_to = "exists") %>%
    dplyr::filter(exists != 0) %>%
    dplyr::filter(! ptm %in% to_exclude_ptms) %>%
    dplyr::select(-c(exists)) %>%
    mutate(ptm_cat = case_when(grepl("_all", ptm) ~ "all", 
                               grepl("_type", ptm) ~ "glycan_type", 
                               grepl("_forms", ptm) ~ "glycan_nforms", 
                               grepl("_hit|_nohit", ptm) ~ "hits",
                               grepl(ref_ptm, ptm) ~ "ref"
    )) %>%
    mutate(ptm_cat = fct_relevel(ptm_cat, c("ref","all", "glycan_type", "glycan_nforms"))) %>%
    mutate(raw_nforms = as.numeric(str_replace_all(ptm, "glycan_forms_\\(|,[0-9]+\\]", ""))) %>%
    mutate(ptm = fct_reorder(ptm, raw_nforms)) %>%
    dplyr::select(-raw_nforms) %>%
    arrange(ptm_cat, ptm)
  
  
  # create output list with different dataframes
  out_list <- list(glycoform_df = glycoform_data, 
                   glycosite_wide = site_level_df, 
                   glycosite_long = site_level_long_df)
  
  return(out_list)

}

# process N Glyco tables
process_nglyco_table <- function(input_file, ref_proteome_df, custom_breaks = NULL) {
  
  # read ptm data
  ptm_data <- read_csv(input_file, name_repair = make.names)
  
  # tidy data
  annotated_data <- ptm_data %>%
    dplyr::select(uniprot_id = Protein.ID, mod_pep = Modified.Peptide, 
                  start_position = Protein.Start, assign_mod = Assigned.Modifications, 
                  obs_mod = Observed.Modifications) %>%
    mutate(n_position = str_extract(assign_mod, "[0-9]+N\\(") %>% 
             str_extract("[0-9]+") %>% as.numeric()) %>%
    mutate(AA = "N") %>%
    mutate(n_position = n_position - 1) %>%
    mutate(n_position = start_position + n_position) %>%
    dplyr::select(-assign_mod, -start_position) %>%
    separate(obs_mod, into = c("glycan_content", "modification_mass"), sep = " % ") 
  
  # parse proteome as a df and subset only to proteins that are detected
  nglyco_proteome <- ref_proteome_df %>%
    subset(uniprot_id %in% unique(annotated_data$uniprot_id))
  
  # find glycosilable asparagines
  glycosilable_df <- get_glycosilable_asparagines(fasta_df = nglyco_proteome) %>%
    dplyr::filter(glycosilable_asparagine == 1)
  
  # annotate glycan types
  glycan_type_data <- annotated_data %>%
    annotate_glycan_type()
  
  # remove glycopeptide level, keeping the glycoform data
  glycoform_data <- glycan_type_data %>%
    dplyr::transmute(uniprot_id, AA, glycan_content, 
                     modification_mass, position = n_position, glycan_all,
                     glycan_type) %>%
    ## MULTIPLE GLYCOPEPTIDES WITH SAME GLYCOSITE AND GLYCOFORM COLLAPSED HERE
    distinct()
  
  # create bins of glycoforms
  per_site_n_forms <- glycoform_data %>%
    group_by(uniprot_id, AA, position) %>%
    summarise(n_forms = n()) %>%
    ungroup()
  
  if(is.null(custom_breaks)) {
    qs <- c(0, quantile(per_site_n_forms$n_forms, seq(0,1,length.out = 4)))
  } else {
    qs <- custom_breaks
  }
  per_site_quant_nforms <- per_site_n_forms %>%
    mutate(n_form_cat = cut(per_site_n_forms$n_forms, breaks = qs)) %>%
    dplyr::transmute(uniprot_id,AA, position, nform_range = paste0("glycan_forms_", n_form_cat), val = 1) %>%
    pivot_wider(names_from = nform_range, values_from = val) %>%
    replace(is.na(.), 0)
  
  # collapse to site level
  site_level_df <- glycoform_data %>%
    dplyr::transmute(uniprot_id, AA, position, 
                     glycan_type = paste0("glycan_type_", glycan_type)) %>%
    pivot_longer(-c(uniprot_id, AA, position)) %>%
    dplyr::select(-name) %>%
    distinct() %>%
    mutate(sub = 1) %>%
    pivot_wider(id_cols = c(uniprot_id, AA, position), names_from = value, values_from = sub, values_fill = 0) %>%
    left_join(., per_site_quant_nforms, by = c("uniprot_id", "AA", "position")) %>%
    mutate(glycan_all = 1) %>%
    full_join(., glycosilable_df, by = c("uniprot_id", "AA", "position")) %>%
    replace(is.na(.), 0) %>%
    mutate(is_asparagine = ifelse(AA == "N", 1, 0)) %>%
    mutate(ng_asparagine = ifelse(glycan_all == 0 & is_asparagine == 1, 1, 0),
           ng_glycosilable_asparagine = ifelse(glycan_all == 0 & glycosilable_asparagine == 1 , 1, 0)) 
  
  # create long format dataframe
  to_exclude_ptms <- c("ng_asparagine", "is_asparagine", "glycosilable_asparagine", "glycosilable_st")
  ref_ptm <- "ng_glycosilable_asparagine"
  site_level_long_df <- site_level_df %>%
    pivot_longer(-c(uniprot_id, AA, position), names_to = "ptm", values_to = "exists") %>%
    dplyr::filter(exists != 0) %>%
    dplyr::filter(! ptm %in% to_exclude_ptms) %>%
    dplyr::select(-c(exists)) %>%
    mutate(ptm_cat = case_when(grepl("_all", ptm) ~ "all", 
                               grepl("_type", ptm) ~ "glycan_type", 
                               grepl("_forms", ptm) ~ "glycan_nforms", 
                               grepl("_hit|_nohit", ptm) ~ "hits",
                               grepl(ref_ptm, ptm) ~ "ref"
    )) %>%
    mutate(ptm_cat = fct_relevel(ptm_cat, c("ref","all", "glycan_type", "glycan_nforms"))) %>%
    mutate(raw_nforms = as.numeric(str_replace_all(ptm, "glycan_forms_\\(|,[0-9]+\\]", ""))) %>%
    mutate(ptm = fct_reorder(ptm, raw_nforms)) %>%
    dplyr::select(-raw_nforms) %>%
    arrange(ptm_cat, ptm)
  
  # create output list with different dataframes
  out_list <- list(glycoform_df = glycoform_data, 
                   glycosite_wide = site_level_df, 
                   glycosite_long = site_level_long_df)
  
  return(out_list)
  
}

# annotate site long DF
annotate_site_long_df <- function(site_long_df, topo_file, 
                                  afold_file, sift_file, anno_package) {
  
  # read annotation files
  topo_df <- read_csv(topo_file) %>%
    dplyr::filter(description %in% c("Extracellular", "Cytoplasmic"))
  afold_df <- read_csv(afold_file)
  sift_df <- read_tsv(sift_file) 
  
  site_unique_df <- site_long_df %>%
    dplyr::select(uniprot_id, AA, position) %>%
    distinct() 
  
  annotated_topo <- site_unique_df %>%
    left_join(., topo_df, by = c("uniprot_id" = "uniprot_id")) %>%
    drop_na() %>%
    mutate(is_in = ifelse(position >= start & position <= end, 1, 0)) %>%
    dplyr::select(-c(start,end)) %>%
    distinct() %>%
    group_by(uniprot_id, position, description) %>%
    summarise(is_in = sum(is_in)) %>%
    ungroup() %>%
    pivot_wider(names_from = description, values_from = is_in, values_fill = 0)
  
  annotated_afold <- site_unique_df %>%
    left_join(., afold_df, by = c("uniprot_id" = "protein_id", "AA" = "AA", "position" = "position")) %>%
    drop_na() %>%
    mutate(not_idr_low_acc_5 = case_when(IDR == 0 & low_acc_5 == 1 ~ 1,
                                         TRUE ~ 0),
           not_idr_high_acc_5 = case_when(IDR == 0 & high_acc_5 == 1 ~ 1,
                                          IDR == 0 & high_acc_5 == 0 ~ 0,
                                          TRUE ~ as.numeric(NA)))
  
  annotated_sift <- site_unique_df %>%
    mutate(gene_id = mapIds(x = anno_package, keys = uniprot_id, keytype = "UNIPROT", column = "ENSEMBL")) %>%
    left_join(., sift_df, by = c("gene_id" = "Gene_id", "position" = "Position_of_amino_acid_substitution")) %>%
    dplyr::select(-gene_id) %>%
    group_by(uniprot_id, AA, position) %>%
    summarise(sift_score = mean(SIFT_score, na.rm = TRUE)) %>%
    ungroup() %>%
    drop_na()
  
  site_long_annotated <- site_long_df %>%
    left_join(., annotated_topo, by = c("uniprot_id", "position")) %>%
    left_join(., annotated_afold, by = c("uniprot_id", "AA", "position")) %>%
    left_join(., annotated_sift, by = c("uniprot_id", "AA", "position"))
  
  return(site_long_annotated)
  
}

# ORA from region and long df
region_ora <- function(site_long_df, int_ptm, ref_ptm, int_roi) {
  
  # create cont mat
  cont_mat <- site_long_df %>%
    dplyr::filter(ptm %in% c(int_ptm, ref_ptm)) %>%
    dplyr::select(ptm, !!sym(int_roi)) %>%
    drop_na() %>%
    group_by(ptm, !!sym(int_roi)) %>%
    dplyr::count() %>%
    pivot_wider(names_from = !!sym(int_roi), values_from = n) %>%
    column_to_rownames("ptm") %>%
    as.matrix()
  
  # substitute NAs with ones
  cont_mat[is.na(cont_mat)] <- 1
  
  # ORA
  fisher_test <- fisher.test(cont_mat)
  odds_ratio <- fisher_test$estimate
  p_value <- fisher_test$p.value
  # output df
  out_df <- data.frame(ptm = int_ptm, motif = int_roi, p_value = p_value,
                       n_ptm_in_motif = cont_mat[int_ptm, "1"],
                       n_ptm_out_motif = cont_mat[int_ptm, "0"],
                       n_bg_in_motif = cont_mat[ref_ptm, "1"],
                       n_bg_out_motif = cont_mat[ref_ptm, "0"],
                       log_odds = log(odds_ratio))
  return(out_df)
}

# perform ora from site long annotated df on a given set of ROIs
ora_from_site_long_df <- function(input_annotated_df, int_rois) {
  
  # create ptm_cat df to keep factor orders
  ptm_cat_df <- input_annotated_df %>%
    dplyr::select(ptm, ptm_cat) %>%
    distinct()
  
  # get ref ptm
  ref_ptm <- unique(as.character(input_annotated_df$ptm[input_annotated_df$ptm_cat == "ref"]))
  
  # get all PTMs that are going to be analyzed
  ptms_to_analyze <- levels(input_annotated_df$ptm) %>%
    .[.!= ref_ptm]
  
  # create all the potential combinations of ptms to analyze and ROIs
  comb_df <- expand_grid(ptms_to_analyze, int_rois)
  
  # iterate over the rows of comb df
  ora_results_df <- pbapply(comb_df, 1, function(x) {
    
    int_ptm <- x["ptms_to_analyze"]
    int_roi <- x["int_rois"]
    
    outdf <- region_ora(site_long_df = input_annotated_df, int_ptm = int_ptm,
                        int_roi = int_roi, ref_ptm = ref_ptm)
    
    return(outdf)
    
  }) %>%
    bind_rows()
  
  # add reference and compute fractions
  ptm_fractions_data <- left_join(ptm_cat_df, ora_results_df, by = "ptm") %>%
    mutate(is_sig = ifelse(p_value <= 0.05, "p < 0.05", "p > 0.05")) %>%
    mutate(fraction_ptm = n_ptm_in_motif / (n_ptm_in_motif + n_ptm_out_motif),
           fraction_bg = n_bg_in_motif / (n_bg_in_motif + n_bg_out_motif)) %>%
    dplyr::filter(!is.na(motif))
  
  fraction_background <- ptm_fractions_data %>%
    dplyr::select(motif, fraction_bg) %>%
    transmute(ptm = ref_ptm, motif = motif, fraction_ptm = fraction_bg, p_value = 1, is_sig = NA, ptm_cat = "ref") %>%
    distinct()
  
  out_df <- bind_rows(ptm_fractions_data, fraction_background) %>%
    mutate(ptm = factor(ptm, 
                        levels = levels(ptm_cat_df$ptm)[levels(ptm_cat_df$ptm) %in% ptm]),
           ptm_cat = factor(ptm_cat, 
                            levels = levels(ptm_cat_df$ptm_cat)[levels(ptm_cat_df$ptm_cat) %in% ptm_cat]))
  
  return(out_df)
  
}

# pairwise wilcox comparison
pairwise_wilcox_comparison <- function(int_df, int_col) {
  
  ref_ptm <- unique(as.character(int_df$ptm[int_df$ptm_cat == "ref"]))
  ref_vector <- int_df %>%
    dplyr::filter(ptm == ref_ptm) %>%
    pull(!!sym(int_col))
  ref_median <- median(ref_vector, na.rm = TRUE)
  comp_df <- int_df %>%
    group_by(ptm, ptm_cat) %>%
    summarise(wilcox_p = wilcox.test(x = !!sym(int_col), y = ref_vector)$p.value) %>%
    mutate(is_sig = ifelse(wilcox_p <= 0.05, "p < 0.05", "p > 0.05"))
  
  return(comp_df)
}