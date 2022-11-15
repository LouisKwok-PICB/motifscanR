get_motif_ix <- function(pwms, seqs, genome, p.cutoff, thread, random.seed, loc.file){
  cutoff_motifs_matrix <- NULL
  if (file_test("-f", loc.file)) {
    message("Found motif score cutoff matrix.\n")
    load(loc.file)
  }else{
    message("Not found motif score cutoff matrix.\n")
    cutoff_motifs_matrix <- motifs_cutoff(pwms, genome, random.seed=random.seed)
    save(cutoff_motifs_matrix, file = loc.file)
  }
  message("Scanning motifs...\n")
  motif_sites <- scan_motif(cutoff_motifs_matrix, seqs, strand = 'both', out = "matches",
                            p_value = p.cutoff, remove_dup = TRUE, thread = thread)
  motif_sites_length <- matrix(sapply(motif_sites, length), ncol = dim(motif_sites)[2])
  return(motif_sites_length)
}

get_motif_ix_plus <- function(pwms, seqs, genome, p.cutoff, thread, random.seed, loc.file){
  cutoff_motifs_matrix <- NULL
  if (file_test("-f", loc.file)) {
    message("Found motif score cutoff matrix.\n")
    load(loc.file)
  }else{
    message("Not found motif score cutoff matrix.\n")
    cutoff_motifs_matrix <- motifs_cutoff(pwms, genome, random.seed=random.seed)
    save(cutoff_motifs_matrix, file = loc.file)
  }
  message("Scanning motifs...\n")

  motif_sites <- scan_motif(cutoff_motifs_matrix, seqs, strand = 'both', out = "scores",
                            p_value = p.cutoff, remove_dup = FALSE, thread = thread)

  motif_sites_score <- matrix(sapply(motif_sites, function(x){
    if (length(x) > 0) {
      max(x)
    }else{
      NA
    }
  }), ncol = dim(motif_sites)[2])
  motif_sites_length <- matrix(sapply(motif_sites, length), ncol = dim(motif_sites)[2])
  return(list(motifScore=motif_sites_score,
              motifScans=motif_sites_length,
              motifCount=motif_sites_length))
}


get_motif_positions <- function(pwms, seqs, genome, p.cutoff, thread, random.seed, loc.file){
  cutoff_motifs_matrix <- NULL
  if (file_test("-f", loc.file)) {
    message("Found motif score cutoff matrix.\n")
    load(loc.file)
  }else{
    message("Not found motif score cutoff matrix.\n")
    cutoff_motifs_matrix <- motifs_cutoff(pwms, genome, random.seed=random.seed)
    save(cutoff_motifs_matrix, file = loc.file)
  }
  message("Scanning motifs...\n")
  motif_sites <- scan_motif(cutoff_motifs_matrix, seqs, strand = 'both', out = "positions",
                            p_value = p.cutoff, remove_dup = FALSE, thread = thread)
  motif_sites_scores <- matrix(sapply(motif_sites, function(x){
    if (length(x$scores) > 0) {
      x$scores
      }else{
        NULL
        }
    }), ncol = dim(motif_sites)[2])
  motif_sites_locations <- matrix(sapply(motif_sites, function(x){
    if (length(x$locations) > 0) {
      x$locations
    }else{
      NULL
    }
  }), ncol = dim(motif_sites)[2])
  motif_sites_strands <- matrix(sapply(motif_sites, function(x){
    if (length(x$strands) > 0) {
      x$strands
    }else{
      NULL
    }
  }), ncol = dim(motif_sites)[2])
  return(list(motifScores=motif_sites_scores,
              motifLocations=motif_sites_locations,
              motifStrands=motif_sites_strands))
}


# --------scan motif for each sequence------------
scan_motif <- function(pwms, seqs, strand = c('+', '-', 'both'),
                       out = c("matches", "scores", "positions"),
                       p_value, remove_dup = TRUE, thread = 1){
  if(any(unlist(lapply(pwms, function(x){is.na(x@score[as.character(as.name(p_value))])})))){
    stop(paste0("PWM has no motif score cutoff set for P-value ", p_value))
  }
  cl <- makeCluster(thread)
  motif_sites <- parallel::parSapply(cl, pwms, function(pwm){scan_pwm(seqs, pwm, strand,
                                                                      p_value, remove_dup, out)})
  stopCluster(cl)
  return(motif_sites)
}

scan_pwm <- function(seqs, pwm, strand, p_value,
                     remove_dup, out){
  if (strand %in% c('both', '+')){
    sites_fwd <- scan_pwm_by_strand(seqs, pwm=pwm, strand='+', p_value,
                                    remove_dup, out)
  }
  if (strand %in% c('both', '-')){
    sites_rev <- scan_pwm_by_strand(seqs, pwm=pwm, strand='-', p_value,
                                    remove_dup, out)
  }
  if (out == "positions") {
    if (strand == '+'){
      sites <- sites_fwd
    }else{
      if(strand == '-'){
        sites <- sites_rev
      }else{
        sites <- lapply(1:length(seqs), function(seqidx){
          scores <- c(sites_fwd[[seqidx]]$scores, sites_rev[[seqidx]]$scores)
          locations <- c(sites_fwd[[seqidx]]$locations, sites_rev[[seqidx]]$locations)
          strands <- c(sites_fwd[[seqidx]]$strands, sites_rev[[seqidx]]$strands)
          return(list(scores=scores,
                      locations=locations,
                      strands=strands))
        })
      }
    }
    return(sites)
  }else{
    if (strand == '+'){
      sites <- sites_fwd
    }else{
      if(strand == '-'){
        sites <- sites_rev
      }else{
        sites <- lapply(1:length(seqs), function(seqidx){
          c(sites_fwd[[seqidx]], sites_rev[[seqidx]])
        })
      }
    }
    return(sites)
  }
}


scan_pwm_by_strand <- function(sequences, pwm, strand,
                               p_value, remove_dup, out){
  if (strand == '+') {
    sequences <- sequences
  }else{
    if (strand == '-') {
      sequences <- as.character(reverseComplement(DNAStringSet(sequences)))
    }else{
      stop(paste0("invalid strand option: ", strand))
    }
  }
  score_cutoff <- pwm@score[as.character(as.name(p_value))]
  sliding_scores <- sliding_motif_score(
    as.matrix(pwm), length(pwm), MaxScore(pwm), sequences)
  sites <- pick_motif_max(
    sliding_scores, score_cutoff, strand=strand, out=out)
  if (remove_dup){
    sites <- lapply(sites, function(x){
      if(length(x) > 0){
        max(x)
      }else{
        x
      }})
  }
  return(sites)
}


pick_motif_max <- function(sliding_scores, cutoff,
                           strand='+', out=c("matches", "scores", "positions")){
  if (!(strand %in% c('+', '-'))) {
    stop(paste0("expect '+' or '-' for strand, got ", strand))
  }
  sites <- lapply(seq_along(sliding_scores), function(region_idx){
    scores <- sliding_scores[[region_idx]]
    if (strand == '-') {
      scores <- rev(scores)
    }
    sites_by_region <- scores[scores >= cutoff]
    if (out == "positions") {
      sites_by_region_loc <- which(scores >= cutoff)
      sites_by_region_strand <- rep(strand, length(sites_by_region_loc))
      return(list(scores=sites_by_region,
                  locations=sites_by_region_loc,
                  strands=sites_by_region_strand))
    }
    return(sites_by_region)
  })
  return(sites)
}


# ----- motif cutoff matrix -------
motifs_cutoff <- function(pwms, genome,
                          random_n=1000000, random.seed=NULL){
  if (random_n < 100) {
    stop("each motif must have at least 100 sampling scores")
  }
  message("Generating motifs score cutoff matrix...\n")
  max_len <- do.call(max, lapply(pwms, length))
  message("Selecting random sequences from genome...")
  random_seqs <- random_genome(random_n, genome, max_len, random.seed)
  rc_random_seqs <- as.character(reverseComplement(DNAStringSet(random_seqs)))
  pwms_cutoff <- PWMCutoffList(pwms)
  message('Begin to scan motif score...')
  start_time <- Sys.time()
  sampling_scores <- t(as.matrix(as.data.frame(lapply(pwms_cutoff@listData,
                                                      function(pwm){
                                                        pwm_name <- paste0(c(pwm@ID, pwm@name), collapse = '_')
                                                        message(paste('Scanning', pwm_name,
                                                                      paste0('[', match(pwm_name, names(pwms_cutoff)), '/', length(pwms_cutoff), ']')))
                                                        pwm_matrix <- as.matrix(pwm)
                                                        scores <- motif_score(pwm_matrix, length(pwm),
                                                                              MaxScore(pwm), random_seqs)
                                                        rc_scores <- motif_score(pwm_matrix, length(pwm),
                                                                                 MaxScore(pwm), rc_random_seqs)
                                                        score_max <- apply(rbind(scores, rc_scores), 2, max)
                                                        return(score_max)
                                                      }))))
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time)
  message(paste("Took", round(duration[[1]], 2),  units(duration), "to run.\n"))
  pwms_add_score <- setCutoff(pwms_cutoff, sampling_scores)
  return(pwms_add_score)
}


# ------ random select n seqs from genome ------
random_genome <- function(n_select, genome, len, random.seed = NULL){
  set.seed(random.seed)

  message("Calculating chromosome weight...")
  start_time <- Sys.time()
  if (class(genome) == "FaFile") {
    chorm_length <- seqlengths(genome)
    chrom_sizes_sum <- sum(chorm_length)
    chrom_weight <- chorm_length / chrom_sizes_sum
  }else{
    chorm_length <- unlist(lapply(as.list(genome), length))
    chrom_sizes_sum <- sum(chorm_length)
    chrom_weight <- chorm_length / chrom_sizes_sum
  }
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time)
  message(paste("Took", round(duration[[1]], 2),  units(duration), "to run.\n"))

  message("Generating random sequences...")
  random_start <- Sys.time()
  random_chroms <- sample(names(chrom_weight), n_select,
                          prob = chrom_weight, replace = TRUE)

  seqs_loc <- data.frame(seqnames=random_chroms, start=0)
  for(chrom in unique(seqs_loc$seqnames)){
    seqs_loc_idx <- rownames(seqs_loc[seqs_loc$seqnames == chrom,])
    seqs_loc[seqs_loc_idx,]$start <- sample(1:(chorm_length[chrom] - len), length(seqs_loc_idx),
                                            replace = TRUE)
  }
  seqs_loc$end <- seqs_loc$start + len - 1
  seqs_gr <- GRanges(seqs_loc)
  seqs <- as.character(BSgenome::getSeq(genome, seqs_gr))
  seqs <- seqs[letterFrequency(DNAStringSet(seqs), 'N', as.prob = TRUE) == 0]

  while(length(seqs) < n_select){
    add_random_chrom_idx <- sample(1:n_select, n_select - length(seqs), replace=TRUE)
    add_seqs_loc <- data.frame(seqnames=random_chroms[add_random_chrom_idx], start=0)
    for(chrom in unique(add_seqs_loc$seqnames)){
      add_seqs_loc_idx <- rownames(add_seqs_loc[add_seqs_loc$seqnames == chrom,])
      add_seqs_loc[add_seqs_loc_idx,]$start <- sample(1:(chorm_length[chrom] - len), length(add_seqs_loc_idx),
                                                      replace = TRUE)
    }
    add_seqs_loc$end <- add_seqs_loc$start + len - 1
    add_seqs_gr <- GRanges(add_seqs_loc)
    add_seqs <- as.character(BSgenome::getSeq(genome, add_seqs_gr))
    add_seqs <- add_seqs[letterFrequency(DNAStringSet(add_seqs), 'N', as.prob = TRUE) == 0]

    seqs <- c(seqs, add_seqs)
  }

  random_end <- Sys.time()
  random_duration <- difftime(random_end, random_start)
  message(paste("Took", round(random_duration[[1]], 2),  units(random_duration), "to run.\n"))
  return(seqs)
}


