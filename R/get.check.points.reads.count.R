

# Write a new get.check.points.reads.count
# Convert checkpoints back to DNA coordinates (with splits on the bins that overlap both)
# Count reads using rsamtools
# - unique fragments (by read name)
# - correct strand (first mate on the opposite strand, second mate on reverse strand)
# - mapping quality/unpaired reads/multimapped reads

# get check point reads count
.get.check.points.reads.count<- function(ibam,anno,bam,check_points,PARAMETERS){
  
  
  # RNA2DNA
  RNA2DNA = anno$left:anno$right
  RNA2DNA = RNA2DNA[anno$DNA2RNA > 0]
  
  # Introns
  introns = GenomicRanges::GRanges(seqnames = anno$chr, IRanges::IRanges(which(anno$DNA2RNA == 0)+anno$left), strand = anno$strand)
  
  # Identifying GRanges
  .which = GenomicRanges::GRanges(seqnames = anno$chr, IRanges::IRanges(check_points, check_points), strand = anno$strand)
  .which = GenomicRanges::resize(.which, width = PARAMETERS$WINDOW_WIDTH, fix = "center")
  end(.which) = ifelse(end(.which) > anno$exome_length, anno$exome_length, end(.which))
  end(.which) = RNA2DNA[end(.which)]
  start(.which) = ifelse(start(.which) < 0, 1, start(.which))
  start(.which) = RNA2DNA[start(.which)]
  
  # Removing Introns
  .which = S4Vectors::split(.which, 1:length(.which))
  split.vec = sort(rep(1:length(.which), length(introns)))
  introns.rep = rep(introns, length(.which))
  introns.rep = S4Vectors::split(introns.rep, split.vec)
  
  .which = GenomicRanges::setdiff(.which, introns.rep)
  .which = unlist(.which)
  
  # Determining which strand to use
  ID.strand = structure(c("+", "-"), names = c("-", "+")) 
  if(PARAMETERS$STRANDED == "forward"){
    ID.strand.m1 = anno$strand
    ID.strand.m2 = ID.strand[anno$strand]
  } else if (PARAMETERS$STRANDED == "reverse"){
    ID.strand.m1 = ID.strand[anno$strand]
    ID.strand.m2 = anno$strand
  } else{
    ID.strand.m1 = ID.strand.m2 = ID.strand
  }
  
  # Importing bam file
  if(PARAMETERS$PAIRED){
    
    # Flags
    flag.m1 = scanBamFlag(isProperPair = T, isFirstMateRead = T)
    flag.m2 = scanBamFlag(isProperPair = T, isSecondMateRead = T)
    
    # What
    what = c("qname", "strand", "pos", "mapq")
    
    # Importing Bams
    param.m1 = ScanBamParam(which=.which, what=what, flag = flag.m1)
    bam.m1 <- scanBam(bam[ibam], param=param.m1)
    
    param.m2 = ScanBamParam(which=.which, what=what, flag = flag.m2)
    bam.m2 <- scanBam(bam[ibam], param=param.m2)
    
    bin.indices = names(.which)
    reads.names = list()
    
    # Counting Unique Read IDs
    for(i in 1:length(.which)){
      
      # First mates
      qname.m1 = bam.m1[[i]]$qname
      strand.m1 = bam.m1[[i]]$strand
      mapq.m1 = bam.m1[[i]]$mapq
      qname.m1 = qname.m1[strand.m1 %in% ID.strand.m1 & mapq.m1 >= 255]
      
      # Second mates
      qname.m2 = bam.m2[[i]]$qname
      strand.m2 = bam.m2[[i]]$strand
      mapq.m2 = bam.m2[[i]]$mapq
      qname.m2 = qname.m2[strand.m2 %in% ID.strand.m2 & mapq.m2 >= 255]
      
      reads.names[[bin.indices[i]]] = unique(c(reads.names[[bin.indices[i]]], qname.m1, qname.m2))
    }
    
    
  } else { # PARAMETERS$PAIRED = FALSE
    
    # Flag
    flag = scanBamFlag()
    
    # What
    what = c("qname", "strand", "pos", "mapq")
    
    # Importing Bams
    param = ScanBamParam(which=.which, what=what, flag = flag)
    bam.m1 <- scanBam(bam[ibam], param=param)
    
    bin.indices = names(.which)
    reads.names = list()
    
    # Counting Unique Read IDs
    for(i in 1:length(.which)){
      
      # First mates
      qname.m1 = bam.m1[[i]]$qname
      strand.m1 = bam.m1[[i]]$strand
      mapq.m1 = bam.m1[[i]]$mapq
      qname.m1 = qname.m1[strand.m1 %in% ID.strand.m1 & mapq.m1 >= 255]
      
      reads.names[[bin.indices[i]]] = unique(c(reads.names[[bin.indices[i]]], qname.m1))
    }
    
  }
  
  check_points_count = unlist(lapply(reads.names, length))
  names(check_points_count) = check_points
  
  if (PARAMETERS$REMOVE_LOCAL_TAG_ANOMALITIES==TRUE) {check_points_count=.remove.local.anomalities(check_points_count)}
  
  # return result
  return(check_points_count)
}