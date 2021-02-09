
# get check point reads count
.get.check.points.reads.count<- function(ibam,anno,bam,check_points,PARAMETERS){
  
  # Determining scanBam Parameters
  if(PARAMETERS$PAIRED){
    flag = scanBamFlag(isProperPair = T, isFirstMateRead = T)
  } else { # PARAMETERS$PAIRED = FALSE
    flag = scanBamFlag()
  }
  
  # prepare bam parameters
  which <- IRangesList(chr_unclear=IRanges(anno$left, anno$right))
  names(which)=anno$chr
  what = c("strand", "pos", "mapq", "qwidth") # , "isize")
  param <- ScanBamParam(which=which, what=what, flag = flag)
  
  # read bam file
  ba <- scanBam(bam[ibam], param=param)
  pos=ba[[1]]$pos-anno$left+1
  strand=ba[[1]]$strand
  mapq=ba[[1]]$mapq
  qwidth = ba[[1]]$qwidth
  # isize = ba[[1]]$isize
  
  # Filtering by Strandedness
  if(PARAMETERS$STRANDED == "forward"){
    ID = which(strand == anno$strand)
  } else if (PARAMETERS$STRANDED == "reverse"){
    ID = which(strand != anno$strand & strand != "*")
  } else {
    ID = 1:length(strand)
  }
  pos = pos[ID]
  strand = strand[ID]
  mapq = mapq[ID]
  qwidth = qwidth[ID]
  # isize = isize[ID]
  
  # Filtering by mapq
  mapq[which(is.na(mapq))]=255
  ID=which(mapq>=PARAMETERS$MINIMAL_MAPQ)
  pos=pos[ID]
  strand=strand[ID]
  qwidth = qwidth[ID]
  # isize = isize[ID]
  
  # Filtering negative pos
  ID=which(pos>0)
  pos=pos[ID]
  strand=strand[ID]
  qwidth = qwidth[ID]
  # isize = isize[ID]
  
  # convert pos into rna
  rna_pos=anno$DNA2RNA[pos]
  on_rna_id=which( ((rna_pos>0) + !is.na(rna_pos)) ==2)
  rna_pos=rna_pos[on_rna_id]
  strand=strand[on_rna_id]
  qwidth = qwidth[on_rna_id]
  # isize = isize[on_rna_id]
  
  if(!PARAMETERS$PAIRED){
    
    # divide into strand
    pos_ID=which(strand=="+")
    neg_ID=which(strand=="-")
    pos_pos=rna_pos[pos_ID]
    neg_pos=rna_pos[neg_ID]
    
    # shift
    pos_pos=pos_pos+round(PARAMETERS$FRAGMENT_LENGTH/2);
    neg_pos=neg_pos+PARAMETERS$READ_LENGTH-round(PARAMETERS$FRAGMENT_LENGTH/2)
    
  } else if (PARAMETERS$PAIRED){
    
    # divide into strand
    pos_ID=which(strand=="+")
    neg_ID=which(strand=="-")
    pos_pos=rna_pos[pos_ID]
    neg_pos=rna_pos[neg_ID]
    qwidth_pos=qwidth[pos_ID]
    qwidth_neg=qwidth[neg_ID]
    # isize_pos=isize[pos_ID]
    # isize_neg=isize[neg_ID]
    
    # shift
    # I really liked the idea of using isize, ut if this is the same as TLEN, in 
    # sam file specifications, this is very confusing. 
    # Also, this doesn't work if we convert from DNA2RNA first
    # pos_pos=pos_pos+qwidth_pos+round(abs(isize_pos)/2);
    # neg_pos=neg_pos-round(abs(isize_neg)/2)
    pos_pos=pos_pos+qwidth_pos+round(abs(PARAMETERS$FRAGMENT_LENGTH)/2);
    neg_pos=neg_pos-round(abs(PARAMETERS$FRAGMENT_LENGTH)/2)
    
  }
  
  # merge the two
  pos=c(pos_pos,neg_pos)
  pos=pos[pos > 0]
  pos=pos[pos < anno$exome_length]
  
  # get direct count
  check_points_count=check_points
  pos_table=table(pos)
  
  # smooth
  if (PARAMETERS$REMOVE_LOCAL_TAG_ANOMALITIES==TRUE) {pos_table=.remove.local.anomalities(pos_table)}
  
  # get count
  pos_mapped = as.numeric(names(pos_table))
  for (i in 1:length(check_points)) { 
    ID=which(abs(check_points[i]-pos_mapped)*2 < PARAMETERS$WINDOW_WIDTH)
    check_points_count[i]=sum(pos_table[ID])
  }
  
  # return result
  return(check_points_count)
}