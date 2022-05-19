
#' Adds transparency to colour
#'
#' @param cols Vector of colours.
#' @param alphas Single value or vector of alphas
#' @param ... Passed to rgb
#' @return rgb colours with transparency set.
colAlpha = function(cols,alphas,...) {
  if(length(alphas)==1)
    alphas = rep(alphas,length(cols))
  tmp = col2rgb(cols)
  sapply(seq_len(ncol(tmp)),function(e) rgb(tmp[1,e],tmp[2,e],tmp[3,e],alphas[e]*255,maxColorValue=255,...))
}



annotateBTB = function(btb.fp,PDID,minSegLen=1e6,subCl.minSegLen=1e7,tgtChrs=c(1:23,'X'),longFormat=T,removeBalancedSegs = FALSE,method = c('totalCN','allelicRatio')){
  chromInfo = read.delim('/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',sep = '\t')
  
  #### Read in CN profile for both major and minor clone from Battenberg ####
  btb = gsub('summary.csv','subclones.txt.gz',btb.fp)
  btb = read.delim(btb,sep = '\t')
  btb$posID = paste0(btb$chr,'_',btb$startpos,'_',btb$endpos)
  btb$Idx = c(1:nrow(btb))
  
  # Keep CN segments on chr of interest only
  btb = btb[btb$chr %in% tgtChrs,]
  btb$chr = ifelse(btb$chr == 'X',23,as.numeric(btb$chr))
  
  # remove segments <= 1Mb
  btb$segLen = btb$endpos - btb$startpos + 1
  if(sum(btb$segLen <= minSegLen) >0 ){
    message(sprintf('Removing %d short segment for sample %s',sum(btb$segLen <= minSegLen),PDID))  
  }
  btb = btb[btb$segLen > minSegLen,]
  
  # Get major clone vs minor clone CN profile
  btb.majCl = tibble()
  btb.subCl = tibble()
  
  for(i in 1:nrow(btb)){
    if(is.na(btb$frac2_A[i])){
      tmp = btb[i,c('Idx','chr','startpos','endpos','nMaj1_A','nMin1_A','frac1_A','posID','segLen')]
      colnames(tmp) = c('Idx','Chr','Start','Stop','matNum','patNum','frac','posID','segLen')
      btb.majCl = rbind(btb.majCl,tmp)
    }else{
      if(btb$frac1_A[i] > btb$frac2_A[i]){
        tmp.maj = btb[i,c('Idx','chr','startpos','endpos','nMaj1_A','nMin1_A','frac1_A','posID','segLen')]
        tmp.min = btb[i,c('Idx','chr','startpos','endpos','nMaj2_A','nMin2_A','frac2_A','posID','segLen')]
      }else{
        tmp.maj = btb[i,c('Idx','chr','startpos','endpos','nMaj2_A','nMin2_A','frac2_A','posID','segLen')]
        tmp.min = btb[i,c('Idx','chr','startpos','endpos','nMaj1_A','nMin1_A','frac1_A','posID','segLen')]
      }
      colnames(tmp.maj) = c('Idx','Chr','Start','Stop','matNum','patNum','frac','posID','segLen')
      colnames(tmp.min) = c('Idx','Chr','Start','Stop','matNum','patNum','frac','posID','segLen')
      btb.majCl = rbind(btb.majCl,tmp.maj)
      btb.subCl = rbind(btb.subCl,tmp.min)
    }
  }
  
  
  
  btb.majCl$type = 'maj'
  btb.subCl$type = 'sub'
  message(sprintf('%s - Lowest frac for Major clone is %f',PDID,min(btb.majCl$frac)))
  message(sprintf('%s - Lowest frac for Sub clone is %f',PDID,min(btb.subCl$frac)))
  
  # Remove subclone CNA with < 0.1% cells and shorter than 50Mb
  message(sprintf('%s - %d subclone CNA fragment due to low frac',PDID,sum(btb.subCl$frac<0.2)))
  message(sprintf('%s - %d subclone CNA fragment due to short length < 20Mb',PDID,sum(btb.subCl$segLen<2e7)))
  message(sprintf('%s - %d subclone CNA fragment due to short length < 50Mb',PDID,sum(btb.subCl$segLen<5e7)))
  btb.subCl = btb.subCl[btb.subCl$frac >=0.1 & btb.subCl$segLen >= subCl.minSegLen,]
  
  
  
  # Merge clonal CNprofile (if any) together
  all = rbind(btb.subCl,btb.majCl)
  all$tumTot = all$matNum+all$patNum
  all$tumFrac = all$matNum/(all$tumTot)
  all$tot2min = paste0(all$tumTot,':',all$patNum)
  all$newIdx = all$Idx
  
  # Process primary CNprofile (major clone)
  data = all[all$type == 'maj',]
  
  # Smoothing step:
  modSegs = data.frame()
  toRemove = c()
  
  for(chr in unique(data$Chr)){
    chromLen = chromInfo[chromInfo$chrom == chr & chromInfo$arm == 'q',]$end*1000
    chr.data = data[data$Chr == chr,]
    
    
    if(method=='totalCN'){ # to plot total copy number changes (eg. copyKat)
      for(config in unique(chr.data$tot2min)){
        cf.len.perc = sum(chr.data[chr.data$tot2min == config,]$segLen)*100/chromLen 
        #config = unique(chr.data[chr.data$tumFrac == tumFrac,]$config)
        message(sprintf('\nConfig %s for chr %s occupies %f perc',config,chr,cf.len.perc))
        
        if(cf.len.perc >= 90){
          message(sprintf('Merging config %s for chr %s',config,chr))
          breakFlag=1
          removeFlag=0
          toRemove = c(toRemove,unique(chr.data$Idx))
          
          tmp = chr.data[chr.data$tot2min == config,][1,]
          tmp$Idx = paste(chr.data[chr.data$tot2min == config,]$Idx,collapse = '_')
          tmp$newIdx = paste0(tmp$newIdx,'.', 1)
          tmp$Stop = chromLen
          tmp$Start = 1
          tmp$tumFrac = paste(unique(chr.data[chr.data$tot2min == config,]$tumFrac),collapse = '_') 
          
        }else if(cf.len.perc <= 10){ # remove if it's <10 % of Chr
          removeFlag=1
          message(sprintf('Removing config %s for chr %s',config, chr))
          tmp=NULL
          #tmp = data.frame(Idx=0,Chr=chr,Start=1,Stop=chromLen,normTot=2,normMin=1,tumTot = 2,tumMin = 1,segLen=chromLen, config='2:1',idx=0)
          breakFlag=0
        }else{
          # Merge fragments
          a = chr.data[chr.data$tot2min == config,]
          rmIdx.a = c()
          rmIdx.a.all = c()
          newSeg.a = tibble()
          
          if(nrow(a) > 1){
            for(f in 2:nrow(a)){
              if((a$Start[f] - a$Stop[f-1]) <= minSegLen){
                message('Merging fragments')
                rmIdx.a = c(rmIdx.a,a$Idx[c(f,(f-1))])
                rmIdx.a.all = c(rmIdx.a.all,a$Idx[c(f,(f-1))])
              }
              
              # Do the removal
              
              # If 2 fragments are far apart, or if this is the last segment --> stop and perform merger
              if(((a$Start[f] - a$Stop[f-1]) > minSegLen) || f == nrow(a)){
                if(length(rmIdx.a) > 1){
                  a.toRemove = a[a$Idx %in% unique(rmIdx.a),]
                  # Write new merged segment
                  segLen = (max(a.toRemove$Stop) - min(a.toRemove$Start)) + 1
                  patNum = as.numeric(sapply(strsplit(config,split = ':'),'[',2))
                  tumTot = as.numeric(sapply(strsplit(config,split = ':'),'[',1))
                  matNum = tumTot-patNum
                  tmp.a = data.frame(Idx=paste(unique(rmIdx.a),collapse = '_'),
                                     Chr=chr,Start=min(a.toRemove$Start),Stop=max(a.toRemove$Stop),
                                     matNum = matNum,tumTot = tumTot,patNum = patNum,segLen=segLen,type='maj',
                                     posID = '0',frac = paste(a.toRemove$frac,collapse = '_'),
                                     tumFrac = matNum/tumTot,tot2min = paste0(tumTot,':',patNum),newIdx = paste0(unique(rmIdx.a)[1],'.1'))
                  newSeg.a=rbind(newSeg.a,tmp.a)
                  rmIdx.a = c()
                }else if(length(rmIdx.a) == 1){
                  message('OOPSS')
                  print(rmIdx.a)
                  stop()
                }
              }
            } 
            if(length(rmIdx.a.all) > 1){
              toRemove = c(toRemove,unique(rmIdx.a.all))
              tmp = newSeg.a    
            }else{
              tmp=NULL
            }
          }else{
            tmp = NULL
          }
          
          
          breakFlag=0
          removeFlag=0
        }
        if(!is.null(tmp)){
          #tmp$tumFrac = tmp$matNum/(tmp$tumTot)  
          #if(sum(is.na(tmp$matNum))>0 & unique(tmp$patNum)==0){
          #  tmp$tumFrac = 1
          #}else{
          #  
          #}
          
          modSegs = rbind(modSegs,tmp)
        }
        if(removeFlag == 1){
          toRemove = c(toRemove,c(chr.data[chr.data$tot2min == config,]$Idx))  
        }
        if(breakFlag == 1){
          break
        }
        
      }  
      
      
      
      
      
    }else if (method == 'allelicRatio'){ 
      for(tumFrac in unique(chr.data$tumFrac)){
        tumFrac.len.perc = sum(chr.data[chr.data$tumFrac == tumFrac,]$segLen)*100/chromLen 
        config = unique(chr.data[chr.data$tumFrac == tumFrac,]$tot2min)
        message(sprintf('\ntumFrac %f for chr %s occupies %f perc, with config %s',tumFrac,chr,tumFrac.len.perc,config))
        
        if(tumFrac.len.perc >= 90){
          message(sprintf('Merging tumFrac %f for chr %s',tumFrac,chr))
          breakFlag=1
          removeFlag=0
          toRemove = c(toRemove,unique(chr.data$Idx))
          if(length(config) > 1){
            message(sprintf('%s - >1 config detected for tumFrac %f',PDID,tumFrac))
            #tumMin = sapply(strsplit(config,split = ':'),'[',2)
            #tumTot = sapply(strsplit(config,split = ':'),'[',1)
            #if(unique(tumMin) == 0){
            #  tmp = data.frame(Idx=0,Chr=chr,Start=1,Stop=chromLen,normTot=2,normMin=1,tumTot = paste(tumTot,collapse=','),tumMin = 0,segLen=chromLen, config=config,idx=0)
            #}else{
            #  stop('WHOOPS')
            #}
          }
          #}else if(length(config) == 1){
          tmp = chr.data[chr.data$tumFrac == tumFrac,][1,]
          tmp$Idx = paste(chr.data[chr.data$tumFrac == tumFrac,]$Idx,collapse = '_')
          tmp$newIdx = paste0(tmp$newIdx,'.', 1)
          tmp$Stop = chromLen
          tmp$Start = 1
          #}
          
        }else if(tumFrac.len.perc <= 10){ # remove if it's <10 % of Chr
          removeFlag=1
          message(sprintf('Removing config %s for chr %s',tumFrac, chr))
          tmp=NULL
          #tmp = data.frame(Idx=0,Chr=chr,Start=1,Stop=chromLen,normTot=2,normMin=1,tumTot = 2,tumMin = 1,segLen=chromLen, config='2:1',idx=0)
          breakFlag=0
        }else{
          # Merge fragments
          a = chr.data[chr.data$tumFrac == tumFrac,]
          rmIdx.a = c()
          rmIdx.a.all = c()
          newSeg.a = tibble()
          
          if(nrow(a) > 1){
            for(f in 2:nrow(a)){
              if((a$Start[f] - a$Stop[f-1]) <= minSegLen){
                message('Merging fragments')
                rmIdx.a = c(rmIdx.a,a$Idx[c(f,(f-1))])
                rmIdx.a.all = c(rmIdx.a.all,a$Idx[c(f,(f-1))])
              }
              if(((a$Start[f] - a$Stop[f-1]) > minSegLen) || f == nrow(a)){
                if(length(rmIdx.a) > 1){
                  a.toRemove = a[a$Idx %in% unique(rmIdx.a),]
                  # Write new merged segment
                  segLen = (max(a.toRemove$Stop) - min(a.toRemove$Start)) + 1
                  patNum = paste(unique(a.toRemove$patNum),collapse = '_')
                  tumTot = paste(unique(a.toRemove$tumTot),collapse = '_')
                  matNum = paste(unique(a.toRemove$matNum),collapse = '_')
                  #patNum = as.numeric(sapply(strsplit(unique(a.toRemove$),split = ':'),'[',2))
                  #tumTot = as.numeric(sapply(strsplit(config,split = ':'),'[',1))
                  #matNum = tumTot-patNum
                  tmp.a = data.frame(Idx=paste(unique(rmIdx.a),collapse = '_'),
                                     Chr=chr,Start=min(a.toRemove$Start),Stop=max(a.toRemove$Stop),
                                     matNum = matNum,tumTot = tumTot,patNum = patNum,segLen=segLen,type='maj',
                                     posID = '0',frac = paste(a.toRemove$frac,collapse = '_'),
                                     tumFrac = tumFrac,tot2min = paste0(tumTot,':',patNum),newIdx = paste0(unique(rmIdx.a)[1],'.1'))
                  newSeg.a=rbind(newSeg.a,tmp.a)
                  rmIdx.a = c()
                  
                }else if(length(rmIdx.a) == 1){
                  message('OOPSS')
                  print(rmIdx.a)
                  stop()
                }
              }
            } 
            if(length(rmIdx.a.all) > 1){
              toRemove = c(toRemove,unique(rmIdx.a.all))
              tmp = newSeg.a    
            }else{
              tmp=NULL
            }
          }else{
            tmp = NULL
          }
          
          
          breakFlag=0
          removeFlag=0
        }
        if(!is.null(tmp)){
          #tmp$patNum = as.numeric(tmp$tumMin)
          #tmp$matNum = as.numeric(tmp$tumTot) - as.numeric(tmp$tumMin)
          #if(sum(is.na(tmp$matNum))>0 & unique(tmp$patNum)==0){
          #  tmp$tumFrac = 1
          #}else{
          #  tmp$tumFrac = tmp$matNum/(tmp$patNum+tmp$matNum)  
          #}
          
          modSegs = rbind(modSegs,tmp)
        }
        if(removeFlag == 1){
          toRemove = c(toRemove,c(chr.data[chr.data$tumFrac == tumFrac,]$Idx))  
        }
        if(breakFlag == 1){
          break
        }
        
      }  
    }
    
  }
  
  
  if((!is.null(toRemove)) & (nrow(modSegs)>0)){
    data = data[!data$Idx %in% toRemove,]
    data = rbind(data, modSegs)
  }
  
  
  data = data[order(data$Chr),]
  
  
  #### Generate final output ####
  new_data = data.frame()
  # Add missing segments
  for(chr in unique(chromInfo$chrom)){
    if(!chr %in% tgtChrs){
      next
    }
    
    chromLen = max(chromInfo[chromInfo$chrom == chr,]$end)*1000
    
    
    if(!chr %in% unique(data$Chr)){
      tmp = data.frame(Idx=0,Chr=chr,Start=1,Stop=chromLen,matNum=1,patNum=1,frac=1,segLen=chromLen,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=paste0('0.',chr))
      tmp$posID = paste0(tmp$Chr,'_',tmp$Start,'_',tmp$Stop)
    }else{
      tmp = data[data$Chr==chr,]
    }
    
    # Add missing segments at the beginning of each Chr
    if(min(tmp$Start) > 1){
      new.tmp = data.frame(Idx=0,Chr=chr,Start=1,Stop=min(tmp$Start)-1,matNum=1,patNum=1,frac=1,segLen=min(tmp$Start)-1,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=paste0('0.',chr))
      new.tmp$posID = paste0(new.tmp$Chr,'_',new.tmp$Start,'_',new.tmp$Stop)
      tmp = rbind(new.tmp,tmp)
    }
    
    # Add missing segments at the end of each Chr
    if(max(tmp$Stop) < chromLen){
      new.tmp = data.frame(Idx=0,Chr=chr,Start=max(tmp$Stop)+1,Stop=chromLen,matNum=1,patNum=1,frac=1,segLen=chromLen,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=paste0('0.1.',chr))
      new.tmp$posID = paste0(new.tmp$Chr,'_',new.tmp$Start,'_',new.tmp$Stop)
      tmp = rbind(new.tmp,tmp)
    }else if(max(tmp$Stop) > chromLen){
      message('Weird STOP...')
      print(chromLen)
      print(tmp[tmp$Stop == max(tmp$Stop),])
      tmp[tmp$Stop == max(tmp$Stop),]$Stop = chromLen
    }
    tmp = tmp[order(tmp$Start),]
    final_tmp=tmp
    for(i in 1:nrow(tmp)){
      if(tmp$Stop[i] == max(tmp$Stop)){
        message(sprintf('%s - Chr %d : DONE!',PDID,chr))
      }else if(tmp$Stop[i] > tmp$Start[i+1]){
        message(sprintf('%s - Chr %d : Overlapping segments: %s',PDID,chr,tmp[i:i+1,]))
      }else if(tmp$Stop[i] < tmp$Start[i+1]-1){
        message(sprintf('%s - Chr %d : Adding segment...',PDID,chr))
        new.tmp = data.frame(Idx=0,Chr=chr,Start=tmp$Stop[i]+1,Stop=tmp$Start[i+1]-1,matNum=1,patNum=1,frac=1,segLen=chromLen,type='maj',tumTot=2,tot2min='2:1',tumFrac=0.5,newIdx=paste0('0.1.',chr))
        new.tmp$posID = paste0(new.tmp$Chr,'_',new.tmp$Start,'_',new.tmp$Stop)
        colnames(new.tmp) = colnames(tmp)
        final_tmp=rbind(new.tmp,final_tmp)
      }else if(tmp$Stop[i] == tmp$Start[i+1]-1){
        message(sprintf('%s - Chr %d : GREAT!',PDID,chr))
      }else{
        message(sprintf('%s: What is happening? chr %s, i = %d',PDID,chr,i))
        print(tmp[i:(i+1),])
        stop()
      }
    }
    
    new_data = rbind(new_data,final_tmp)
  }
  
  
  
  
  ### Add subclone CNAs ###
  subCl.CNA = all[all$type == 'sub',]
  for(chr in unique(subCl.CNA$Chr)){
    maxChromLen = as.vector((chromInfo[(chromInfo$chrom == chr) & (chromInfo$arm == 'q'),]$end)*1000)
    idxs = subCl.CNA[(subCl.CNA$Chr == chr) & (subCl.CNA$Stop > maxChromLen),]$Idx
    if(max(subCl.CNA[subCl.CNA$Chr == chr,]$Stop) > maxChromLen){
      subCl.CNA$Stop[subCl.CNA$Idx %in% idxs] <- rep(maxChromLen,length(idxs))
    }
  }
  new_data = rbind(new_data,subCl.CNA)
  
  ### Process output further, transforming to longer format, add abspos
  if(!'X' %in% tgtChrs){
    new_data$Chr = as.numeric(new_data$Chr)  
  }
  
  
  if(longFormat == T){
    out = pivot_longer(new_data,cols = c('Start','Stop'),names_to = 'posType',values_to = 'pos')
    
    out$idx = c(1:nrow(out))
    # Match max length for each chromosome
    for(chrom in unique(out$Chr)){
      #current_maxChromLen = max(data[data$Chr == chrom,]$pos)
      maxChromLen = as.vector((chromInfo[(chromInfo$chrom == chrom) & (chromInfo$arm == 'q'),]$end)*1000)
      row = out[(out$Chr == chrom) & (out$pos > maxChromLen),]$idx
      if(length(row) > 0){
        print(sprintf('Heey! %d',chrom))
        #out$pos[row] <- rep(maxChromLen,length(row))
      }
    }
    
    # Get absolute genomic position
    out$abspos_kb = out$pos/1000 # if chromosome 1, abspos = pos
    for(r in 1:nrow(out)){
      chrom = out$Chr[r]
      if (chrom > 1){
        out$abspos_kb[r] = out$abspos_kb[r] + (chromInfo[(chromInfo$chrom == (chrom-1)) & (chromInfo$arm == 'q'),]$abspos_kb)
      }
    }
    
  }else{
    out = new_data
  }
  
  # Remove all lines with segLen < 1e6
  #out = out[out$segLen >= 1e6,]
  
  # Remove balanced segment
  if(removeBalancedSegs==T){
    out = out[out$tumFrac != 0.5,]
  }
  
  out$idx = c(1:nrow(out))
  #out = pivot_wider(out,names_from = 'posType',values_from = 'pos')
  if(longFormat == T){
    out=out[order(out$Chr,out$abspos_kb),]  
  }else{
    out=out[order(out$Chr,out$Start),]  
  }
  
  
  # Write BTB output
  #out2 = out[,c("Chr","Start","Stop","matNum","patNum","type" )]
  #out2$subcloneCN = ifelse(out2$type == 'sub',T,F)
  #out2=out2[,-c(6)]
  #out2$PDID = PDID
  #write_delim(out2,paste0('../../BTB_output_used/',PDID,'_btbCN.txt'),delim = '\t',col_names = T)
  
  return(out)
}









processBTB = function(btb.fp,PDID,tgtChrs=c(1:23,'X'),longFormat=T,removeBalancedSegs = FALSE,method = c('totalCN','allelicRatio'),minSegLen=1e6,subCl.minSegLen=1e7,keepClonalSegs=F){
  require(GenomicRanges)
  chromInfo = read.delim('/lustre/scratch117/casm/team274/mt22/chrom_abspos_kb.txt',sep = '\t')
  
  #### Read in CN profile for both major and minor clone from Battenberg ####
  btb = gsub('summary.csv','subclones.txt.gz',btb.fp)
  btb = read.delim(btb,sep = '\t')
  btb$posID = paste0(btb$chr,':',btb$startpos,'_',btb$endpos)
  btb$idx = c(1:nrow(btb))
  
  ### Remove clonality if subclone fraction is only <= 5%
  w = which(!is.na(btb$frac2_A) & btb$frac2_A < 0.05)
  btb[w,c("nMaj2_A","nMin2_A","frac2_A")] = NA
  if(sum(btb$frac1_A < 0.05) > 0){
    stop('Weird!')
  }
  
  # Keep CN segments on chr of interest only
  btb = btb[btb$chr %in% tgtChrs,]
  btb$chr = ifelse(btb$chr == 'X',23,as.numeric(btb$chr))
  #btb$segLen = btb$endpos - btb$startpos + 1
  
  btb = GRanges(btb$chr,IRanges(btb$startpos,btb$endpos),
                chr=btb$chr,idx = btb$idx,posID=btb$posID,
                nMaj1_A=btb$nMaj1_A,nMin1_A=btb$nMin1_A, frac1_A=btb$frac1_A,
                nMaj2_A=btb$nMaj2_A,nMin2_A=btb$nMin2_A, frac2_A=btb$frac2_A)
  
  finalSegs_noClonalSegs = GRanges()
  finalClonalSegs = GRanges()
  
  for(chr in tgtChrs){
    #message(sprintf('Chromosome %s',chr))
    chr_segs = GRanges()
    data = btb[seqnames(btb) == chr]
    
    # Fill in the gaps - they should all be diploids
    gap_segs = GRanges(gaps(data),chr=chr,patNum=1, matNum=1,totCN=2,tumFrac = 0.5)
    # Fill in the end of the chromosome too (if needed)
    chromLen = max(chromInfo[chromInfo$chrom == chr,]$end)*1000
    if(max(end(data)) < chromLen){
      tmp = GRanges(chr,IRanges(max(end(data))+1,chromLen),chr=chr,patNum=1, matNum=1,totCN=2,tumFrac = 0.5)
      gap_segs = append(gap_segs,tmp)
    }else if (max(end(data)) > chromLen){
      if(start(data[end(data) == max(end(data))]) < chromLen){
        end(data[end(data) == max(end(data))]) = chromLen  
      }else{
        #remove that seg
        data = data[end(data) != max(end(data))]
      }
      
    }
    # If max(end(data)) is longer than chromLen (just leave it for now...)
    chr_segs = append(chr_segs,gap_segs)
    
    
    # Extract all "noClonal" segments
    no_clonalSegs = data[which(is.na(data$frac2_A))]
    if(length(no_clonalSegs)>0){
      mcols(no_clonalSegs) = data.frame(chr=chr,idx = no_clonalSegs$idx, patNum=no_clonalSegs$nMin1_A, matNum=no_clonalSegs$nMaj1_A,
                                        totCN=no_clonalSegs$nMin1_A+no_clonalSegs$nMaj1_A,tumFrac = no_clonalSegs$nMaj1_A/(no_clonalSegs$nMin1_A+no_clonalSegs$nMaj1_A))  
    }
    
    # Extract all "clonal" segments
    clonal = data[which((!is.na(data$frac2_A)))] 
    #clonal = btb[which((!is.na(btb$frac2_A)) & (pmin(btb$frac1_A,btb$frac2_A) >= 0.1)),c('Idx','chr','startpos','endpos','nMaj1_A','nMin1_A','frac1_A','nMaj2_A','nMin2_A','frac2_A','posID','segLen')] 
    if(length(clonal) > 0){
      message(sprintf('%s - chr %s: clonal segments found!',PDID,chr))
      clonal$segType = 'clonal'  
      
      ### Sort out overlapping segments
      # There should be no overlaps between no_change and altered / clonal.
      # If there are overlaps between altered and clonal - remove the overlapping segments from altered
      o = findOverlapPairs(no_clonalSegs,clonal)
      # Remove overlapping segments from "no_clonalSegs" group (to be replaced by disjoin segments only)
      no_clonalSegs = no_clonalSegs[!no_clonalSegs$idx %in% o@first$idx]
      # If there are any overlaping segments between "no_clonalSegs" and "clonal"
      if(length(o)>0){
        message(sprintf('%s - chr %s: OVERLAPPING clonal segments found!',PDID,chr))
        no_clonalSegs_o = append(o@first,o@second)
        disjoin_ranges = disjoin(no_clonalSegs_o)
        disjoin_ranges = disjoin_ranges[!overlapsAny(disjoin_ranges,clonal)]
        disjoin_ranges_no_clonalSegs = findOverlapPairs(disjoin_ranges,o@first)
        # Check that disjoin_ranges_no_clonalSegs@second should always be of type "no_clonalSegs"
        if(sum(disjoin_ranges_no_clonalSegs@second$type != 'no_clonalSegs') > 0){
          stop(sprintf('%s - chr %s: Something strange in disjoin happend...',PDID,chr))
        }
        # Each new no_clonalSegs segment should match to a unique segment from o3@first
        for(i in 1:length(disjoin_ranges_no_clonalSegs)){
          pair = disjoin_ranges_no_clonalSegs[i]
          start(pair@second) = start(pair@first)
          end(pair@second) = end(pair@first)
          no_clonalSegs = append(no_clonalSegs,pair@second)
        }
      }
      
      # Restructure the clonal segments
      
      for(i in 1:length(clonal)){
        seg = clonal[i]
        # Define maj vs sub clones
        if(seg$frac1_A > seg$frac2_A){
          maj_seg = GRanges(seg$chr,IRanges(start(seg),end(seg)),chr=seg$chr,clonalType='maj',
                            patNum=seg$nMin1_A, matNum=seg$nMaj1_A,totCN=seg$nMin1_A+seg$nMaj1_A,tumFrac = seg$nMaj1_A/(seg$nMaj1_A+seg$nMin1_A))
          sub_seg = GRanges(seg$chr,IRanges(start(seg),end(seg)),chr=seg$chr,clonalType='sub',
                            patNum=seg$nMin2_A, matNum=seg$nMaj2_A,totCN=seg$nMin2_A+seg$nMaj2_A,tumFrac = seg$nMaj2_A/(seg$nMaj2_A+seg$nMin2_A))
        }else{
          sub_seg = GRanges(seg$chr,IRanges(start(seg),end(seg)),chr=seg$chr,clonalType='sub',
                            patNum=seg$nMin1_A, matNum=seg$nMaj1_A,totCN=seg$nMin1_A+seg$nMaj1_A,tumFrac = seg$nMaj1_A/(seg$nMaj1_A+seg$nMin1_A))
          maj_seg = GRanges(seg$chr,IRanges(start(seg),end(seg)),chr=seg$chr,clonalType='maj',
                            patNum=seg$nMin2_A, matNum=seg$nMaj2_A,totCN=seg$nMin2_A+seg$nMaj2_A,tumFrac = seg$nMaj2_A/(seg$nMaj2_A+seg$nMin2_A))
        }
        finalClonalSegs = append(finalClonalSegs,sub_seg)
        finalClonalSegs = append(finalClonalSegs,maj_seg)
      }
      
    }
    
    
    chr_segs = append(chr_segs,no_clonalSegs)
    
    ### Merging segments of the same type and matNum:patNum configuration
    final_chr_segs = GRanges()
    chr_segs$tot2min = paste0(chr_segs$totCN,":",chr_segs$patNum)
    
    if(method == 'allelicRatio'){
      # Special scenarios for "PD46693" chr10 and 11
      if((PDID == "PD46693" & chr %in% c(10,11)) |(PDID == "PD43255" & chr %in% c(18)) ){
        message(sprintf('%s chr%s: Replacing cn segments',PDID, chr))
        # Add artificial piece in
        final_chr_segs = GRanges(chr,IRanges(1,chromLen),chr=chr,tot2min='x:0',tumFrac = 1)
      }
      for(tumFrac in unique(chr_segs$tumFrac)){
        ### - MAYBE: If >= 95% of the chromosome has the same tumFrac (config) --> make the entire chromosome having the same config
        
        #print(config)
        chr_segs.sub = chr_segs[chr_segs$tumFrac == tumFrac]
        chr_segs.sub_ranges = GenomicRanges::reduce(chr_segs.sub)
        mcols(chr_segs.sub_ranges) = data.frame(chr=chr,tot2min=paste(unique(chr_segs.sub$tot2min),collapse = '_'),tumFrac = unique(chr_segs.sub$tumFrac))
        final_chr_segs = append(final_chr_segs,chr_segs.sub_ranges)
      }
    }else if(method == 'totalCN'){
      # Special scenarios for "PD46693" chr10 and 11
      #if((PDID == "PD46693" & chr %in% c(10,11)) |(PDID == "PD43255" & chr %in% c(18)) ){
      #  message(sprintf('%s chr%s: Replacing cn segments',PDID, chr))
      #  # Add artificial piece in
      #  final_chr_segs = GRanges(chr,IRanges(1,chromLen),chr=chr,tot2min='x:0',tumFrac = 1)
      #}
      for(totCN in unique(chr_segs$totCN)){
        ### - MAYBE: If >= 95% of the chromosome has the same tumFrac (config) --> make the entire chromosome having the same config
        
        #print(config)
        chr_segs.sub = chr_segs[chr_segs$totCN == totCN]
        chr_segs.sub_ranges = GenomicRanges::reduce(chr_segs.sub)
        mcols(chr_segs.sub_ranges) = data.frame(chr=chr,totCN=unique(chr_segs.sub$totCN),tumFrac = paste(unique(chr_segs.sub$tumFrac),collapse = '_'))
        final_chr_segs = append(final_chr_segs,chr_segs.sub_ranges)
      }
    }
      
    
    
    
    
    
    #final_chr_segs$tot2min = paste0(final_chr_segs$totCN,":",final_chr_segs$patNum)
    #final_chr_segs$posID = paste0(seqnames(final_chr_segs),':',start(final_chr_segs),'_',end(final_chr_segs))
    
    finalSegs_noClonalSegs = append(finalSegs_noClonalSegs,final_chr_segs)
  }
  
  
  
  ### Add subclone CNAs ###
  if(keepClonalSegs == T){
    finalSegs_noClonalSegs$clonalType = 'none'
    finalSegs = append(finalSegs_noClonalSegs,finalClonalSegs)
  }else{
    finalSegs = finalSegs_noClonalSegs
  }
  finalSegs$Start = start(finalSegs)
  finalSegs$Stop = end(finalSegs)
  #finalSegs$tot2min = paste0(finalSegs$totCN,":",finalSegs$patNum)
  finalSegs$posID = paste0(seqnames(finalSegs),':',start(finalSegs),'_',end(finalSegs))
  finalSegs = finalSegs[width(finalSegs) > minSegLen]
  
  
  
  
  if(longFormat == T){
    data = as.data.frame(mcols(finalSegs))
    out = pivot_longer(data,cols = c('Start','Stop'),names_to = 'posType',values_to = 'pos')
    out$idx = c(1:nrow(out))
    # Match max length for each chromosome
    for(chrom in unique(out$chr)){
      #current_maxChromLen = max(data[data$Chr == chrom,]$pos)
      maxChromLen = as.vector((chromInfo[(chromInfo$chrom == chrom) & (chromInfo$arm == 'q'),]$end)*1000)
      row = out[(out$chr == chrom) & (out$pos > maxChromLen),]$idx
      if(length(row) > 0){
        print(sprintf('Heey! %d %d',chrom,length(row)))
      }
    }
    
    # Get absolute genomic position
    out$abspos_kb = out$pos/1000 # if chromosome 1, abspos = pos
    for(r in 1:nrow(out)){
      chrom = out$chr[r]
      if (chrom > 1){
        out$abspos_kb[r] = out$abspos_kb[r] + (chromInfo[(chromInfo$chrom == (chrom-1)) & (chromInfo$arm == 'q'),]$abspos_kb)
      }
    }
    out$idx = c(1:nrow(out))
    out=out[order(out$chr,out$abspos_kb),]  
    if(removeBalancedSegs==T){
      out = out[out$tumFrac != 0.5,]
    }
  }else{
    out = finalSegs
    # Remove balanced segment
    if(removeBalancedSegs==T){
      out = out[out$tumFrac != 0.5]
    }
    out=out[order(out$chr,out$Start)]  
    out$idx = c(1:length(out))
    names(out) = seqnames(out)
  }
  
  # Remove all lines with segLen < 1e6
  #out = out[out$segLen >= 1e6,]
  return(out)
}






#' Save raster and non-raster versions of plots
#'
#' Saves a pdf and png version of everything.
#'
#' @param baseName Name of file to save.  File name without extension.
#' @param plotFun Function to generate the plot.  Must take two arguments, noFrame and noPlot which control if the plot part / frame part are drawn.
#' @param width Width in inches.
#' @param heights Height in inches.
#' @param res Resolution in pixels per inch for rasterised image.
#' @param rawData If provided, save the raw data used to make the plot.
#' @param row.names Passed to write.table
#' @param col.names Passed to write.table
#' @param ... Passed to plotFun
saveFig = function(baseName,plotFun,width=4,height=3,res=300,rawData=NULL,row.names=FALSE,col.names=TRUE,...){
  #Save the base versions
  pdf(paste0(baseName,'.pdf'),width=width,height=height,useDingbats=FALSE)
  plotFun(noFrame=FALSE,noPlot=FALSE,...)
  dev.off()
  png(paste0(baseName,'.png'),width=width,height=height,res=300,units='in')
  plotFun(...)
  dev.off()
  #And the various bits of frame in vector format
  pdf(paste0(baseName,'_frame.pdf'),width=width,height=height,useDingbats=FALSE)
  plotFun(noPlot=TRUE,noFrame=FALSE,...)
  dev.off()
  #And the various bits of the actual plot in raster format
  png(paste0(baseName,'_plot.png'),width=width,height=height,res=300,units='in')
  plotFun(noPlot=FALSE,noFrame=TRUE,...)
  dev.off()
  #Save raw data if provided
  if(!is.null(rawData)){
    write.table(rawData,paste0(baseName,'_rawData.tsv'),
                quote=FALSE,
                sep='\t',
                row.names=row.names,
                col.names=col.names)
  }
}



#' Generate coverage (logR) and BAF from BAM file
#'
#' The most "raw" way to look at copy number status of a BAM file is to calculate the coverage and alternate allele frequency at heterozygous SNPs.  This function uses alleleCounter to generate this from a BAM file, for a set of heterozygous SNPs.  You can specify the SNPs however you like, but usually this would either be using a pre-specified panel of SNPs, or using the \code{\link{findHetSNPs}} function to identify potential heterozygous SNPs directly.
#'
#' Will optionally run ASCAT to infer CN segments (in unmatched mode) and produces a basic plot to visualise the raw data.  If ASCAT is run, there is some extra information that is needed.  These values almost all relate to guessing if a given SNP is heterozygous in the germline.  The allele frequency \code{AF} is used as a prior and posterior probabilities of heterozygosity and homozygosity (in the germline) is determined using the provided error rate (for homozygosity) and maximal abbarent cell fraction (\code{abbCellFrac} for heterozygosity).  Any SNP with posterior probability of homozygosity greater than \code{minHomProb} or observed BAF less than \code{minHetAF} will be coded as homozygous.  This errs on the side of marking things as homozygous if the data is not clear, which is the correct thing to do.
#'
#' @inheritParams phaseSNPsFromCN
#' @param BAM The BAM file to calculate coverage and BAF on.
#' @param hSNPs The heterozygous SNPs to investigate.  Defaults to 1000 genomes SNPs.  Must be a GRanges object with REF and ALT columns.
#' @param AF The population allele frequency of the alternate allele.  If NULL, defaults to 0.05
#' @param gcFile ASCAT GC correction file for the SNPs provided.  If given, used to GC correct logR.
#' @param chrs Which chromosomes to plot and/or run ASCAT on.
#' @param errRate Assumed sequencing error rate.
#' @param abbCellFrac Assumed fraction of abberant cells for determining heterozygous SNPs.  Higher values are more stringent.
#' @param minHomProb Any SNP with posterior probability of homozygosity in the germline greater than this will be code as homozygous.
#' @param minHetAF Any SNP with BAF below this will be coded as homozygous.
#' @param runASCAT Should we run ASCAT?
#' @return The hSNPs object with extra columns, particularly total coverage and BAF.  If ASCAT has been run, columns associated with the fit are also added and the full ASCAT output is stored at \code{hSNPs@metadata$ascat.output}.
#' @importFrom stats median
#' @export
generateCoverageAndBAF = function(BAM,refGenome,hSNPs=SNPs1k,AF=hSNPs$AF,outPath=NULL,gcFile=NULL,
                                  alleleCounterParams=list(f=3,F=3852,m=20,q=35),chrs=c(1:22,'X'),errRate=0.01,abbCellFrac=0.95,minHomProb=0.05,minHetAF=0.1,runASCAT=FALSE,nParallel=1){
  ###############################
  # Preliminary checks and setup
  if(is.null(AF))
    AF=0.05
  hSNPs$AF = AF
  ##############################
  # Run alleleCounter as needed
  #Set the non-optional allelecounter parameters
  alleleCounterParams$x=FALSE
  alleleCounterParams$bams = BAM
  alleleCounterParams$tgtLoci = hSNPs
  alleleCounterParams$refGenome = refGenome
  alleleCounterParams$outputs = outPath
  alleleCounterParams$nParallel = nParallel
  rawCnts = do.call(alleleCounter,alleleCounterParams)[[1]]
  ###################################
  # Extract the values we care about
  rawCnts$refCnts = as.matrix(mcols(rawCnts)[,c('A','C','G','T')])[cbind(seq_along(rawCnts$REF),match(rawCnts$REF,c('A','C','G','T')))]
  rawCnts$altCnts = as.matrix(mcols(rawCnts)[,c('A','C','G','T')])[cbind(seq_along(rawCnts$ALT),match(rawCnts$ALT,c('A','C','G','T')))]
  rawCnts$sourceBAM = BAM
  rawCnts$coverage = rawCnts$refCnts + rawCnts$altCnts
  rawCnts$BAF = rawCnts$altCnts/rawCnts$coverage
  chrs = intersect(chrs,unique(as.character(seqnames(rawCnts[rawCnts$coverage>0]))))
  #Determine which ones are likely heterozygous in the germline
  #The model is that we assume an error rate, that by default is conservatively chosen to be a high value of 0.01 (see this paper https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1659-6).  Then we need to specify what the highest purity we expect is.  Obviously if it's pure tumour, then LOH will result in complete loss of one allele and we cannot tell which are heterozygous and which are homozygous.  The default value is 95% purity and the model assumes a loss of one allele in the aberant cell (i.e. haploid).
  #We also incoporate the allele frequency, if known, as a prior.  Otherwise a uniform prior of 5% allele frequency is used.
  #Convert everything to the less than 0.5 case
  pp = ifelse(rawCnts$AF>0.5,1-rawCnts$AF,rawCnts$AF)
  x = ifelse(rawCnts$altCnts>rawCnts$altCnts/rawCnts$coverage,rawCnts$coverage-rawCnts$altCnts,rawCnts$altCnts)
  pHom = (1-pp*(1-pp))*dbinom(x,rawCnts$coverage,errRate)
  pHet = pp*(1-pp)*dbinom(x,rawCnts$coverage,(1-abbCellFrac)/(2-abbCellFrac))
  #Workaround for numerical issues
  pHom[is.na(pHom)] = 0
  pHet[is.na(pHet)] = 0
  normFac = pHet + pHom
  pHet = pHet/normFac
  pHom = pHom/normFac
  rawCnts$isHet = !(pHom>minHomProb | abs(rawCnts$BAF-0.5)>(0.5-minHetAF))
  #Now calculate logR.  Because the definition is normalised to the median, dividing by any constant here will change nothing.
  rawCnts$logR = log(rawCnts$coverage) - median(log(rawCnts$coverage))
  #Run ascat if required
  if(runASCAT){
    message("Fitting CN states using unmatched ASCAT")
    #Make temporary object to manipulate, dropping zeros
    tmp = rawCnts[rawCnts$coverage>0,]
    nom = paste0('BAM',1)
    #Drop anything that is not a target chromosome
    chrsLocal = intersect(chrs,unique(as.character(seqnames(tmp))))
    tmp = tmp[as.character(seqnames(tmp)) %in% chrsLocal,]
    #Reorder as required
    tmp = tmp[order(factor(as.character(seqnames(tmp)),levels=chrsLocal),start(tmp)),]
    #Make position object
    SNPpos = data.frame(row.names=as.character(tmp),chrs=as.character(seqnames(tmp)),pos=start(tmp))
    ch = split(seq_along(tmp),as.character(seqnames(tmp)))[chrs]
    names(ch) = NULL
    #Construct in ASCAT format
    inDat = list(Tumor_LogR = matrix(tmp$logR,ncol=1,dimnames=list(as.character(tmp),nom)),
                 Tumor_BAF = matrix(tmp$BAF,ncol=1,dimnames=list(as.character(tmp),nom)), 
                 Tumor_LogR_segmented = NULL, 
                 Tumor_BAF_segmented = NULL, 
                 Germline_LogR = NULL, 
                 Germline_BAF = NULL, 
                 SNPpos = SNPpos,
                 ch = ch,
                 chr = ASCAT:::split_genome(SNPpos),
                 chrs = chrsLocal,
                 samples = nom,
                 gender = 'XX',
                 sexchromosomes = 'XX',
                 failedarrays = NULL
    )
    #Do we have a GC content file.  If so, use it to correct
    if(!is.null(gcFile))
      inDat = ASCAT::ascat.GCcorrect(inDat,gcFile)
    #ascat.plotRawData(inDat,'~')
    isHom = list(germlinegenotypes = matrix(!tmp$isHet,ncol=1,dimnames=dimnames(inDat$Tumor_BAF)),failedarrays=NULL)
    inDat = ASCAT::ascat.aspcf(inDat, ascat.gg=isHom)
    #ascat.plotSegmentedData(inDat,'~')
    ascat.output = ASCAT::ascat.runAscat(inDat,gamma=1)
    #ascat.plotSunrise(ascat.output$distance_matrix[[1]],ascat.output$ploidy,ascat.output$aberrantcellfraction)
    #Include 
    rawCnts@metadata$ascat.output = ascat.output
    rawCnts$nA = NA
    rawCnts$nB = NA
    m = match(rownames(ascat.output$nA),as.character(rawCnts))
    rawCnts$nA[m] = ascat.output$nA[,1]
    rawCnts$nB[m] = ascat.output$nB[,1]
    if(!is.null(gcFile)){
      rawCnts$correctedLogR = NA
      rawCnts$correctedLogR[m] = inDat$Tumor_LogR[,1]
    }
  }
  return(rawCnts)
}







CNA.MCMC_inhouse = function (clu, fttmat, bins, cut.cor, n.cores, seed) {
  CON <- NULL
  for (i in min(clu):max(clu)) {
    data.c <- apply(fttmat[, which(clu == i)], 1, median)
    CON <- cbind(CON, data.c)
    i <- i + 1
  }
  norm.mat.sm <- exp(CON)
  n <- nrow(norm.mat.sm)
  BR <- NULL
  for (c in 1:ncol(norm.mat.sm)) {
    breks <- c(seq(1, as.integer(n/bins - 1) * bins, bins), n)
    bre <- NULL
    for (i in 1:(length(breks) - 2)) {
      a1 <- max(mean(norm.mat.sm[breks[i]:breks[i + 1], 
                                 c]), 0.001)
      set.seed(seed)
      posterior1 <- MCMCpack::MCpoissongamma(norm.mat.sm[breks[i]:breks[i + 
                                                                          1], c], a1, 1, mc = 1000)
      a2 <- max(mean(norm.mat.sm[(breks[i + 1] + 1):breks[i + 
                                                            2], c]), 0.001)
      set.seed(seed)
      posterior2 <- MCMCpack::MCpoissongamma(norm.mat.sm[(breks[i + 1] + 1):breks[i + 2], c], a2, 1, mc = 1000)
      if (ks.test(posterior1, posterior2)$statistic[[1]] > cut.cor) {
        bre <- c(bre, breks[i + 1])
      }
      i <- i + 1
    }
    breks <- sort(unique(c(1, bre, n)))
    BR <- sort(unique(c(BR, breks)))
    c <- c + 1
  }
  norm.mat.sm <- exp(fttmat)
  seg <- function(z) {
    x <- numeric(n)
    for (i in 1:(length(BR) - 1)) {
      a <- max(mean(norm.mat.sm[BR[i]:BR[i + 1], z]), 0.001)
      set.seed(seed)
      posterior1 <- MCMCpack::MCpoissongamma(norm.mat.sm[BR[i]:BR[i + 1], z], a, 1, mc = 1000)
      x[BR[i]:BR[i + 1]] <- mean(posterior1)
      i <- i + 1
    }
    x <- log(x)
  }
  seg.test <- parallel::mclapply(1:ncol(norm.mat.sm), seg, mc.cores = n.cores)
  logCNA <- matrix(unlist(seg.test), ncol = ncol(norm.mat.sm), byrow = FALSE)
  res <- list(logCNA, BR)
  names(res) <- c("logCNA", "breaks")
  return(res)
}








copykat_inhouse = function (rawmat = rawdata, id.type = "S", cell.line = "no", 
                       ngene.chr = 5, LOW.DR = 0.05, UP.DR = 0.1, win.size = 25, 
                       norm.cell.names = "", KS.cut = 0.1, sam.name = "", distance = "euclidean", 
                       output.seg = "FALSE", plot.genes = "TRUE", genome = "hg20", 
                       n.cores = 1, seed=2397) 
{
  start_time <- Sys.time()
  set.seed(seed)
  sample.name <- paste(sam.name, "_copykat_", sep = "")
  print("running copykat v1.0.8 updated 02/25/2022 introduced mm10 module, fixed typos")
  print("step1: read and filter data ...")
  print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data", sep = ""))
  genes.raw <- apply(rawmat, 2, function(x) (sum(x > 0)))
  if (sum(genes.raw > 200) == 0) 
    stop("none cells have more than 200 genes")
  if (sum(genes.raw < 100) > 1) {
    rawmat <- rawmat[, -which(genes.raw < 200)]
    print(paste("filtered out ", sum(genes.raw <= 200), 
                " cells with less than 200 genes; remaining ", ncol(rawmat), 
                " cells", sep = ""))
  }
  der <- apply(rawmat, 1, function(x) (sum(x > 0)))/ncol(rawmat)
  if (sum(der > LOW.DR) >= 1) {
    rawmat <- rawmat[which(der > LOW.DR), ]
    print(paste(nrow(rawmat), " genes past LOW.DR filtering", sep = ""))
  }
  WNS1 <- "data quality is ok"
  if (nrow(rawmat) < 7000) {
    WNS1 <- "low data quality"
    UP.DR <- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }
  print("step 2: annotations gene coordinates ...")
  if (genome == "hg20") {
    anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type)
  }else if (genome == "mm10") {
    anno.mat <- annotateGenes.mm10(mat = rawmat, ID.type = id.type)
    dim(rawmat)
  }
  anno.mat <- anno.mat[order(as.numeric(anno.mat$abspos), decreasing = FALSE), ]
  
  if (genome == "hg20") {
    HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
    toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]), HLAs))
    if (length(toRev) > 0) {
      anno.mat <- anno.mat[-toRev, ]
    }
  }                     
  
  
  ToRemov2 <- NULL
  for (i in 8:ncol(anno.mat)) {
    cell <- cbind(anno.mat$chromosome_name, anno.mat[, i])
    cell <- cell[cell[, 2] != 0, ]
    if (length(as.numeric(cell)) < 5) {
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    else if (length(rle(cell[, 1])$length) < length(unique((anno.mat$chromosome_name))) | 
             min(rle(cell[, 1])$length) < ngene.chr) {
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    i <- i + 1
  }
  
  
  if (length(ToRemov2) == (ncol(anno.mat) - 7)) 
    stop("all cells are filtered")
  if (length(ToRemov2) > 0) {
    anno.mat <- anno.mat[, -which(colnames(anno.mat) %in% 
                                    ToRemov2)]
  }
  rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
  norm.mat <- log(sqrt(rawmat3) + sqrt(rawmat3 + 1))
  norm.mat <- apply(norm.mat, 2, function(x) (x <- x - mean(x)))
  colnames(norm.mat) <- colnames(rawmat3)
  print("step 3: smoothing data with dlm ...")
  dlm.sm <- function(c) {
    model <- dlm::dlmModPoly(order = 1, dV = 0.16, dW = 0.001)
    x <- dlm::dlmSmooth(norm.mat[, c], model)$s
    x <- x[2:length(x)]
    x <- x - mean(x)
  }
  test.mc <- parallel::mclapply(1:ncol(norm.mat), dlm.sm, 
                                mc.cores = n.cores)
  norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), 
                            byrow = FALSE)
  colnames(norm.mat.smooth) <- colnames(norm.mat)
  print("step 4: measuring baselines ...")
  if (cell.line == "yes") {
    print("running pure cell line mode")
    relt <- baseline.synthetic(norm.mat = norm.mat.smooth, 
                               min.cells = 10, n.cores = n.cores)
    norm.mat.relat <- relt$expr.relat
    CL <- relt$cl
    WNS <- "run with cell line mode"
    preN <- NULL
  }
  else if (length(norm.cell.names) > 1) {
    NNN <- length(colnames(norm.mat.smooth)[which(colnames(norm.mat.smooth) %in% 
                                                    norm.cell.names)])
    print(paste(NNN, " known normal cells found in dataset", 
                sep = ""))
    if (NNN == 0) 
      stop("known normal cells provided; however none existing in testing dataset")
    print("run with known normal...")
    basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% 
                                             norm.cell.names)], 1, median)
    print("baseline is from known input")
    d <- parallelDist::parDist(t(norm.mat.smooth), threads = n.cores, 
                               method = "euclidean")
    km <- 6
    fit <- hclust(d, method = "ward.D2")
    CL <- cutree(fit, km)
    while (!all(table(CL) > 5)) {
      km <- km - 1
      CL <- cutree(fit, k = km)
      if (km == 2) {
        break
      }
    }
    WNS <- "run with known normal"
    preN <- norm.cell.names
    norm.mat.relat <- norm.mat.smooth - basel
  }
  else {
    basa <- baseline.norm.cl(norm.mat.smooth = norm.mat.smooth, 
                             min.cells = 5, n.cores = n.cores)
    basel <- basa$basel
    WNS <- basa$WNS
    preN <- basa$preN
    CL <- basa$cl
    if (WNS == "unclassified.prediction") {
      basa <- baseline.GMM(CNA.mat = norm.mat.smooth, 
                           max.normal = 5, mu.cut = 0.05, Nfraq.cut = 0.99, 
                           RE.before = basa, n.cores = n.cores)
      basel <- basa$basel
      WNS <- basa$WNS
      preN <- basa$preN
    }
    norm.mat.relat <- norm.mat.smooth - basel
  }
  DR2 <- apply(rawmat3, 1, function(x) (sum(x > 0)))/ncol(rawmat3)
  norm.mat.relat <- norm.mat.relat[which(DR2 >= UP.DR), ]
  anno.mat2 <- anno.mat[which(DR2 >= UP.DR), ]
  ToRemov3 <- NULL
  for (i in 8:ncol(anno.mat2)) {
    cell <- cbind(anno.mat2$chromosome_name, anno.mat2[, i])
    cell <- cell[cell[, 2] != 0, ]
    if (length(as.numeric(cell)) < 5) {
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }
    else if (length(rle(cell[, 1])$length) < length(unique((anno.mat$chromosome_name))) | 
             min(rle(cell[, 1])$length) < ngene.chr) {
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }
    i <- i + 1
  }
  if (length(ToRemov3) == ncol(norm.mat.relat)) 
    stop("all cells are filtered")
  if (length(ToRemov3) > 0) {
    norm.mat.relat <- norm.mat.relat[, -which(colnames(norm.mat.relat) %in% ToRemov3)]
  }
  CL <- CL[which(names(CL) %in% colnames(norm.mat.relat))]
  CL <- CL[order(match(names(CL), colnames(norm.mat.relat)))]
  
  print("step 5: segmentation...")
  results <- CNA.MCMC_inhouse(clu = CL, fttmat = norm.mat.relat, bins = win.size, 
                         cut.cor = KS.cut, n.cores = n.cores, seed = seed)
  if (length(results$breaks) < 25) {
    print("too few breakpoints detected; decreased KS.cut to 50%")
    results <- CNA.MCMC_inhouse(clu = CL, fttmat = norm.mat.relat, 
                           bins = win.size, cut.cor = 0.5 * KS.cut, n.cores = n.cores, seed = seed)
  }
  if (length(results$breaks) < 25) {
    print("too few breakpoints detected; decreased KS.cut to 75%")
    results <- CNA.MCMC_inhouse(clu = CL, fttmat = norm.mat.relat, 
                           bins = win.size, cut.cor = 0.5 * 0.5 * KS.cut, n.cores = n.cores, seed = seed)
  }
  if (length(results$breaks) < 25) 
    stop("too few segments; try to decrease KS.cut; or improve data")
  colnames(results$logCNA) <- colnames(norm.mat.relat)
  results.com <- apply(results$logCNA, 2, function(x) (x <- x - mean(x)))
  RNA.copycat <- cbind(anno.mat2[, 1:7], results.com)
  write.table(RNA.copycat, paste(sample.name, "CNA_raw_results_gene_by_cell.txt", sep = ""), sep = "\t", row.names = FALSE, quote = F)
  if (genome == "hg20") {
    print("step 6: convert to genomic bins...")
    Aj <- convert.all.bins.hg20(DNA.mat = DNA.hg20, RNA.mat = RNA.copycat, 
                                n.cores = n.cores)
    uber.mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
    
    print("step 7: adjust baseline ...")
    if (cell.line == "yes") {
      mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
      write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.txt", sep = ""), sep = "\t", row.names = FALSE, quote = F)
      if (distance == "euclidean") {
        hcc <- hclust(parallelDist::parDist(t(mat.adj), threads = n.cores, method = distance), method = "ward.D")
      } else {
        hcc <- hclust(as.dist(1 - cor(mat.adj, method = distance)), method = "ward.D")
      }
      saveRDS(hcc, file = paste(sample.name, "clustering_results.rds", sep = ""))
      
      print("step 8: ploting heatmap ...")
      my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
      chr <- as.numeric(Aj$DNA.adj$chrom)%%2 + 1
      rbPal1 <- colorRampPalette(c("black", "grey"))
      CHR <- rbPal1(2)[as.numeric(chr)]
      chr1 <- cbind(CHR, CHR)
      if (ncol(mat.adj) < 3000) {
        h <- 10
      } else {
        h <- 15
      }
      col_breaks = c(seq(-1, -0.4, length = 50), seq(-0.4, -0.2, length = 150), seq(-0.2, 0.2, length = 600), 
                     seq(0.2, 0.4, length = 150), seq(0.4, 1, length = 50))
      if (distance == "euclidean") {
        jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), 
             height = h * 250, width = 4000, res = 100)
        heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) parallelDist::parDist(x, 
                                                                                            threads = n.cores, method = distance), 
                  hclustfun = function(x) hclust(x, 
                  method = "ward.D"), ColSideColors = chr1, 
                  Colv = NA, Rowv = TRUE, notecol = "black", 
                  col = my_palette, breaks = col_breaks, key = TRUE, 
                  keysize = 1, density.info = "none", trace = "none", 
                  cexRow = 0.1, cexCol = 0.1, cex.main = 1, 
                  cex.lab = 0.1, symm = F, symkey = F, symbreaks = T, 
                  cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), 
                  cex.main = 4, margins = c(10, 10))
        dev.off()
        if (plot.genes == "TRUE") {
          dim(results.com)
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom)%%2 + 1
          rbPal1g <- colorRampPalette(c("black", "grey"))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg, CHRg)
          pdf(paste(sample.name, "with_genes_heatmap.pdf", sep = ""), height = h * 2.5, width = 40)
          p=heatmap.3(t(results.com), dendrogram = "r", 
                    distfun = function(x) parallelDist::parDist(x, threads = n.cores, method = distance), 
                    hclustfun = function(x) hclust(x, method = "ward.D"), 
                    ColSideColors = chr1g, Colv = NA, Rowv = TRUE, 
                    notecol = "black", col = my_palette, breaks = col_breaks, 
                    key = TRUE, keysize = 1, density.info = "none", 
                    trace = "none", cexRow = 0.1, cexCol = 0.1, 
                    cex.main = 1, cex.lab = 0.1, symm = F, symkey = F, 
                    symbreaks = T, cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 10))
          print(p)
          dev.off()
        }
      }
      else {
        jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), height = h * 250, width = 4000, res = 100)
        heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) as.dist(1 - cor(t(x), method = distance)), 
                  hclustfun = function(x) hclust(x, method = "ward.D"), ColSideColors = chr1, 
                  Colv = NA, Rowv = TRUE, notecol = "black", 
                  col = my_palette, breaks = col_breaks, key = TRUE, 
                  keysize = 1, density.info = "none", trace = "none", 
                  cexRow = 0.1, cexCol = 0.1, cex.main = 1, 
                  cex.lab = 0.1, symm = F, symkey = F, symbreaks = T, 
                  cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), 
                  cex.main = 4, margins = c(10, 10))
        dev.off()
        if (plot.genes == "TRUE") {
          dim(results.com)
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom)%%2 + 1
          rbPal1g <- colorRampPalette(c("black", "grey"))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg, CHRg)
          pdf(paste(sample.name, "with_genes_heatmap.pdf", sep = ""), height = h * 2.5, width = 40)
          p=heatmap.3(t(results.com), dendrogram = "r", 
                    distfun = function(x) as.dist(1 - cor(t(x), method = distance)), 
                    hclustfun = function(x) hclust(x, method = "ward.D"), ColSideColors = chr1g, 
                    Colv = NA, Rowv = TRUE, notecol = "black", 
                    col = my_palette, breaks = col_breaks, key = TRUE, 
                    keysize = 1, density.info = "none", trace = "none", 
                    cexRow = 0.1, cexCol = 0.1, cex.main = 1, 
                    cex.lab = 0.1, symm = F, symkey = F, symbreaks = T, 
                    cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), 
                    cex.main = 4, margins = c(10, 10))
          print(p)
          dev.off()
        }
      }
      end_time <- Sys.time()
      print(end_time - start_time)
      reslts <- list(cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
      names(reslts) <- c("CNAmat", "hclustering")
      return(reslts)
    }
    else {
      if (distance == "euclidean") {
        hcc <- hclust(parallelDist::parDist(t(uber.mat.adj), 
                                            threads = n.cores, method = distance), method = "ward.D")
      }
      else {
        hcc <- hclust(as.dist(1 - cor(uber.mat.adj, 
                                      method = distance)), method = "ward.D")
      }
      hc.umap <- cutree(hcc, 2)
      names(hc.umap) <- colnames(results.com)
      cl.ID <- NULL
      for (i in 1:max(hc.umap)) {
        cli <- names(hc.umap)[which(hc.umap == i)]
        pid <- length(intersect(cli, preN))/length(cli)
        cl.ID <- c(cl.ID, pid)
        i <- i + 1
      }
      com.pred <- names(hc.umap)
      com.pred[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
      com.pred[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "aneuploid"
      names(com.pred) <- names(hc.umap)
      results.com.rat <- uber.mat.adj - apply(uber.mat.adj[, 
                                                           which(com.pred == "diploid")], 1, mean)
      results.com.rat <- apply(results.com.rat, 2, function(x) (x <- x - 
                                                                  mean(x)))
      results.com.rat.norm <- results.com.rat[, which(com.pred == 
                                                        "diploid")]
      dim(results.com.rat.norm)
      cf.h <- apply(results.com.rat.norm, 1, sd)
      base <- apply(results.com.rat.norm, 1, mean)
      adjN <- function(j) {
        a <- results.com.rat[, j]
        a[abs(a - base) <= 0.25 * cf.h] <- mean(a)
        a
      }
      mc.adjN <- parallel::mclapply(1:ncol(results.com.rat), 
                                    adjN, mc.cores = n.cores)
      adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), byrow = FALSE)
      colnames(adj.results) <- colnames(results.com.rat)
      mat.adj <- t(t(adj.results) - apply(adj.results, 2, mean))
      
      
      print("step 8: final prediction ...")
      if (distance == "euclidean") {
        hcc <- hclust(parallelDist::parDist(t(mat.adj), 
                                            threads = n.cores, method = distance), method = "ward.D")
      }
      else {
        hcc <- hclust(as.dist(1 - cor(mat.adj, method = distance)), 
                      method = "ward.D")
      }
      hc.umap <- cutree(hcc, 2)
      names(hc.umap) <- colnames(results.com)
      saveRDS(hcc, file = paste(sample.name, "clustering_results.rds", 
                                sep = ""))
      cl.ID <- NULL
      for (i in 1:max(hc.umap)) {
        cli <- names(hc.umap)[which(hc.umap == i)]
        pid <- length(intersect(cli, preN))/length(cli)
        cl.ID <- c(cl.ID, pid)
        i <- i + 1
      }
      com.preN <- names(hc.umap)
      com.preN[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
      com.preN[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "aneuploid"
      names(com.preN) <- names(hc.umap)
      if (WNS == "unclassified.prediction") {
        com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
        com.preN[which(com.preN == "aneuploid")] <- "c2:aneuploid:low.conf"
      }
      print("step 9: saving results...")
      "%!in%" <- function(x, y) !(x %in% y)
      ndef <- colnames(rawmat)[which(colnames(rawmat) %!in% 
                                       names(com.preN)[1:290])]
      if (length(ndef) > 0) {
        res <- data.frame(cbind(c(names(com.preN), ndef), 
                                c(com.preN, rep("not.defined", length(ndef)))))
        colnames(res) <- c("cell.names", "copykat.pred")
      }
      write.table(res, paste(sample.name, "prediction.txt", 
                             sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
      write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, 
                                                           "CNA_results.txt", sep = ""), sep = "\t", row.names = FALSE, 
                  quote = F)
      print("step 10: ploting heatmap ...")
      my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, 
                                                                  name = "RdBu")))(n = 999)
      chr <- as.numeric(Aj$DNA.adj$chrom)%%2 + 1
      rbPal1 <- colorRampPalette(c("black", "grey"))
      CHR <- rbPal1(2)[as.numeric(chr)]
      chr1 <- cbind(CHR, CHR)
      rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, 
                                                          name = "Dark2")[2:1])
      compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]
      cells <- rbind(compreN_pred, compreN_pred)
      if (ncol(mat.adj) < 3000) {
        h <- 10
      }
      else {
        h <- 15
      }
      col_breaks = c(seq(-1, -0.4, length = 50), seq(-0.4, 
                                                     -0.2, length = 150), seq(-0.2, 0.2, length = 600), 
                     seq(0.2, 0.4, length = 150), seq(0.4, 1, length = 50))
      if (distance == "euclidean") {
        jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), 
             height = h * 250, width = 4000, res = 100)
        heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) parallelDist::parDist(x, 
                                                                                            threads = n.cores, method = distance), hclustfun = function(x) hclust(x, 
                                                                                                                                                                  method = "ward.D"), ColSideColors = chr1, 
                  RowSideColors = cells, Colv = NA, Rowv = TRUE, 
                  notecol = "black", col = my_palette, breaks = col_breaks, 
                  key = TRUE, keysize = 1, density.info = "none", 
                  trace = "none", cexRow = 0.1, cexCol = 0.1, 
                  cex.main = 1, cex.lab = 0.1, symm = F, symkey = F, 
                  symbreaks = T, cex = 1, main = paste(WNS1, 
                                                       "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 
                                                                                                       10))
        legend("topright", paste("pred.", names(table(com.preN)), 
                                 sep = ""), pch = 15, col = RColorBrewer::brewer.pal(n = 8, 
                                                                                     name = "Dark2")[2:1], cex = 1)
        dev.off()
        if (plot.genes == "TRUE") {
          dim(results.com)
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom)%%2 + 1
          rbPal1g <- colorRampPalette(c("black", "grey"))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg, CHRg)
          pdf(paste(sample.name, "with_genes_heatmap.pdf", 
                    sep = ""), height = h * 2.5, width = 40)
          p = heatmap.3(t(results.com), dendrogram = "r", 
                    distfun = function(x) parallelDist::parDist(x, 
                                                                threads = n.cores, method = distance), 
                    hclustfun = function(x) hclust(x, method = "ward.D"), 
                    ColSideColors = chr1g, RowSideColors = cells, 
                    Colv = NA, Rowv = TRUE, notecol = "black", 
                    col = my_palette, breaks = col_breaks, key = TRUE, 
                    keysize = 1, density.info = "none", trace = "none", 
                    cexRow = 0.1, cexCol = 0.1, cex.main = 1, 
                    cex.lab = 0.1, symm = F, symkey = F, symbreaks = T, 
                    cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), 
                    cex.main = 4, margins = c(10, 10))
          print(p)
          dev.off()
        }
      }
      else {
        jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), 
             height = h * 250, width = 4000, res = 100)
        heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) as.dist(1 - 
                                                                                cor(t(x), method = distance)), hclustfun = function(x) hclust(x, 
                                                                                                                                              method = "ward.D"), ColSideColors = chr1, 
                  RowSideColors = cells, Colv = NA, Rowv = TRUE, 
                  notecol = "black", col = my_palette, breaks = col_breaks, 
                  key = TRUE, keysize = 1, density.info = "none", 
                  trace = "none", cexRow = 0.1, cexCol = 0.1, 
                  cex.main = 1, cex.lab = 0.1, symm = F, symkey = F, 
                  symbreaks = T, cex = 1, main = paste(WNS1, 
                                                       "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 
                                                                                                       10))
        legend("topright", paste("pred.", names(table(com.preN)), 
                                 sep = ""), pch = 15, col = RColorBrewer::brewer.pal(n = 8, 
                                                                                     name = "Dark2")[2:1], cex = 1)
        dev.off()
        if (plot.genes == "TRUE") {
          dim(results.com)
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom)%%2 + 1
          rbPal1g <- colorRampPalette(c("black", "grey"))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg, CHRg)
          pdf(paste(sample.name, "with_genes_heatmap.pdf", 
                    sep = ""), height = h * 2.5, width = 40)
          p= heatmap.3(t(results.com), dendrogram = "r", 
                    distfun = function(x) as.dist(1 - cor(t(x), 
                                                          method = distance)), hclustfun = function(x) hclust(x, 
                                                                                                              method = "ward.D"), ColSideColors = chr1g, 
                    RowSideColors = cells, Colv = NA, Rowv = TRUE, 
                    notecol = "black", col = my_palette, breaks = col_breaks, 
                    key = TRUE, keysize = 1, density.info = "none", 
                    trace = "none", cexRow = 0.1, cexCol = 0.1, 
                    cex.main = 1, cex.lab = 0.1, symm = F, symkey = F, 
                    symbreaks = T, cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), cex.main = 4, margins = c(10,10))
          print(p)
          dev.off()
        }
      }
      if (output.seg == "TRUE") {
        print("generating seg files for IGV viewer")
        thisRatio <- cbind(Aj$RNA.adj[, 1:3], mat.adj)
        Short <- NULL
        chr <- rle(thisRatio$chrom)[[2]]
        for (c in 4:ncol(thisRatio)) {
          for (x in 1:length(chr)) {
            thisRatio.sub <- thisRatio[which(thisRatio$chrom == 
                                               chr[x]), ]
            seg.mean.sub <- rle(thisRatio.sub[, c])[[2]]
            rle.length.sub <- rle(thisRatio.sub[, c])[[1]]
            num.mark.sub <- seq(1, length(rle.length.sub), 
                                1)
            loc.start.sub <- seq(1, length(rle.length.sub), 
                                 1)
            loc.end.sub <- seq(1, length(rle.length.sub), 
                               1)
            len <- 0
            j <- 1
            for (j in 1:length(rle.length.sub)) {
              num.mark.sub[j] <- rle.length.sub[j]
              loc.start.sub[j] <- thisRatio.sub$chrompos[len + 
                                                           1]
              len <- num.mark.sub[j] + len
              loc.end.sub[j] <- thisRatio.sub$chrompos[len]
              j <- j + 1
            }
            ID <- rep(colnames(thisRatio[c]), times = length(rle.length.sub))
            chrom <- rep(chr[x], times = length(rle.length.sub))
            Short.sub <- cbind(ID, chrom, loc.start.sub, 
                               loc.end.sub, num.mark.sub, seg.mean.sub)
            Short <- rbind(Short, Short.sub)
            x <- x + 1
          }
          c <- c + 1
        }
        colnames(Short) <- c("ID", "chrom", "loc.start", 
                             "loc.end", "num.mark", "seg.mean")
        head(Short)
        write.table(Short, paste(sample.name, "CNA_results.seg", 
                                 sep = ""), row.names = FALSE, quote = FALSE, 
                    sep = "\t")
      }
      end_time <- Sys.time()
      print(end_time - start_time)
      reslts <- list(res, cbind(Aj$RNA.adj[, 1:3], mat.adj), 
                     hcc)
      names(reslts) <- c("prediction", "CNAmat", "hclustering")
      return(reslts)
    }
  }
  if (genome == "mm10") {
    uber.mat.adj <- data.matrix(results.com)
    dim(uber.mat.adj)
    if (distance == "euclidean") {
      hcc <- hclust(parallelDist::parDist(t(uber.mat.adj), 
                                          threads = n.cores, method = distance), method = "ward.D")
    }
    else {
      hcc <- hclust(as.dist(1 - cor(uber.mat.adj, method = distance)), 
                    method = "ward.D")
    }
    hc.umap <- cutree(hcc, 2)
    names(hc.umap) <- colnames(results.com)
    cl.ID <- NULL
    for (i in 1:max(hc.umap)) {
      cli <- names(hc.umap)[which(hc.umap == i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i <- i + 1
    }
    com.pred <- names(hc.umap)
    com.pred[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
    com.pred[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "aneuploid"
    names(com.pred) <- names(hc.umap)
    results.com.rat <- uber.mat.adj - apply(uber.mat.adj[, 
                                                         which(com.pred == "diploid")], 1, mean)
    results.com.rat <- apply(results.com.rat, 2, function(x) (x <- x - 
                                                                mean(x)))
    results.com.rat.norm <- results.com.rat[, which(com.pred == 
                                                      "diploid")]
    dim(results.com.rat.norm)
    cf.h <- apply(results.com.rat.norm, 1, sd)
    base <- apply(results.com.rat.norm, 1, mean)
    adjN <- function(j) {
      a <- results.com.rat[, j]
      a[abs(a - base) <= 0.25 * cf.h] <- mean(a)
      a
    }
    mc.adjN <- parallel::mclapply(1:ncol(results.com.rat), 
                                  adjN, mc.cores = n.cores)
    adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), 
                          byrow = FALSE)
    colnames(adj.results) <- colnames(results.com.rat)
    mat.adj <- t(t(adj.results) - apply(adj.results, 2, 
                                        mean))
    print("step 8: final prediction ...")
    if (distance == "euclidean") {
      hcc <- hclust(parallelDist::parDist(t(mat.adj), 
                                          threads = n.cores, method = distance), method = "ward.D")
    }
    else {
      hcc <- hclust(as.dist(1 - cor(mat.adj, method = distance)), 
                    method = "ward.D")
    }
    hc.umap <- cutree(hcc, 2)
    names(hc.umap) <- colnames(results.com)
    saveRDS(hcc, file = paste(sample.name, "clustering_results.rds", 
                              sep = ""))
    cl.ID <- NULL
    for (i in 1:max(hc.umap)) {
      cli <- names(hc.umap)[which(hc.umap == i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i <- i + 1
    }
    com.preN <- names(hc.umap)
    com.preN[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
    com.preN[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "aneuploid"
    names(com.preN) <- names(hc.umap)
    if (WNS == "unclassified.prediction") {
      com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
      com.preN[which(com.preN == "aneuploid")] <- "c2:aneuploid:low.conf"
    }
    print("step 9: saving results...")
    "%!in%" <- function(x, y) !(x %in% y)
    ndef <- colnames(rawmat)[which(colnames(rawmat) %!in% 
                                     names(com.preN))]
    if (length(ndef) > 0) {
      res <- data.frame(cbind(c(names(com.preN), ndef), 
                              c(com.preN, rep("not.defined", length(ndef)))))
      colnames(res) <- c("cell.names", "copykat.pred")
    }
    write.table(res, paste(sample.name, "prediction.txt", 
                           sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(cbind(anno.mat2[, 1:7], mat.adj), paste(sample.name, 
                                                        "CNA_results.txt", sep = ""), sep = "\t", row.names = FALSE, 
                quote = F)
    print("step 10: ploting heatmap ...")
    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, 
                                                                name = "RdBu")))(n = 999)
    rownames(mat.adj) <- anno.mat2$mgi_symbol
    chrg <- as.numeric(anno.mat2$chromosome_name)%%2 + 1
    rle(as.numeric(anno.mat2$chromosome_name))
    rbPal1g <- colorRampPalette(c("black", "grey"))
    CHRg <- rbPal1g(2)[as.numeric(chrg)]
    chr1g <- cbind(CHRg, CHRg)
    rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, 
                                                        name = "Dark2")[2:1])
    compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]
    cells <- rbind(compreN_pred, compreN_pred)
    if (ncol(mat.adj) < 3000) {
      h <- 10
    }
    else {
      h <- 15
    }
    col_breaks = c(seq(-1, -0.4, length = 50), seq(-0.4, 
                                                   -0.2, length = 150), seq(-0.2, 0.2, length = 600), 
                   seq(0.2, 0.4, length = 150), seq(0.4, 1, length = 50))
    if (distance == "euclidean") {
      pdf(paste(sample.name, "with_genes_heatmap.pdf", 
                sep = ""), height = h * 2.5, width = 40)
      p= heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) parallelDist::parDist(x, 
                                                                                          threads = n.cores, method = distance), hclustfun = function(x) hclust(x, 
                                                                                                                                                                method = "ward.D"), ColSideColors = chr1g, RowSideColors = cells, 
                Colv = NA, Rowv = TRUE, notecol = "black", col = my_palette, 
                breaks = col_breaks, key = TRUE, keysize = 1, 
                density.info = "none", trace = "none", cexRow = 0.1, 
                cexCol = 0.1, cex.main = 1, cex.lab = 0.1, symm = F, 
                symkey = F, symbreaks = T, cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 10))
      print(p)
      dev.off()
    }
    else {
      pdf(paste(sample.name, "with_genes_heatmap1.pdf", 
                sep = ""), height = h * 2.5, width = 40)
      p=heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) as.dist(1 - 
                                                                              cor(t(x), method = distance)), hclustfun = function(x) hclust(x, 
                                                                                                                                            method = "ward.D"), ColSideColors = chr1g, RowSideColors = cells, 
                Colv = NA, Rowv = TRUE, notecol = "black", col = my_palette, 
                breaks = col_breaks, key = TRUE, keysize = 1, 
                density.info = "none", trace = "none", cexRow = 0.1, 
                cexCol = 0.1, cex.main = 1, cex.lab = 0.1, symm = F, 
                symkey = F, symbreaks = T, cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 10))
      print(p)
      dev.off()
    }
    if (output.seg == "TRUE") {
      print("generating seg files for IGV viewer")
      thisRatio <- cbind(anno.mat2[, c(2, 3, 1)], mat.adj)
      Short <- NULL
      chr <- rle(thisRatio$chromosome_name)[[2]]
      for (c in 4:ncol(thisRatio)) {
        for (x in 1:length(chr)) {
          thisRatio.sub <- thisRatio[which(thisRatio$chromosome_name == 
                                             chr[x]), ]
          seg.mean.sub <- rle(thisRatio.sub[, c])[[2]]
          rle.length.sub <- rle(thisRatio.sub[, c])[[1]]
          num.mark.sub <- seq(1, length(rle.length.sub), 
                              1)
          loc.start.sub <- seq(1, length(rle.length.sub), 
                               1)
          loc.end.sub <- seq(1, length(rle.length.sub), 
                             1)
          len <- 0
          j <- 1
          for (j in 1:length(rle.length.sub)) {
            num.mark.sub[j] <- rle.length.sub[j]
            loc.start.sub[j] <- thisRatio.sub$start_position[len + 
                                                               1]
            len <- num.mark.sub[j] + len
            loc.end.sub[j] <- thisRatio.sub$start_position[len]
            j <- j + 1
          }
          ID <- rep(colnames(thisRatio[c]), times = length(rle.length.sub))
          chrom <- rep(chr[x], times = length(rle.length.sub))
          Short.sub <- cbind(ID, chrom, loc.start.sub, 
                             loc.end.sub, num.mark.sub, seg.mean.sub)
          Short <- rbind(Short, Short.sub)
          x <- x + 1
        }
        c <- c + 1
      }
      colnames(Short) <- c("ID", "chrom", "loc.start", 
                           "loc.end", "num.mark", "seg.mean")
      write.table(Short, paste(sample.name, "CNA_results.seg", 
                               sep = ""), row.names = FALSE, quote = FALSE, 
                  sep = "\t")
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    reslts <- list(res, cbind(anno.mat2[, 1:7], mat.adj), 
                   hcc)
    names(reslts) <- c("prediction", "CNAmat", "hclustering")
    return(reslts)
  }
}

















copyKat_inhouse_v2 = function (rawmat = rawdata, id.type = "S", cell.line = "no", 
          ngene.chr = 5, LOW.DR = 0.05, UP.DR = 0.1, win.size = 25, 
          norm.cell.names = "", KS.cut = 0.1, sam.name = "", distance = "euclidean", 
          output.seg = "FALSE", plot.genes = "TRUE", genome = "hg20", 
          n.cores = 1) 
{
  start_time <- Sys.time()
  set.seed(1234)
  sample.name <- paste(sam.name, "_copykat_", sep = "")
  print("running copykat v1.0.8 updated 02/25/2022 introduced mm10 module, fixed typos")
  print("step1: read and filter data ...")
  print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data", 
              sep = ""))
  genes.raw <- apply(rawmat, 2, function(x) (sum(x > 0)))
  if (sum(genes.raw > 200) == 0) 
    stop("none cells have more than 200 genes")
  if (sum(genes.raw < 100) > 1) {
    rawmat <- rawmat[, -which(genes.raw < 200)]
    print(paste("filtered out ", sum(genes.raw <= 200), 
                " cells with less than 200 genes; remaining ", ncol(rawmat), 
                " cells", sep = ""))
  }
  der <- apply(rawmat, 1, function(x) (sum(x > 0)))/ncol(rawmat)
  if (sum(der > LOW.DR) >= 1) {
    rawmat <- rawmat[which(der > LOW.DR), ]
    print(paste(nrow(rawmat), " genes past LOW.DR filtering", 
                sep = ""))
  }
  WNS1 <- "data quality is ok"
  if (nrow(rawmat) < 7000) {
    WNS1 <- "low data quality"
    UP.DR <- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }
  print("step 2: annotations gene coordinates ...")
  if (genome == "hg20") {
    anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type)
  }else if (genome == "mm10") {
    anno.mat <- annotateGenes.mm10(mat = rawmat, ID.type = id.type)
    dim(rawmat)
  }
  anno.mat <- anno.mat[order(as.numeric(anno.mat$abspos), 
                             decreasing = FALSE), ]
  if (genome == "hg20") {
    HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
    toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]), 
                                               HLAs))
    if (length(toRev) > 0) {
      anno.mat <- anno.mat[-toRev, ]
    }
  }
  ToRemov2 <- NULL
  for (i in 8:ncol(anno.mat)) {
    cell <- cbind(anno.mat$chromosome_name, anno.mat[, i])
    cell <- cell[cell[, 2] != 0, ]
    if (length(as.numeric(cell)) < 5) {
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    else if (length(rle(cell[, 1])$length) < length(unique((anno.mat$chromosome_name))) | 
             min(rle(cell[, 1])$length) < ngene.chr) {
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    i <- i + 1
  }
  if (length(ToRemov2) == (ncol(anno.mat) - 7)) 
    stop("all cells are filtered")
  if (length(ToRemov2) > 0) {
    anno.mat <- anno.mat[, -which(colnames(anno.mat) %in% 
                                    ToRemov2)]
  }
  rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
  norm.mat <- log(sqrt(rawmat3) + sqrt(rawmat3 + 1))
  norm.mat <- apply(norm.mat, 2, function(x) (x <- x - mean(x)))
  colnames(norm.mat) <- colnames(rawmat3)
  print("step 3: smoothing data with dlm ...")
  dlm.sm <- function(c) {
    model <- dlm::dlmModPoly(order = 1, dV = 0.16, dW = 0.001)
    x <- dlm::dlmSmooth(norm.mat[, c], model)$s
    x <- x[2:length(x)]
    x <- x - mean(x)
  }
  test.mc <- parallel::mclapply(1:ncol(norm.mat), dlm.sm, 
                                mc.cores = n.cores)
  norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), 
                            byrow = FALSE)
  colnames(norm.mat.smooth) <- colnames(norm.mat)
  print("step 4: measuring baselines ...")
  if (cell.line == "yes") {
    print("running pure cell line mode")
    relt <- baseline.synthetic(norm.mat = norm.mat.smooth, 
                               min.cells = 10, n.cores = n.cores)
    norm.mat.relat <- relt$expr.relat
    CL <- relt$cl
    WNS <- "run with cell line mode"
    preN <- NULL
  }else if (length(norm.cell.names) > 1) {
    NNN <- length(colnames(norm.mat.smooth)[which(colnames(norm.mat.smooth) %in% 
                                                    norm.cell.names)])
    print(paste(NNN, " known normal cells found in dataset", 
                sep = ""))
    if (NNN == 0) 
      stop("known normal cells provided; however none existing in testing dataset")
    print("run with known normal...")
    basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% 
                                             norm.cell.names)], 1, median)
    print("baseline is from known input")
    d <- parallelDist::parDist(t(norm.mat.smooth), threads = n.cores, 
                               method = "euclidean")
    km <- 6
    fit <- hclust(d, method = "ward.D2")
    CL <- cutree(fit, km)
    while (!all(table(CL) > 5)) {
      km <- km - 1
      CL <- cutree(fit, k = km)
      if (km == 2) {
        break
      }
    }
    WNS <- "run with known normal"
    preN <- norm.cell.names
    norm.mat.relat <- norm.mat.smooth - basel
  }else {
    basa <- baseline.norm.cl(norm.mat.smooth = norm.mat.smooth, 
                             min.cells = 5, n.cores = n.cores)
    basel <- basa$basel
    WNS <- basa$WNS
    preN <- basa$preN
    CL <- basa$cl
    if (WNS == "unclassified.prediction") {
      basa <- baseline.GMM(CNA.mat = norm.mat.smooth, 
                           max.normal = 5, mu.cut = 0.05, Nfraq.cut = 0.99, 
                           RE.before = basa, n.cores = n.cores)
      basel <- basa$basel
      WNS <- basa$WNS
      preN <- basa$preN
    }
    norm.mat.relat <- norm.mat.smooth - basel
  }
  DR2 <- apply(rawmat3, 1, function(x) (sum(x > 0)))/ncol(rawmat3)
  norm.mat.relat <- norm.mat.relat[which(DR2 >= UP.DR), ]
  anno.mat2 <- anno.mat[which(DR2 >= UP.DR), ]
  ToRemov3 <- NULL
  for (i in 8:ncol(anno.mat2)) {
    cell <- cbind(anno.mat2$chromosome_name, anno.mat2[, 
                                                       i])
    cell <- cell[cell[, 2] != 0, ]
    if (length(as.numeric(cell)) < 5) {
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }else if (length(rle(cell[, 1])$length) < length(unique((anno.mat$chromosome_name))) | 
              min(rle(cell[, 1])$length) < ngene.chr) {
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }
    i <- i + 1
  }
  if (length(ToRemov3) == ncol(norm.mat.relat)) 
    stop("all cells are filtered")
  if (length(ToRemov3) > 0) {
    norm.mat.relat <- norm.mat.relat[, -which(colnames(norm.mat.relat) %in% 
                                                ToRemov3)]
  }
  CL <- CL[which(names(CL) %in% colnames(norm.mat.relat))]
  CL <- CL[order(match(names(CL), colnames(norm.mat.relat)))]
  print("step 5: segmentation...")
  results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat, bins = win.size, 
                      cut.cor = KS.cut, n.cores = n.cores)
  if (length(results$breaks) < 25) {
    print("too few breakpoints detected; decreased KS.cut to 50%")
    results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat, 
                        bins = win.size, cut.cor = 0.5 * KS.cut, n.cores = n.cores)
  }
  if (length(results$breaks) < 25) {
    print("too few breakpoints detected; decreased KS.cut to 75%")
    results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat, 
                        bins = win.size, cut.cor = 0.5 * 0.5 * KS.cut, n.cores = n.cores)
  }
  if (length(results$breaks) < 25) 
    stop("too few segments; try to decrease KS.cut; or improve data")
  colnames(results$logCNA) <- colnames(norm.mat.relat)
  results.com <- apply(results$logCNA, 2, function(x) (x <- x - 
                                                         mean(x)))
  RNA.copycat <- cbind(anno.mat2[, 1:7], results.com)
  write.table(RNA.copycat, paste(sample.name, "CNA_raw_results_gene_by_cell.txt", 
                                 sep = ""), sep = "\t", row.names = FALSE, quote = F)
  if (genome == "hg20") {
    print("step 6: convert to genomic bins...")
    Aj <- convert.all.bins.hg20(DNA.mat = DNA.hg20, RNA.mat = RNA.copycat, 
                                n.cores = n.cores)
    uber.mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
    print("step 7: adjust baseline ...")
    if (cell.line == "yes") {
      mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
      write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, 
                                                           "CNA_results.txt", sep = ""), sep = "\t", row.names = FALSE, 
                  quote = F)
      if (distance == "euclidean") {
        hcc <- hclust(parallelDist::parDist(t(mat.adj), 
                                            threads = n.cores, method = distance), method = "ward.D")
      }
      else {
        hcc <- hclust(as.dist(1 - cor(mat.adj, method = distance)), 
                      method = "ward.D")
      }
      saveRDS(hcc, file = paste(sample.name, "clustering_results.rds", 
                                sep = ""))
      print("step 8: ploting heatmap ...")
      my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, 
                                                                  name = "RdBu")))(n = 999)
      chr <- as.numeric(Aj$DNA.adj$chrom)%%2 + 1
      rbPal1 <- colorRampPalette(c("black", "grey"))
      CHR <- rbPal1(2)[as.numeric(chr)]
      chr1 <- cbind(CHR, CHR)
      if (ncol(mat.adj) < 3000) {
        h <- 10
      }else {
        h <- 15
      }
      col_breaks = c(seq(-1, -0.4, length = 50), seq(-0.4, 
                                                     -0.2, length = 150), seq(-0.2, 0.2, length = 600), 
                     seq(0.2, 0.4, length = 150), seq(0.4, 1, length = 50))
      if (distance == "euclidean") {
        jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), 
             height = h * 250, width = 4000, res = 100)
        heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) parallelDist::parDist(x, 
                                                                                            threads = n.cores, method = distance), hclustfun = function(x) hclust(x, 
                                                                                                                                                                  method = "ward.D"), ColSideColors = chr1, 
                  Colv = NA, Rowv = TRUE, notecol = "black", 
                  col = my_palette, breaks = col_breaks, key = TRUE, 
                  keysize = 1, density.info = "none", trace = "none", 
                  cexRow = 0.1, cexCol = 0.1, cex.main = 1, 
                  cex.lab = 0.1, symm = F, symkey = F, symbreaks = T, 
                  cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), 
                  cex.main = 4, margins = c(10, 10))
        dev.off()
        if (plot.genes == "TRUE") {
          dim(results.com)
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom)%%2 + 1
          rbPal1g <- colorRampPalette(c("black", "grey"))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg, CHRg)
          pdf(paste(sample.name, "with_genes_heatmap.pdf", 
                    sep = ""), height = h * 2.5, width = 40)
          heatmap.3(t(results.com), dendrogram = "r", 
                    distfun = function(x) parallelDist::parDist(x, 
                                                                threads = n.cores, method = distance), 
                    hclustfun = function(x) hclust(x, method = "ward.D"), 
                    ColSideColors = chr1g, Colv = NA, Rowv = TRUE, 
                    notecol = "black", col = my_palette, breaks = col_breaks, 
                    key = TRUE, keysize = 1, density.info = "none", 
                    trace = "none", cexRow = 0.1, cexCol = 0.1, 
                    cex.main = 1, cex.lab = 0.1, symm = F, symkey = F, 
                    symbreaks = T, cex = 1, main = paste(WNS1, 
                                                         "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 
                                                                                                         10))
          dev.off()
        }
      }else {
        jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), 
             height = h * 250, width = 4000, res = 100)
        heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) as.dist(1 - 
                                                                                cor(t(x), method = distance)), hclustfun = function(x) hclust(x, 
                                                                                                                                              method = "ward.D"), ColSideColors = chr1, 
                  Colv = NA, Rowv = TRUE, notecol = "black", 
                  col = my_palette, breaks = col_breaks, key = TRUE, 
                  keysize = 1, density.info = "none", trace = "none", 
                  cexRow = 0.1, cexCol = 0.1, cex.main = 1, 
                  cex.lab = 0.1, symm = F, symkey = F, symbreaks = T, 
                  cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), 
                  cex.main = 4, margins = c(10, 10))
        dev.off()
        if (plot.genes == "TRUE") {
          dim(results.com)
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom)%%2 + 1
          rbPal1g <- colorRampPalette(c("black", "grey"))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg, CHRg)
          pdf(paste(sample.name, "with_genes_heatmap.pdf", 
                    sep = ""), height = h * 2.5, width = 40)
          heatmap.3(t(results.com), dendrogram = "r", 
                    distfun = function(x) as.dist(1 - cor(t(x), 
                                                          method = distance)), hclustfun = function(x) hclust(x, 
                                                                                                              method = "ward.D"), ColSideColors = chr1g, 
                    Colv = NA, Rowv = TRUE, notecol = "black", 
                    col = my_palette, breaks = col_breaks, key = TRUE, 
                    keysize = 1, density.info = "none", trace = "none", 
                    cexRow = 0.1, cexCol = 0.1, cex.main = 1, 
                    cex.lab = 0.1, symm = F, symkey = F, symbreaks = T, 
                    cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), 
                    cex.main = 4, margins = c(10, 10))
          dev.off()
        }
      }
      end_time <- Sys.time()
      print(end_time - start_time)
      reslts <- list(cbind(Aj$RNA.adj[, 1:3], mat.adj), 
                     hcc)
      names(reslts) <- c("CNAmat", "hclustering")
      return(reslts)
    }else {
      if (distance == "euclidean") {
        hcc <- hclust(parallelDist::parDist(t(uber.mat.adj), 
                                            threads = n.cores, method = distance), method = "ward.D")
      }else {
        hcc <- hclust(as.dist(1 - cor(uber.mat.adj, 
                                      method = distance)), method = "ward.D")
      }
      hc.umap <- cutree(hcc, 2)
      names(hc.umap) <- colnames(results.com)
      cl.ID <- NULL
      for (i in 1:max(hc.umap)) {
        cli <- names(hc.umap)[which(hc.umap == i)]
        pid <- length(intersect(cli, preN))/length(cli)
        cl.ID <- c(cl.ID, pid)
        i <- i + 1
      }
      com.pred <- names(hc.umap)
      com.pred[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
      com.pred[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "aneuploid"
      names(com.pred) <- names(hc.umap)
      results.com.rat <- uber.mat.adj - apply(uber.mat.adj[, 
                                                           which(com.pred == "diploid")], 1, mean)
      results.com.rat <- apply(results.com.rat, 2, function(x) (x <- x - 
                                                                  mean(x)))
      results.com.rat.norm <- results.com.rat[, which(com.pred == 
                                                        "diploid")]
      dim(results.com.rat.norm)
      cf.h <- apply(results.com.rat.norm, 1, sd)
      base <- apply(results.com.rat.norm, 1, mean)
      adjN <- function(j) {
        a <- results.com.rat[, j]
        a[abs(a - base) <= 0.25 * cf.h] <- mean(a)
        a
      }
      mc.adjN <- parallel::mclapply(1:ncol(results.com.rat), 
                                    adjN, mc.cores = n.cores)
      adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), 
                            byrow = FALSE)
      colnames(adj.results) <- colnames(results.com.rat)
      mat.adj <- t(t(adj.results) - apply(adj.results, 
                                          2, mean))
      print("step 8: final prediction ...")
      if (distance == "euclidean") {
        hcc <- hclust(parallelDist::parDist(t(mat.adj), 
                                            threads = n.cores, method = distance), method = "ward.D")
      }else {
        hcc <- hclust(as.dist(1 - cor(mat.adj, method = distance)), 
                      method = "ward.D")
      }
      hc.umap <- cutree(hcc, 2)
      names(hc.umap) <- colnames(results.com)
      saveRDS(hcc, file = paste(sample.name, "clustering_results.rds", 
                                sep = ""))
      cl.ID <- NULL
      for (i in 1:max(hc.umap)) {
        cli <- names(hc.umap)[which(hc.umap == i)]
        pid <- length(intersect(cli, preN))/length(cli)
        cl.ID <- c(cl.ID, pid)
        i <- i + 1
      }
      com.preN <- names(hc.umap)
      com.preN[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
      com.preN[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "aneuploid"
      names(com.preN) <- names(hc.umap)
      if (WNS == "unclassified.prediction") {
        com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
        com.preN[which(com.preN == "aneuploid")] <- "c2:aneuploid:low.conf"
      }
      print("step 9: saving results...")
      "%!in%" <- function(x, y) !(x %in% y)
      ndef <- colnames(rawmat)[which(colnames(rawmat) %!in% 
                                       names(com.preN))]
      if (length(ndef) > 0) {
        res <- data.frame(cbind(c(names(com.preN), ndef), 
                                c(com.preN, rep("not.defined", length(ndef)))))
        colnames(res) <- c("cell.names", "copykat.pred")
      }
      write.table(res, paste(sample.name, "prediction.txt", 
                             sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
      write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, 
                                                           "CNA_results.txt", sep = ""), sep = "\t", row.names = FALSE, 
                  quote = F)
      print("step 10: ploting heatmap ...")
      my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, 
                                                                  name = "RdBu")))(n = 999)
      chr <- as.numeric(Aj$DNA.adj$chrom)%%2 + 1
      rbPal1 <- colorRampPalette(c("black", "grey"))
      CHR <- rbPal1(2)[as.numeric(chr)]
      chr1 <- cbind(CHR, CHR)
      rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, 
                                                          name = "Dark2")[2:1])
      compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]
      cells <- rbind(compreN_pred, compreN_pred)
      if (ncol(mat.adj) < 3000) {
        h <- 10
      }else {
        h <- 15
      }
      col_breaks = c(seq(-1, -0.4, length = 50), seq(-0.4, 
                                                     -0.2, length = 150), seq(-0.2, 0.2, length = 600), 
                     seq(0.2, 0.4, length = 150), seq(0.4, 1, length = 50))
      if (distance == "euclidean") {
        jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), 
             height = h * 250, width = 4000, res = 100)
        heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) parallelDist::parDist(x, 
                                                                                            threads = n.cores, method = distance), hclustfun = function(x) hclust(x, 
                                                                                                                                                                  method = "ward.D"), ColSideColors = chr1, 
                  RowSideColors = cells, Colv = NA, Rowv = TRUE, 
                  notecol = "black", col = my_palette, breaks = col_breaks, 
                  key = TRUE, keysize = 1, density.info = "none", 
                  trace = "none", cexRow = 0.1, cexCol = 0.1, 
                  cex.main = 1, cex.lab = 0.1, symm = F, symkey = F, 
                  symbreaks = T, cex = 1, main = paste(WNS1, 
                                                       "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 
                                                                                                       10))
        legend("topright", paste("pred.", names(table(com.preN)), 
                                 sep = ""), pch = 15, col = RColorBrewer::brewer.pal(n = 8, 
                                                                                     name = "Dark2")[2:1], cex = 1)
        dev.off()
        if (plot.genes == "TRUE") {
          dim(results.com)
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom)%%2 + 1
          rbPal1g <- colorRampPalette(c("black", "grey"))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg, CHRg)
          pdf(paste(sample.name, "with_genes_heatmap.pdf", 
                    sep = ""), height = h * 2.5, width = 40)
          
          p=heatmap.3(t(results.com), dendrogram = "r", 
                      distfun = function(x) parallelDist::parDist(x, 
                                                                  threads = n.cores, method = distance), 
                      hclustfun = function(x) hclust(x, method = "ward.D"), 
                      ColSideColors = chr1g, RowSideColors = cells, 
                      Colv = NA, Rowv = TRUE, notecol = "black", 
                      col = my_palette, breaks = col_breaks, key = TRUE, 
                      keysize = 1, density.info = "none", trace = "none", 
                      cexRow = 0.1, cexCol = 0.1, cex.main = 1, 
                      cex.lab = 0.1, symm = F, symkey = F, symbreaks = T, 
                      cex = 1, main = paste(WNS1, "; ", WNS, sep = ""), 
                      cex.main = 4, margins = c(10, 10))
          print(p)
          
          dev.off()
        }
      }else {
        jpeg(paste(sample.name, "heatmap.jpeg", sep = ""), 
             height = h * 250, width = 4000, res = 100)
        heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) as.dist(1 - 
                                                                                cor(t(x), method = distance)), hclustfun = function(x) hclust(x, 
                                                                                                                                              method = "ward.D"), ColSideColors = chr1, 
                  RowSideColors = cells, Colv = NA, Rowv = TRUE, 
                  notecol = "black", col = my_palette, breaks = col_breaks, 
                  key = TRUE, keysize = 1, density.info = "none", 
                  trace = "none", cexRow = 0.1, cexCol = 0.1, 
                  cex.main = 1, cex.lab = 0.1, symm = F, symkey = F, 
                  symbreaks = T, cex = 1, main = paste(WNS1, 
                                                       "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 
                                                                                                       10))
        legend("topright", paste("pred.", names(table(com.preN)), 
                                 sep = ""), pch = 15, col = RColorBrewer::brewer.pal(n = 8, 
                                                                                     name = "Dark2")[2:1], cex = 1)
        dev.off()
        if (plot.genes == "TRUE") {
          dim(results.com)
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom)%%2 + 1
          rbPal1g <- colorRampPalette(c("black", "grey"))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg, CHRg)
          pdf(paste(sample.name, "with_genes_heatmap.pdf", 
                    sep = ""), height = h * 2.5, width = 40)
          heatmap.3(t(results.com), dendrogram = "r", 
                    distfun = function(x) as.dist(1 - cor(t(x), 
                                                          method = distance)), hclustfun = function(x) hclust(x, 
                                                                                                              method = "ward.D"), ColSideColors = chr1g, 
                    RowSideColors = cells, Colv = NA, Rowv = TRUE, 
                    notecol = "black", col = my_palette, breaks = col_breaks, 
                    key = TRUE, keysize = 1, density.info = "none", 
                    trace = "none", cexRow = 0.1, cexCol = 0.1, 
                    cex.main = 1, cex.lab = 0.1, symm = F, symkey = F, 
                    symbreaks = T, cex = 1, main = paste(WNS1, 
                                                         "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 
                                                                                                         10))
          dev.off()
        }
      }
      if (output.seg == "TRUE") {
        print("generating seg files for IGV viewer")
        thisRatio <- cbind(Aj$RNA.adj[, 1:3], mat.adj)
        Short <- NULL
        chr <- rle(thisRatio$chrom)[[2]]
        for (c in 4:ncol(thisRatio)) {
          for (x in 1:length(chr)) {
            thisRatio.sub <- thisRatio[which(thisRatio$chrom == 
                                               chr[x]), ]
            seg.mean.sub <- rle(thisRatio.sub[, c])[[2]]
            rle.length.sub <- rle(thisRatio.sub[, c])[[1]]
            num.mark.sub <- seq(1, length(rle.length.sub), 
                                1)
            loc.start.sub <- seq(1, length(rle.length.sub), 
                                 1)
            loc.end.sub <- seq(1, length(rle.length.sub), 
                               1)
            len <- 0
            j <- 1
            for (j in 1:length(rle.length.sub)) {
              num.mark.sub[j] <- rle.length.sub[j]
              loc.start.sub[j] <- thisRatio.sub$chrompos[len + 
                                                           1]
              len <- num.mark.sub[j] + len
              loc.end.sub[j] <- thisRatio.sub$chrompos[len]
              j <- j + 1
            }
            ID <- rep(colnames(thisRatio[c]), times = length(rle.length.sub))
            chrom <- rep(chr[x], times = length(rle.length.sub))
            Short.sub <- cbind(ID, chrom, loc.start.sub, 
                               loc.end.sub, num.mark.sub, seg.mean.sub)
            Short <- rbind(Short, Short.sub)
            x <- x + 1
          }
          c <- c + 1
        }
        colnames(Short) <- c("ID", "chrom", "loc.start", 
                             "loc.end", "num.mark", "seg.mean")
        head(Short)
        write.table(Short, paste(sample.name, "CNA_results.seg", 
                                 sep = ""), row.names = FALSE, quote = FALSE, 
                    sep = "\t")
      }
      end_time <- Sys.time()
      print(end_time - start_time)
      reslts <- list(res, cbind(Aj$RNA.adj[, 1:3], mat.adj), 
                     hcc)
      names(reslts) <- c("prediction", "CNAmat", "hclustering")
      return(reslts)
    }
  }
  if (genome == "mm10") {
    uber.mat.adj <- data.matrix(results.com)
    dim(uber.mat.adj)
    if (distance == "euclidean") {
      hcc <- hclust(parallelDist::parDist(t(uber.mat.adj), 
                                          threads = n.cores, method = distance), method = "ward.D")
    }
    else {
      hcc <- hclust(as.dist(1 - cor(uber.mat.adj, method = distance)), 
                    method = "ward.D")
    }
    hc.umap <- cutree(hcc, 2)
    names(hc.umap) <- colnames(results.com)
    cl.ID <- NULL
    for (i in 1:max(hc.umap)) {
      cli <- names(hc.umap)[which(hc.umap == i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i <- i + 1
    }
    com.pred <- names(hc.umap)
    com.pred[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
    com.pred[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "aneuploid"
    names(com.pred) <- names(hc.umap)
    results.com.rat <- uber.mat.adj - apply(uber.mat.adj[, 
                                                         which(com.pred == "diploid")], 1, mean)
    results.com.rat <- apply(results.com.rat, 2, function(x) (x <- x - 
                                                                mean(x)))
    results.com.rat.norm <- results.com.rat[, which(com.pred == 
                                                      "diploid")]
    dim(results.com.rat.norm)
    cf.h <- apply(results.com.rat.norm, 1, sd)
    base <- apply(results.com.rat.norm, 1, mean)
    adjN <- function(j) {
      a <- results.com.rat[, j]
      a[abs(a - base) <= 0.25 * cf.h] <- mean(a)
      a
    }
    mc.adjN <- parallel::mclapply(1:ncol(results.com.rat), 
                                  adjN, mc.cores = n.cores)
    adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), 
                          byrow = FALSE)
    colnames(adj.results) <- colnames(results.com.rat)
    mat.adj <- t(t(adj.results) - apply(adj.results, 2, 
                                        mean))
    print("step 8: final prediction ...")
    if (distance == "euclidean") {
      hcc <- hclust(parallelDist::parDist(t(mat.adj), 
                                          threads = n.cores, method = distance), method = "ward.D")
    }
    else {
      hcc <- hclust(as.dist(1 - cor(mat.adj, method = distance)), 
                    method = "ward.D")
    }
    hc.umap <- cutree(hcc, 2)
    names(hc.umap) <- colnames(results.com)
    saveRDS(hcc, file = paste(sample.name, "clustering_results.rds", 
                              sep = ""))
    cl.ID <- NULL
    for (i in 1:max(hc.umap)) {
      cli <- names(hc.umap)[which(hc.umap == i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i <- i + 1
    }
    com.preN <- names(hc.umap)
    com.preN[which(hc.umap == which(cl.ID == max(cl.ID)))] <- "diploid"
    com.preN[which(hc.umap == which(cl.ID == min(cl.ID)))] <- "aneuploid"
    names(com.preN) <- names(hc.umap)
    if (WNS == "unclassified.prediction") {
      com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
      com.preN[which(com.preN == "aneuploid")] <- "c2:aneuploid:low.conf"
    }
    print("step 9: saving results...")
    "%!in%" <- function(x, y) !(x %in% y)
    ndef <- colnames(rawmat)[which(colnames(rawmat) %!in% 
                                     names(com.preN))]
    if (length(ndef) > 0) {
      res <- data.frame(cbind(c(names(com.preN), ndef), 
                              c(com.preN, rep("not.defined", length(ndef)))))
      colnames(res) <- c("cell.names", "copykat.pred")
    }
    write.table(res, paste(sample.name, "prediction.txt", 
                           sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(cbind(anno.mat2[, 1:7], mat.adj), paste(sample.name, 
                                                        "CNA_results.txt", sep = ""), sep = "\t", row.names = FALSE, 
                quote = F)
    print("step 10: ploting heatmap ...")
    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, 
                                                                name = "RdBu")))(n = 999)
    rownames(mat.adj) <- anno.mat2$mgi_symbol
    chrg <- as.numeric(anno.mat2$chromosome_name)%%2 + 1
    rle(as.numeric(anno.mat2$chromosome_name))
    rbPal1g <- colorRampPalette(c("black", "grey"))
    CHRg <- rbPal1g(2)[as.numeric(chrg)]
    chr1g <- cbind(CHRg, CHRg)
    rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, 
                                                        name = "Dark2")[2:1])
    compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]
    cells <- rbind(compreN_pred, compreN_pred)
    if (ncol(mat.adj) < 3000) {
      h <- 10
    }
    else {
      h <- 15
    }
    col_breaks = c(seq(-1, -0.4, length = 50), seq(-0.4, 
                                                   -0.2, length = 150), seq(-0.2, 0.2, length = 600), 
                   seq(0.2, 0.4, length = 150), seq(0.4, 1, length = 50))
    if (distance == "euclidean") {
      pdf(paste(sample.name, "with_genes_heatmap.pdf", 
                sep = ""), height = h * 2.5, width = 40)
      heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) parallelDist::parDist(x, 
                                                                                          threads = n.cores, method = distance), hclustfun = function(x) hclust(x, 
                                                                                                                                                                method = "ward.D"), ColSideColors = chr1g, RowSideColors = cells, 
                Colv = NA, Rowv = TRUE, notecol = "black", col = my_palette, 
                breaks = col_breaks, key = TRUE, keysize = 1, 
                density.info = "none", trace = "none", cexRow = 0.1, 
                cexCol = 0.1, cex.main = 1, cex.lab = 0.1, symm = F, 
                symkey = F, symbreaks = T, cex = 1, main = paste(WNS1, 
                                                                 "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 
                                                                                                                 10))
      dev.off()
    }
    else {
      pdf(paste(sample.name, "with_genes_heatmap1.pdf", 
                sep = ""), height = h * 2.5, width = 40)
      heatmap.3(t(mat.adj), dendrogram = "r", distfun = function(x) as.dist(1 - 
                                                                              cor(t(x), method = distance)), hclustfun = function(x) hclust(x, 
                                                                                                                                            method = "ward.D"), ColSideColors = chr1g, RowSideColors = cells, 
                Colv = NA, Rowv = TRUE, notecol = "black", col = my_palette, 
                breaks = col_breaks, key = TRUE, keysize = 1, 
                density.info = "none", trace = "none", cexRow = 0.1, 
                cexCol = 0.1, cex.main = 1, cex.lab = 0.1, symm = F, 
                symkey = F, symbreaks = T, cex = 1, main = paste(WNS1, 
                                                                 "; ", WNS, sep = ""), cex.main = 4, margins = c(10, 
                                                                                                                 10))
      dev.off()
    }
    if (output.seg == "TRUE") {
      print("generating seg files for IGV viewer")
      thisRatio <- cbind(anno.mat2[, c(2, 3, 1)], mat.adj)
      Short <- NULL
      chr <- rle(thisRatio$chromosome_name)[[2]]
      for (c in 4:ncol(thisRatio)) {
        for (x in 1:length(chr)) {
          thisRatio.sub <- thisRatio[which(thisRatio$chromosome_name == 
                                             chr[x]), ]
          seg.mean.sub <- rle(thisRatio.sub[, c])[[2]]
          rle.length.sub <- rle(thisRatio.sub[, c])[[1]]
          num.mark.sub <- seq(1, length(rle.length.sub), 
                              1)
          loc.start.sub <- seq(1, length(rle.length.sub), 
                               1)
          loc.end.sub <- seq(1, length(rle.length.sub), 
                             1)
          len <- 0
          j <- 1
          for (j in 1:length(rle.length.sub)) {
            num.mark.sub[j] <- rle.length.sub[j]
            loc.start.sub[j] <- thisRatio.sub$start_position[len + 
                                                               1]
            len <- num.mark.sub[j] + len
            loc.end.sub[j] <- thisRatio.sub$start_position[len]
            j <- j + 1
          }
          ID <- rep(colnames(thisRatio[c]), times = length(rle.length.sub))
          chrom <- rep(chr[x], times = length(rle.length.sub))
          Short.sub <- cbind(ID, chrom, loc.start.sub, 
                             loc.end.sub, num.mark.sub, seg.mean.sub)
          Short <- rbind(Short, Short.sub)
          x <- x + 1
        }
        c <- c + 1
      }
      colnames(Short) <- c("ID", "chrom", "loc.start", 
                           "loc.end", "num.mark", "seg.mean")
      write.table(Short, paste(sample.name, "CNA_results.seg", 
                               sep = ""), row.names = FALSE, quote = FALSE, 
                  sep = "\t")
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    reslts <- list(res, cbind(anno.mat2[, 1:7], mat.adj), 
                   hcc)
    names(reslts) <- c("prediction", "CNAmat", "hclustering")
    return(reslts)
  }
}


