library(BSgenome.Mmusculus.UCSC.mm9)
library(ShortRead)

getChIPSeqPeaks<-function(chip.bam.file,contr.bam.file,fragment.size=200,chromosome,output_peak_file){
	
		chr<-chromosome
		print(paste("Reading in :",chr,sep=""))
		fragment.size=fragment.size
		chr.size<-seqlengths(Mmusculus)[chr]
		which <- GRanges(seqnames=chr,IRanges(1,chr.size))
		what <- c("rname", "strand", "pos", "qwidth","seq","qual")
		param <- ScanBamParam(which = which,what=what)
	
		##Reading in the BAM files## 	
		chip.bam <- scanBam(chip.bam.file, param=param)
		contr.bam <- scanBam(contr.bam.file, param=param)
		name <- names(chip.bam)
		###Making GRanges for ChIP##
		chip.qual    <-QualityScaledBStringSet(chip.bam[[name]]$seq, chip.bam[[name]]$qual)
		chip.read.length <-width(chip.qual)[1]
		chip.qa   <- FastqQuality(quality(chip.qual))
		chip.qa   <- as(chip.qa, "matrix")
		chip.qa.avg   <-as.integer(rowMeans(chip.qa))	
		chip.GRanges<-GRanges(seqnames=as.vector(chip.bam[[name]]$rname),IRanges(start=chip.bam[[name]]$pos,width=chip.read.length),strand=chip.bam[[name]]$strand,qual=chip.qa.avg)	
		chip.GRanges.filtered <-chip.GRanges[elementMetadata(chip.GRanges)$qual >=10]
	
		chip.GRanges.dup <- chip.GRanges.filtered
		chip.GRanges.dup <- chip.GRanges.dup[end(chip.GRanges.dup) < chr.size-fragment.size,]
		chip.GRanges.dup <- chip.GRanges.dup[start(chip.GRanges.dup) >fragment.size,]
		chip.GRanges.dup <- resize(chip.GRanges.dup, width = fragment.size)	
		seqlengths(chip.GRanges.dup)<-chr.size
		chip.lib.size <- length(chip.GRanges.dup)
	
		###Making GRanges for Control##
	
		contr.qual <-QualityScaledBStringSet(contr.bam[[name]]$seq, contr.bam[[name]]$qual)
		contr.read.length <-width(contr.qual)[1]
		contr.qa   <- FastqQuality(quality(contr.qual))
		contr.qa   <- as(contr.qa, "matrix")
		contr.qa.avg   <-as.integer(rowMeans(contr.qa))	
		contr.GRanges<-GRanges(seqnames=as.vector(contr.bam[[name]]$rname),IRanges(start=contr.bam[[name]]$pos,width=contr.read.length),strand=contr.bam[[name]]$strand,qual=contr.qa.avg)	
		contr.GRanges.filtered <-contr.GRanges[elementMetadata(contr.GRanges)$qual >=10]
	
		contr.GRanges.dup <- contr.GRanges.filtered
		contr.GRanges.dup <- contr.GRanges.dup[end(contr.GRanges.dup) < chr.size-fragment.size,]
		contr.GRanges.dup <- contr.GRanges.dup[start(contr.GRanges.dup) >fragment.size,]
		contr.GRanges.dup <- resize(contr.GRanges.dup, width = fragment.size)	
		seqlengths(contr.GRanges.dup)<-chr.size
		control.lib.size <- length(contr.GRanges.dup)
	
		###Select GRanges from the same size of sequencing depth##
		chip.GRanges.used =GRanges()
		contr.GRanges.used = GRanges()
	
		if(control.lib.size >= chip.lib.size){
			chip.GRanges.used  <-chip.GRanges.dup
			contr.GRanges.used <-sample(contr.GRanges.dup,chip.lib.size,replace=FALSE)
		}
		if(chip.lib.size > control.lib.size){
			chip.GRanges.used  <-chip.GRanges.dup
			contr.GRanges.used <-contr.GRanges.dup
			#print(paste(chr,"has number of read in ChIP less than a control=",control.lib.size-chip.lib.size,"from",control.lib.size,":",chip.lib.size))
		}
	
		###Making coverage for ChIP and Control for unique Mapping#####
		chip.GRanges.uniq <- chip.GRanges.used[!duplicated(paste(start(chip.GRanges.used),strand(chip.GRanges.used),sep=":")),]	
		chip.cov <-coverage(chip.GRanges.uniq)
		contr.GRanges.uniq <-contr.GRanges.used[!duplicated(paste(start(contr.GRanges.used),strand(contr.GRanges.used),sep=":")),]
		contr.cov<-coverage(contr.GRanges.uniq)
	
		##Making Islands###
		chip.islands<-slice(chip.cov[[chr]],lower=1)
		contr.islands<-slice(contr.cov[[chr]],lower=1)
		############Peak Spliting#####################
		peak.islands.chip<-chip.islands
		peak.islands.contr<-Views(contr.cov[[chr]],start(peak.islands.chip),end(peak.islands.chip))
	
		enriched <-viewMaxs(peak.islands.chip) > viewMaxs(peak.islands.contr)
		chip.enriched <- peak.islands.chip[enriched]
	
		regions <-data.frame()
		chip.enriched.sub1 <-chip.enriched[width(chip.enriched) <=1000]
	
		if(length(chip.enriched.sub1) > 0){	
			regions<-data.frame(start=start(chip.enriched.sub1),end=end(chip.enriched.sub1))
		}
			chip.enriched.sub2 <-chip.enriched[width(chip.enriched) > 1000]
	
		if(length(chip.enriched.sub2) > 0){
		
			for(i in 1:length(chip.enriched.sub2)){
			
			data.rle = chip.enriched.sub2[[i]]
			y=as.vector(data.rle)
			x=1:length(y)
			
			ss=smooth.spline(x, y, df=2,spar=0.60)
			
			my_vector<-rep(0,length(y))
			my_vector[ss$x]<-ss$y
			
			my_vector[ss$x]<-ss$y
			rve = rle(my_vector)
			pos = cumsum(rve$lengths)
			val.left = c(rve$values[1], rve$values[1:(length(rve$values)-1)])
			val.right = c(rve$values[2:(length(rve$values))], rve$values[length(rve$values)])
			
			pos.max = pos[rve$values>val.left & rve$values>val.right] - floor(rve$lengths[rve$values>val.left & rve$values>val.right]/2-0.5)
			pos.min = pos[rve$values<=val.left & rve$values<=val.right] - floor(rve$lengths[rve$values<=val.left & rve$values<=val.right]/2-0.5)
			value.max<-rve$values[rve$values>val.left & rve$values>val.right]
			value.min<-rve$values[rve$values<=val.left & rve$values<=val.right]
			
			value.min =value.min[-1]
			value.min =value.min[-length(value.min)]
			
			pos.min =pos.min[-1]
			pos.min =pos.min[-length(pos.min)]
			
			value.min.f=c(value.min,0)		
			z<-c(value.min.f > 0.8*value.max)
			pos.min<-c(pos.min,pos[length(pos)])
			f.pos<-pos.min[z==FALSE]
			f.pos.uniq<-unique(f.pos)
			
				if(length(f.pos.uniq)>1){
					f.pos.f<-c(f.pos.uniq[2:length(f.pos.uniq)],f.pos.uniq[length(f.pos.uniq)])
					f<-c(f.pos.f-f.pos.uniq >400)
					final.position<-f.pos.uniq[f==TRUE]
					final.position<-c(pos[1],final.position,pos[length(pos)])
					start = c(final.position[1:(length(final.position)-1)])
					end = c(final.position[2:(length(final.position))])
					frame<-data.frame(start=start+start(chip.enriched.sub2[i]),end=end+start(chip.enriched.sub2[i]))
					regions <-rbind(regions,frame)
				}else{
					final.position<-c(pos[1],f.pos.uniq,pos[length(pos)])
					final.position.uniq<-unique(final.position)
					start = c(final.position.uniq[1])
					end = c(final.position.uniq[2])
					frame<-data.frame(start=start+start(chip.enriched.sub2[i]),end=end+start(chip.enriched.sub2[i]))
					regions <-rbind(regions,frame)
				}	
			}
		}
		############Peak Detection####################
		chip.islands <-Views(chip.cov[[chr]],regions$start,regions$end)
		contr.islands<-Views(contr.cov[[chr]],regions$start,regions$end)
	
		enriched <-viewMaxs(chip.islands) > viewMaxs(contr.islands)
		chip.enriched <- chip.islands[enriched]
	
		############Get Number of reads##############
		chr.peaks<-data.frame(chromosome=chr,start=start(chip.enriched),end=end(chip.enriched))
		region.iranges <-RangedData(space=as.character(chr),IRanges(start=start(chip.enriched),end=end(chip.enriched)),idx=1:length(chip.enriched))
		####before normalized##
		bchip <- countOverlaps(region.iranges,subject=chip.GRanges.uniq,minoverlap = 100L)
		bcontr <- countOverlaps(region.iranges,subject=contr.GRanges.uniq,minoverlap = 100L)
		###############################################
		chip.peak.loc <- viewWhichMaxs(chip.enriched) 
		## Choosing center on flat peaks 
		flat <- width(viewRangeMaxs(chip.enriched))>1
		chip.peak.loc[flat] <- chip.peak.loc[flat]+floor((width(viewRangeMaxs(chip.enriched))[flat])/2)
	
		chip.peak.val  <- viewMaxs(chip.enriched)
		contr.peak.val <- as.numeric(contr.cov[[chr]][chip.peak.loc])
	
		chr.peaks$peak_center<-chip.peak.loc
		##############################################
		chr.peaks$chip.Ntags <-bchip
		chr.peaks$contr.Ntags <-bcontr
		##############################################
		chr.peaks$chip.max<-chip.peak.val
		chr.peaks$control.max<-contr.peak.val
		chr.peaks<-subset(chr.peaks,chip.Ntags>contr.Ntags)
	
		p.values.Ntags <- dnbinom(chr.peaks$contr.Ntags,chr.peaks$chip.Ntags,0.5,log=F)
		p.adjusted.Ntags <- p.adjust(p.values.Ntags, method="BH")
	
		chr.peaks$padj.Ntags<-p.adjusted.Ntags
	
		p.values.maxCov <- dnbinom(chr.peaks$control.max,chr.peaks$chip.max,0.5,log=F)
		p.adjusted.maxCov <- p.adjust(p.values.maxCov , method="BH")
	
		chr.peaks$padj.maxCov<-p.adjusted.maxCov 
		chr.peaks$fold.Ntags<-(chr.peaks$chip.Ntags+1)/(chr.peaks$contr.Ntags+1)
		chr.peaks$fold.max<-(chr.peaks$chip.max+1)/(chr.peaks$control.max+1)
	
		write.table(chr.peaks,file=output_peak_file,append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
		print(paste("Finish! :",chr,sep=""))
}

