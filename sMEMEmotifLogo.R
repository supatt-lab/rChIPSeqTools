library(seqLogo)
library(Biostrings)

getMemeMotifPssm = function( memeOut) {
	######get input sites#####################
	infoLine<-grep("^data",memeOutput,perl=T)
	info<- unlist( strsplit( memeOutput[infoLine],"N=") )
	TotalNsite<-info[2]
	TotalNsite<-gsub("\\s+","",TotalNsite, perl=TRUE)
	##########################################
	lines <- grep( "^MOTIF\\s+\\d", memeOutput, perl=T)
	pssms <- list()
	
	for (i in 1:length(lines)) {
		m.start <- paste("Motif ", i, " position-specific probability matrix", sep="")
		m.line1 <- grep(m.start, memeOut)
		
		line <- memeOutput[ lines[i] ]
		splitted <- unlist( strsplit( line, "[\\t\\s]", perl=T ) )
		splitted <- splitted[ splitted != "" ]
		sites <- as.integer( splitted[ 8 ] )
		
		
		if ( length( m.line1 ) > 0 ) {
			m.desc <- unlist( strsplit( memeOut[m.line1 + 2], " " ) )
			winLen <- as.numeric(m.desc[6])
			e.val  <- as.numeric(m.desc[10])
			pssm <- matrix(, nrow=winLen, ncol=4)
			for (j in 1:winLen) {
				pssm[j,] <- as.numeric( unlist( strsplit(memeOut[m.line1 + 2 + j], " ") )[c(2,4,6,8)] )
			}
			tmp.p <- list()
			tmp.p$pssm <- pssm
			tmp.p$e.val <- e.val
			tmp.p$site<-sites
			tmp.p$total<-TotalNsite
			pssms[[i]] <- tmp.p
			rm ( tmp.p )
		} else {
			tmp.p <- list()
			tmp.p$pssm <- NULL
			tmp.p$e.val <- 99999
			tmp.p$e.site <- NULL
			tmp.p$total<-TotalNsite
			pssms[[i]] <- tmp.p
			
		}
	}
	return ( pssms )
}

draw_logo<-function(memeOut){
	pssms<-getMemeMotifPssm (memeOutput)
	pdf(file=logo.out.file, width=12, height=8, onefile=T, bg="transparent",family = "Helvetica",fonts = NULL)
	
	#modify the seqLogo function
	mySeqLogo = seqLogo::seqLogo
	
	bad = (sapply( body(mySeqLogo), "==", "grid.newpage()") |
				sapply( body(mySeqLogo), "==", "par(ask = FALSE)"))
	body(mySeqLogo)[bad] = NULL
	#define dimensions
	ncol <- 3
	nrow <- 4
	
	#create plot
	grid.newpage()
	
	i<-0
	for(row in 1:nrow){
		for(col in 1:ncol){
			i<-i+1
			vp <- viewport(x = (col-1)/ncol, y = 1-(row-1)/nrow, w = 1/ncol, h = 
							1/nrow, just = c("left", "top"))
			pushViewport(vp)
			if(i <=length(pssms)){
				pssm<-t(pssms[[i]]$pssm)
				mySeqLogo(pssm,xfontsize=6,yfontsize=6)
				text<-paste("Sites=",pssms[[i]]$site,"/","Total inputs=",pssms[[i]]$total)
				grid.text(text, y = unit(1, "npc") -unit(0.2,"inches"),just="center",gp=gpar(col="red", fontsize=8))
				upViewport()
			}
		}
	}
	
	#create plot
	grid.newpage()
	
	i<-0
	for(row in 1:nrow){
		for(col in 1:ncol){
			i<-i+1
			vp <- viewport(x = (col-1)/ncol, y = 1-(row-1)/nrow, w = 1/ncol, h = 
							1/nrow, just = c("left", "top"))
			pushViewport(vp)
			if(i <=length(pssms)){
				pssm<-t(pssms[[i]]$pssm)
				mySeqLogo(reverseComplement(pssm),xfontsize=6,yfontsize=6)
				text<-paste("Sites=",pssms[[i]]$site,"/","Total inputs=",pssms[[i]]$total)
				grid.text(text, y = unit(1, "npc") -unit(0.2,"inches"),just="center",gp=gpar(col="red", fontsize=8))
				upViewport()
			}
		}
	}
	dev.off()
}


