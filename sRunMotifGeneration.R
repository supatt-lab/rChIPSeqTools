source("sMEMEmotifLogo.R")


meme.files<-c("example.data/Demo.meme.result.txt")

meme.results <- scan(file=meme.files, what='character', sep="\n")
memeOutput<-meme.results
logo.out.file<-paste(meme.files,".logo.pdf",sep="")
draw_logo(memeOutput)


