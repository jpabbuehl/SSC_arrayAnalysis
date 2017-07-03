analysis <- function(gse,chipType,file=gse) {
  gse='GSE51068'
  chipType='HGU219'
	library(aroma.affymetrix)
	library(gtools)

	path<-getwd()

	cdf <- AffymetrixCdfFile$byChipType(chipType)
	cs <- AffymetrixCelSet$byName(gse, cdf=cdf)
  
	arrayname<- getFullNames(cs)
	bc <- RmaBackgroundCorrection(cs)
	csBC <- process(bc,verbose=verbose)
	qn <- QuantileNormalization(csBC, typesToUpdate="pm")
	#qn <- ScaleNormalization(csBC,targetAvg=4400)
	#qn<- ProbeLevelTransform(csBC)
	csN <- process(qn, verbose=verbose)
	plm <- RmaPlm(csN)
	fit(plm, verbose=verbose)

	## Check quality RLE and Nuse ###
	qam <- QualityAssessmentModel(plm)

	filename <- sprintf("%s,plotRle.png", getName(qam));
	arrays <- 1:length(arrayname)
	png(filename, width=150+30*length(arrays), height=500);
	plotRle(qam,show.names=FALSE);
  axis(side=1, at=seq(along=arrays), labels=1:length(arrays), las=2);
	dev.off();

	filename <- sprintf("%s,plotNuse.png", getName(qam));
	arrays <- 1:length(arrayname)
	png(filename, width=150+30*length(arrays), height=500);
  plotNuse(qam);
  axis(side=3, at=seq(along=arrays), labels=1:length(arrays), las=2);
	dev.off();
	
	ces <- getChipEffectSet(plm)
	#Exclude bad arrays
	plotRle(qam,show.names=FALSE);
	axis(side=1, at=seq(along=arrays), labels=1:length(arrays), las=2);
	ces <- arrayexclusion(ces)

	## load phenodata ###
	data1 <- extractExpressionSet(ces, units=NULL, addNames=TRUE)
	samples <- read.table(paste("metadata\\",gse,".csv",sep=""),row.names=1, header = TRUE, sep=";")
	coordinates <- match(row.names(samples),as.factor(sampleNames(data)))
	samples<-samples[coordinates[!is.na(coordinates)],]
	pd <- new("AnnotatedDataFrame", data=samples)
	phenoData(data) <- pd

	## Save expressionSet for further analysis ##
	save(data,file=paste("eset/",file,".Rda",sep=""))
	cat('finish\n')
}

arrayexclusion <- function(ces){
   input <- ask("Which array to exclude ? if several, use comma for each ")
   input<-as.numeric(as.vector(unlist(strsplit(input,","))))
   if(length(input)<1) { return (ces) }
   else { return (extract(ces, -input))}
}