library(Basic4Cseq)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm9)
library(GenomicAlignments)
library(caTools)

args<-commandArgs(T)

firstcutter = args[1]
secondcutter = args[2]
species = args[4]
if(species=='Mus_musculus'){
createVirtualFragmentLibrary(chosenGenome = Mmusculus,
firstCutter = firstcutter,
secondCutter = secondcutter,
readLength = 50,
libraryName = args[3],
useAllData = F)} else {
	createVirtualFragmentLibrary(chosenGenome = Hsapiens,
firstCutter = firstcutter,
secondCutter = secondcutter,
readLength = 50,
libraryName = args[3],
useAllData = F)
}