#extra line to include channelrhodopsin stimulations on heatmap
#from v14, the fwd/back angles thing works
#change from v11, color scale for heatmap is correct (previous was reversed)
#change from v12, this one corrects the startframes in forward and backward durations/initiations

### set moving average number (the number of frames you want to average)
spanformedian <- 3
### set threshold number, fwd n>threshold, back n < -threshold, pause in between
fwdthreshold <- 1
backthreshold <- 1
### count only durations that are longer than x
durationlength <- 2

################################################################################## 
extractdata<- function(x, ...) {
	#prep sublist strings
	sublist<-NULL
	if(length(c(...))!=0)
	{
		sublist<-paste("$",paste(..., sep="$"),sep="")
	}
	#below line is old slow code
	#returnobjct<-NULL
	#this is faster
	returnobjct<-as.list(rep(NA,length(x)))
	for(i in 1:length(x))
	{
		#append is slow
		#returnobjct <-?(returnobjct,list(eval(parse(text=paste("x[[i]]",sublist,sep="")))))

		returnobjct[i] <-list(eval(parse(text=paste("x[[i]]",sublist,sep=""))))
		#x[[i]]$manualdata$crudegenotype
	}
	returnobjct
	#sublist
}

##########################################################


samplenameslist<-system("ls | grep '.txt'",intern = TRUE)
#eliminate memo.txt
for(i in 1:length(samplenameslist))
{
	if(samplenameslist[i]=="memo.txt")
	{
		samplenameslist[i]<-list(NULL)
	}
}
samplenamesvec<-unlist(samplenameslist)

##########################################################
#put each data into list
# This part can be used for common application
data_list <- NULL
for(i in 1:length(samplenameslist))
{
	samplename<-strsplit(samplenamesvec[i],".txt")[[1]]
	temp <- read.csv(samplenameslist[i], header = T, sep = ",")
	tempvel <- temp$velocity
	
		#deal with outliner. 
	#big negative values are caused by failure of plate shift detection
	temp$velocity<-ifelse(abs(temp$velocity)>50,NA,temp$velocity)
	
	templist <- list(genotype=strsplit(samplename,"_")[[1]][1], samplenumber=strsplit(samplename,"_")[[1]][2], rawdata= temp)
	data_list<-append(data_list, list(templist))
	
	
}


# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "velocity")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


velocitymatrix<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedmatrix<-NULL
for(i in 1:length(data_list))
{

	correctedmatrix<-extractdata(data_list, "rawdata","velocity")[[i]]
	velocitymatrix[1:length(correctedmatrix),i] <- as.matrix(correctedmatrix)
}

write.csv(velocitymatrix, "velocitymatrix.csv")

newvelocitymatrix <- data.matrix(velocitymatrix, rownames.force = NA)

##############################draw histogram
quartz()

hist(newvelocitymatrix, breaks=seq(-70,70,by=1), freq=F, xlim=c(-70,70), ylim=c(0,0.15))
dev.copy(png, file=paste("hist.png", sep=""),width=800, height=300)
dev.off()

################plot heatmap
library(RColorBrewer)
cols<-colorRampPalette(brewer.pal(11,"RdBu"))(100)
quartz()
#par(mfrow=c(1,2))
#plot.new()
image(0:nrow(velocitymatrix), 1:ncol(velocitymatrix), velocitymatrix,zlim=c(-70,70), col=cols, axes=FALSE, xlab="", ylab="") 
#the next rows if included will plot lines showing channelrhodopsin stim periods
#stims <- c(200, 210, 510, 520, 820, 830)
#stims <- c(200, 300, 600, 700, 1000, 1100)
#abline(v=(stims), lty=2, col="blue")
axis(side=1, las=1, at=100*0:1500)
dev.copy(png, file=paste("heatmap3.png", sep=""),width=800, height=300)
dev.off()

#plot color scale bar for heatmap


quartz()
my.colors = colorRampPalette(c(cols))
z=matrix(1:100,nrow=1)
x=1
y=seq(-70, 70, len=100) 
image(x,y,z,col=my.colors(100),axes=FALSE,xlab="",ylab="")
axis(4, at=c(-60,-40, -20, 0, 20, 40,60))
dev.copy(png, file=paste("heatmapcolorbar.png", sep=""), width=100, height=400)
dev.off()

####remove NAs from original data

velocitymatrixNoNA <- apply(velocitymatrix, 2, na.omit)

nrowcorrected <- NULL
for(i in 1:length(data_list))
{
	noNAcorrectedmatrix<-velocitymatrixNoNA[[i]]
	nrowcorrected <-cbind(nrowcorrected, length(noNAcorrectedmatrix))
}

maxrow<-max(nrowcorrected)
maxrow

which(nrowcorrected==maxrow)

#make matrix from uneven length data
noNAvelocitymatrix <-matrix(NA, ncol=length(data_list), nrow=maxrow)
noNAcorrectedmatrix <- NULL
for(i in 1:length(data_list))
{
	noNAcorrectedmatrix<-velocitymatrixNoNA[[i]]
	noNAvelocitymatrix[1:length(noNAcorrectedmatrix),i] <- as.matrix(noNAcorrectedmatrix)
	
}


############################################################
#calculate runningmedian and write it into a matrix
avgnoNAvelocitymatrix<-matrix(NA, ncol=length(data_list),nrow= maxrow)
noNAcorrectedmatrix<-NULL
for(i in 1:length(data_list))
{ 
	noNAcorrectedmatrix<-velocitymatrixNoNA[[i]]
	median <- function(x) {runmed(x, spanformedian)}
	movavgmatrix <- median(noNAcorrectedmatrix)
	avgnoNAvelocitymatrix[1:length(movavgmatrix),i] <- as.matrix(movavgmatrix)
}

#######################for forward initiations
rlefwd <- function(x) {
	rle(x>1)
}

f <-apply(avgnoNAvelocitymatrix,2,rlefwd)

newfwd <- function(x) {data.frame(number=x$values, duration=x$lengths, startframe=cumsum(x$lengths), row.names=NULL)}
newnewfwd <- lapply(f, newfwd)

#previous start frame is off (number is end of period not start)so first subset all of third column minus the last element
functionB <- function(x){head(x[,3], -1)} 
trythis <- lapply(newnewfwd, functionB)

#add 1 to the start of the vector
functiontry <- function(x){new <- c(1, x)}
trynewfwd <- lapply(trythis, functiontry)

#add the corrected vector to the data frame
trypleaseworkfwd <- mapply(cbind, newnewfwd, "newstart"=trynewfwd, SIMPLIFY=FALSE)

remove3 <- function(x){x[,-3]}
fwdframe <- lapply(trypleaseworkfwd, remove3)

#sapply(newnewfwd, write.table, file="fwdlogical.csv", append=TRUE, sep=",")

#this part filters out short non-fwd eg(if FWD=5frames, non-FWD=1-2frames, FWD=8 frames, will count as FWD=15)
tryfunction <- function(x) {
	transformer <- transform(x, newcolumn=ifelse((x[,1]	==FALSE & x[,2]<4), TRUE, x[,1]))
	}
pleasework <- lapply(fwdframe,tryfunction)


#convert TRUE/FALSE into 1/0
integerize <- function(x){
	x*1}
nological <- lapply(pleasework,integerize)

#drop first column
columndrop <- function(x){(x[,-1])}
newnological <- lapply(nological, columndrop)

nnl <- newnological

#omitNAs
omit <- function(x){
	x[complete.cases(x),]	
}

omitNAs <- lapply(nnl,omit)

#rebuilddataframe
#rebuild <- function(x){
#	data.frame(x$n, x$startframe, x$newcolumn)
#}
#
#newdf <- lapply(omitNAs,rebuild)

#make an index column for 1s and 0s
func.rle <- function(x){
rle(x$newcolumn)
}

rle.idx <- lapply(omitNAs, func.rle)

new <-function(x){rep(seq(x$lengths), x$lengths)}
newnew <- lapply(rle.idx, new)

#add index column to data frame
newnewdf <- mapply(cbind, omitNAs, "new"=newnew, SIMPLIFY=FALSE)
library(plyr)

#add rows when index value is the same
function.ans <- function(x){
	ddply(x, .(new), colwise(sum))
	}
ans <- lapply(newnewdf, function.ans)

#this works to extract the first value of the TTT sequence (ie the correct startframe)
functiontry <- function(x){
	x[!duplicated(x[,4]),]
}
correctstarts<-lapply(newnewdf, functiontry)

#extract startframe as a vector
extractsf<-function(x){x[,2]}
sfvector <-lapply(correctstarts, extractsf)

newcolumn <- mapply(cbind, ans, "startframe"=sfvector, SIMPLIFY=FALSE)


#keep rows where newcolumn is not zero
subset.func <- function(x){
subset(x, x[,4]!=0)
}
newfwdframe <- lapply(newcolumn, subset.func)

#now getrid of rows with duration < 2
#select lengths with velocities above threshold
duration <- function(x) {
	subset(x, duration > durationlength)
	}

fwdstarts <- lapply(newfwdframe, duration)

#keep n(durations) and newsf(start frames)
extractor <- function(x) {
	x[,c(5,2)]
	
}

forwardstats <- lapply(fwdstarts, extractor)

sapply(forwardstats, write.table, file="forwardstats.csv", append=TRUE, sep=",", row.names=FALSE)

#calculate number of fwdinititations
numberfwdinitiations <- lapply(forwardstats, nrow)

#count mean duration of intitation per animal
fwddurationfunc <- function(x){
	mean(x[,2])
}

meanfwdduration <- lapply(fwdstarts, fwddurationfunc)

#count rows in velocitymatrix minus NA values
newfunction <- colSums(!is.na(avgnoNAvelocitymatrix))
newnew <- as.list(newfunction)

#calculate number of minutes by dividing rows by 600
numberminutes <- function(x){x/600}
minutes<- lapply(newnew, numberminutes)

numberofminutes <- as.numeric(minutes)

fwdinitiationsperminute <- as.numeric(numberfwdinitiations)/numberofminutes

fwdsummary <-cbind(numberfwdinitiations, fwdinitiationsperminute, meanfwdduration)

write.csv(fwdsummary, "fwdsummary.csv", row.names=FALSE)

###############################################################calculate backwards initiations 

rleback <- function(x) {
	rle(x<(-1))
}

b <-apply(avgnoNAvelocitymatrix,2,rleback)

newback <- function(x) {data.frame(number=x$values, duration=x$lengths, startframe=cumsum(x$lengths), row.names=NULL)}
newnewback <- lapply(b, newback)

#previous start frame is off (number is end of period not start)so first subset all of third column minus the last element
functionB <- function(x){head(x[,3], -1)} 
trynow <- lapply(newnewback, functionB)

#add 1 to the start of the vector
functiontry <- function(x){new <- c(1, x)}
trynewback <- lapply(trynow, functiontry)

#add the corrected vector to the data frame
trypleaseworkback <- mapply(cbind, newnewback, "newstart"=trynewback, SIMPLIFY=FALSE)

#remove column 3 <- function(x){x[,-3]}
backframe <- lapply(trypleaseworkback, remove3)

#this part filters out short non-fwd eg(if FWD=5frames, non-FWD=1-2frames, FWD=8 frames, will count as FWD=15)
tryfunction <- function(x) {
	transformer <- transform(x, newcolumn=ifelse((x[,1]	==FALSE & x[,2]<4), TRUE, x[,1]))
	}
pleaseworkback <- lapply(backframe,tryfunction)


#convert TRUE/FALSE into 1/0
integerize <- function(x){
	x*1}
nologicalback <- lapply(pleaseworkback,integerize)

#drop first column
columndrop <- function(x){(x[,-1])}
newnologicalback <- lapply(nologicalback, columndrop)

nnlback <- newnologicalback

#omitNAs
omit <- function(x){
	x[complete.cases(x),]	
}

omitNAsback <- lapply(nnlback,omit)

#rebuilddataframe
#rebuild <- function(x){
#	data.frame(x$n, x$startframe, x$newcolumn)
#}
#
#newdf <- lapply(omitNAs,rebuild)

#make an index column for 1s and 0s
func.rle <- function(x){
rle(x$newcolumn)
}

rle.idx.back <- lapply(omitNAsback, func.rle)

#new <-function(x){rep(seq(x$lengths), x$lengths)}
newnewback <- lapply(rle.idx.back, new)

#add index column to data frame
newnewdfback <- mapply(cbind, omitNAsback, "new"=newnewback, SIMPLIFY=FALSE)
library(plyr)

#add rows when index value is the same
function.ans <- function(x){
	ddply(x, .(new), colwise(sum))
	}
ansback <- lapply(newnewdfback, function.ans)

#this works to extract the first value of the TTT sequence (ie the correct startframe)
functiontry <- function(x){
	x[!duplicated(x[,4]),]
}
correctstartsback<-lapply(newnewdfback, functiontry)

#extract startframe as a vector
extractsf<-function(x){x[,2]}
sfvectorback <-lapply(correctstartsback, extractsf)

newcolumnback <- mapply(cbind, ansback, "startframe"=sfvectorback, SIMPLIFY=FALSE)


#keep rows where newcolumn is not zero
subset.func <- function(x){
subset(x, x[,4]!=0)
}
newback <- lapply(newcolumnback, subset.func)

#now getrid of rows with duration < 2
#select lengths with velocities above threshold
duration <- function(x) {
	subset(x, duration > durationlength)
	}

backstarts <- lapply(newback, duration)

#keep n(durations) and newsf(start frames)
extractor <- function(x) {
	x[,c(5,2)]
	
}

backstats <- lapply(backstarts, extractor)

sapply(backstats, write.table, file="backwardstats.csv", append=TRUE, sep=",", row.names=FALSE)

#calculate number of fwdinititations
numberbackinitiations <- lapply(backstats, nrow)

#count mean duration of intitation per animal
fwddurationfunc <- function(x){
	mean(x[,2])
}

meanbackduration <- lapply(backstarts, fwddurationfunc)

#count rows in velocitymatrix minus NA values
newfunction <- colSums(!is.na(avgnoNAvelocitymatrix))
newnew <- as.list(newfunction)

#calculate number of minutes by dividing rows by 600
numberminutes <- function(x){x/600}
minutes<- lapply(newnew, numberminutes)

numberofminutes <- as.numeric(minutes)

backinitiationsperminute <- as.numeric(numberbackinitiations)/numberofminutes

backsummary <-cbind(numberbackinitiations, backinitiationsperminute, meanbackduration)

write.csv(backsummary, "backsummary.csv")



#######################group data to make stacked bar graph
#add all fwd durations(number of frames)
totalframes <- sum(as.numeric(newnew))
fwdframes <- sum(avgnoNAvelocitymatrix > 1, na.rm=TRUE)
backframes <- sum(avgnoNAvelocitymatrix < -1, na.rm=TRUE)
pauseframes <- totalframes-fwdframes-backframes
timetable <- data.frame(forward=fwdframes/totalframes, backward=backframes/totalframes, pause=pauseframes/totalframes)
write.csv(timetable, "directionpercentage.csv")

##############################################################calculate all angles

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_2")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix2<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_2")[[i]]
	anglematrix2[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

#write.csv(anglematrix1, "anglematrix1.csv")

absangle2 <-apply(anglematrix2,2,abs)

fun <- function(x){
	if(is.numeric(x)) {mean(x, na.rm=T)}
	
}

angle2<-apply(absangle2,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_3")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix3<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_3")[[i]]
	anglematrix3[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle3 <-apply(anglematrix3,2,abs)

angle3<-apply(absangle3,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_4")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix4<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_4")[[i]]
	anglematrix4[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle4 <-apply(anglematrix4,2,abs)

angle4<-apply(absangle4,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_5")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix5<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_5")[[i]]
	anglematrix5[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle5 <-apply(anglematrix5,2,abs)

angle5<-apply(absangle5,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_6")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix6<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_6")[[i]]
	anglematrix6[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle6 <-apply(anglematrix6,2,abs)

angle6<-apply(absangle6,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_7")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix7<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_7")[[i]]
	anglematrix7[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle7 <-apply(anglematrix7,2,abs)

angle7<-apply(absangle7,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_8")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix8<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_8")[[i]]
	anglematrix8[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle8 <-apply(anglematrix8,2,abs)

angle8<-apply(absangle8,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_9")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix9<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_9")[[i]]
	anglematrix9[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle9 <-apply(anglematrix9,2,abs)

angle9<-apply(absangle9,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_10")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix10<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_10")[[i]]
	anglematrix10[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle10 <-apply(anglematrix10,2,abs)

angle10<-apply(absangle10,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_11")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix11<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_11")[[i]]
	anglematrix11[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle11 <-apply(anglematrix11,2,abs)

angle11<-apply(absangle11,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_12")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix12<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_12")[[i]]
	anglematrix12[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle12 <-apply(anglematrix12,2,abs)

angle12<-apply(absangle12,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_13")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix13<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_13")[[i]]
	anglematrix13[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle13 <-apply(anglematrix13,2,abs)

angle13<-apply(absangle13,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_14")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix14<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_14")[[i]]
	anglematrix14[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle14 <-apply(anglematrix14,2,abs)

angle14<-apply(absangle14,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_15")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix15<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_15")[[i]]
	anglematrix15[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle15 <-apply(anglematrix15,2,abs)

angle15<-apply(absangle15,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_16")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix16<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_16")[[i]]
	anglematrix16[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle16 <-apply(anglematrix16,2,abs)

angle16<-apply(absangle16,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_17")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix17<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_17")[[i]]
	anglematrix17[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle17 <-apply(anglematrix17,2,abs)

angle17<-apply(absangle17,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_18")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix18<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_18")[[i]]
	anglematrix18[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle18 <-apply(anglematrix18,2,abs)

angle18<-apply(absangle18,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_19")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix19<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_19")[[i]]
	anglematrix19[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle19 <-apply(anglematrix19,2,abs)

angle19<-apply(absangle19,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_20")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix20<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_20")[[i]]
	anglematrix20[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle20 <-apply(anglematrix20,2,abs)

angle20<-apply(absangle20,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_21")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix21<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_21")[[i]]
	anglematrix21[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle21 <-apply(anglematrix21,2,abs)

angle21<-apply(absangle21,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_22")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix22<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_22")[[i]]
	anglematrix22[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle22 <-apply(anglematrix22,2,abs)

angle22<-apply(absangle22,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_23")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix23<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_23")[[i]]
	anglematrix23[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle23 <-apply(anglematrix23,2,abs)

angle23<-apply(absangle23,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_24")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix24<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_24")[[i]]
	anglematrix24[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle24 <-apply(anglematrix24,2,abs)

angle24<-apply(absangle24,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_25")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix25<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_25")[[i]]
	anglematrix25[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle25 <-apply(anglematrix25,2,abs)

angle25<-apply(absangle25,2,fun)


# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_26")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix26<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_26")[[i]]
	anglematrix26[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle26 <-apply(anglematrix26,2,abs)

angle26<-apply(absangle26,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_27")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix27<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_27")[[i]]
	anglematrix27[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle27 <-apply(anglematrix27,2,abs)

angle27<-apply(absangle27,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_28")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix28<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_28")[[i]]
	anglematrix28[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle28 <-apply(anglematrix28,2,abs)

angle28<-apply(absangle28,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_29")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix29<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_29")[[i]]
	anglematrix29[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle29 <-apply(anglematrix29,2,abs)

angle29<-apply(absangle29,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_30")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix30<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_30")[[i]]
	anglematrix30[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle30 <-apply(anglematrix30,2,abs)

angle30<-apply(absangle30,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_31")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix31<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_31")[[i]]
	anglematrix31[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle31 <-apply(anglematrix31,2,abs)

angle31<-apply(absangle31,2,fun)

# in the case of acquisition duration is valiable...
nrowofeach<-NULL
for(i in 1:length(data_list))
{
	correctedmatrix<-extractdata(data_list, "rawdata", "angle_32")[[i]]
	nrowofeach<-cbind(nrowofeach ,length(correctedmatrix))
}
nrowofeach
maxrow<-max(nrowofeach)
maxrow

which(nrowofeach == maxrow)


anglematrix32<-matrix(NA, ncol=length(data_list),nrow= maxrow)
correctedanglematrix<-NULL
for(i in 1:length(data_list))
{

	correctedanglematrix<-extractdata(data_list, "rawdata","angle_32")[[i]]
	anglematrix32[1:length(correctedanglematrix),i] <- as.matrix(correctedanglematrix)
}

absangle32 <-apply(anglematrix32,2,abs)

angle32<-apply(absangle32,2,fun)

allangles <- data.frame(angle2, angle3, angle4, angle5, angle6, angle7, angle8, angle9, angle10, angle11, angle12, angle13, angle14, angle15, angle16, angle17, angle18, angle19, angle20, angle21, angle22, angle23, angle24, angle25, angle26, angle27, angle28, angle29, angle30, angle31, angle32)

write.csv(allangles, "allangles.csv")



##################################this part returns abs mean angle for forward and backward periods
#note that forward and backward are taken from original data where velocity>1 and velocity <-1. This will be slightly different than the frames calculated for forward/backward initiations/durations becuase those are from the runmean smoothing

list.df <- lapply(data_list, data.frame)

runmedCol <- function(x){runmed(x[c("rawdata.velocity")], k=3)}

subsetfwdangles <- function(x){
	subset(x, rawdata.velocity > 1, select=rawdata.angle_1:rawdata.angle_33)}
	
fwdangle <- lapply(list.df, subsetfwdangles)


absfwdangle <- lapply(fwdangle, abs)

meanabsfwdangle <- lapply(absfwdangle, colMeans)

fwdangles.df <- data.frame(t(sapply(meanabsfwdangle,c)))

subsetbackangles <- function(x){
	subset(x, rawdata.velocity < -1, select=rawdata.angle_1:rawdata.angle_33)}

backangle <- lapply(list.df, subsetbackangles)


absbackangle <- lapply(backangle, abs)

meanabsbackangle <- lapply(absbackangle, colMeans)

backangles.df <- data.frame(t(sapply(meanabsbackangle,c)))

write.csv(backangles.df, "backangles.csv")

write.csv(fwdangles.df, "forwardangles.csv")

################################replot curvature with increased range
##########################################################
# rainbow color

# 16x obj um/pix with  4x binning
#objscale<-0.39375*4
# 10x
#objscale<-0.63*4
#4x
objscale<-1.575*4

##########################################################

getInnerproduct<-function(ax,ay, bx,by)
{
	ax*bx+ay*by
}
#getInnerproduct(0,1,1,2)

getOuterproduct<-function(ax,ay, bx,by)
{
	ax*by-ay*bx
}
#getOuterproduct(0,1,1,2)


##########################################################

#get data files name
filenameslist<-system("ls | grep '.txt'",intern = TRUE)
#eliminate memo.txt
for(i in 1:length(filenameslist))
{
	if(filenameslist[i]=="memo.txt")
	{
		filenameslist[i]<-list(NULL)
	}
}
filenamesvec<-unlist(filenameslist)
filenamesvec
##########################################################
#put each data into list
# This part can be used for common application
data_list<-NULL
for(i in 1:length(filenamesvec))
{	
	
	data_list<-append(data_list, list(read.csv(filenamesvec[i])))
}
save(data_list, file=paste(format(Sys.time(), "%Y%m%d_%H_%M"),"data_list.Rdata",sep=""))



jointvecmatrixlist<-NULL

for(i in 1:length(data_list))
{
	workingdf<-data_list[[i]]
	header<-colnames(workingdf)
	
	#stage pos
	tempxpos<-workingdf$xpos
	tempypos<-workingdf$ypos
	
	#joint pos and direction
	jointsdf<-as.matrix(workingdf[,grep("joint",header)])
	#widthsdf<-as.matrix(workingdf[,grep("width",header)])
	jointsnum<-length(jointsdf[1,])/2

	jointx<-jointsdf[,(1:ncol(jointsdf))%%2-1==0]
	jointy<-jointsdf[,(1:ncol(jointsdf))%%2==0]

	#absolute position of joints
	absx<- jointx* objscale-tempxpos
	absy<- jointy* objscale-tempypos
	
	#motion x y
	motionx<-absx[2:nrow(absx),]-absx[1:(nrow(absx)-1),]
	motiony<-absy[2:nrow(absy),]-absy[1:(nrow(absy)-1),]
	#add NA at first frame
	motionx<-rbind(rep(NA,ncol(motionx)),motionx)
	motiony<-rbind(rep(NA,ncol(motiony)),motiony)
	
	#direction vector
	dirx<-absx[,3:ncol(absx)]-absx[,1:(ncol(absx)-2)]
	diry<-absy[,3:ncol(absy)]-absy[,1:(ncol(absy)-2)]
	
	#direction vec length
	dirlength<-sqrt(dirx^2+diry^2)
	
	#velocity
	jointvel<- -getInnerproduct(dirx/dirlength, diry/dirlength, motionx[,2:(ncol(motionx)-1)], motiony[,2:(ncol(motiony)-1)])
	
	#put NA at the frame can not extracted worm posture
	jointvel[is.na(workingdf$velocity),]<-NA
	jointvecmatrixlist <-append(jointvecmatrixlist,list(jointvel))

}




matrixlist<-NULL
for(i in 1:length(data_list))
{
	workingdf<-data_list[[i]]
	header<-colnames(workingdf)
	
	jointsdf<-as.matrix(workingdf[,grep("joint",header)])
	widthsdf<-as.matrix(workingdf[,grep("width",header)])
	jointsnum<-length(jointsdf[1,])/2
	anglemat<-matrix(NA,nrow=nrow(jointsdf),ncol=jointsnum-2)
	jointx<-jointsdf[,(1:ncol(jointsdf))%%2-1==0]
	jointy<-jointsdf[,(1:ncol(jointsdf))%%2==0]
	
	vecax<-jointx[,2:(ncol(jointx)-1)]-jointx[,1:(ncol(jointx)-2)]
	vecay<-jointy[,2:(ncol(jointx)-1)]-jointy[,1:(ncol(jointx)-2)]
	vecbx<-jointx[,3:(ncol(jointx))]-jointx[,2:(ncol(jointx)-1)]
	vecby<-jointy[,3:(ncol(jointx))]-jointy[,2:(ncol(jointx)-1)]
	outer<-getOuterproduct(vecax, vecay, vecbx, vecby)
	inner<-getInnerproduct(vecax, vecay, vecbx, vecby)
	anglemat<-atan2(outer, inner)
	#put NA at the frame can not extructed worm posture
	anglemat[is.na(workingdf$velocity),]<-NA
	

	matrixlist<-append(matrixlist,list(anglemat))

}

samplenames<-unlist(lapply(strsplit(filenamesvec,"\\."),function(x){x[1]}))
for(i in 1:length(matrixlist))
{
	tempmat<-matrixlist[[i]]
	quartz()
	par(mar=c(2,0,0,0))
	# This color code emphasize curvature but makes bit harder to see diff. between deeper and sharrower bend? 
	#use other color code depend on what you want to show
cols2 <- rgb(c(0:255,rep(255,255),255:0,rep(0,255*3))/255,c(rep(0,255),0:255,rep(255,255*2),255:0,rep(0,255))/255,c(rep(0,255*3),0:255,rep(255,255),255:0)/255)	

	image(1:nrow(tempmat),1:ncol(tempmat), tempmat,col=cols2,zlim=c(-0.9,0.9))
	dev.copy(png, file=paste(samplenames[i],"curvatureplot.png",sep=""), width=nrow(anglemat), height=300)
	dev.off()
}

##########################################################
#plot color bar for new curvature plots

quartz()
my.colors = colorRampPalette(c(cols2))
z=matrix(1:100,nrow=1)
x=1
y=seq(-0.9, 0.9, len=100) 
image(x,y,z,col=my.colors(100),axes=FALSE,xlab="",ylab="")
axis(4, at=c(-.8, -.4, 0, .4, .8))
dev.copy(png, file=paste("curvaturecolorbar.png", sep=""), width=100, height=400)