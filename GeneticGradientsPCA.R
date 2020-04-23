### This script computes the principal components (PCA) for a user-specified genetic dataset (SNPs or STR) in Arlequin format (.arp) and a file with the latitude and longitude coordinates for all the individuals of the genetic dataset. It also plot the first PC maps (by default it plosts the first three PCs but it could print the first 100 PCs if implemented)
### This script is distributed with copyright under the GNU GENERAL PUBLIC LICENSE. The GNU General Public License is a free, copyleft license for software and other kinds of works. A copy of the license is distributed with the script.
### Citation: Catarina Branco, Nicolas Ray, Mathias Currat, Miguel Arenas. 2020. Influence of Paleolithic Range Contraction, Admixture and Long-Distance Dispersal on Genetic Gradients of Modern Humans in Asia. Molecular Ecology.
### Miguel Arenas, 2020 (c).
###
### Information to run this script
### This script has been designed to analyze the genetic data of the study cited above (70 populations, 20 individuals per population and 500 SNPs per individual). If the genetic dataset is different, the script must be adapted.
### 1. This script runs on R, https://cran.r-project.org/
### 2. The libraries "maps" and "fields" must be installed
### 3. Introduce the path of the working directory (where the files with the genetic and geographic information are placed) in myplotsFolder
### 4. Introduce the number of genetic datsets that will be analyzed in nCoalSims
### 5. Introduce the name of genetic datsets that will be analyzed in OriginalNameInput. Excuple the .arp extension. Example, for 100 datasets (so nCoalSims <- 100) with names Data_1.arp to Data_100.arp, introduce here only Data_
### 6. The file with the coordinates for every individual of the sample (in the same order of the .arp file) must be named "coord.txt"
### 7. Ready. Under some operative systems some warnings may be displayed when running the script.
### 8. For every genetic dataset, the script prints two images with the first three PC gradients. For PC1, one of the images includes the connection between the positive and negative PC centroids. Black points show the sample locations. For all the genetic datasets an image is printed with a black line (connection between the positive and negative PC centroids) of every datset and a green line for the median (slope and ordinate in the origin) of the black lines.
### Check the provided illustrative example that includes input and output files.


library(maps)
library(fields)

### YOU MUST CHECK THESE THINGS BEFORE RUNNING
myplotsFolder <- paste("/Users/Miguel/Desktop/SoftwareGeneticGradientsAsia",sep="")       # path
nCoalSims <- 100
# number of genetic datasets(.arp files)
OriginalNameInput<-paste("GeneticData_ARP_",sep="") 				# basic name of inputs: Settings_GeneSamples_small_1.arp
setwd(myplotsFolder)


### Analyze all sims for this scenario
matrix1_results <- matrix(data=0, ncol=6, nrow=nCoalSims) # make a results matrix
matrix1_resultsSlIn <- matrix(data=0, ncol=6, nrow=nCoalSims) # make a results matrix
matrix1_resultsCoords <- matrix(data=0, ncol=12, nrow=nCoalSims) # make a results matrix
matrix2_C_1 <- matrix(data=0, ncol=2, nrow=48735) # make a results matrix
matrix2_PC_1 <- matrix(data=0, ncol=1, nrow=48735) # make a results matrix
matrix2_C_2 <- matrix(data=0, ncol=2, nrow=48735) # make a results matrix
matrix2_PC_2 <- matrix(data=0, ncol=1, nrow=48735) # make a results matrix
matrix2_C_3 <- matrix(data=0, ncol=2, nrow=48735) # make a results matrix
matrix2_PC_3 <- matrix(data=0, ncol=1, nrow=48735) # make a results matrix
matrix1_resultsSlIn_Cent <- matrix(data=0, ncol=6, nrow=nCoalSims) # make a results matrix
matrix1_resultsCoords_Cent <- matrix(data=0, ncol=12, nrow=nCoalSims) # make a results matrix

n <- 1
while (n <= nCoalSims) # for each genetic dataset
	{
	## processing ARP file function
	myprocess = function(
	input= string1, 
	Nbsample = 70,
	Nbloci = 500,
	Nbindiv = rep(20, 70))
		{
		X <- scan(file = input, what = character(), sep = "\t", quiet = TRUE, skip = 0, nlines = 0, comment.char = "#")
		Totindiv <- sum(Nbindiv)
		lx <- length(X)
		X <- X[30:lx]
		M <- matrix(NA, nrow = Totindiv, ncol = (2 + Nbloci)) 
		for(k in 1:Nbsample)
			{
			shift <- 1
			if (k == 1) k.s <- 0 else k.s <- sum(Nbindiv[1:(k-1)]) 
			for (i in 1:(Nbindiv[k]))
				{
				X[shift + (i-1)*(2+Nbloci)] <- k
				M[k.s + i, ] <- as.numeric( X[ ( shift + (i-1)*(2+Nbloci) ):(shift -1 +i*(2+Nbloci) )] )
				shift <- shift + 1
				}
			lx <- length(X)
			X<-X[(10 + shift+i*(2+Nbloci)):lx]
			}
		return(M[,-2])
		}

	colrp <- colorRampPalette(c("yellow","orange","red","brown"))(20)
    colrp_white <- colorRampPalette(c("white","white","white","white"))(20)

	# Name of file
	NameInput<-paste(OriginalNameInput,n,".arp",sep="")
	string1 <- paste("",NameInput,sep="")
	coord = read.table("coord.txt")

	# get the genotypes from the ARP file 
	genotype <- myprocess(Nbsample= 70, Nbloci = 500, Nbindiv = rep(20, 70))
	genotype <- genotype[,-1]

	# compute PCA
	objInit <- prcomp(genotype, scale = T)

	
	#### Cumulatives ####
	## PCA components
	# $sdev	
	# $rotation
	# $center
	# $scale
	# $x

	if (n == 1) # first replicate
		{
		cumObj <- objInit
		print("> Working on file:")
		print(string1)

		# print(objInit$sdev)
		# print(cumObj$sdev)
		
        TotalCumValue <- vector()
        
        for (b in 1:3)
        {
            if (b == 1) # PC1
            {
                
                for (kkk in 1:70)
                {
                    TotalCumValue[kkk] <- 0
                }
                
                for (a in 1:1400)
                {
                    if (a >= 1 && a < 21)
                    {
                        TotalCumValue[1] <- TotalCumValue[1] + objInit$x[a,b]
                    }
                    if (a >= 21 && a < 41)
                    {
                        TotalCumValue[2] <- TotalCumValue[2] + objInit$x[a,b]
                    }
                    if (a >= 41 && a < 61)
                    {
                        TotalCumValue[3] <- TotalCumValue[3] + objInit$x[a,b]
                    }
                    if (a >= 61 && a < 81)
                    {
                        TotalCumValue[4] <- TotalCumValue[4] + objInit$x[a,b]
                    }
                    if (a >= 81 && a < 101)
                    {
                        TotalCumValue[5] <- TotalCumValue[5] + objInit$x[a,b]
                    }
                    if (a >= 101 && a < 121)
                    {
                        TotalCumValue[6] <- TotalCumValue[6] + objInit$x[a,b]
                    }
                    if (a >= 121 && a < 141)
                    {
                        TotalCumValue[7] <- TotalCumValue[7] + objInit$x[a,b]
                    }
                    if (a >= 141 && a < 161)
                    {
                        TotalCumValue[8] <- TotalCumValue[8] + objInit$x[a,b]
                    }
                    if (a >= 161 && a < 181)
                    {
                        TotalCumValue[9] <- TotalCumValue[9] + objInit$x[a,b]
                    }
                    if (a >= 181 && a < 201)
                    {
                        TotalCumValue[10] <- TotalCumValue[10] + objInit$x[a,b]
                    }
                    if (a >= 201 && a < 221)
                    {
                        TotalCumValue[11] <- TotalCumValue[11] + objInit$x[a,b]
                    }
                    if (a >= 221 && a < 241)
                    {
                        TotalCumValue[12] <- TotalCumValue[12] + objInit$x[a,b]
                    }
                    if (a >= 241 && a < 261)
                    {
                        TotalCumValue[13] <- TotalCumValue[13] + objInit$x[a,b]
                    }
                    if (a >= 261 && a < 281)
                    {
                        TotalCumValue[14] <- TotalCumValue[14] + objInit$x[a,b]
                    }
                    if (a >= 281 && a < 301)
                    {
                        TotalCumValue[15] <- TotalCumValue[15] + objInit$x[a,b]
                    }
                    if (a >= 301 && a < 321)
                    {
                        TotalCumValue[16] <- TotalCumValue[16] + objInit$x[a,b]
                    }
                    if (a >= 321 && a < 341)
                    {
                        TotalCumValue[17] <- TotalCumValue[17] + objInit$x[a,b]
                    }
                    if (a >= 341 && a < 361)
                    {
                        TotalCumValue[18] <- TotalCumValue[18] + objInit$x[a,b]
                    }
                    if (a >= 361 && a < 381)
                    {
                        TotalCumValue[19] <- TotalCumValue[19] + objInit$x[a,b]
                    }
                    if (a >= 381 && a < 401)
                    {
                        TotalCumValue[20] <- TotalCumValue[20] + objInit$x[a,b]
                    }
                    if (a >= 401 && a < 421)
                    {
                        TotalCumValue[21] <- TotalCumValue[21] + objInit$x[a,b]
                    }
                    if (a >= 421 && a < 441)
                    {
                        TotalCumValue[22] <- TotalCumValue[22] + objInit$x[a,b]
                    }
                    if (a >= 441 && a < 461)
                    {
                        TotalCumValue[23] <- TotalCumValue[23] + objInit$x[a,b]
                    }
                    if (a >= 461 && a < 481)
                    {
                        TotalCumValue[24] <- TotalCumValue[24] + objInit$x[a,b]
                    }
                    if (a >= 481 && a < 501)
                    {
                        TotalCumValue[25] <- TotalCumValue[25] + objInit$x[a,b]
                    }
                    if (a >= 501 && a < 521)
                    {
                        TotalCumValue[26] <- TotalCumValue[26] + objInit$x[a,b]
                    }
                    if (a >= 521 && a < 541)
                    {
                        TotalCumValue[27] <- TotalCumValue[27] + objInit$x[a,b]
                    }
                    if (a >= 541 && a < 561)
                    {
                        TotalCumValue[28] <- TotalCumValue[28] + objInit$x[a,b]
                    }
                    if (a >= 561 && a < 581)
                    {
                        TotalCumValue[29] <- TotalCumValue[29] + objInit$x[a,b]
                    }
                    if (a >= 581 && a < 601)
                    {
                        TotalCumValue[30] <- TotalCumValue[30] + objInit$x[a,b]
                    }
                    if (a >= 601 && a < 621)
                    {
                        TotalCumValue[31] <- TotalCumValue[31] + objInit$x[a,b]
                    }
                    if (a >= 621 && a < 641)
                    {
                        TotalCumValue[32] <- TotalCumValue[32] + objInit$x[a,b]
                    }
                    if (a >= 641 && a < 661)
                    {
                        TotalCumValue[33] <- TotalCumValue[33] + objInit$x[a,b]
                    }
                    if (a >= 661 && a < 681)
                    {
                        TotalCumValue[34] <- TotalCumValue[34] + objInit$x[a,b]
                    }
                    if (a >= 681 && a < 701)
                    {
                        TotalCumValue[35] <- TotalCumValue[35] + objInit$x[a,b]
                    }
                    if (a >= 701 && a < 721)
                    {
                        TotalCumValue[36] <- TotalCumValue[36] + objInit$x[a,b]
                    }
                    if (a >= 721 && a < 741)
                    {
                        TotalCumValue[37] <- TotalCumValue[37] + objInit$x[a,b]
                    }
                    if (a >= 741 && a < 761)
                    {
                        TotalCumValue[38] <- TotalCumValue[38] + objInit$x[a,b]
                    }
                    if (a >= 761 && a < 781)
                    {
                        TotalCumValue[39] <- TotalCumValue[39] + objInit$x[a,b]
                    }
                    if (a >= 781 && a < 801)
                    {
                        TotalCumValue[40] <- TotalCumValue[40] + objInit$x[a,b]
                    }
                    if (a >= 801 && a < 821)
                    {
                        TotalCumValue[41] <- TotalCumValue[41] + objInit$x[a,b]
                    }
                    if (a >= 821 && a < 841)
                    {
                        TotalCumValue[42] <- TotalCumValue[42] + objInit$x[a,b]
                    }
                    if (a >= 841 && a < 861)
                    {
                        TotalCumValue[43] <- TotalCumValue[43] + objInit$x[a,b]
                    }
                    if (a >= 861 && a < 881)
                    {
                        TotalCumValue[44] <- TotalCumValue[44] + objInit$x[a,b]
                    }
                    if (a >= 881 && a < 901)
                    {
                        TotalCumValue[45] <- TotalCumValue[45] + objInit$x[a,b]
                    }
                    if (a >= 901 && a < 921)
                    {
                        TotalCumValue[46] <- TotalCumValue[46] + objInit$x[a,b]
                    }
                    if (a >= 921 && a < 941)
                    {
                        TotalCumValue[47] <- TotalCumValue[47] + objInit$x[a,b]
                    }
                    if (a >= 941 && a < 961)
                    {
                        TotalCumValue[48] <- TotalCumValue[48] + objInit$x[a,b]
                    }
                    if (a >= 961 && a < 981)
                    {
                        TotalCumValue[49] <- TotalCumValue[49] + objInit$x[a,b]
                    }
                    if (a >= 981 && a < 1001)
                    {
                        TotalCumValue[50] <- TotalCumValue[50] + objInit$x[a,b]
                    }
                    if (a >= 1001 && a < 1021)
                    {
                        TotalCumValue[51] <- TotalCumValue[51] + objInit$x[a,b]
                    }
                    if (a >= 1021 && a < 1041)
                    {
                        TotalCumValue[52] <- TotalCumValue[52] + objInit$x[a,b]
                    }
                    if (a >= 1041 && a < 1061)
                    {
                        TotalCumValue[53] <- TotalCumValue[53] + objInit$x[a,b]
                    }
                    if (a >= 1061 && a < 1081)
                    {
                        TotalCumValue[54] <- TotalCumValue[54] + objInit$x[a,b]
                    }
                    if (a >= 1081 && a < 1101)
                    {
                        TotalCumValue[55] <- TotalCumValue[55] + objInit$x[a,b]
                    }
                    if (a >= 1101 && a < 1121)
                    {
                        TotalCumValue[56] <- TotalCumValue[56] + objInit$x[a,b]
                    }
                    if (a >= 1121 && a < 1141)
                    {
                        TotalCumValue[57] <- TotalCumValue[57] + objInit$x[a,b]
                    }
                    if (a >= 1141 && a < 1161)
                    {
                        TotalCumValue[58] <- TotalCumValue[58] + objInit$x[a,b]
                    }
                    if (a >= 1161 && a < 1181)
                    {
                        TotalCumValue[59] <- TotalCumValue[59] + objInit$x[a,b]
                    }
                    if (a >= 1181 && a < 1201)
                    {
                        TotalCumValue[60] <- TotalCumValue[60] + objInit$x[a,b]
                    }
                    if (a >= 1201 && a < 1221)
                    {
                        TotalCumValue[61] <- TotalCumValue[61] + objInit$x[a,b]
                    }
                    if (a >= 1221 && a < 1241)
                    {
                        TotalCumValue[62] <- TotalCumValue[62] + objInit$x[a,b]
                    }
                    if (a >= 1241 && a < 1261)
                    {
                        TotalCumValue[63] <- TotalCumValue[63] + objInit$x[a,b]
                    }
                    if (a >= 1261 && a < 1281)
                    {
                        TotalCumValue[64] <- TotalCumValue[64] + objInit$x[a,b]
                    }
                    if (a >= 1281 && a < 1301)
                    {
                        TotalCumValue[65] <- TotalCumValue[65] + objInit$x[a,b]
                    }
                    if (a >= 1301 && a < 1321)
                    {
                        TotalCumValue[66] <- TotalCumValue[66] + objInit$x[a,b]
                    }
                    if (a >= 1321 && a < 1341)
                    {
                        TotalCumValue[67] <- TotalCumValue[67] + objInit$x[a,b]
                    }
                    if (a >= 1341 && a < 1361)
                    {
                        TotalCumValue[68] <- TotalCumValue[68] + objInit$x[a,b]
                    }
                    if (a >= 1361 && a < 1381)
                    {
                        TotalCumValue[69] <- TotalCumValue[69] + objInit$x[a,b]
                    }
                    if (a >= 1381 && a < 1401)
                    {
                        TotalCumValue[70] <- TotalCumValue[70] + objInit$x[a,b]
                    }                    
                } # end for a
                
                minValuePC1 <- min(TotalCumValue)
                maxValuePC1 <- max(TotalCumValue)
                xxmin <- -1000
                xxmax <- -1000
                
                for (kkk in 1:70)
                {
                    if (minValuePC1 == TotalCumValue[kkk])
                    {
                        xxmin <- kkk
                    }
                    if (maxValuePC1 == TotalCumValue[kkk])
                    {
                        xxmax <- kkk
                    }
                }
                
                minLinePC1 <- (xxmin * 20) - 20 + 1
                maxLinePC1 <- (xxmax * 20) - 20 + 1
                
                #kkk <- 20 # highest/lowest 20 points
                #ndx <- order(objInit$x[,b], decreasing = T)[1:kkk]
                #maxValuePC1_vector <- objInit$x[,b][ndx]
                #ndx <- order(objInit$x[,b])[1:kkk]
                #minValuePC1_vector <- objInit$x[,b][ndx]
            } # end b == 1
            
            if (b == 2) # PC2
            {
                for (kkk in 1:70)
                {
                    TotalCumValue[kkk] <- 0
                }
                
                for (a in 1:1400)
                {
                    if (a >= 1 && a < 21)
                    {
                        TotalCumValue[1] <- TotalCumValue[1] + objInit$x[a,b]
                    }
                    if (a >= 21 && a < 41)
                    {
                        TotalCumValue[2] <- TotalCumValue[2] + objInit$x[a,b]
                    }
                    if (a >= 41 && a < 61)
                    {
                        TotalCumValue[3] <- TotalCumValue[3] + objInit$x[a,b]
                    }
                    if (a >= 61 && a < 81)
                    {
                        TotalCumValue[4] <- TotalCumValue[4] + objInit$x[a,b]
                    }
                    if (a >= 81 && a < 101)
                    {
                        TotalCumValue[5] <- TotalCumValue[5] + objInit$x[a,b]
                    }
                    if (a >= 101 && a < 121)
                    {
                        TotalCumValue[6] <- TotalCumValue[6] + objInit$x[a,b]
                    }
                    if (a >= 121 && a < 141)
                    {
                        TotalCumValue[7] <- TotalCumValue[7] + objInit$x[a,b]
                    }
                    if (a >= 141 && a < 161)
                    {
                        TotalCumValue[8] <- TotalCumValue[8] + objInit$x[a,b]
                    }
                    if (a >= 161 && a < 181)
                    {
                        TotalCumValue[9] <- TotalCumValue[9] + objInit$x[a,b]
                    }
                    if (a >= 181 && a < 201)
                    {
                        TotalCumValue[10] <- TotalCumValue[10] + objInit$x[a,b]
                    }
                    if (a >= 201 && a < 221)
                    {
                        TotalCumValue[11] <- TotalCumValue[11] + objInit$x[a,b]
                    }
                    if (a >= 221 && a < 241)
                    {
                        TotalCumValue[12] <- TotalCumValue[12] + objInit$x[a,b]
                    }
                    if (a >= 241 && a < 261)
                    {
                        TotalCumValue[13] <- TotalCumValue[13] + objInit$x[a,b]
                    }
                    if (a >= 261 && a < 281)
                    {
                        TotalCumValue[14] <- TotalCumValue[14] + objInit$x[a,b]
                    }
                    if (a >= 281 && a < 301)
                    {
                        TotalCumValue[15] <- TotalCumValue[15] + objInit$x[a,b]
                    }
                    if (a >= 301 && a < 321)
                    {
                        TotalCumValue[16] <- TotalCumValue[16] + objInit$x[a,b]
                    }
                    if (a >= 321 && a < 341)
                    {
                        TotalCumValue[17] <- TotalCumValue[17] + objInit$x[a,b]
                    }
                    if (a >= 341 && a < 361)
                    {
                        TotalCumValue[18] <- TotalCumValue[18] + objInit$x[a,b]
                    }
                    if (a >= 361 && a < 381)
                    {
                        TotalCumValue[19] <- TotalCumValue[19] + objInit$x[a,b]
                    }
                    if (a >= 381 && a < 401)
                    {
                        TotalCumValue[20] <- TotalCumValue[20] + objInit$x[a,b]
                    }
                    if (a >= 401 && a < 421)
                    {
                        TotalCumValue[21] <- TotalCumValue[21] + objInit$x[a,b]
                    }
                    if (a >= 421 && a < 441)
                    {
                        TotalCumValue[22] <- TotalCumValue[22] + objInit$x[a,b]
                    }
                    if (a >= 441 && a < 461)
                    {
                        TotalCumValue[23] <- TotalCumValue[23] + objInit$x[a,b]
                    }
                    if (a >= 461 && a < 481)
                    {
                        TotalCumValue[24] <- TotalCumValue[24] + objInit$x[a,b]
                    }
                    if (a >= 481 && a < 501)
                    {
                        TotalCumValue[25] <- TotalCumValue[25] + objInit$x[a,b]
                    }
                    if (a >= 501 && a < 521)
                    {
                        TotalCumValue[26] <- TotalCumValue[26] + objInit$x[a,b]
                    }
                    if (a >= 521 && a < 541)
                    {
                        TotalCumValue[27] <- TotalCumValue[27] + objInit$x[a,b]
                    }
                    if (a >= 541 && a < 561)
                    {
                        TotalCumValue[28] <- TotalCumValue[28] + objInit$x[a,b]
                    }
                    if (a >= 561 && a < 581)
                    {
                        TotalCumValue[29] <- TotalCumValue[29] + objInit$x[a,b]
                    }
                    if (a >= 581 && a < 601)
                    {
                        TotalCumValue[30] <- TotalCumValue[30] + objInit$x[a,b]
                    }
                    if (a >= 601 && a < 621)
                    {
                        TotalCumValue[31] <- TotalCumValue[31] + objInit$x[a,b]
                    }
                    if (a >= 621 && a < 641)
                    {
                        TotalCumValue[32] <- TotalCumValue[32] + objInit$x[a,b]
                    }
                    if (a >= 641 && a < 661)
                    {
                        TotalCumValue[33] <- TotalCumValue[33] + objInit$x[a,b]
                    }
                    if (a >= 661 && a < 681)
                    {
                        TotalCumValue[34] <- TotalCumValue[34] + objInit$x[a,b]
                    }
                    if (a >= 681 && a < 701)
                    {
                        TotalCumValue[35] <- TotalCumValue[35] + objInit$x[a,b]
                    }
                    if (a >= 701 && a < 721)
                    {
                        TotalCumValue[36] <- TotalCumValue[36] + objInit$x[a,b]
                    }
                    if (a >= 721 && a < 741)
                    {
                        TotalCumValue[37] <- TotalCumValue[37] + objInit$x[a,b]
                    }
                    if (a >= 741 && a < 761)
                    {
                        TotalCumValue[38] <- TotalCumValue[38] + objInit$x[a,b]
                    }
                    if (a >= 761 && a < 781)
                    {
                        TotalCumValue[39] <- TotalCumValue[39] + objInit$x[a,b]
                    }
                    if (a >= 781 && a < 801)
                    {
                        TotalCumValue[40] <- TotalCumValue[40] + objInit$x[a,b]
                    }
                    if (a >= 801 && a < 821)
                    {
                        TotalCumValue[41] <- TotalCumValue[41] + objInit$x[a,b]
                    }
                    if (a >= 821 && a < 841)
                    {
                        TotalCumValue[42] <- TotalCumValue[42] + objInit$x[a,b]
                    }
                    if (a >= 841 && a < 861)
                    {
                        TotalCumValue[43] <- TotalCumValue[43] + objInit$x[a,b]
                    }
                    if (a >= 861 && a < 881)
                    {
                        TotalCumValue[44] <- TotalCumValue[44] + objInit$x[a,b]
                    }
                    if (a >= 881 && a < 901)
                    {
                        TotalCumValue[45] <- TotalCumValue[45] + objInit$x[a,b]
                    }
                    if (a >= 901 && a < 921)
                    {
                        TotalCumValue[46] <- TotalCumValue[46] + objInit$x[a,b]
                    }
                    if (a >= 921 && a < 941)
                    {
                        TotalCumValue[47] <- TotalCumValue[47] + objInit$x[a,b]
                    }
                    if (a >= 941 && a < 961)
                    {
                        TotalCumValue[48] <- TotalCumValue[48] + objInit$x[a,b]
                    }
                    if (a >= 961 && a < 981)
                    {
                        TotalCumValue[49] <- TotalCumValue[49] + objInit$x[a,b]
                    }
                    if (a >= 981 && a < 1001)
                    {
                        TotalCumValue[50] <- TotalCumValue[50] + objInit$x[a,b]
                    }
                    if (a >= 1001 && a < 1021)
                    {
                        TotalCumValue[51] <- TotalCumValue[51] + objInit$x[a,b]
                    }
                    if (a >= 1021 && a < 1041)
                    {
                        TotalCumValue[52] <- TotalCumValue[52] + objInit$x[a,b]
                    }
                    if (a >= 1041 && a < 1061)
                    {
                        TotalCumValue[53] <- TotalCumValue[53] + objInit$x[a,b]
                    }
                    if (a >= 1061 && a < 1081)
                    {
                        TotalCumValue[54] <- TotalCumValue[54] + objInit$x[a,b]
                    }
                    if (a >= 1081 && a < 1101)
                    {
                        TotalCumValue[55] <- TotalCumValue[55] + objInit$x[a,b]
                    }
                    if (a >= 1101 && a < 1121)
                    {
                        TotalCumValue[56] <- TotalCumValue[56] + objInit$x[a,b]
                    }
                    if (a >= 1121 && a < 1141)
                    {
                        TotalCumValue[57] <- TotalCumValue[57] + objInit$x[a,b]
                    }
                    if (a >= 1141 && a < 1161)
                    {
                        TotalCumValue[58] <- TotalCumValue[58] + objInit$x[a,b]
                    }
                    if (a >= 1161 && a < 1181)
                    {
                        TotalCumValue[59] <- TotalCumValue[59] + objInit$x[a,b]
                    }
                    if (a >= 1181 && a < 1201)
                    {
                        TotalCumValue[60] <- TotalCumValue[60] + objInit$x[a,b]
                    }
                    if (a >= 1201 && a < 1221)
                    {
                        TotalCumValue[61] <- TotalCumValue[61] + objInit$x[a,b]
                    }
                    if (a >= 1221 && a < 1241)
                    {
                        TotalCumValue[62] <- TotalCumValue[62] + objInit$x[a,b]
                    }
                    if (a >= 1241 && a < 1261)
                    {
                        TotalCumValue[63] <- TotalCumValue[63] + objInit$x[a,b]
                    }
                    if (a >= 1261 && a < 1281)
                    {
                        TotalCumValue[64] <- TotalCumValue[64] + objInit$x[a,b]
                    }
                    if (a >= 1281 && a < 1301)
                    {
                        TotalCumValue[65] <- TotalCumValue[65] + objInit$x[a,b]
                    }
                    if (a >= 1301 && a < 1321)
                    {
                        TotalCumValue[66] <- TotalCumValue[66] + objInit$x[a,b]
                    }
                    if (a >= 1321 && a < 1341)
                    {
                        TotalCumValue[67] <- TotalCumValue[67] + objInit$x[a,b]
                    }
                    if (a >= 1341 && a < 1361)
                    {
                        TotalCumValue[68] <- TotalCumValue[68] + objInit$x[a,b]
                    }
                    if (a >= 1361 && a < 1381)
                    {
                        TotalCumValue[69] <- TotalCumValue[69] + objInit$x[a,b]
                    }
                    if (a >= 1381 && a < 1401)
                    {
                        TotalCumValue[70] <- TotalCumValue[70] + objInit$x[a,b]
                    }                          
                   
                } # end for a
                
                minValuePC2 <- min(TotalCumValue)
                maxValuePC2 <- max(TotalCumValue)
                xxmin <- -1000
                xxmax <- -1000
                
                for (kkk in 1:70)
                {
                    if (minValuePC2 == TotalCumValue[kkk])
                    {
                        xxmin <- kkk
                    }
                    if (maxValuePC2 == TotalCumValue[kkk])
                    {
                        xxmax <- kkk
                    }
                }
                
                minLinePC2 <- (xxmin * 20) - 20 + 1
                maxLinePC2 <- (xxmax * 20) - 20 + 1
                
                #kkk <- 20 # highest/lowest 20 points
                #ndx <- order(objInit$x[,b], decreasing = T)[1:kkk]
                #maxValuePC2_vector <- objInit$x[,b][ndx]
                #ndx <- order(objInit$x[,b])[1:kkk]
                #minValuePC2_vector <- objInit$x[,b][ndx]
            } # end b == 2
            
            if (b == 3) # PC3
            {
                for (kkk in 1:70)
                {
                    TotalCumValue[kkk] <- 0
                }
                
                for (a in 1:1400)
                {
                    if (a >= 1 && a < 21)
                    {
                        TotalCumValue[1] <- TotalCumValue[1] + objInit$x[a,b]
                    }
                    if (a >= 21 && a < 41)
                    {
                        TotalCumValue[2] <- TotalCumValue[2] + objInit$x[a,b]
                    }
                    if (a >= 41 && a < 61)
                    {
                        TotalCumValue[3] <- TotalCumValue[3] + objInit$x[a,b]
                    }
                    if (a >= 61 && a < 81)
                    {
                        TotalCumValue[4] <- TotalCumValue[4] + objInit$x[a,b]
                    }
                    if (a >= 81 && a < 101)
                    {
                        TotalCumValue[5] <- TotalCumValue[5] + objInit$x[a,b]
                    }
                    if (a >= 101 && a < 121)
                    {
                        TotalCumValue[6] <- TotalCumValue[6] + objInit$x[a,b]
                    }
                    if (a >= 121 && a < 141)
                    {
                        TotalCumValue[7] <- TotalCumValue[7] + objInit$x[a,b]
                    }
                    if (a >= 141 && a < 161)
                    {
                        TotalCumValue[8] <- TotalCumValue[8] + objInit$x[a,b]
                    }
                    if (a >= 161 && a < 181)
                    {
                        TotalCumValue[9] <- TotalCumValue[9] + objInit$x[a,b]
                    }
                    if (a >= 181 && a < 201)
                    {
                        TotalCumValue[10] <- TotalCumValue[10] + objInit$x[a,b]
                    }
                    if (a >= 201 && a < 221)
                    {
                        TotalCumValue[11] <- TotalCumValue[11] + objInit$x[a,b]
                    }
                    if (a >= 221 && a < 241)
                    {
                        TotalCumValue[12] <- TotalCumValue[12] + objInit$x[a,b]
                    }
                    if (a >= 241 && a < 261)
                    {
                        TotalCumValue[13] <- TotalCumValue[13] + objInit$x[a,b]
                    }
                    if (a >= 261 && a < 281)
                    {
                        TotalCumValue[14] <- TotalCumValue[14] + objInit$x[a,b]
                    }
                    if (a >= 281 && a < 301)
                    {
                        TotalCumValue[15] <- TotalCumValue[15] + objInit$x[a,b]
                    }
                    if (a >= 301 && a < 321)
                    {
                        TotalCumValue[16] <- TotalCumValue[16] + objInit$x[a,b]
                    }
                    if (a >= 321 && a < 341)
                    {
                        TotalCumValue[17] <- TotalCumValue[17] + objInit$x[a,b]
                    }
                    if (a >= 341 && a < 361)
                    {
                        TotalCumValue[18] <- TotalCumValue[18] + objInit$x[a,b]
                    }
                    if (a >= 361 && a < 381)
                    {
                        TotalCumValue[19] <- TotalCumValue[19] + objInit$x[a,b]
                    }
                    if (a >= 381 && a < 401)
                    {
                        TotalCumValue[20] <- TotalCumValue[20] + objInit$x[a,b]
                    }
                    if (a >= 401 && a < 421)
                    {
                        TotalCumValue[21] <- TotalCumValue[21] + objInit$x[a,b]
                    }
                    if (a >= 421 && a < 441)
                    {
                        TotalCumValue[22] <- TotalCumValue[22] + objInit$x[a,b]
                    }
                    if (a >= 441 && a < 461)
                    {
                        TotalCumValue[23] <- TotalCumValue[23] + objInit$x[a,b]
                    }
                    if (a >= 461 && a < 481)
                    {
                        TotalCumValue[24] <- TotalCumValue[24] + objInit$x[a,b]
                    }
                    if (a >= 481 && a < 501)
                    {
                        TotalCumValue[25] <- TotalCumValue[25] + objInit$x[a,b]
                    }
                    if (a >= 501 && a < 521)
                    {
                        TotalCumValue[26] <- TotalCumValue[26] + objInit$x[a,b]
                    }
                    if (a >= 521 && a < 541)
                    {
                        TotalCumValue[27] <- TotalCumValue[27] + objInit$x[a,b]
                    }
                    if (a >= 541 && a < 561)
                    {
                        TotalCumValue[28] <- TotalCumValue[28] + objInit$x[a,b]
                    }
                    if (a >= 561 && a < 581)
                    {
                        TotalCumValue[29] <- TotalCumValue[29] + objInit$x[a,b]
                    }
                    if (a >= 581 && a < 601)
                    {
                        TotalCumValue[30] <- TotalCumValue[30] + objInit$x[a,b]
                    }
                    if (a >= 601 && a < 621)
                    {
                        TotalCumValue[31] <- TotalCumValue[31] + objInit$x[a,b]
                    }
                    if (a >= 621 && a < 641)
                    {
                        TotalCumValue[32] <- TotalCumValue[32] + objInit$x[a,b]
                    }
                    if (a >= 641 && a < 661)
                    {
                        TotalCumValue[33] <- TotalCumValue[33] + objInit$x[a,b]
                    }
                    if (a >= 661 && a < 681)
                    {
                        TotalCumValue[34] <- TotalCumValue[34] + objInit$x[a,b]
                    }
                    if (a >= 681 && a < 701)
                    {
                        TotalCumValue[35] <- TotalCumValue[35] + objInit$x[a,b]
                    }
                    if (a >= 701 && a < 721)
                    {
                        TotalCumValue[36] <- TotalCumValue[36] + objInit$x[a,b]
                    }
                    if (a >= 721 && a < 741)
                    {
                        TotalCumValue[37] <- TotalCumValue[37] + objInit$x[a,b]
                    }
                    if (a >= 741 && a < 761)
                    {
                        TotalCumValue[38] <- TotalCumValue[38] + objInit$x[a,b]
                    }
                    if (a >= 761 && a < 781)
                    {
                        TotalCumValue[39] <- TotalCumValue[39] + objInit$x[a,b]
                    }
                    if (a >= 781 && a < 801)
                    {
                        TotalCumValue[40] <- TotalCumValue[40] + objInit$x[a,b]
                    }
                    if (a >= 801 && a < 821)
                    {
                        TotalCumValue[41] <- TotalCumValue[41] + objInit$x[a,b]
                    }
                    if (a >= 821 && a < 841)
                    {
                        TotalCumValue[42] <- TotalCumValue[42] + objInit$x[a,b]
                    }
                    if (a >= 841 && a < 861)
                    {
                        TotalCumValue[43] <- TotalCumValue[43] + objInit$x[a,b]
                    }
                    if (a >= 861 && a < 881)
                    {
                        TotalCumValue[44] <- TotalCumValue[44] + objInit$x[a,b]
                    }
                    if (a >= 881 && a < 901)
                    {
                        TotalCumValue[45] <- TotalCumValue[45] + objInit$x[a,b]
                    }
                    if (a >= 901 && a < 921)
                    {
                        TotalCumValue[46] <- TotalCumValue[46] + objInit$x[a,b]
                    }
                    if (a >= 921 && a < 941)
                    {
                        TotalCumValue[47] <- TotalCumValue[47] + objInit$x[a,b]
                    }
                    if (a >= 941 && a < 961)
                    {
                        TotalCumValue[48] <- TotalCumValue[48] + objInit$x[a,b]
                    }
                    if (a >= 961 && a < 981)
                    {
                        TotalCumValue[49] <- TotalCumValue[49] + objInit$x[a,b]
                    }
                    if (a >= 981 && a < 1001)
                    {
                        TotalCumValue[50] <- TotalCumValue[50] + objInit$x[a,b]
                    }
                    if (a >= 1001 && a < 1021)
                    {
                        TotalCumValue[51] <- TotalCumValue[51] + objInit$x[a,b]
                    }
                    if (a >= 1021 && a < 1041)
                    {
                        TotalCumValue[52] <- TotalCumValue[52] + objInit$x[a,b]
                    }
                    if (a >= 1041 && a < 1061)
                    {
                        TotalCumValue[53] <- TotalCumValue[53] + objInit$x[a,b]
                    }
                    if (a >= 1061 && a < 1081)
                    {
                        TotalCumValue[54] <- TotalCumValue[54] + objInit$x[a,b]
                    }
                    if (a >= 1081 && a < 1101)
                    {
                        TotalCumValue[55] <- TotalCumValue[55] + objInit$x[a,b]
                    }
                    if (a >= 1101 && a < 1121)
                    {
                        TotalCumValue[56] <- TotalCumValue[56] + objInit$x[a,b]
                    }
                    if (a >= 1121 && a < 1141)
                    {
                        TotalCumValue[57] <- TotalCumValue[57] + objInit$x[a,b]
                    }
                    if (a >= 1141 && a < 1161)
                    {
                        TotalCumValue[58] <- TotalCumValue[58] + objInit$x[a,b]
                    }
                    if (a >= 1161 && a < 1181)
                    {
                        TotalCumValue[59] <- TotalCumValue[59] + objInit$x[a,b]
                    }
                    if (a >= 1181 && a < 1201)
                    {
                        TotalCumValue[60] <- TotalCumValue[60] + objInit$x[a,b]
                    }
                    if (a >= 1201 && a < 1221)
                    {
                        TotalCumValue[61] <- TotalCumValue[61] + objInit$x[a,b]
                    }
                    if (a >= 1221 && a < 1241)
                    {
                        TotalCumValue[62] <- TotalCumValue[62] + objInit$x[a,b]
                    }
                    if (a >= 1241 && a < 1261)
                    {
                        TotalCumValue[63] <- TotalCumValue[63] + objInit$x[a,b]
                    }
                    if (a >= 1261 && a < 1281)
                    {
                        TotalCumValue[64] <- TotalCumValue[64] + objInit$x[a,b]
                    }
                    if (a >= 1281 && a < 1301)
                    {
                        TotalCumValue[65] <- TotalCumValue[65] + objInit$x[a,b]
                    }
                    if (a >= 1301 && a < 1321)
                    {
                        TotalCumValue[66] <- TotalCumValue[66] + objInit$x[a,b]
                    }
                    if (a >= 1321 && a < 1341)
                    {
                        TotalCumValue[67] <- TotalCumValue[67] + objInit$x[a,b]
                    }
                    if (a >= 1341 && a < 1361)
                    {
                        TotalCumValue[68] <- TotalCumValue[68] + objInit$x[a,b]
                    }
                    if (a >= 1361 && a < 1381)
                    {
                        TotalCumValue[69] <- TotalCumValue[69] + objInit$x[a,b]
                    }
                    if (a >= 1381 && a < 1401)
                    {
                        TotalCumValue[70] <- TotalCumValue[70] + objInit$x[a,b]
                    }                                              
                } # end for a
                
                minValuePC3 <- min(TotalCumValue)
                maxValuePC3 <- max(TotalCumValue)
                xxmin <- -1000
                xxmax <- -1000
                
                for (kkk in 1:70)
                {
                    if (minValuePC3 == TotalCumValue[kkk])
                    {
                        xxmin <- kkk
                    }
                    if (maxValuePC3 == TotalCumValue[kkk])
                    {
                        xxmax <- kkk
                    }
                }
                
                minLinePC3 <- (xxmin * 20) - 20 + 1
                maxLinePC3 <- (xxmax * 20) - 20 + 1
                
                #kkk <- 20 # highest/lowest 20 points
                #ndx <- order(objInit$x[,b], decreasing = T)[1:kkk]
                #maxValuePC3_vector <- objInit$x[,b][ndx]
                #ndx <- order(objInit$x[,b])[1:kkk]
                #minValuePC3_vector <- objInit$x[,b][ndx]
                
            } # end b == 3
        } # end for b

		matrix1_results[1,1]<-minLinePC1
		matrix1_results[1,2]<-maxLinePC1
		matrix1_results[1,3]<-minLinePC2
		matrix1_results[1,4]<-maxLinePC2
		matrix1_results[1,5]<-minLinePC3
		matrix1_results[1,6]<-maxLinePC3

		} else # the rest of genetic datasets
		{
		print("> Working on file:")
		print(string1)

		for (a in 1:100) #sdev
			{
			cumObj$sdev[a] <- cumObj$sdev[a] + objInit$sdev[a]
			}
		for (a in 1:100) #rotation
			{
			for (b in 1:100)
				{
				cumObj$rotation[a,b] <- cumObj$rotation[a,b] + objInit$rotation[a,b]
				}
			}
		for (a in 1:100) #center
			{
			cumObj$center[a] <- cumObj$center[a] + objInit$center[a]
			}
		for (a in 1:100) #scale
			{
			cumObj$scale[a] <- cumObj$scale[a] + objInit$scale[a]
			}
		for (a in 1:1400) #x
			{
			for (b in 1:100)
				{
				cumObj$x[a,b] <- cumObj$x[a,b] + objInit$x[a,b]
				}
			}

        TotalCumValue <- vector()
        for (b in 1:3)
        {
            if (b == 1) # PC1
            {
                
                for (kkk in 1:70)
                {
                    TotalCumValue[kkk] <- 0
                }
                
                for (a in 1:1400)
                {
                    if (a >= 1 && a < 21)
                    {
                        TotalCumValue[1] <- TotalCumValue[1] + objInit$x[a,b]
                    }
                    if (a >= 21 && a < 41)
                    {
                        TotalCumValue[2] <- TotalCumValue[2] + objInit$x[a,b]
                    }
                    if (a >= 41 && a < 61)
                    {
                        TotalCumValue[3] <- TotalCumValue[3] + objInit$x[a,b]
                    }
                    if (a >= 61 && a < 81)
                    {
                        TotalCumValue[4] <- TotalCumValue[4] + objInit$x[a,b]
                    }
                    if (a >= 81 && a < 101)
                    {
                        TotalCumValue[5] <- TotalCumValue[5] + objInit$x[a,b]
                    }
                    if (a >= 101 && a < 121)
                    {
                        TotalCumValue[6] <- TotalCumValue[6] + objInit$x[a,b]
                    }
                    if (a >= 121 && a < 141)
                    {
                        TotalCumValue[7] <- TotalCumValue[7] + objInit$x[a,b]
                    }
                    if (a >= 141 && a < 161)
                    {
                        TotalCumValue[8] <- TotalCumValue[8] + objInit$x[a,b]
                    }
                    if (a >= 161 && a < 181)
                    {
                        TotalCumValue[9] <- TotalCumValue[9] + objInit$x[a,b]
                    }
                    if (a >= 181 && a < 201)
                    {
                        TotalCumValue[10] <- TotalCumValue[10] + objInit$x[a,b]
                    }
                    if (a >= 201 && a < 221)
                    {
                        TotalCumValue[11] <- TotalCumValue[11] + objInit$x[a,b]
                    }
                    if (a >= 221 && a < 241)
                    {
                        TotalCumValue[12] <- TotalCumValue[12] + objInit$x[a,b]
                    }
                    if (a >= 241 && a < 261)
                    {
                        TotalCumValue[13] <- TotalCumValue[13] + objInit$x[a,b]
                    }
                    if (a >= 261 && a < 281)
                    {
                        TotalCumValue[14] <- TotalCumValue[14] + objInit$x[a,b]
                    }
                    if (a >= 281 && a < 301)
                    {
                        TotalCumValue[15] <- TotalCumValue[15] + objInit$x[a,b]
                    }
                    if (a >= 301 && a < 321)
                    {
                        TotalCumValue[16] <- TotalCumValue[16] + objInit$x[a,b]
                    }
                    if (a >= 321 && a < 341)
                    {
                        TotalCumValue[17] <- TotalCumValue[17] + objInit$x[a,b]
                    }
                    if (a >= 341 && a < 361)
                    {
                        TotalCumValue[18] <- TotalCumValue[18] + objInit$x[a,b]
                    }
                    if (a >= 361 && a < 381)
                    {
                        TotalCumValue[19] <- TotalCumValue[19] + objInit$x[a,b]
                    }
                    if (a >= 381 && a < 401)
                    {
                        TotalCumValue[20] <- TotalCumValue[20] + objInit$x[a,b]
                    }
                    if (a >= 401 && a < 421)
                    {
                        TotalCumValue[21] <- TotalCumValue[21] + objInit$x[a,b]
                    }
                    if (a >= 421 && a < 441)
                    {
                        TotalCumValue[22] <- TotalCumValue[22] + objInit$x[a,b]
                    }
                    if (a >= 441 && a < 461)
                    {
                        TotalCumValue[23] <- TotalCumValue[23] + objInit$x[a,b]
                    }
                    if (a >= 461 && a < 481)
                    {
                        TotalCumValue[24] <- TotalCumValue[24] + objInit$x[a,b]
                    }
                    if (a >= 481 && a < 501)
                    {
                        TotalCumValue[25] <- TotalCumValue[25] + objInit$x[a,b]
                    }
                    if (a >= 501 && a < 521)
                    {
                        TotalCumValue[26] <- TotalCumValue[26] + objInit$x[a,b]
                    }
                    if (a >= 521 && a < 541)
                    {
                        TotalCumValue[27] <- TotalCumValue[27] + objInit$x[a,b]
                    }
                    if (a >= 541 && a < 561)
                    {
                        TotalCumValue[28] <- TotalCumValue[28] + objInit$x[a,b]
                    }
                    if (a >= 561 && a < 581)
                    {
                        TotalCumValue[29] <- TotalCumValue[29] + objInit$x[a,b]
                    }
                    if (a >= 581 && a < 601)
                    {
                        TotalCumValue[30] <- TotalCumValue[30] + objInit$x[a,b]
                    }
                    if (a >= 601 && a < 621)
                    {
                        TotalCumValue[31] <- TotalCumValue[31] + objInit$x[a,b]
                    }
                    if (a >= 621 && a < 641)
                    {
                        TotalCumValue[32] <- TotalCumValue[32] + objInit$x[a,b]
                    }
                    if (a >= 641 && a < 661)
                    {
                        TotalCumValue[33] <- TotalCumValue[33] + objInit$x[a,b]
                    }
                    if (a >= 661 && a < 681)
                    {
                        TotalCumValue[34] <- TotalCumValue[34] + objInit$x[a,b]
                    }
                    if (a >= 681 && a < 701)
                    {
                        TotalCumValue[35] <- TotalCumValue[35] + objInit$x[a,b]
                    }
                    if (a >= 701 && a < 721)
                    {
                        TotalCumValue[36] <- TotalCumValue[36] + objInit$x[a,b]
                    }
                    if (a >= 721 && a < 741)
                    {
                        TotalCumValue[37] <- TotalCumValue[37] + objInit$x[a,b]
                    }
                    if (a >= 741 && a < 761)
                    {
                        TotalCumValue[38] <- TotalCumValue[38] + objInit$x[a,b]
                    }
                    if (a >= 761 && a < 781)
                    {
                        TotalCumValue[39] <- TotalCumValue[39] + objInit$x[a,b]
                    }
                    if (a >= 781 && a < 801)
                    {
                        TotalCumValue[40] <- TotalCumValue[40] + objInit$x[a,b]
                    }
                    if (a >= 801 && a < 821)
                    {
                        TotalCumValue[41] <- TotalCumValue[41] + objInit$x[a,b]
                    }
                    if (a >= 821 && a < 841)
                    {
                        TotalCumValue[42] <- TotalCumValue[42] + objInit$x[a,b]
                    }
                    if (a >= 841 && a < 861)
                    {
                        TotalCumValue[43] <- TotalCumValue[43] + objInit$x[a,b]
                    }
                    if (a >= 861 && a < 881)
                    {
                        TotalCumValue[44] <- TotalCumValue[44] + objInit$x[a,b]
                    }
                    if (a >= 881 && a < 901)
                    {
                        TotalCumValue[45] <- TotalCumValue[45] + objInit$x[a,b]
                    }
                    if (a >= 901 && a < 921)
                    {
                        TotalCumValue[46] <- TotalCumValue[46] + objInit$x[a,b]
                    }
                    if (a >= 921 && a < 941)
                    {
                        TotalCumValue[47] <- TotalCumValue[47] + objInit$x[a,b]
                    }
                    if (a >= 941 && a < 961)
                    {
                        TotalCumValue[48] <- TotalCumValue[48] + objInit$x[a,b]
                    }
                    if (a >= 961 && a < 981)
                    {
                        TotalCumValue[49] <- TotalCumValue[49] + objInit$x[a,b]
                    }
                    if (a >= 981 && a < 1001)
                    {
                        TotalCumValue[50] <- TotalCumValue[50] + objInit$x[a,b]
                    }
                    if (a >= 1001 && a < 1021)
                    {
                        TotalCumValue[51] <- TotalCumValue[51] + objInit$x[a,b]
                    }
                    if (a >= 1021 && a < 1041)
                    {
                        TotalCumValue[52] <- TotalCumValue[52] + objInit$x[a,b]
                    }
                    if (a >= 1041 && a < 1061)
                    {
                        TotalCumValue[53] <- TotalCumValue[53] + objInit$x[a,b]
                    }
                    if (a >= 1061 && a < 1081)
                    {
                        TotalCumValue[54] <- TotalCumValue[54] + objInit$x[a,b]
                    }
                    if (a >= 1081 && a < 1101)
                    {
                        TotalCumValue[55] <- TotalCumValue[55] + objInit$x[a,b]
                    }
                    if (a >= 1101 && a < 1121)
                    {
                        TotalCumValue[56] <- TotalCumValue[56] + objInit$x[a,b]
                    }
                    if (a >= 1121 && a < 1141)
                    {
                        TotalCumValue[57] <- TotalCumValue[57] + objInit$x[a,b]
                    }
                    if (a >= 1141 && a < 1161)
                    {
                        TotalCumValue[58] <- TotalCumValue[58] + objInit$x[a,b]
                    }
                    if (a >= 1161 && a < 1181)
                    {
                        TotalCumValue[59] <- TotalCumValue[59] + objInit$x[a,b]
                    }
                    if (a >= 1181 && a < 1201)
                    {
                        TotalCumValue[60] <- TotalCumValue[60] + objInit$x[a,b]
                    }
                    if (a >= 1201 && a < 1221)
                    {
                        TotalCumValue[61] <- TotalCumValue[61] + objInit$x[a,b]
                    }
                    if (a >= 1221 && a < 1241)
                    {
                        TotalCumValue[62] <- TotalCumValue[62] + objInit$x[a,b]
                    }
                    if (a >= 1241 && a < 1261)
                    {
                        TotalCumValue[63] <- TotalCumValue[63] + objInit$x[a,b]
                    }
                    if (a >= 1261 && a < 1281)
                    {
                        TotalCumValue[64] <- TotalCumValue[64] + objInit$x[a,b]
                    }
                    if (a >= 1281 && a < 1301)
                    {
                        TotalCumValue[65] <- TotalCumValue[65] + objInit$x[a,b]
                    }
                    if (a >= 1301 && a < 1321)
                    {
                        TotalCumValue[66] <- TotalCumValue[66] + objInit$x[a,b]
                    }
                    if (a >= 1321 && a < 1341)
                    {
                        TotalCumValue[67] <- TotalCumValue[67] + objInit$x[a,b]
                    }
                    if (a >= 1341 && a < 1361)
                    {
                        TotalCumValue[68] <- TotalCumValue[68] + objInit$x[a,b]
                    }
                    if (a >= 1361 && a < 1381)
                    {
                        TotalCumValue[69] <- TotalCumValue[69] + objInit$x[a,b]
                    }
                    if (a >= 1381 && a < 1401)
                    {
                        TotalCumValue[70] <- TotalCumValue[70] + objInit$x[a,b]
                    }                                             
                } # end for a
                
                minValuePC1 <- min(TotalCumValue)
                maxValuePC1 <- max(TotalCumValue)
                xxmin <- -1000
                xxmax <- -1000
                
                for (kkk in 1:70)
                {
                    if (minValuePC1 == TotalCumValue[kkk])
                    {
                        xxmin <- kkk
                    }
                    if (maxValuePC1 == TotalCumValue[kkk])
                    {
                        xxmax <- kkk
                    }
                }
                
                minLinePC1 <- (xxmin * 20) - 20 + 1
                maxLinePC1 <- (xxmax * 20) - 20 + 1
                
                #kkk <- 20 # highest/lowest 20 points
                #ndx <- order(objInit$x[,b], decreasing = T)[1:kkk]
                #maxValuePC1_vector <- objInit$x[,b][ndx]
                #ndx <- order(objInit$x[,b])[1:kkk]
                #minValuePC1_vector <- objInit$x[,b][ndx]
            } # end b == 1
            
            if (b == 2) # PC2
            {
                for (kkk in 1:70)
                {
                    TotalCumValue[kkk] <- 0
                }
                
                for (a in 1:1400)
                {
                    if (a >= 1 && a < 21)
                    {
                        TotalCumValue[1] <- TotalCumValue[1] + objInit$x[a,b]
                    }
                    if (a >= 21 && a < 41)
                    {
                        TotalCumValue[2] <- TotalCumValue[2] + objInit$x[a,b]
                    }
                    if (a >= 41 && a < 61)
                    {
                        TotalCumValue[3] <- TotalCumValue[3] + objInit$x[a,b]
                    }
                    if (a >= 61 && a < 81)
                    {
                        TotalCumValue[4] <- TotalCumValue[4] + objInit$x[a,b]
                    }
                    if (a >= 81 && a < 101)
                    {
                        TotalCumValue[5] <- TotalCumValue[5] + objInit$x[a,b]
                    }
                    if (a >= 101 && a < 121)
                    {
                        TotalCumValue[6] <- TotalCumValue[6] + objInit$x[a,b]
                    }
                    if (a >= 121 && a < 141)
                    {
                        TotalCumValue[7] <- TotalCumValue[7] + objInit$x[a,b]
                    }
                    if (a >= 141 && a < 161)
                    {
                        TotalCumValue[8] <- TotalCumValue[8] + objInit$x[a,b]
                    }
                    if (a >= 161 && a < 181)
                    {
                        TotalCumValue[9] <- TotalCumValue[9] + objInit$x[a,b]
                    }
                    if (a >= 181 && a < 201)
                    {
                        TotalCumValue[10] <- TotalCumValue[10] + objInit$x[a,b]
                    }
                    if (a >= 201 && a < 221)
                    {
                        TotalCumValue[11] <- TotalCumValue[11] + objInit$x[a,b]
                    }
                    if (a >= 221 && a < 241)
                    {
                        TotalCumValue[12] <- TotalCumValue[12] + objInit$x[a,b]
                    }
                    if (a >= 241 && a < 261)
                    {
                        TotalCumValue[13] <- TotalCumValue[13] + objInit$x[a,b]
                    }
                    if (a >= 261 && a < 281)
                    {
                        TotalCumValue[14] <- TotalCumValue[14] + objInit$x[a,b]
                    }
                    if (a >= 281 && a < 301)
                    {
                        TotalCumValue[15] <- TotalCumValue[15] + objInit$x[a,b]
                    }
                    if (a >= 301 && a < 321)
                    {
                        TotalCumValue[16] <- TotalCumValue[16] + objInit$x[a,b]
                    }
                    if (a >= 321 && a < 341)
                    {
                        TotalCumValue[17] <- TotalCumValue[17] + objInit$x[a,b]
                    }
                    if (a >= 341 && a < 361)
                    {
                        TotalCumValue[18] <- TotalCumValue[18] + objInit$x[a,b]
                    }
                    if (a >= 361 && a < 381)
                    {
                        TotalCumValue[19] <- TotalCumValue[19] + objInit$x[a,b]
                    }
                    if (a >= 381 && a < 401)
                    {
                        TotalCumValue[20] <- TotalCumValue[20] + objInit$x[a,b]
                    }
                    if (a >= 401 && a < 421)
                    {
                        TotalCumValue[21] <- TotalCumValue[21] + objInit$x[a,b]
                    }
                    if (a >= 421 && a < 441)
                    {
                        TotalCumValue[22] <- TotalCumValue[22] + objInit$x[a,b]
                    }
                    if (a >= 441 && a < 461)
                    {
                        TotalCumValue[23] <- TotalCumValue[23] + objInit$x[a,b]
                    }
                    if (a >= 461 && a < 481)
                    {
                        TotalCumValue[24] <- TotalCumValue[24] + objInit$x[a,b]
                    }
                    if (a >= 481 && a < 501)
                    {
                        TotalCumValue[25] <- TotalCumValue[25] + objInit$x[a,b]
                    }
                    if (a >= 501 && a < 521)
                    {
                        TotalCumValue[26] <- TotalCumValue[26] + objInit$x[a,b]
                    }
                    if (a >= 521 && a < 541)
                    {
                        TotalCumValue[27] <- TotalCumValue[27] + objInit$x[a,b]
                    }
                    if (a >= 541 && a < 561)
                    {
                        TotalCumValue[28] <- TotalCumValue[28] + objInit$x[a,b]
                    }
                    if (a >= 561 && a < 581)
                    {
                        TotalCumValue[29] <- TotalCumValue[29] + objInit$x[a,b]
                    }
                    if (a >= 581 && a < 601)
                    {
                        TotalCumValue[30] <- TotalCumValue[30] + objInit$x[a,b]
                    }
                    if (a >= 601 && a < 621)
                    {
                        TotalCumValue[31] <- TotalCumValue[31] + objInit$x[a,b]
                    }
                    if (a >= 621 && a < 641)
                    {
                        TotalCumValue[32] <- TotalCumValue[32] + objInit$x[a,b]
                    }
                    if (a >= 641 && a < 661)
                    {
                        TotalCumValue[33] <- TotalCumValue[33] + objInit$x[a,b]
                    }
                    if (a >= 661 && a < 681)
                    {
                        TotalCumValue[34] <- TotalCumValue[34] + objInit$x[a,b]
                    }
                    if (a >= 681 && a < 701)
                    {
                        TotalCumValue[35] <- TotalCumValue[35] + objInit$x[a,b]
                    }
                    if (a >= 701 && a < 721)
                    {
                        TotalCumValue[36] <- TotalCumValue[36] + objInit$x[a,b]
                    }
                    if (a >= 721 && a < 741)
                    {
                        TotalCumValue[37] <- TotalCumValue[37] + objInit$x[a,b]
                    }
                    if (a >= 741 && a < 761)
                    {
                        TotalCumValue[38] <- TotalCumValue[38] + objInit$x[a,b]
                    }
                    if (a >= 761 && a < 781)
                    {
                        TotalCumValue[39] <- TotalCumValue[39] + objInit$x[a,b]
                    }
                    if (a >= 781 && a < 801)
                    {
                        TotalCumValue[40] <- TotalCumValue[40] + objInit$x[a,b]
                    }
                    if (a >= 801 && a < 821)
                    {
                        TotalCumValue[41] <- TotalCumValue[41] + objInit$x[a,b]
                    }
                    if (a >= 821 && a < 841)
                    {
                        TotalCumValue[42] <- TotalCumValue[42] + objInit$x[a,b]
                    }
                    if (a >= 841 && a < 861)
                    {
                        TotalCumValue[43] <- TotalCumValue[43] + objInit$x[a,b]
                    }
                    if (a >= 861 && a < 881)
                    {
                        TotalCumValue[44] <- TotalCumValue[44] + objInit$x[a,b]
                    }
                    if (a >= 881 && a < 901)
                    {
                        TotalCumValue[45] <- TotalCumValue[45] + objInit$x[a,b]
                    }
                    if (a >= 901 && a < 921)
                    {
                        TotalCumValue[46] <- TotalCumValue[46] + objInit$x[a,b]
                    }
                    if (a >= 921 && a < 941)
                    {
                        TotalCumValue[47] <- TotalCumValue[47] + objInit$x[a,b]
                    }
                    if (a >= 941 && a < 961)
                    {
                        TotalCumValue[48] <- TotalCumValue[48] + objInit$x[a,b]
                    }
                    if (a >= 961 && a < 981)
                    {
                        TotalCumValue[49] <- TotalCumValue[49] + objInit$x[a,b]
                    }
                    if (a >= 981 && a < 1001)
                    {
                        TotalCumValue[50] <- TotalCumValue[50] + objInit$x[a,b]
                    }
                    if (a >= 1001 && a < 1021)
                    {
                        TotalCumValue[51] <- TotalCumValue[51] + objInit$x[a,b]
                    }
                    if (a >= 1021 && a < 1041)
                    {
                        TotalCumValue[52] <- TotalCumValue[52] + objInit$x[a,b]
                    }
                    if (a >= 1041 && a < 1061)
                    {
                        TotalCumValue[53] <- TotalCumValue[53] + objInit$x[a,b]
                    }
                    if (a >= 1061 && a < 1081)
                    {
                        TotalCumValue[54] <- TotalCumValue[54] + objInit$x[a,b]
                    }
                    if (a >= 1081 && a < 1101)
                    {
                        TotalCumValue[55] <- TotalCumValue[55] + objInit$x[a,b]
                    }
                    if (a >= 1101 && a < 1121)
                    {
                        TotalCumValue[56] <- TotalCumValue[56] + objInit$x[a,b]
                    }
                    if (a >= 1121 && a < 1141)
                    {
                        TotalCumValue[57] <- TotalCumValue[57] + objInit$x[a,b]
                    }
                    if (a >= 1141 && a < 1161)
                    {
                        TotalCumValue[58] <- TotalCumValue[58] + objInit$x[a,b]
                    }
                    if (a >= 1161 && a < 1181)
                    {
                        TotalCumValue[59] <- TotalCumValue[59] + objInit$x[a,b]
                    }
                    if (a >= 1181 && a < 1201)
                    {
                        TotalCumValue[60] <- TotalCumValue[60] + objInit$x[a,b]
                    }
                    if (a >= 1201 && a < 1221)
                    {
                        TotalCumValue[61] <- TotalCumValue[61] + objInit$x[a,b]
                    }
                    if (a >= 1221 && a < 1241)
                    {
                        TotalCumValue[62] <- TotalCumValue[62] + objInit$x[a,b]
                    }
                    if (a >= 1241 && a < 1261)
                    {
                        TotalCumValue[63] <- TotalCumValue[63] + objInit$x[a,b]
                    }
                    if (a >= 1261 && a < 1281)
                    {
                        TotalCumValue[64] <- TotalCumValue[64] + objInit$x[a,b]
                    }
                    if (a >= 1281 && a < 1301)
                    {
                        TotalCumValue[65] <- TotalCumValue[65] + objInit$x[a,b]
                    }
                    if (a >= 1301 && a < 1321)
                    {
                        TotalCumValue[66] <- TotalCumValue[66] + objInit$x[a,b]
                    }
                    if (a >= 1321 && a < 1341)
                    {
                        TotalCumValue[67] <- TotalCumValue[67] + objInit$x[a,b]
                    }
                    if (a >= 1341 && a < 1361)
                    {
                        TotalCumValue[68] <- TotalCumValue[68] + objInit$x[a,b]
                    }
                    if (a >= 1361 && a < 1381)
                    {
                        TotalCumValue[69] <- TotalCumValue[69] + objInit$x[a,b]
                    }
                    if (a >= 1381 && a < 1401)
                    {
                        TotalCumValue[70] <- TotalCumValue[70] + objInit$x[a,b]
                    }                                              
                } # end for a
                
                minValuePC2 <- min(TotalCumValue)
                maxValuePC2 <- max(TotalCumValue)
                xxmin <- -1000
                xxmax <- -1000
                
                for (kkk in 1:70)
                {
                    if (minValuePC2 == TotalCumValue[kkk])
                    {
                        xxmin <- kkk
                    }
                    if (maxValuePC2 == TotalCumValue[kkk])
                    {
                        xxmax <- kkk
                    }
                }
                
                minLinePC2 <- (xxmin * 20) - 20 + 1
                maxLinePC2 <- (xxmax * 20) - 20 + 1
                
                #kkk <- 20 # highest/lowest 20 points
                #ndx <- order(objInit$x[,b], decreasing = T)[1:kkk]
                #maxValuePC2_vector <- objInit$x[,b][ndx]
                #ndx <- order(objInit$x[,b])[1:kkk]
                #minValuePC2_vector <- objInit$x[,b][ndx]
            } # end b == 2
            
            if (b == 3) # PC3
            {
                for (kkk in 1:70)
                {
                    TotalCumValue[kkk] <- 0
                }
                
                for (a in 1:1400)
                {
                    if (a >= 1 && a < 21)
                    {
                        TotalCumValue[1] <- TotalCumValue[1] + objInit$x[a,b]
                    }
                    if (a >= 21 && a < 41)
                    {
                        TotalCumValue[2] <- TotalCumValue[2] + objInit$x[a,b]
                    }
                    if (a >= 41 && a < 61)
                    {
                        TotalCumValue[3] <- TotalCumValue[3] + objInit$x[a,b]
                    }
                    if (a >= 61 && a < 81)
                    {
                        TotalCumValue[4] <- TotalCumValue[4] + objInit$x[a,b]
                    }
                    if (a >= 81 && a < 101)
                    {
                        TotalCumValue[5] <- TotalCumValue[5] + objInit$x[a,b]
                    }
                    if (a >= 101 && a < 121)
                    {
                        TotalCumValue[6] <- TotalCumValue[6] + objInit$x[a,b]
                    }
                    if (a >= 121 && a < 141)
                    {
                        TotalCumValue[7] <- TotalCumValue[7] + objInit$x[a,b]
                    }
                    if (a >= 141 && a < 161)
                    {
                        TotalCumValue[8] <- TotalCumValue[8] + objInit$x[a,b]
                    }
                    if (a >= 161 && a < 181)
                    {
                        TotalCumValue[9] <- TotalCumValue[9] + objInit$x[a,b]
                    }
                    if (a >= 181 && a < 201)
                    {
                        TotalCumValue[10] <- TotalCumValue[10] + objInit$x[a,b]
                    }
                    if (a >= 201 && a < 221)
                    {
                        TotalCumValue[11] <- TotalCumValue[11] + objInit$x[a,b]
                    }
                    if (a >= 221 && a < 241)
                    {
                        TotalCumValue[12] <- TotalCumValue[12] + objInit$x[a,b]
                    }
                    if (a >= 241 && a < 261)
                    {
                        TotalCumValue[13] <- TotalCumValue[13] + objInit$x[a,b]
                    }
                    if (a >= 261 && a < 281)
                    {
                        TotalCumValue[14] <- TotalCumValue[14] + objInit$x[a,b]
                    }
                    if (a >= 281 && a < 301)
                    {
                        TotalCumValue[15] <- TotalCumValue[15] + objInit$x[a,b]
                    }
                    if (a >= 301 && a < 321)
                    {
                        TotalCumValue[16] <- TotalCumValue[16] + objInit$x[a,b]
                    }
                    if (a >= 321 && a < 341)
                    {
                        TotalCumValue[17] <- TotalCumValue[17] + objInit$x[a,b]
                    }
                    if (a >= 341 && a < 361)
                    {
                        TotalCumValue[18] <- TotalCumValue[18] + objInit$x[a,b]
                    }
                    if (a >= 361 && a < 381)
                    {
                        TotalCumValue[19] <- TotalCumValue[19] + objInit$x[a,b]
                    }
                    if (a >= 381 && a < 401)
                    {
                        TotalCumValue[20] <- TotalCumValue[20] + objInit$x[a,b]
                    }
                    if (a >= 401 && a < 421)
                    {
                        TotalCumValue[21] <- TotalCumValue[21] + objInit$x[a,b]
                    }
                    if (a >= 421 && a < 441)
                    {
                        TotalCumValue[22] <- TotalCumValue[22] + objInit$x[a,b]
                    }
                    if (a >= 441 && a < 461)
                    {
                        TotalCumValue[23] <- TotalCumValue[23] + objInit$x[a,b]
                    }
                    if (a >= 461 && a < 481)
                    {
                        TotalCumValue[24] <- TotalCumValue[24] + objInit$x[a,b]
                    }
                    if (a >= 481 && a < 501)
                    {
                        TotalCumValue[25] <- TotalCumValue[25] + objInit$x[a,b]
                    }
                    if (a >= 501 && a < 521)
                    {
                        TotalCumValue[26] <- TotalCumValue[26] + objInit$x[a,b]
                    }
                    if (a >= 521 && a < 541)
                    {
                        TotalCumValue[27] <- TotalCumValue[27] + objInit$x[a,b]
                    }
                    if (a >= 541 && a < 561)
                    {
                        TotalCumValue[28] <- TotalCumValue[28] + objInit$x[a,b]
                    }
                    if (a >= 561 && a < 581)
                    {
                        TotalCumValue[29] <- TotalCumValue[29] + objInit$x[a,b]
                    }
                    if (a >= 581 && a < 601)
                    {
                        TotalCumValue[30] <- TotalCumValue[30] + objInit$x[a,b]
                    }
                    if (a >= 601 && a < 621)
                    {
                        TotalCumValue[31] <- TotalCumValue[31] + objInit$x[a,b]
                    }
                    if (a >= 621 && a < 641)
                    {
                        TotalCumValue[32] <- TotalCumValue[32] + objInit$x[a,b]
                    }
                    if (a >= 641 && a < 661)
                    {
                        TotalCumValue[33] <- TotalCumValue[33] + objInit$x[a,b]
                    }
                    if (a >= 661 && a < 681)
                    {
                        TotalCumValue[34] <- TotalCumValue[34] + objInit$x[a,b]
                    }
                    if (a >= 681 && a < 701)
                    {
                        TotalCumValue[35] <- TotalCumValue[35] + objInit$x[a,b]
                    }
                    if (a >= 701 && a < 721)
                    {
                        TotalCumValue[36] <- TotalCumValue[36] + objInit$x[a,b]
                    }
                    if (a >= 721 && a < 741)
                    {
                        TotalCumValue[37] <- TotalCumValue[37] + objInit$x[a,b]
                    }
                    if (a >= 741 && a < 761)
                    {
                        TotalCumValue[38] <- TotalCumValue[38] + objInit$x[a,b]
                    }
                    if (a >= 761 && a < 781)
                    {
                        TotalCumValue[39] <- TotalCumValue[39] + objInit$x[a,b]
                    }
                    if (a >= 781 && a < 801)
                    {
                        TotalCumValue[40] <- TotalCumValue[40] + objInit$x[a,b]
                    }
                    if (a >= 801 && a < 821)
                    {
                        TotalCumValue[41] <- TotalCumValue[41] + objInit$x[a,b]
                    }
                    if (a >= 821 && a < 841)
                    {
                        TotalCumValue[42] <- TotalCumValue[42] + objInit$x[a,b]
                    }
                    if (a >= 841 && a < 861)
                    {
                        TotalCumValue[43] <- TotalCumValue[43] + objInit$x[a,b]
                    }
                    if (a >= 861 && a < 881)
                    {
                        TotalCumValue[44] <- TotalCumValue[44] + objInit$x[a,b]
                    }
                    if (a >= 881 && a < 901)
                    {
                        TotalCumValue[45] <- TotalCumValue[45] + objInit$x[a,b]
                    }
                    if (a >= 901 && a < 921)
                    {
                        TotalCumValue[46] <- TotalCumValue[46] + objInit$x[a,b]
                    }
                    if (a >= 921 && a < 941)
                    {
                        TotalCumValue[47] <- TotalCumValue[47] + objInit$x[a,b]
                    }
                    if (a >= 941 && a < 961)
                    {
                        TotalCumValue[48] <- TotalCumValue[48] + objInit$x[a,b]
                    }
                    if (a >= 961 && a < 981)
                    {
                        TotalCumValue[49] <- TotalCumValue[49] + objInit$x[a,b]
                    }
                    if (a >= 981 && a < 1001)
                    {
                        TotalCumValue[50] <- TotalCumValue[50] + objInit$x[a,b]
                    }
                    if (a >= 1001 && a < 1021)
                    {
                        TotalCumValue[51] <- TotalCumValue[51] + objInit$x[a,b]
                    }
                    if (a >= 1021 && a < 1041)
                    {
                        TotalCumValue[52] <- TotalCumValue[52] + objInit$x[a,b]
                    }
                    if (a >= 1041 && a < 1061)
                    {
                        TotalCumValue[53] <- TotalCumValue[53] + objInit$x[a,b]
                    }
                    if (a >= 1061 && a < 1081)
                    {
                        TotalCumValue[54] <- TotalCumValue[54] + objInit$x[a,b]
                    }
                    if (a >= 1081 && a < 1101)
                    {
                        TotalCumValue[55] <- TotalCumValue[55] + objInit$x[a,b]
                    }
                    if (a >= 1101 && a < 1121)
                    {
                        TotalCumValue[56] <- TotalCumValue[56] + objInit$x[a,b]
                    }
                    if (a >= 1121 && a < 1141)
                    {
                        TotalCumValue[57] <- TotalCumValue[57] + objInit$x[a,b]
                    }
                    if (a >= 1141 && a < 1161)
                    {
                        TotalCumValue[58] <- TotalCumValue[58] + objInit$x[a,b]
                    }
                    if (a >= 1161 && a < 1181)
                    {
                        TotalCumValue[59] <- TotalCumValue[59] + objInit$x[a,b]
                    }
                    if (a >= 1181 && a < 1201)
                    {
                        TotalCumValue[60] <- TotalCumValue[60] + objInit$x[a,b]
                    }
                    if (a >= 1201 && a < 1221)
                    {
                        TotalCumValue[61] <- TotalCumValue[61] + objInit$x[a,b]
                    }
                    if (a >= 1221 && a < 1241)
                    {
                        TotalCumValue[62] <- TotalCumValue[62] + objInit$x[a,b]
                    }
                    if (a >= 1241 && a < 1261)
                    {
                        TotalCumValue[63] <- TotalCumValue[63] + objInit$x[a,b]
                    }
                    if (a >= 1261 && a < 1281)
                    {
                        TotalCumValue[64] <- TotalCumValue[64] + objInit$x[a,b]
                    }
                    if (a >= 1281 && a < 1301)
                    {
                        TotalCumValue[65] <- TotalCumValue[65] + objInit$x[a,b]
                    }
                    if (a >= 1301 && a < 1321)
                    {
                        TotalCumValue[66] <- TotalCumValue[66] + objInit$x[a,b]
                    }
                    if (a >= 1321 && a < 1341)
                    {
                        TotalCumValue[67] <- TotalCumValue[67] + objInit$x[a,b]
                    }
                    if (a >= 1341 && a < 1361)
                    {
                        TotalCumValue[68] <- TotalCumValue[68] + objInit$x[a,b]
                    }
                    if (a >= 1361 && a < 1381)
                    {
                        TotalCumValue[69] <- TotalCumValue[69] + objInit$x[a,b]
                    }
                    if (a >= 1381 && a < 1401)
                    {
                        TotalCumValue[70] <- TotalCumValue[70] + objInit$x[a,b]
                    }                                             
                } # end for a
                
                minValuePC3 <- min(TotalCumValue)
                maxValuePC3 <- max(TotalCumValue)
                xxmin <- -1000
                xxmax <- -1000
                
                for (kkk in 1:70)
                {
                    if (minValuePC3 == TotalCumValue[kkk])
                    {
                        xxmin <- kkk
                    }
                    if (maxValuePC3 == TotalCumValue[kkk])
                    {
                        xxmax <- kkk
                    }
                }
                
                minLinePC3 <- (xxmin * 20) - 20 + 1
                maxLinePC3 <- (xxmax * 20) - 20 + 1
                
                #kkk <- 20 # highest/lowest 20 points
                #ndx <- order(objInit$x[,b], decreasing = T)[1:kkk]
                #maxValuePC3_vector <- objInit$x[,b][ndx]
                #ndx <- order(objInit$x[,b])[1:kkk]
                #minValuePC3_vector <- objInit$x[,b][ndx]
            } # end b == 3
        }

        # print(objInit$sdev)
		# print(cumObj$sdev)

		matrix1_results[n,1]<-minLinePC1
		matrix1_results[n,2]<-maxLinePC1
		matrix1_results[n,3]<-minLinePC2
		matrix1_results[n,4]<-maxLinePC2
		matrix1_results[n,3]<-minLinePC3
		matrix1_results[n,4]<-maxLinePC3
		}
		
	### Calculation of slopes and intercepts for PC1, PC2 and PC3
	
	xMinPC1 <- coord[minLinePC1,1]
	yMinPC1 <- coord[minLinePC1,2]
	xMaxPC1 <- coord[maxLinePC1,1]
	yMaxPC1 <- coord[maxLinePC1,2]
	xMinPC2 <- coord[minLinePC2,1]
	yMinPC2 <- coord[minLinePC2,2]
	xMaxPC2 <- coord[maxLinePC2,1]
	yMaxPC2 <- coord[maxLinePC2,2]
	xMinPC3 <- coord[minLinePC3,1]
	yMinPC3 <- coord[minLinePC3,2]
	xMaxPC3 <- coord[maxLinePC3,1]
	yMaxPC3 <- coord[maxLinePC3,2]
	SlopePC1 <- (yMaxPC1 - yMinPC1)/(xMaxPC1 - xMinPC1)
	InterceptPC1 <- yMinPC1 -SlopePC1*xMinPC1 
	SlopePC2 <- (yMaxPC2 - yMinPC2)/(xMaxPC2 - xMinPC2)
	InterceptPC2 <- yMinPC2 -SlopePC2*xMinPC2
	SlopePC3 <- (yMaxPC3 - yMinPC3)/(xMaxPC3 - xMinPC3)
	InterceptPC3 <- yMinPC3 -SlopePC3*xMinPC3
	matrix1_resultsSlIn[n,1] <- SlopePC1
	matrix1_resultsSlIn[n,2] <- InterceptPC1
	matrix1_resultsSlIn[n,3] <- SlopePC2
	matrix1_resultsSlIn[n,4] <- InterceptPC2
	matrix1_resultsSlIn[n,5] <- SlopePC3
	matrix1_resultsSlIn[n,6] <- InterceptPC3
	matrix1_resultsCoords[n,1] <- xMinPC1
	matrix1_resultsCoords[n,2] <- yMinPC1
	matrix1_resultsCoords[n,3] <- xMaxPC1
	matrix1_resultsCoords[n,4] <- yMaxPC1
	matrix1_resultsCoords[n,5] <- xMinPC2
	matrix1_resultsCoords[n,6] <- yMinPC2
	matrix1_resultsCoords[n,7] <- xMaxPC2
	matrix1_resultsCoords[n,8] <- yMaxPC2
	matrix1_resultsCoords[n,9] <- xMinPC3
	matrix1_resultsCoords[n,10] <- yMinPC3
	matrix1_resultsCoords[n,11] <- xMaxPC3
	matrix1_resultsCoords[n,12] <- yMaxPC3
	
	######################################################################################
	# Interpolated PC values, do this for PC1, 2 and 3, now below is for fit = 3....

    for (r in 1:1400) # size of counterN in the loop below
		{
		matrix2_C_1[r,1] <- 0
		matrix2_C_1[r,2] <- 0
		matrix2_PC_1[r] <- 0
		matrix2_C_2[r,1] <- 0
		matrix2_C_2[r,2] <- 0
		matrix2_PC_2[r] <- 0
		matrix2_C_3[r,1] <- 0
		matrix2_C_3[r,2] <- 0
		matrix2_PC_3[r] <- 0
		}
        
	for (i in 1:3)
		{
		fit<- Krig(coord, objInit$x[,i], theta= 1500, m = 1, Distance="rdist.earth",na.rm=TRUE)
		
		q <- -10.0 # x min
		w <- 11.0 # y min
		counterN <- 0
		while (q <= 170)  #x max
			{
			while (w <= 78)  #y max
				{
				counterN <- counterN + 1
				xnew <- rbind( c( q,w))
				PCvalueNew <- predict( fit, xnew)	
					
				if (i == 1)
					{
                        matrix2_C_1[counterN,1] <- q
                        matrix2_C_1[counterN,2] <- w
                        matrix2_PC_1[counterN] <- PCvalueNew
					}
				if (i == 2)
					{
                        matrix2_C_2[counterN,1] <- q
                        matrix2_C_2[counterN,2] <- w
                        matrix2_PC_2[counterN] <- PCvalueNew
					}
				if (i == 3)
					{
                        matrix2_C_3[counterN,1] <- q
                        matrix2_C_3[counterN,2] <- w
                        matrix2_PC_3[counterN] <- PCvalueNew
					}
				w<-w+0.5
				}
				
			w <- 11.0
			q <- q+0.5
			}
		}
	
	#############
	# calculation of centroids in Asia
	
	## PC1
    q <- -10.0
    w <- 11.0
	counterN <- 0
	PositSum <- 0
	NegatSum <- 0
	PositPCandC_SumLong <- 0
	PositPCandC_SumLat <- 0
	NegatPCandC_SumLong <- 0
	NegatPCandC_SumLat <- 0
	while (q <= 170)  #x
		{
		while (w <= 78)  #y
			{
			counterN <- counterN + 1
			if (matrix2_PC_1[counterN] > 0) # positive value
				{
				PositSum <- PositSum + matrix2_PC_1[counterN]
				PositPCandC_SumLong <- PositPCandC_SumLong + (matrix2_PC_1[counterN]*matrix2_C_1[counterN,1])
				PositPCandC_SumLat <- PositPCandC_SumLat + (matrix2_PC_1[counterN]*matrix2_C_1[counterN,2])
				} else # negative value
				{
				NegatSum <- NegatSum + matrix2_PC_1[counterN]
				NegatPCandC_SumLong <- NegatPCandC_SumLong + (matrix2_PC_1[counterN]*matrix2_C_1[counterN,1])
				NegatPCandC_SumLat <- NegatPCandC_SumLat + (matrix2_PC_1[counterN]*matrix2_C_1[counterN,2])
				}
            w<-w+0.5
			}
		w <- 11.0
		q <- q+0.5
		}
	
	mean_LongPosit_1 <- PositPCandC_SumLong / PositSum
	mean_LatPosit_1 <- PositPCandC_SumLat / PositSum
	mean_LongNegativ_1 <- NegatPCandC_SumLong / NegatSum
	mean_LatNegativ_1 <- NegatPCandC_SumLat / NegatSum
	
	# Gradient line
	SlopeCentroid_1 <- (mean_LatNegativ_1 - mean_LatPosit_1)/(mean_LongNegativ_1 - mean_LongPosit_1)
	InterceptCentroid_1 <- mean_LatPosit_1 -SlopeCentroid_1*mean_LongPosit_1 
	matrix1_resultsSlIn_Cent[n,1] <- SlopeCentroid_1
	matrix1_resultsSlIn_Cent[n,2] <- InterceptCentroid_1
	matrix1_resultsCoords_Cent[n,1] <- mean_LongPosit_1
	matrix1_resultsCoords_Cent[n,2] <- mean_LatPosit_1
	matrix1_resultsCoords_Cent[n,3] <- mean_LongNegativ_1
	matrix1_resultsCoords_Cent[n,4] <- mean_LatNegativ_1
	
	## PC2
    q <- -10.0
    w <- 11.0
	counterN <- 0
	PositSum <- 0
	NegatSum <- 0
	PositPCandC_SumLong <- 0
	PositPCandC_SumLat <- 0
	NegatPCandC_SumLong <- 0
	NegatPCandC_SumLat <- 0
	while (q <= 170)  #x
		{
		while (w <= 78)  #y
			{
			counterN <- counterN + 1
			
			if (matrix2_PC_2[counterN] > 0) # positive value
				{
				PositSum <- PositSum + matrix2_PC_2[counterN]
				PositPCandC_SumLong <- PositPCandC_SumLong + (matrix2_PC_2[counterN]*matrix2_C_2[counterN,1])
				PositPCandC_SumLat <- PositPCandC_SumLat + (matrix2_PC_2[counterN]*matrix2_C_2[counterN,2])
				} else # negative value
				{
				NegatSum <- NegatSum + matrix2_PC_2[counterN]
				NegatPCandC_SumLong <- NegatPCandC_SumLong + (matrix2_PC_2[counterN]*matrix2_C_2[counterN,1])
				NegatPCandC_SumLat <- NegatPCandC_SumLat + (matrix2_PC_2[counterN]*matrix2_C_2[counterN,2])
				}
			w<-w+0.5
			}
		w <- 11.0
		q <- q+0.5
		}

    mean_LongPosit_2 <- PositPCandC_SumLong / PositSum
	mean_LatPosit_2 <- PositPCandC_SumLat / PositSum
	mean_LongNegativ_2 <- NegatPCandC_SumLong / NegatSum
	mean_LatNegativ_2 <- NegatPCandC_SumLat / NegatSum
	
	# Gradient line
	SlopeCentroid_2 <- (mean_LatNegativ_2 - mean_LatPosit_2)/(mean_LongNegativ_2 - mean_LongPosit_2)
	InterceptCentroid_2 <- mean_LatPosit_2 -SlopeCentroid_2*mean_LongPosit_2 
	matrix1_resultsSlIn_Cent[n,3] <- SlopeCentroid_2
	matrix1_resultsSlIn_Cent[n,4] <- InterceptCentroid_2
	matrix1_resultsCoords_Cent[n,5] <- mean_LongPosit_2
	matrix1_resultsCoords_Cent[n,6] <- mean_LatPosit_2
	matrix1_resultsCoords_Cent[n,7] <- mean_LongNegativ_2
	matrix1_resultsCoords_Cent[n,8] <- mean_LatNegativ_2
	
	## PC3
        q <- -10.0
        w <- 11.0
		counterN <- 0
		PositSum <- 0
		NegatSum <- 0
		PositPCandC_SumLong <- 0
		PositPCandC_SumLat <- 0
		NegatPCandC_SumLong <- 0
		NegatPCandC_SumLat <- 0
		while (q <= 170)  #x
			{
			while (w <= 78)  #y
				{
				counterN <- counterN + 1
				if (matrix2_PC_3[counterN] > 0) # positive value
					{
					PositSum <- PositSum + matrix2_PC_3[counterN]
					PositPCandC_SumLong <- PositPCandC_SumLong + (matrix2_PC_3[counterN]*matrix2_C_3[counterN,1])
					PositPCandC_SumLat <- PositPCandC_SumLat + (matrix2_PC_3[counterN]*matrix2_C_3[counterN,2])
					} else # negative value
					{
					NegatSum <- NegatSum + matrix2_PC_3[counterN]
					NegatPCandC_SumLong <- NegatPCandC_SumLong + (matrix2_PC_3[counterN]*matrix2_C_3[counterN,1])
					NegatPCandC_SumLat <- NegatPCandC_SumLat + (matrix2_PC_3[counterN]*matrix2_C_3[counterN,2])
					}
				w<-w+0.5
				}
			w <- 11.0
			q <- q+0.5
			}
		mean_LongPosit_3 <- PositPCandC_SumLong / PositSum
		mean_LatPosit_3 <- PositPCandC_SumLat / PositSum
		mean_LongNegativ_3 <- NegatPCandC_SumLong / NegatSum
		mean_LatNegativ_3 <- NegatPCandC_SumLat / NegatSum
		
		# Gradient line
		SlopeCentroid_3 <- (mean_LatNegativ_3 - mean_LatPosit_3)/(mean_LongNegativ_3 - mean_LongPosit_3)
		InterceptCentroid_3 <- mean_LatPosit_3 -SlopeCentroid_3*mean_LongPosit_3 
		matrix1_resultsSlIn_Cent[n,5] <- SlopeCentroid_3
		matrix1_resultsSlIn_Cent[n,6] <- InterceptCentroid_3
		matrix1_resultsCoords_Cent[n,9] <- mean_LongPosit_3
		matrix1_resultsCoords_Cent[n,10] <- mean_LatPosit_3
		matrix1_resultsCoords_Cent[n,11] <- mean_LongNegativ_3
		matrix1_resultsCoords_Cent[n,12] <- mean_LatNegativ_3
		
	#
	# plot centroids on the map "ResultsCentroid.png"
            # plot centroids on the map "ResultsCentroid0.png"
            
            figureName<-paste(NameInput,"_Centroid.png",sep="")
            png(figureName, width = 1500, height = 300,pointsize = 12)
            par(mfrow = c(1,3))
            
            for (i in 1:3)
            {
                fit<- Krig(coord, objInit$x[,i], theta= 1500, m = 1, Distance="rdist.earth",na.rm=TRUE)
                surface(fit, type ="I", col = colrp, xlab ="", ylab="", las = 1, main = paste("pc", i ), extrap = T)
                map(add=T, interior=F)
                points(coord, pch = 20, cex = 1.5)
                
                for (w in 1:nCoalSims)
                {
                    if (i == 1)
                    {
                        segments(mean_LongPosit_1, mean_LatPosit_1, mean_LongNegativ_1, mean_LatNegativ_1, col= 'pink', lwd = 3.0) # coord x y, x y
                        #abline(a=InterceptCentroid_1, b=SlopeCentroid_1, col = "black") # intercept, slope
                    }
                    if (i == 2)
                    {
                        #segments(mean_LongPosit_2, mean_LatPosit_2, mean_LongNegativ_2, mean_LatNegativ_2, col= 'pink', lwd = 3.0) # coord x y, x y
                        #abline(a=InterceptCentroid_2, b=SlopeCentroid_2, col = "black") # intercept, slope
                    }
                    if (i == 3)
                    {
                        #segments(mean_LongPosit_3, mean_LatPosit_3, mean_LongNegativ_3, mean_LatNegativ_3, col= 'pink', lwd = 3.0) # coord x y, x y
                        #abline(a=InterceptCentroid_3, b=SlopeCentroid_3, col = "black") # intercept, slope
                    }
                }
            }
            dev.off()
		
        #############
        # plot centroids on the map "ResultsCentroid.png"
        #	figureName<-paste(NameInput,"Centroid.png",sep="")
        figureName<-paste(NameInput,"_Basic.png",sep="")
        png(figureName, width = 1500, height = 300, pointsize = 12)
        
        par(mfrow = c(1,3))
        
        for (i in 1:3)
        {
            fit<- Krig(coord, objInit$x[,i], theta= 1500, m = 1, Distance="rdist.earth",na.rm=TRUE)
            surface(fit, type ="I", col = colrp, xlab ="", ylab="", las = 1, main = paste("pc", i ), extrap = T)
            map(add=T, interior=F)
            points(coord, pch = 20, cex = 1.5)
            #		for (w in 1:nCoalSims)
            #			{
            #			if (i == 1)
            #				{
            #				 segments(mean_LongPosit_1, mean_LatPosit_1, mean_LongNegativ_1, mean_LatNegativ_1, col= 'pink') # coord x y, x y
            #				#abline(a=InterceptCentroid_1, b=SlopeCentroid_1, col = "black") # intercept, slope
            #				}
            #			if (i == 2)
            #				{
            #				segments(mean_LongPosit_2, mean_LatPosit_2, mean_LongNegativ_2, mean_LatNegativ_2, col= 'pink') # coord x y, x y
            #				#abline(a=InterceptCentroid_2, b=SlopeCentroid_2, col = "black") # intercept, slope
            #				}
            #			if (i == 3)
            #				{
            #				segments(mean_LongPosit_3, mean_LatPosit_3, mean_LongNegativ_3, mean_LatNegativ_3, col= 'pink') # coord x y, x y
            #				#abline(a=InterceptCentroid_3, b=SlopeCentroid_3, col = "black") # intercept, slope			
            #				}
            #			}
        }
        dev.off()

    ###
	n<-n+1
	}	# end of genetic datasets
n <- 1	
	

### Calculating PC stats among datasets ###
print("\n> Calculating stats..")
if (nCoalSims == 1)
	{
	obj<-cumObj
	print("Only 1 replicate")
	} else 	# the rest of replicates
	{
	obj<-cumObj

	for (a in 1:100) #sdev
			{
			obj$sdev[a] <- cumObj$sdev[a] / nCoalSims
			}
	for (a in 1:100) #rotation
		{
		for (b in 1:100)
			{
			obj$rotation[a,b] <- cumObj$rotation[a,b] / nCoalSims
			}
		}
	for (a in 1:100) #center
		{
		obj$center[a] <- cumObj$center[a] / nCoalSims
		}
	for (a in 1:100) #scale
		{
		obj$scale[a] <- cumObj$scale[a] / nCoalSims
		}
	for (a in 1:1400) #x
		{
		for (b in 1:100)
			{
			obj$x[a,b] <- cumObj$x[a,b] / nCoalSims
			}
		}
	}	
# print(obj$sdev)

# median
meanSlopePC1_MM <- median(matrix1_resultsSlIn[,1], na.rm = TRUE)
meanSlopePC2_MM <- median(matrix1_resultsSlIn[,3], na.rm = TRUE)
meanSlopePC3_MM <- median(matrix1_resultsSlIn[,5], na.rm = TRUE)
meanInterceptPC1_MM <- median(matrix1_resultsSlIn[,2], na.rm = TRUE)
meanInterceptPC2_MM <- median(matrix1_resultsSlIn[,4], na.rm = TRUE)
meanInterceptPC3_MM <- median(matrix1_resultsSlIn[,6], na.rm = TRUE)


############### Average of centroid lines ############
medianSlopePC1 <- median(matrix1_resultsSlIn_Cent[,1], na.rm = TRUE)
medianSlopePC2 <- median(matrix1_resultsSlIn_Cent[,3], na.rm = TRUE)
medianSlopePC3 <- median(matrix1_resultsSlIn_Cent[,5], na.rm = TRUE)
medianInterceptPC1 <- median(matrix1_resultsSlIn_Cent[,2], na.rm = TRUE)
medianInterceptPC2 <- median(matrix1_resultsSlIn_Cent[,4], na.rm = TRUE)
medianInterceptPC3 <- median(matrix1_resultsSlIn_Cent[,6], na.rm = TRUE)
Q1_SlopePC1 <- quantile(matrix1_resultsSlIn_Cent[,1], 0.999, na.rm = TRUE)
Q1_SlopePC2 <- quantile(matrix1_resultsSlIn_Cent[,3], 0.999, na.rm = TRUE)
Q1_SlopePC3 <- quantile(matrix1_resultsSlIn_Cent[,5], 0.999, na.rm = TRUE)
Q1_InterceptPC1 <- quantile(matrix1_resultsSlIn_Cent[,2], 0.999, na.rm = TRUE)
Q1_InterceptPC2 <- quantile(matrix1_resultsSlIn_Cent[,4], 0.999, na.rm = TRUE)
Q1_InterceptPC3 <- quantile(matrix1_resultsSlIn_Cent[,6], 0.999, na.rm = TRUE)
Q2_SlopePC1 <- quantile(matrix1_resultsSlIn_Cent[,1], 0.001, na.rm = TRUE)
Q2_SlopePC2 <- quantile(matrix1_resultsSlIn_Cent[,3], 0.001, na.rm = TRUE)
Q2_SlopePC3 <- quantile(matrix1_resultsSlIn_Cent[,5], 0.001, na.rm = TRUE)
Q2_InterceptPC1 <- quantile(matrix1_resultsSlIn_Cent[,2], 0.001, na.rm = TRUE)
Q2_InterceptPC2 <- quantile(matrix1_resultsSlIn_Cent[,4], 0.001, na.rm = TRUE)
Q2_InterceptPC3 <- quantile(matrix1_resultsSlIn_Cent[,6], 0.001, na.rm = TRUE)
Q3_SlopePC1 <- quantile(matrix1_resultsSlIn_Cent[,1], 0.5, na.rm = TRUE)
Q3_SlopePC2 <- quantile(matrix1_resultsSlIn_Cent[,3], 0.5, na.rm = TRUE)
Q3_SlopePC3 <- quantile(matrix1_resultsSlIn_Cent[,5], 0.5, na.rm = TRUE)
Q3_InterceptPC1 <- quantile(matrix1_resultsSlIn_Cent[,2], 0.5, na.rm = TRUE)
Q3_InterceptPC2 <- quantile(matrix1_resultsSlIn_Cent[,4], 0.5, na.rm = TRUE)
Q3_InterceptPC3 <- quantile(matrix1_resultsSlIn_Cent[,6], 0.5, na.rm = TRUE)

squareroot<-sqrt(nCoalSims) # new
cocient<-1/squareroot # new

sd_SlopePC1 <- sd(matrix1_resultsSlIn_Cent[,1], na.rm = TRUE)
sd_SlopePC2 <- sd(matrix1_resultsSlIn_Cent[,3], na.rm = TRUE)
sd_SlopePC3 <- sd(matrix1_resultsSlIn_Cent[,5], na.rm = TRUE)
sd_InterceptPC1 <- sd(matrix1_resultsSlIn_Cent[,2], na.rm = TRUE)
sd_InterceptPC2 <- sd(matrix1_resultsSlIn_Cent[,4], na.rm = TRUE)
sd_InterceptPC3 <- sd(matrix1_resultsSlIn_Cent[,6], na.rm = TRUE)
se_SlopePC1 <- sd_SlopePC1*cocient
se_SlopePC2 <- sd_SlopePC2*cocient
se_SlopePC3 <- sd_SlopePC3*cocient
se_InterceptPC1 <- sd_InterceptPC1*cocient
se_InterceptPC2 <- sd_InterceptPC2*cocient
se_InterceptPC3 <- sd_InterceptPC3*cocient
se95_SlopePC1 <- se_SlopePC1*1.96 # 95% significance
se95_SlopePC2 <- se_SlopePC2*1.96
se95_SlopePC3 <- se_SlopePC3*1.96
se95_InterceptPC1 <- se_InterceptPC1*1.96
se95_InterceptPC2 <- se_InterceptPC2*1.96
se95_InterceptPC3 <- se_InterceptPC3*1.96
PmedianSlopePC1 <- medianSlopePC1 + se95_SlopePC1
PmedianSlopePC2 <- medianSlopePC2 + se95_SlopePC2
PmedianSlopePC3 <- medianSlopePC3 + se95_SlopePC3
PmedianInterceptPC1 <- medianInterceptPC1 + se95_InterceptPC1
PmedianInterceptPC2 <- medianInterceptPC2 + se95_InterceptPC2
PmedianInterceptPC3 <- medianInterceptPC3 + se95_InterceptPC3
NmedianSlopePC1 <- medianSlopePC1 - se95_SlopePC1
NmedianSlopePC2 <- medianSlopePC2 - se95_SlopePC2
NmedianSlopePC3 <- medianSlopePC3 - se95_SlopePC3
NmedianInterceptPC1 <- medianInterceptPC1 - se95_InterceptPC1
NmedianInterceptPC2 <- medianInterceptPC2 - se95_InterceptPC2
NmedianInterceptPC3 <- medianInterceptPC3 - se95_InterceptPC3

# Plot empty plot for lines (using variable obj)
figureName<-paste("ResultsLinesCentroids",".png",sep="")
png(figureName, width = 1500, height = 300,pointsize = 12)
par(mfrow = c(1,3))

for (i in 1:3)
{
    fit<- Krig(coord, obj$x[,i], theta= 1500, m = 1, Distance="rdist.earth",na.rm=TRUE)
    #surface(fit, type ="I", col = "white", xlab ="", ylab="", las = 1, main = paste("pc", i ), extrap = T)
    surface(fit, type ="I", col = colrp_white, xlab ="", ylab="", las = 1, main = paste("pc", i ), extrap = T)
    map(add=T, interior=F)
    points(coord, pch = 20, cex = 1.5)
    
    for (w in 1:nCoalSims) # all coal sims
    {
        if (i == 1)
        {
            # segments(matrix1_resultsCoords_Cent[w,3], matrix1_resultsCoords_Cent[w,4], matrix1_resultsCoords_Cent[w,1], matrix1_resultsCoords_Cent[w,2], col= 'pink') # coord - OK
            abline(a=matrix1_resultsSlIn_Cent[w,2], b=matrix1_resultsSlIn_Cent[w,1], col = "black", lwd = 2.5) # intercept, slope
        }
        if (i == 2)
        {
            # segments(matrix1_resultsCoords_Cent[w,7], matrix1_resultsCoords_Cent[w,8], matrix1_resultsCoords_Cent[w,5], matrix1_resultsCoords_Cent[w,6], col= 'pink') # coord - OK
            #abline(a=matrix1_resultsSlIn_Cent[w,4], b=matrix1_resultsSlIn_Cent[w,3], col = "black", lwd = 2.5)
        }
        if (i == 3)
        {
            # segments(matrix1_resultsCoords_Cent[w,11], matrix1_resultsCoords_Cent[w,12], matrix1_resultsCoords_Cent[w,9], matrix1_resultsCoords_Cent[w,10], col= 'pink') # coord - OK
            #abline(a=matrix1_resultsSlIn_Cent[w,6], b=matrix1_resultsSlIn_Cent[w,5], col = "black", lwd = 2.5)
        }
    }
    
    # print median
    if (i == 1)
    {
        abline(a=medianInterceptPC1, b=medianSlopePC1, col = "green", lty = "solid", lwd = 4.5) # intercept, slope
        #abline(a=Q1_InterceptPC1, b=Q1_SlopePC1, col = "blue", lty="dashed") # intercept, slope
        #abline(a=Q2_InterceptPC1, b=Q2_SlopePC1, col = "blue", lty="dashed") # intercept, slope
        ##abline(a=Q3_InterceptPC1, b=Q3_SlopePC1, col = "blue", lty="dashed") # intercept, slope
        #abline(a=PmedianInterceptPC1, b=PmedianSlopePC1, col = "red", lty="dashed") # intercept, slope
        #abline(a=NmedianInterceptPC1, b=NmedianSlopePC1, col = "red", lty="dashed") # intercept, slope
    }
    if (i == 2)
    {
        #abline(a=medianInterceptPC2, b=medianSlopePC2, col = "green", lty = "solid", lwd = 4.5) # intercept, slope
        #abline(a=Q1_InterceptPC2, b=Q1_SlopePC2, col = "blue", lty="dashed") # intercept, slope
        #abline(a=Q2_InterceptPC2, b=Q2_SlopePC2, col = "blue", lty="dashed") # intercept, slope
        ##abline(a=Q3_InterceptPC2, b=Q3_SlopePC2, col = "blue", lty="dashed") # intercept, slope
        #abline(a=PmedianInterceptPC2, b=PmedianSlopePC2, col = "red", lty="dashed") # intercept, slope
        #abline(a=NmedianInterceptPC2, b=NmedianSlopePC2, col = "red", lty="dashed") # intercept, slope
    }
    if (i == 3)
    {
        #abline(a=medianInterceptPC3, b=medianSlopePC3, col = "green", lty = "solid", lwd = 4.5) # intercept, slope
        #abline(a=Q1_InterceptPC3, b=Q1_SlopePC3, col = "blue", lty="dashed") # intercept, slope
        #abline(a=Q2_InterceptPC3, b=Q2_SlopePC3, col = "blue", lty="dashed") # intercept, slope
        ##abline(a=Q3_InterceptPC3, b=Q3_SlopePC3, col = "blue", lty="dashed") # intercept, slope
        #abline(a=PmedianInterceptPC3, b=PmedianSlopePC3, col = "red", lty="dashed") # intercept, slope
        #abline(a=NmedianInterceptPC3, b=NmedianSlopePC3, col = "red", lty="dashed") # intercept, slope
    }
}
dev.off()






