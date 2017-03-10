### R code from vignette source 'rasp.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: rasp.Rnw:60-61
###################################################
library(rasp)


###################################################
### code chunk number 2: rasp.Rnw:66-68
###################################################
data(YRI)
dim(YRI)


###################################################
### code chunk number 3: rasp.Rnw:73-76
###################################################
genderYRI <- attr(YRI, "gender")
head(genderYRI)
table(genderYRI)


###################################################
### code chunk number 4: rasp.Rnw:93-94
###################################################
plotTernary(YRI, "ENSG00000160741", transcripts=1:3)


###################################################
### code chunk number 5: rasp.Rnw:105-106
###################################################
plotTernary(YRI, "ENSG00000160741", transcripts=1:3)


###################################################
### code chunk number 6: rasp.Rnw:120-122
###################################################
plotAllIso("ENSG00000005448", data = YRI, group = genderYRI,
                      inds = 1:10)


###################################################
### code chunk number 7: rasp.Rnw:148-151
###################################################
gene <- YRI[YRI$gene_id=="ENSG00000160741", 3:ncol(YRI)]
mod <- testRasp(gene, genderYRI)
mod


###################################################
### code chunk number 8: rasp.Rnw:158-159
###################################################
plotTernary(YRI, "ENSG00000160741", transcripts=1:3, group=genderYRI)


###################################################
### code chunk number 9: rasp.Rnw:172-176
###################################################
ngenes <- 5
sel.genes <-  names(table(YRI[,1]))[1:ngenes]
res <- rasp(YRI[YRI[,1]%in%sel.genes, ], genderYRI, mc.cores = 1, 
            geneidCol = 1, expressionCols = 3:ncol(YRI))


###################################################
### code chunk number 10: rasp.Rnw:183-184
###################################################
print(res)


###################################################
### code chunk number 11: rasp.Rnw:191-196
###################################################
set.seed(1234)
g <- as.character(genderYRI)
g[g == "female"] <- ifelse(runif(sum(g == "female"))<.5, "female1", "female2")
resG <- rasp(YRI[YRI[,1]%in%sel.genes, ], factor(g), mc.cores = 1,
             geneidCol = 1, expressionCols = 3:ncol(YRI), testGroup = TRUE)


###################################################
### code chunk number 12: rasp.Rnw:199-200
###################################################
print(resG)


###################################################
### code chunk number 13: rasp.Rnw:209-210
###################################################
sessionInfo()


