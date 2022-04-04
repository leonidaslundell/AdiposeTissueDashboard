library(data.table)
#create the plots so that the render time speeds up
load("Adipose gene expression analysis/data/results.Rdata")|> try(silent = T)
resSig <- lapply(resSig, as.data.table)
res <- lapply(res, as.data.table)

res <- res[c("NGT_PREvsPOST", "NGT_PREvsREC", "T2D_PREvsPOST", "T2D_PREvsREC")]
resSig <- resSig[c("NGT_PREvsPOST", "NGT_PREvsREC", "T2D_PREvsPOST", "T2D_PREvsREC")]

#create average expression

cpmsAll <- removeBatchEffect(x = cpm(dataDGE), batch = dataDGE$samples$sex)
rownames(cpmsAll) <- res$NGT_PREvsPOST$Symbol
cpmsAll <- cpmsAll[!is.na(res$NGT_PREvsPOST$Symbol),]

cpms <- sapply(list(NGT_PREvsPOST = c("NGT_PRE", "NGT_POST"),
                    NGT_PREvsREC = c("NGT_PRE", "NGT_REC"),
                    T2D_PREvsPOST = c("T2D_PRE", "T2D_POST"),
                    T2D_PREvsREC = c("T2D_PRE", "T2D_REC")), \(i){
  rowMeans(cpmsAll[,dataDGE$samples$super_group %in% i])
}) |> as.data.table(keep.rownames = T)

colnames(cpms) <- c("Symbol", "NGT_PREvsPOST", "NGT_PREvsREC", "T2D_PREvsPOST", "T2D_PREvsREC")
setkey(cpms, "Symbol")

res <- lapply(res, \(x){
  x[!is.na(Symbol),c("Symbol", "logFC", "adj.P.Val")]
})
datMD <- cpms[match(res$NGT_PREvsPOST$Symbol, cpms$Symbol),]

for(x in names(res)){
  datMD[,paste0(x, ".logFC")] <- res[[x]]$logFC
  datMD[,paste0(x, ".sig")] <- symnum(res[[x]]$adj.P.Val, 
                                     cutpoints = c(1,0.05,0), 
                                     c("*", "ns")) |> as.factor()
  # datMD[abs(res[[x]]$logFC)<log2(1.1), c(x, paste0(x, ".logFC"), paste0(x, ".sig"))] <- NA
  datMD[res[[x]]$adj.P.Val>0.2, paste0(x, ".sig")] <- NA
}
keep <- (datMD[,-1] |> is.na() |> rowSums()) < 12
datMD <- datMD[keep,]

####################################################
#significance

datRes <- melt(lapply(res, \(x) x[,c("Symbol", "adj.P.Val")]))[,-2]

datRes$group1 <- gsub("vs.*","", datRes$L1)
datRes$group2 <- gsub("_.*vs","", datRes$L1)
datRes$group2 <- gsub("NGT","NGT_", datRes$group2)
datRes$group2 <- gsub("T2D","T2D_", datRes$group2)

datRes$Diagnosis <- gsub("_.*","", datRes$L1)
datRes$value <- symnum(datRes$value, cutpoints = c(1,0.05,0.01,0.001,0), c("***", "**", "*", "ns"))
datRes <- as.data.table(datRes)

###################################################
# points for overall plot
datBox <- t(cpmsAll) |> as.data.frame()
rownames(datBox) == dataDGE$samples$seq.id
datBox$id <- as.factor(dataDGE$samples$id)
datBox$Exercise <- as.factor(dataDGE$samples$exercise)
datBox$Diagnosis <- as.factor(dataDGE$samples$status)

datBox <- reshape2::melt(datBox) |> as.data.table()
datBox$Exercise <- factor(datBox$Exercise, levels = c("PRE", "POST", "REC"))
datBox$Group <- paste0(datBox$Diagnosis,"_",datBox$Exercise)
datBox$Group <- factor(datBox$Group, levels = c( "NGT_PRE", "NGT_POST", "NGT_REC", "T2D_PRE", "T2D_POST", "T2D_REC"))
colnames(datBox)[4] <- "Symbol"

#artifially low values from the remove batch effect
datBox[value<0, value := 0.01] 

####################################################
#save sig for overwriting
save(list = c("datMD", "datBox", "datRes"), file = "~/OneDrive - University of Copenhagen/SharedData/Adipose/Adipose gene expression analysis/interactiveDash/dat.Rdata")
