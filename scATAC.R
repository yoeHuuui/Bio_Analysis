library(SnapATAC)
library(GenomicRanges)

x.sp <- createSnap(file = "SRR12340627.snap", sample="SRR12340627")
barcodes = read.csv("singlecell.csv",head=TRUE)
barcodes = barcodes[2:nrow(barcodes),]
promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1)
UMI = log(barcodes$passed_filters+1, 10)
data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio)
barcodes$promoter_ratio = promoter_ratio
barcodes$UMI <- UMI

barcodes.sel <- barcodes[(barcodes$UMI >=3 & barcodes$UMI <= 5),]
rownames(barcodes.sel) = barcodes.sel$barcode
x.sp <- x.sp[which(x.sp@barcode %in% barcodes.sel$barcode),]

x.sp = addBmatToSnap(x.sp, bin.size=10000, num.cores=5)
x.sp = makeBinary(x.sp, mat="bmat")

chr.exclude = seqlevels(x.sp@feature)[grep("Mt|Pt", seqlevels(x.sp@feature))]
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature)
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]}

bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)

hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
)
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
x.sp = x.sp[, idy, mat="bmat"]

x.sp = runDiffusionMaps(obj=x.sp, input.mat="bmat", num.eigs=30)
plotDimReductPW(obj=x.sp, 
                eigs.dims=1:30, 
                point.size=0.3, 
                point.color="grey", 
                point.shape=19, 
                point.alpha=0.6,
                down.sample=5000, 
                pdf.file.name=NULL, 
                pdf.height=7, 
                pdf.width=7)

x.sp = runKNN(obj=x.sp, eigs.dims=1:20, k=15)
x.sp=runCluster(obj=x.sp, tmp.folder=tempdir(), louvain.lib="R-igraph", seed.use=10)
x.sp@metaData$cluster = x.sp@cluster
x.sp = runViz(obj=x.sp, tmp.folder=tempdir(), dims=2, 
  eigs.dims=1:20, method="Rtsne", seed.use=10)

plotViz(
  obj=x.sp,
  method="tsne", 
  main="10X ATAC",
  point.color=x.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=5000,
  legend.add=FALSE
)

system("which snaptools")
system("which macs2")

library(parallel)
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 200)]
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("SRR1240627.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/Bioinfo/yhsun/miniconda3/envs/py38/bin/snaptools",
    path.to.macs="/Bioinfo/yhsun/miniconda3/envs/py38/bin/macs2",
    gsize=1.0e8,
    buffer.size=500, 
    num.cores=3,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir())
}, mc.cores=5)

peaks.names = system("ls | grep narrowPeak", intern=TRUE)
peak.gr.ls = lapply(peaks.names, function(x){
    peak.df = read.table(x)
    GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
  })
peak.gr = reduce(Reduce(c, peak.gr.ls))
peaks.df = as.data.frame(peak.gr)[,1:3]

write.table(peaks.df,file = "peaks.combined.bed",append=FALSE,
        quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
        row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
        fileEncoding = "")
