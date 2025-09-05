suppressMessages(library(DESeq2, quietly=TRUE))
suppressMessages(library(stringr))
#suppressMessages(library("BiocParallel"))
suppressMessages(library("Hmisc"))
#register(MulticoreParam(16))

calc_dds_nobatch <- function(cts, coldata ){
dds <- DESeqDataSetFromMatrix(countData = round(cts),
							  colData = coldata,
							  design= ~condition)
dds <- suppressMessages(DESeq(dds))
	return(dds)
}

calc_dds_batch <- function(cts, coldata, r_batch ) {
dds <- DESeqDataSetFromMatrix(countData = round(cts), # エラーが出るためroundを追加した
							  colData = coldata,
							  design= ~batch + condition)
dds <- suppressMessages(DESeq(dds))
 return(dds)
}

make_coldata <- function(cts, condition, batch, ref_group){
coldata <- data.frame(factor(condition) )# conditionのデータフレームを作る
rownames(coldata) <- colnames(cts)
colnames(coldata) <-'condition'
# cellをつける
coldata['cell'] <- colnames(cts)
if (batch[1] != 'No batch') {
	coldata['batch'] <- batch
	}
coldata$condition <- relevel(coldata$condition, ref = ref_group)
return(coldata)
}



calc_deseq <- function(dds, coldata, independentFiltering=TRUE, deseq2=FALSE) {
res <- results(dds, independentFiltering = independentFiltering)

cat('Computing rlog...\n')
rld <- rlog(dds, blind=FALSE)

res_table <- data.frame(assay(rld))
unshrunk_table <- data.frame(assay(rld))

row_name <-rownames(res_table)

# combinationsの組み合わせを計算
comb <- combn(x=levels(coldata[["condition"]]),m=2)
ncomb <-  dim(comb)[2]

sink("DESeq2_output.txt")
for (i in c(1: ncomb)) {

		cat('=======================================================\n')
		var1 <- comb[,i][2]
		var2 <- comb[,i][1]

			if (deseq2) { # deseq2のオリジナル表記
				comp_name <- str_replace(i, pattern="condition_", replacement="")
				} else { # homer的に逆転
				comp_name <- paste(var2, var1, sep="_vs_")
				}

			cat(comp_name)
			res <- results(dds, contrast=c("condition", var1, var2), alpha=0.05)
			print(summary(res))
			res <- results(dds, contrast=c("condition", var1, var2), independentFiltering = independentFiltering)

			cat('Filter Threshold: ')
			cat(metadata(res)$filterThreshold)
			# contrastを使うとapeglmが使えなくなる。
			# 改善法はhttp://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators

			resLFC <- suppressMessages(lfcShrink(dds, contrast=c("condition", var1, var2), res=res, type = 'ashr'))
			res_table[paste(comp_name, 'log2FC',sep='.')] <- resLFC['log2FoldChange'][[1]]
			unshrunk_table[paste(comp_name, 'log2FC',sep='.')] <- res['log2FoldChange'][[1]]
			res_table[paste(comp_name, 'pvalue',sep='.')] <- resLFC['pvalue'][[1]]
			res_table[paste(comp_name, 'adj.pvalue',sep='.')] <- resLFC['padj'][[1]]
			unshrunk_table[paste(comp_name, 'pvalue',sep='.')] <- res['pvalue'][[1]]
			unshrunk_table[paste(comp_name, 'adj.pvalue',sep='.')] <- res['padj'][[1]]
		}

sink()

res_table['Gene'] <- row.names(res_table)
meta <- metadata(res)
return_res <- Hmisc::llist(res_table, unshrunk_table, meta)

return(return_res)
}


calc_sva_n <- function(dds, condition) {
suppressMessages(library(sva))
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
cat('=======================================================\n')
cat('Compute svseq to estimate the number of latent factors.\n')


svseq <- try(svaseq(dat, mod, mod0))

if (class(svseq) == "try-error") {
	cat("!!!!svseq error. nSV is set as 2.\n")
	svn <- 2
} else {
svn <- svseq$n.sv # SVAでの変数の数)}
}

cat("\n")
cat('=======================================================\n')
sva_pre_res <- Hmisc::llist(svn, dat, mod, mod0)
}


calc_svseq <- function(sva_pre_res, coldata, dds, return_res, type) {
n.sv <- sva_pre_res$svn
svseq <- try(svaseq(sva_pre_res$dat, sva_pre_res$mod, sva_pre_res$mod0, n.sv))
# エラーがないときだけ処理
if (class(svseq) == "try-error") {
	cat('SVA error. Skip this step.\n\n\n')
	devoff()
} else {
comb <- combn(x=levels(coldata[["condition"]]),m=2)
ncomb <-  dim(comb)[2]
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (j in colnames(coldata)) {
	for (i in 1:1) {
  stripchart(svseq$sv[, i] ~ coldata[[j]], vertical = TRUE, main = paste(j, paste0("SV", i)))
  abline(h = 0)
 }}

ddssva <- dds


factors = c()
for (i in 1:n.sv){
	print(paste0('SV',as.character(i)))
	factors <- c(factors, paste0('SV',as.character(i)))
	ddssva[[paste0('SV',as.character(i))]] <-  svseq$sv[,i]
}

design(ddssva) <- as.formula(paste(paste("~", paste(factors, collapse=" + ")), "+ condition", collapse = ''))

ddssva <- suppressMessages(DESeq(ddssva))
# 全ての条件について表示する
res_table = return_res$res_table
unshrunk_table = return_res$unshrunk_table
sink("SVA_output.txt")
for (i in c(1: ncomb)) {

		cat('=======================================================\n')
		var1 <- comb[,i][2]
		var2 <- comb[,i][1]

			if (deseq2) { # deseq2のオリジナル表記
				comp_name <- str_replace(i, pattern="condition_", replacement="")
				} else { # homer的に逆転
				comp_name <- paste(var2, var1, sep="_vs_")
				comp_name <- paste('SVA1', comp_name, sep = '_')
				}

			cat(comp_name)
			res <- results(ddssva, contrast=c("condition", var1, var2), alpha=0.05)

			print(summary(res))
			res <- results(ddssva, contrast=c("condition", var1, var2), independentFiltering = independentFiltering)
			# contrastを使うとapeglmが使えなくなる。
			# 改善法はhttp://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
			resLFC <- suppressMessages(lfcShrink(ddssva, contrast=c("condition", var1, var2), res=res, type = type))
			res_table[paste(comp_name, 'log2FC',sep='.')] <- resLFC['log2FoldChange'][[1]]
			unshrunk_table[paste(comp_name, 'log2FC',sep='.')] <- res['log2FoldChange'][[1]]
			res_table[paste(comp_name, 'pvalue',sep='.')] <- resLFC['pvalue'][[1]]
			res_table[paste(comp_name, 'adj.pvalue',sep='.')] <- resLFC['padj'][[1]]
			unshrunk_table[paste(comp_name, 'pvalue',sep='.')] <- res['pvalue'][[1]]
			unshrunk_table[paste(comp_name, 'adj.pvalue',sep='.')] <- res['padj'][[1]]
		}
sink()
return_res <- Hmisc::llist(res_table, unshrunk_table)

}
return(return_res)
}