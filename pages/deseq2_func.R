suppressMessages(library(DESeq2, quietly=TRUE))
suppressMessages(library(stringr))
library(apeglm)
#suppressMessages(library("BiocParallel"))
suppressMessages(library("Hmisc"))
#register(MulticoreParam(16))

calc_dds_nobatch <- function(){
  try(dev.off())
  coldata$condition <- relevel(coldata$condition, ref = ref_group)
  dds <<- DESeqDataSetFromMatrix(countData = round(cts),
                                 colData = coldata,
                                 design= ~condition)
  
  if (use_custom_size_factors && !is.null(custom_size_factors)) {
  	print("use cutsom size factors")
    sizeFactors(dds) <<- custom_size_factors
  } else {
    dds <<- estimateSizeFactors(dds)
  }
  
  dds <<- suppressMessages(DESeq(dds))
  png(paste0(res_dir,'/DispersionEstimates.png'))
  plotDispEsts(dds)
  dev.off()
}

calc_dds_batch <- function() {
  try(dev.off())
  coldata$condition <- relevel(coldata$condition, ref = ref_group)
  dds <<- DESeqDataSetFromMatrix(countData = round(cts),
                                 colData = coldata,
                                 design= ~batch + condition)
  
  if (use_custom_size_factors && !is.null(custom_size_factors)) {
  	print("use cutsom size factors")
    sizeFactors(dds) <<- custom_size_factors
  } else {
    dds <<- estimateSizeFactors(dds)
  }
  
  dds <<- suppressMessages(DESeq(dds))
  png(paste0(res_dir,'/DispersionEstimates.png'))
  plotDispEsts(dds)
  dev.off()
}


calc_dds_LRT <- function() {
  sink()
  sink(paste0(res_dir, "/DESeq2_output.txt"))



  # Make sure full_model and reduced_model are formulas, not strings
  if (is.character(full_model)) {
    full_model <<- as.formula(full_model)
    cat("Converted full_model from string to formula\n")
  }
  
  if (is.character(reduced_model)) {
    reduced_model <<- as.formula(reduced_model)
    cat("Converted reduced_model from string to formula\n")
  }

  
  for (i in c(1:dim(coldata)[2])) {
    coldata[,i] <- factor(coldata[,i])
  }
  
  cat('coldata')
  print(coldata)
  
  # 連続変数を処理
  if (exists("continuous_vars") && length(continuous_vars) > 0) {
    cat("Processing continuous variables...\n")
    
    for (col_name in continuous_vars) {
      if (col_name %in% colnames(coldata)) {
        # 連続変数として処理
        cat(paste0("Treating '", col_name, "' as continuous variable\n"))
        
        # 数値に変換
        if (!is.numeric(coldata[[col_name]])) {
          original_values <- coldata[[col_name]]
          # 数字だけを抽出して数値に変換
          numeric_values <- as.numeric(gsub("[^0-9.]", "", as.character(original_values)))
          
          if (any(is.na(numeric_values))) {
            warning(paste0("Cannot convert '", col_name, "' to numeric. Using as factor."))
            coldata[[col_name]] <- factor(coldata[[col_name]])
          } else {
            coldata[[col_name]] <- numeric_values
            cat(paste0("Converted '", col_name, "' to numeric values: ", 
                     paste(head(numeric_values), collapse=", "), "...\n"))
          }
        } else {
          cat(paste0("'", col_name, "' is already numeric\n"))
        }
      }
    }
  }
  
  # 残りの変数を因子型に変換（連続変数以外）
  for (i in c(1:dim(coldata)[2])) {
    col_name <- colnames(coldata)[i]
    if (!exists("continuous_vars") || !(col_name %in% continuous_vars)) {
      cat(paste0("Treating '", col_name, "' as categorical variable\n"))
      coldata[,i] <- factor(coldata[,i])
    }
  }
  
  cat('coldata\n')
  print(coldata)
  
  # ポリノミアル設定の情報表示（情報表示のみで式の変更はしない）
  if (exists("add_polynomial") && add_polynomial && exists("polynomial_variable")) {
    if (exists("use_poly_function")) {
      use_poly_func <- use_poly_function
    } else {
      use_poly_func <- TRUE
    }
    
    if (use_poly_func) {
      cat("Polynomial analysis enabled using poly() function\n")
    } else {
      cat("Polynomial analysis enabled using I() function (explicit powers)\n")
    }
    
    cat("Time variable: ", polynomial_variable, "\n")
    cat("Polynomial degree: ", polynomial_degree, "\n")
    
    # poly()関数使用時の追加設定
    raw_param <- ""
    if (use_poly_func && exists("use_raw") && use_raw) {
      cat("Using raw polynomials (raw=TRUE)\n")
      raw_param <- ", raw=TRUE"
    } else if (use_poly_func) {
      cat("Using orthogonal polynomials (raw=FALSE)\n")
    }
    
    # 変数の型変換を行う（式の変更はしない）
    if (polynomial_variable %in% colnames(coldata) && !is.numeric(coldata[[polynomial_variable]])) {
      original_values <- coldata[[polynomial_variable]]
      numeric_values <- as.numeric(gsub("[^0-9.]", "", as.character(original_values)))
      
      if (!any(is.na(numeric_values))) {
        coldata[[polynomial_variable]] <- numeric_values
        cat("Converted", polynomial_variable, "to numeric values\n")
      }
    }
    
    # poly()関数使用時のユニークポイント数チェック
    if (use_poly_func && polynomial_variable %in% colnames(coldata)) {
      unique_points <- length(unique(coldata[[polynomial_variable]]))
      cat("Number of unique time points:", unique_points, "\n")
      cat("Requested polynomial degree:", polynomial_degree, "\n")
      
      if (polynomial_degree >= unique_points) {
        cat("ERROR: Polynomial degree (", polynomial_degree, ") must be less than unique points (", unique_points, ")\n")
        cat("Maximum polynomial degree for your data: ", unique_points - 1, "\n")
        stop(paste0("Polynomial degree (", polynomial_degree, ") is too high. Maximum allowed: ", unique_points - 1))
      }
    }
  }
  
  # 最終的なモデル式の表示
  cat("Full model: \n")
  print(full_model)
  cat('\nReduced model: \n')
  print(reduced_model)
  
  cat('cts\n')
  print(head(cts))
  cat('coldata\n')
  print(coldata)
  
  cat('coldata final\n')
  print(coldata)
  
  sink()
  
  # DESeq2の実行
  dds <<- DESeqDataSetFromMatrix(countData = round(cts),
                               colData = coldata,
                               design = full_model)
  dds <<- DESeq(dds, test="LRT", reduced = reduced_model)
  png(paste0(res_dir,'/DispersionEstimates.png'))
  plotDispEsts(dds)
  dev.off()
  res <- results(dds, independentFiltering = independentFiltering)
  cat('FDR<0.05: \n')
  print(table(res$padj<0.05))
  cat('\nResults:')
  df <- DataFrame(res@listData, row.names = row.names(cts))
  write.table(df, paste0(res_dir, "/DESeq2_LRT_res.tsv"), quote=FALSE, row.names = TRUE, col.names = NA, sep = '\t')
  print(head(res))
  sink()
}


make_coldata <- function(){
coldata <<- data.frame(factor(condition) )# conditionのデータフレームを作る global変数へ

#cts <<- cts # これらはglobal変数
#condition <<- condition
#batch <<- r_batch
#ref_group <<- ref_group

rownames(coldata) <<- colnames(cts)
colnames(coldata) <<-'condition'
# cellをつける
coldata['cell'] <<- colnames(cts)
if (batch[1] != 'No batch') {
	coldata['batch'] <<- batch
	}
coldata$condition <<- relevel(coldata$condition, ref = ref_group)
#return(coldata)
}


make_coldata2 <- function(){ # rlog計算用
coldata <<- data.frame(factor(condition) )# conditionのデータフレームを作る global変数へ
rownames(coldata) <<- colnames(cts)
colnames(coldata) <<-'condition'
}


calc_deseq <- function() {
try(dev.off())
res <- results(dds, independentFiltering = independentFiltering)

if (rld_calc){
	rld <<- rlog(dds, blind=FALSE)
	res_table <<- data.frame(assay(rld)) #これらはglobal変数
	unshrunk_table <<- data.frame(assay(rld))
} else {

	res_table <<- data.frame(row.names = rownames(cts)) #これらはglobal変数
	unshrunk_table <<- data.frame(row.names = rownames(cts))
}
row_name <<-rownames(res_table)

# combinationsの組み合わせを計算
comb <<- combn(x=levels(coldata[["condition"]]),m=2) # global変数
ncomb <<-  dim(comb)[2]

sink(paste0(res_dir,"/DESeq2_output.txt"))
		#	cat('Filtering Threshold: ')
		#	cat(attr(res,"filterThreshold"))
cat('\n')
cat('independentFiltering: ')
cat(independentFiltering)
cat('\n')


if (independentFiltering)
{
	cat('% of genes filtered out\nFiltering threshold for mean count:\n')
	print( metadata(res)$filterThreshold)
	cat('\n')

}

# apeglmのとき
if (type == 'apeglm') {
	combinations <- resultsNames(dds)
	combinations <- combinations[str_detect(combinations, 'condition')]
		for ( i in combinations) {
			comp_name <- sub('condition_','', i)
				cat(comp_name )
				cat('\n')

			if (deseq2) { # deseq2のオリジナル表記
				} else { # homer的に逆転
				first <- sub('_vs_', '', str_extract(comp_name, '.*_vs_')) # vsの前
				second <- sub('_vs_', '', str_extract(comp_name, '_vs_.*')) # vsの後
				comp_name <- paste(second, first, sep="_vs_")
				}
			cat(comp_name)
			res <- results(dds, name=i, alpha = results_alpha,independentFiltering = independentFiltering)
			print(summary(res))
			resLFC <- suppressMessages(lfcShrink(dds,coef=i, type="apeglm"))
			res_table[paste(comp_name, 'log2FC',sep='.')] <<- resLFC['log2FoldChange'][[1]]
			unshrunk_table[paste(comp_name, 'log2FC',sep='.')] <<- res['log2FoldChange'][[1]]
			res_table[paste(comp_name, 'pvalue',sep='.')] <<- resLFC['pvalue'][[1]]
			res_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- resLFC['padj'][[1]]
			res_table[paste(comp_name, 'stat',sep='.')] <<- resLFC['stat'][[1]]
			unshrunk_table[paste(comp_name, 'pvalue',sep='.')] <<- res['pvalue'][[1]]
			unshrunk_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- res['padj'][[1]]
			unshrunk_table[paste(comp_name, 'stat',sep='.')] <<- res['stat'][[1]]
			png(paste0(res_dir,'/MA.', comp_name,'.shrunken.png'))
			plotMA(resLFC, ylim=c(-2,2))
			dev.off()

		}

	} else {
	for (i in c(1: ncomb)) {

			cat('---------------------------------\n')
			var1 <- comb[,i][2]
			var2 <- comb[,i][1]

			if ( ref_in & !( ref_group %in% c(var1, var2)) ) { # refereceが指定されていて、それが含まれないとき
					next
				}

				if (deseq2) { # deseq2のオリジナル表記
					comp_name <- paste(var1, var2, sep="_vs_")
					} else { # homer的に逆転
					comp_name <- paste(var2, var1, sep="_vs_")
					}
				cat(comp_name )
				cat('\n')
				res <- results(dds, contrast=c("condition", var1, var2), alpha=results_alpha)
				print(summary(res))
				res <- results(dds, contrast=c("condition", var1, var2), independentFiltering = independentFiltering)
				# contrastを使うとapeglmが使えなくなる。
				# 改善法はhttp://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
				#write.table(res,'rest.csv')

				resLFC <- suppressMessages(lfcShrink(dds, contrast=c("condition", var1, var2), res=res, type = type))
				res_table[paste(comp_name, 'log2FC',sep='.')] <<- resLFC['log2FoldChange'][[1]]
				unshrunk_table[paste(comp_name, 'log2FC',sep='.')] <<- res['log2FoldChange'][[1]]
				res_table[paste(comp_name, 'pvalue',sep='.')] <<- resLFC['pvalue'][[1]]
				res_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- resLFC['padj'][[1]]
				unshrunk_table[paste(comp_name, 'pvalue',sep='.')] <<- res['pvalue'][[1]]
				unshrunk_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- res['padj'][[1]]
				unshrunk_table[paste(comp_name, 'stat', sep='.')] <<- res['stat'][[1]]
				png(paste0(res_dir,'/MA.', comp_name,'.png'))
				plotMA(res, ylim=c(-2,2))
				dev.off()
				try(dev.off())
				png(paste0(res_dir,'/MA.', comp_name,'.shrunken.png'))
				plotMA(resLFC, ylim=c(-2,2))
				dev.off()
			}
}
sink()

meta <- metadata(res)
return_res <- Hmisc::llist(res_table, unshrunk_table, meta)
saveRDS(res_table, paste0(temp_dir, '/res_table.RDS'))
write.table(res_table, paste0(res_dir,"/DESeq2_res_LFC-shrunk.tsv"),  quote=FALSE, row.names = TRUE, col.names = NA, sep = '\t')
write.table(unshrunk_table, paste0(res_dir, "/DESeq2_res.tsv"),  quote=FALSE,  row.names = TRUE, col.names = NA, sep = '\t')
#return(return_res)
}


calc_sva_n <- function() {
suppressMessages(library(sva))
dat  <<- counts(dds, normalized = TRUE) #すべてglobal変数へ
idx  <<- rowMeans(dat) > 1
dat  <<- dat[idx, ]
mod  <<- model.matrix(~ condition, colData(dds))
mod0 <<- model.matrix(~   1, colData(dds))

svseq <- try(svaseq(dat, mod, mod0))

if (class(svseq) == "try-error") {
	cat("!!!!svseq error. nSV is set as 2.\n")
	svn <<- 2
} else {
svn <<- svseq$n.sv # SVAでの変数の数)} # global 変数
}
return(svseq$n.sv)
}

calc_sva_n_both <- function() {
  suppressMessages(library(sva))
  
  # データは既にcalc_sva_nで準備されているものを使用
  if (!exists("dat") || !exists("mod")) {
    # calc_sva_nが実行されていない場合は実行
    calc_sva_n()
  }
  
  # BE法（既存の結果を使用）
  n.sv.be <- svn
  
  # Leek法で再計算
  n.sv.leek <- NA
  tryCatch({
    # svaseqでLeek法を使用
    svseq_leek <- svaseq(dat, mod, mod0, n.sv = NULL, numSVmethod = "leek")
    n.sv.leek <- svseq_leek$n.sv
  }, error = function(e) {
    # エラーの場合はnum.sv関数を使用
    tryCatch({
      n.sv.leek <- num.sv(dat, mod, method = "leek")
    }, error = function(e2) {
      cat("Leek method failed:", e2$message, "\n")
      n.sv.leek <- NA
    })
  })
  
  # グローバル変数に保存（Python側で取得用）
  svn_be <<- n.sv.be
  svn_leek <<- n.sv.leek
  
  # 結果を表示
  cat("\n=== SV Estimation Results ===\n")
  cat("BE method (default):", n.sv.be, "\n")
  cat("Leek method:", ifelse(is.na(n.sv.leek), "Failed", n.sv.leek), "\n")
  
  if (!is.na(n.sv.leek) && n.sv.be > n.sv.leek * 2) {
    cat("\nWARNING: BE estimate is much higher than Leek estimate.\n")
    cat("Consider using a more conservative approach.\n")
  }
  
  # リストで返す
  return(list(be = n.sv.be, leek = n.sv.leek))
}

calc_svseq <- function() {
try(dev.off())
tryCatch({

pdf(paste0(res_dir,"/SVA_graph.pdf"))
sva_res_table <<- res_table
sva_unshrunk_table <<- unshrunk_table
#	comb <- combn(x=levels(coldata[["condition"]]),m=2) # global変数にする
#	ncomb <-  dim(comb)[2]
sink(paste0(res_dir,"/SVA_output.txt"))
#　ここからループ

if (sva_calc) {
 svn <<-2
 cat("Calc only 2 surrogate vaiables.")
 } #defaultは2まで
cat('\n\n')

for (x in c(1:svn)) {

	svseq <- suppressMessages(try(svaseq(dat, mod, mod0, x)))
	# エラーがないときだけ処理
	if (class(svseq) == "try-error") {
		cat('SVA error. Skip this step.\n')
		cat('SVA_n :')
		cat(x)
		cat('\n')
		try( dev.off() )
	} else {

	par(mfrow = c(2, 1), mar = c(3,5,3,1))
	for (j in colnames(coldata)) {
		for (i in 1:x) {
	  stripchart(svseq$sv[, i] ~ coldata[[j]], vertical = TRUE, main = paste(j, paste0("SV", i)))
	  abline(h = 0)
	 }}

	ddssva <- dds

	factors = c()
	for (i in 1:x){
		print(paste0('SV',as.character(i)))
		factors <- c(factors, paste0('SV',as.character(i)))
		ddssva[[paste0('SV',as.character(i))]] <-  svseq$sv[,i]
	}

	design(ddssva) <- as.formula(paste(paste("~", paste(factors, collapse=" + ")), "+ condition", collapse = ''))
	cat("design = ")
	cat(paste(paste("~", paste(factors, collapse=" + ")), "+ condition", collapse = ''))
	cat('\n')
	ddssva <- suppressMessages(DESeq(ddssva))

	png(paste0(res_dir,'/SV', as.character(x), '_DispersionEstimates.png'))
	plotDispEsts(ddssva)
	try(dev.off())

	# apeglmのとき
	if (type == 'apeglm') {
		combinations <- resultsNames(ddssva)
		combinations <- combinations[str_detect(combinations, 'condition')]
			for ( i in combinations) {
				comp_name <- sub('condition_','', i)

				if (deseq2) { # deseq2のオリジナル表記
					} else { # homer的に逆転
					first <- sub('_vs_', '', str_extract(comp_name, '.*_vs_')) # vsの前
					second <- sub('_vs_', '', str_extract(comp_name, '_vs_.*')) # vsの後
					comp_name <- paste(second, first, sep="_vs_")
					}
				comp_name <- paste(paste0('SVA', as.character(x)), comp_name, sep = '_')
				cat(comp_name)
				res <- results(ddssva, name=i, alpha = results_alpha, independentFiltering = independentFiltering)
				print(summary(res))
				resLFC <- suppressMessages(lfcShrink(ddssva,coef=i, type="apeglm"))
				sva_res_table[paste(comp_name, 'log2FC',sep='.')] <<- resLFC['log2FoldChange'][[1]]
				sva_unshrunk_table[paste(comp_name, 'log2FC',sep='.')] <<- res['log2FoldChange'][[1]]
				sva_res_table[paste(comp_name, 'pvalue',sep='.')] <<- resLFC['pvalue'][[1]]
				sva_res_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- resLFC['padj'][[1]]
				sva_unshrunk_table[paste(comp_name, 'pvalue',sep='.')] <<- res['pvalue'][[1]]
				sva_unshrunk_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- res['padj'][[1]]
				sva_unshrunk_table[paste(comp_name, 'stat',sep='.')] <<- res['stat'][[1]]
				png(paste0(res_dir,'/MA.', comp_name,'.png'))
				plotMA(res, ylim=c(-2,2))
				png(paste0(res_dir,'/MA.', comp_name,'.shrunken.png'))
				plotMA(resLFC, ylim=c(-2,2))
				try(dev.off())
			}

	} else {

	# 全ての条件について表示する
	for (i in c(1: ncomb)) {
			cat('---------------------------------\n')
			var1 <- comb[,i][2]
			var2 <- comb[,i][1]

				if (deseq2) { # deseq2のオリジナル表記
					comp_name <- paste(var1, var2, sep="_vs_")
					} else { # homer的に逆転
					comp_name <- paste(var2, var1, sep="_vs_")
					}
				comp_name <- paste(paste0('SVA', as.character(x)), comp_name, sep = '_')

				cat(comp_name)
				res <- results(ddssva, contrast=c("condition", var1, var2), alpha=results_alpha)

				cat(summary(res))
				res <- results(ddssva, contrast=c("condition", var1, var2), independentFiltering = independentFiltering)
				# write.table(res,paste0(comp_name,'.csv'))
				# contrastを使うとapeglmが使えなくなる。
				# 改善法はhttp://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
				resLFC <- suppressMessages(lfcShrink(ddssva, contrast=c("condition", var1, var2), res=res, type = type))
				sva_res_table[paste(comp_name, 'log2FC',sep='.')] <<- resLFC['log2FoldChange'][[1]]
				sva_unshrunk_table[paste(comp_name, 'log2FC',sep='.')] <<- res['log2FoldChange'][[1]]
				sva_res_table[paste(comp_name, 'pvalue',sep='.')] <<- resLFC['pvalue'][[1]]
				sva_res_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- resLFC['padj'][[1]]
				sva_unshrunk_table[paste(comp_name, 'pvalue',sep='.')] <<- res['pvalue'][[1]]
				sva_unshrunk_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- res['padj'][[1]]
				sva_unshrunk_table[paste(comp_name, 'stat',sep='.')] <<- res['stat'][[1]]
				tryCatch({
				png(paste0(res_dir,'/MA.', comp_name,'.png'))
				plotMA(res, ylim=c(-2,2))

				png(paste0(res_dir,'/MA.', comp_name,'.shrunken.png'))
				plotMA(resLFC, ylim=c(-2,2))
				try(dev.off())
			}, error = function(e) {
				  print(paste("SVA MA plot error :  ", e))
				})


			}
	}
		if (x < svn) {
			cat('===============================================\n')
		}
}	}#ループの終わり
sink()
#SVA_return_res <- Hmisc::llist(sva_res_table, sva_unshrunk_table)
write.table(sva_res_table, paste0(res_dir,"/SVA_res_LFC-shrunk.tsv"),  quote=FALSE, row.names = TRUE, col.names = NA, sep = '\t')
write.table(sva_unshrunk_table, paste0(res_dir,"/SVA_res.tsv"),  quote=FALSE,  row.names = TRUE, col.names = NA, sep = '\t')
try(dev.off())

}, error = function(e) {
  print(paste("SVA calc error :  ", e))
})


}


calc_ruvseq <- function()
{
suppressMessages(library("RUVSeq"))

ruv_res_table <<- res_table
ruv_unshrunk_table <<- unshrunk_table
#ruv_res_table = return_res$res_table
#ruv_unshrunk_table = return_res$unshrunk_table

# non significant genes の決定
not.sig <- c()
counter <- 1
for (i in resultsNames(dds) ) {
	if (counter ==1){
		counter <- 2 # 最初はパスする
	} else {
		if (str_count(i, ref_group) > 0) { # ref_groupを含む組み合わせについて
			res <- results(dds, name=i)
			NS <- rownames(res)[which(res$pvalue > RUV.alpha)]
			if (length(not.sig) == 0) {
				not.sig <- NS
				} else {
				not.sig <- intersect(not.sig, NS)
				}
		}
	}
	}
cat('-------------\n')
cat('The number of control genes for RUVg: ')
cat(length(not.sig))
cat('\n-------------\n')

# small sample size (n=6)では、k=1
if (length(colnames(cts)) < 7) {
	RUV_n <- 1
} else { RUV_n <- 3}
pdf(paste0(res_dir,"/RUV_graph.pdf"))
sink(paste0(res_dir,'/RUVseq.txt'))
#comb <- combn(x=levels(coldata[["condition"]]),m=2)
#ncomb <-  dim(comb)[2]

for (x in c(1:RUV_n)) {
		set <- newSeqExpressionSet(counts(dds))
	idx  <- rowSums(counts(set) > 5) >= 2
	set  <- set[idx, ]
	set <- betweenLaneNormalization(set, which="upper")
	#not.sig <- rownames(res)[which(res$pvalue > .1)]
	empirical <- rownames(set)[ rownames(set) %in% not.sig ]
	set <- try(RUVg(set, empirical, k=x))
	# エラーがないときだけ処理
	if (class(set) == "try-error") {
		cat('RUV error. Skip this step.\n')
		cat('RUV_n :')
		cat(x)
		cat('\n')
		try(dev.off())
	} else {
	pData(set)

	par(mfrow = c(2, 1), mar = c(3,5,3,1))
	for (j in colnames(coldata)) {
		for (i in 1:x) {
	  stripchart(pData(set)[, i] ~ coldata[[j]], vertical = TRUE, main = paste(j, paste0("RUV", i)))
	  abline(h = 0)
	 }}


	ddsruv <- dds
	ddsruv$W1 <- set$W_1

	if (x>1) {
		ddsruv$W2 <- set$W_2
	}
	if (x>2) {
		ddsruv$W3 <- set$W_3
	}

	factors = c()
	for (i in 1:x){
		factors <- c(factors, paste0('W',as.character(i)))
	}

	design(ddsruv) <- as.formula(paste(paste("~", paste(factors, collapse=" + ")), "+ condition", collapse = ''))
	cat("design = ")
	cat(paste(paste("~", paste(factors, collapse=" + ")), "+ condition", collapse = ''))
	cat('\n')

	ddsruv <- DESeq(ddsruv)

	png(paste0(res_dir,'/RUV', as.character(x), '_DispersionEstimates.png'))
	plotDispEsts(ddsruv)
	dev.off()

	resruv <- results(ddsruv, independentFiltering = independentFiltering)

	# apeglmのとき
	if (type == 'apeglm') {
		combinations <- resultsNames(ddsruv)
		combinations <- combinations[str_detect(combinations, 'condition')]
			for ( i in combinations) {
				comp_name <- sub('condition_','', i)

				if (deseq2) { # deseq2のオリジナル表記
					} else { # homer的に逆転
					first <- sub('_vs_', '', str_extract(comp_name, '.*_vs_')) # vsの前
					second <- sub('_vs_', '', str_extract(comp_name, '_vs_.*')) # vsの後
					comp_name <- paste(second, first, sep="_vs_")
					}
				comp_name <- paste(paste0('RUV', as.character(x)), comp_name, sep = '_')

				cat(comp_name)
				res <- results(ddsruv, name=i, alpha = results_alpha,independentFiltering = independentFiltering)
				print(summary(res))
				resLFC <- suppressMessages(lfcShrink(ddsruv,coef=i, type="apeglm"))
				res_table[paste(comp_name, 'log2FC',sep='.')] <<- resLFC['log2FoldChange'][[1]]
				unshrunk_table[paste(comp_name, 'log2FC',sep='.')] <<- res['log2FoldChange'][[1]]
				res_table[paste(comp_name, 'pvalue',sep='.')] <<- resLFC['pvalue'][[1]]
				res_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- resLFC['padj'][[1]]
				unshrunk_table[paste(comp_name, 'pvalue',sep='.')] <<- res['pvalue'][[1]]
				unshrunk_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- res['padj'][[1]]
				unshrunk_table[paste(comp_name, 'stat',sep='.')] <<- res['stat'][[1]]
				png(paste0(res_dir,'/MA.', comp_name,'.shrunken.png'))
				plotMA(resLFC, ylim=c(-2,2))
				dev.off()
			}

	} else {

	# 全ての条件について表示する

	for (i in c(1: ncomb)) {

			cat('---------------------------------\n')
			var1 <- comb[,i][2]
			var2 <- comb[,i][1]

				if (deseq2) { # deseq2のオリジナル表記
					comp_name <- paste(var1, var2, sep="_vs_")
					} else { # homer的に逆転
					comp_name <- paste(var2, var1, sep="_vs_")
					}

				comp_name <- paste(paste0('RUV', as.character(x)), comp_name, sep = '_')

				cat(comp_name)
				res <- results(ddsruv, contrast=c("condition", var1, var2), alpha=results_alpha)
				print(summary(res))
				res <- results(ddsruv, contrast=c("condition", var1, var2), independentFiltering = independentFiltering)
				# contrastを使うとapeglmが使えなくなる。
				# 改善法はhttp://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
				resLFC <- suppressMessages(lfcShrink(ddsruv, contrast=c("condition", var1, var2), res=res, type = type))
				res_table[paste(comp_name, 'log2FC',sep='.')] <<- resLFC['log2FoldChange'][[1]]
				unshrunk_table[paste(comp_name, 'log2FC',sep='.')] <<- res['log2FoldChange'][[1]]
				res_table[paste(comp_name, 'pvalue',sep='.')] <<- resLFC['pvalue'][[1]]
				res_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- resLFC['padj'][[1]]
				unshrunk_table[paste(comp_name, 'pvalue',sep='.')] <<- res['pvalue'][[1]]
				unshrunk_table[paste(comp_name, 'adj.pvalue',sep='.')] <<- res['padj'][[1]]
				unshrunk_table[paste(comp_name, 'stat',sep='.')] <<- res['stat'][[1]]
				png(paste0(res_dir,'/MA.', comp_name,'.shrunken.png'))
				plotMA(resLFC, ylim=c(-2,2))
				dev.off()
				}
	}
			if (x < RUV_n) {
			cat('===============================================\n')
	}	}
} # loopの最後

#RUV_return_res <- Hmisc::llist(ruv_res_table, ruv_unshrunk_table)
write.table(res_table, paste0(res_dir,"/RUV_res_LFC-shrunk.tsv"),  quote=FALSE, row.names = TRUE, col.names = NA, sep = '\t')
write.table(unshrunk_table, paste0(res_dir,"/RUV_res.tsv"),  quote=FALSE,  row.names = TRUE, col.names = NA, sep = '\t')

dev.off()
sink()
}




calc_rlog_no_group <- function() {
dds <<- DESeqDataSetFromMatrix(countData = round(cts),
							  colData = coldata,
							  design= ~condition) # global変数
rld <<- rlog(dds, blind=TRUE)
write.table(data.frame(assay(rld)), paste0(temp_dir, "/rld.tsv"),  quote=FALSE, row.names = TRUE, col.names = NA, sep = '\t')
}


calc_rlog <- function() {
dds <<- DESeqDataSetFromMatrix(countData = round(cts),
							  colData = coldata,
							  design= ~condition) # global変数
rld <<- rlog(dds, blind=FALSE)
write.table(data.frame(assay(rld)), paste0(temp_dir, "/rld.tsv"),  quote=FALSE, row.names = TRUE, col.names = NA, sep = '\t')
}


calc_vst_no_group <- function() {
dds <<- DESeqDataSetFromMatrix(countData = round(cts),
							  colData = coldata,
							  design= ~condition) # global変数
rld <<- vst(dds, blind=TRUE)
write.table(data.frame(assay(rld)), paste0(temp_dir, "/rld.tsv"),  quote=FALSE, row.names = TRUE, col.names = NA, sep = '\t')
}


calc_vst <- function() {
dds <<- DESeqDataSetFromMatrix(countData = round(cts),
							  colData = coldata,
							  design= ~condition) # global変数
rld <<- vst(dds, blind=FALSE)
write.table(data.frame(assay(rld)), paste0(temp_dir, "/rld.tsv"),  quote=FALSE, row.names = TRUE, col.names = NA, sep = '\t')
}



calc_svseq_LRT <- function() {
sink()

# 関数の最初に追加するコード
if (!dir.exists(res_dir)) {
  dir.create(res_dir, recursive = TRUE)
  cat("Created results directory:", res_dir, "\n")
}

sink(paste0(res_dir,'/SVA_output.txt'))
for (i in c(1:dim(coldata)[2])) {
coldata[,i] <- factor(coldata[,i])
                      }



  # Check if polynomial terms are requested
  if (exists("add_polynomial") && add_polynomial && exists("polynomial_variable")) {
    # 実装方法の確認
    use_poly_func <- TRUE
    if (exists("use_poly_function")) {
      use_poly_func <- use_poly_function
    }
    
    if (use_poly_func) {
      cat("Polynomial analysis enabled using poly() function\n")
    } else {
      cat("Polynomial analysis enabled using I() function (explicit powers)\n")
    }
    
    cat("Time variable: ", polynomial_variable, "\n")
    cat("Polynomial degree: ", polynomial_degree, "\n")
    
    # poly()関数使用時の追加設定
    raw_param <- ""
    if (use_poly_func && exists("use_raw") && use_raw) {
      cat("Using raw polynomials (raw=TRUE)\n")
      raw_param <- ", raw=TRUE"
    } else if (use_poly_func) {
      cat("Using orthogonal polynomials (raw=FALSE)\n")
    }
    
    # Check if the time variable exists in coldata
    if (polynomial_variable %in% colnames(coldata)) {
      # Convert time variable to numeric if it's not already
      if (!is.numeric(coldata[[polynomial_variable]])) {
        original_values <- coldata[[polynomial_variable]]
        # Try to extract numeric part
        numeric_values <- as.numeric(gsub("[^0-9.]", "", as.character(original_values)))
        
        # Check if conversion was successful
        if (any(is.na(numeric_values))) {
          warning("Cannot convert time variable to numeric. Using original values.")
        } else {
          # Store original values and set numeric values
          coldata[[paste0(polynomial_variable, "_original")]] <- original_values
          coldata[[polynomial_variable]] <- numeric_values
          cat("Converted", polynomial_variable, "to numeric:\n")
          for (i in 1:length(original_values)) {
            cat(paste0(original_values[i], " -> ", numeric_values[i], "\n"))
          }
        }
      }
    }



     # ----- 式の処理方法を修正 -----
formula_str <- deparse(full_model)
cat("Original formula string: ", formula_str, "\n")

if (use_poly_func) {
  # poly()関数を使った実装
  if (polynomial_degree >= 2) {
    # 時間変数を検索するより単純な方法
    # 単に末尾に項を追加する方法を採用
    if (raw_param == "") {
      poly_term <- paste0(" + poly(", polynomial_variable, ", degree=", polynomial_degree, ")")
    } else {
      poly_term <- paste0(" + poly(", polynomial_variable, ", degree=", polynomial_degree, ", raw=TRUE)")
    }
    
    formula_str <- paste0(formula_str, poly_term)
    cat("Modified formula string with poly(): ", formula_str, "\n")
  }
} else {
  # I()関数を使った実装（明示的なべき乗）
  if (polynomial_degree >= 2) {
    # 線形項が含まれているかどうかを確認
    linear_pattern <- paste0("([^a-zA-Z0-9_]|^)", polynomial_variable, "([^a-zA-Z0-9_]|$)")
    if (!grepl(linear_pattern, formula_str)) {
      # 線形項がない場合は追加
      formula_str <- paste0(formula_str, " + ", polynomial_variable)
    }
    
    # 2次の項を追加
    quadratic_term <- paste0(" + I(", polynomial_variable, "^2)")
    formula_str <- paste0(formula_str, quadratic_term)
    
    if (polynomial_degree >= 3) {
      # 3次の項を追加
      cubic_term <- paste0(" + I(", polynomial_variable, "^3)")
      formula_str <- paste0(formula_str, cubic_term)
    }
    
    cat("Modified formula string with I(): ", formula_str, "\n")
  }
}


# デバッグ用に式の文字列を出力
cat("Final formula string: ", formula_str, "\n")

# 文字列からフォーミュラに変換
tryCatch({
  full_model <- as.formula(formula_str)
  cat("Successfully converted to formula\n")
}, error = function(e) {
  cat("Error converting string to formula: ", e$message, "\n")
  cat("Using original formula\n")
})
} # polyのとき

full_model <- as.formula(full_model)
reduced_model <- as.formula(reduced_model)
         cat('full model: ')
         print(full_model)
         cat('\nreduced_model: ')
         print(reduced_model)

cat('coldata\n')
print(coldata)
pdf(paste0(res_dir,"/SVA_graph.pdf"))


if (sva_calc) {
 svn <<-2
 cat("Calc only 2 surrogate vaiables.")
 } #defaultは2まで
cat('\n\n')
#　ここからループ
for (x in c(1:svn)) {

cat(paste0('=============================================SVA',str(x),'\n'))
	svseq <- suppressMessages(try(svaseq(dat, mod, mod0, x)))
	# エラーがないときだけ処理
	if (class(svseq) == "try-error") {
		cat('SVA error. Skip this step.\n')
		cat('SVA_n :')
		cat(x)
		cat('\n')
		dev.off()
	} else {

	par(mfrow = c(2, 1), mar = c(3,5,3,1))
	for (j in colnames(coldata)) {
		for (i in 1:x) {
	  stripchart(svseq$sv[, i] ~ coldata[[j]], vertical = TRUE, main = paste(j, paste0("SV", i)))
	  abline(h = 0)
	 }}

	sva_coldata <- coldata

	factors = c()
	for (i in 1:x){
		print(paste0('SV',as.character(i)))
		factors <- c(factors, paste0('SV',as.character(i)))
		sva_coldata[[paste0('SV',as.character(i))]] <-  svseq$sv[,i]
	}

	library(formula.tools)


	sva_full_model <- as.formula(paste0(as.character(full_model), ' + ', paste(factors, collapse=" + ")))
	sva_reduced_model <- as.formula(paste0(as.character(reduced_model), ' + ', paste(factors, collapse=" + ")))
	cat('sva_full_model: \n')
	print(sva_full_model)
	cat('\nsva_reduced_model: \n')
	print(sva_reduced_model)


dds <<- DESeqDataSetFromMatrix(countData = round(cts), # エラーが出るためroundを追加した
							  colData = sva_coldata,
							  design= sva_full_model)
dds <<- DESeq(dds, test="LRT", reduced = sva_reduced_model )
png(paste0(res_dir,'/SVA',as.character(x), '_DispersionEstimates.png'))
plotDispEsts(dds)
dev.off()
res <- results(dds, independentFiltering = independentFiltering)
cat('FDR<0.05: \n')
print(	table(res$padj<0.05) )
cat('\nResults:')
df <- DataFrame(res@listData, row.names = row.names(cts))
write.table(df, paste0(res_dir,"/SVA", as.character(x), "_DESeq2_LRT_res.tsv"),  quote=FALSE,  row.names = TRUE, col.names = NA, sep = '\t')
write.table(df, paste0(res_dir,"/SVA_LRT_res.tsv"), quote=FALSE, row.names = TRUE, col.names = NA, sep = '\t')

print(head(res))

		if (x < svn) {
			cat('===============================================\n')
		}

}
}	#ループの終わり
sink()
dev.off()
}



#==============追加関数
# limma eBayesの実装
calc_limma <- function() {
  suppressMessages(library(limma))
  
  if (limma_count) {
    suppressMessages(library(edgeR))
  }
  
  sink(paste0(res_dir, "/limma_output.txt"))
  cat("Running limma eBayes analysis\n")
  
  # データと実験デザインの準備
  cts_matrix <- as.matrix(cts)
  
  # ファクターとしての条件を設定
  coldata$condition <- factor(coldata$condition)
  coldata$condition <- relevel(coldata$condition, ref = ref_group)
  
  # デザイン行列の作成（~0 + conditionスタイル）
  if (batch[1] != 'No batch') {
    coldata$batch <- factor(coldata$batch)
    design <- model.matrix(~0 + condition + batch, data=coldata)
    colnames(design) <- gsub("^condition", "", colnames(design))
  } else {
    design <- model.matrix(~0 + condition, data=coldata)
    colnames(design) <- gsub("^condition", "", colnames(design))
  }
  
  cat("Design matrix column names:\n")
  print(colnames(design))
  
  # データタイプに応じた処理
  if (limma_count) {
    # RNA-seqカウント
    dge <- DGEList(counts = cts_matrix)
    dge <- calcNormFactors(dge)
    
    tryCatch({
      png(paste0(res_dir, "/voom_plot.png"))
      v <- voom(dge, design, plot=TRUE)
      dev.off()
      fit <- lmFit(v, design)
    }, error = function(e) {
      cat("Error in voom:", e$message, "\n")
      cat("Using log2 counts instead\n")
      v <- log2(dge$counts + 0.5)
      fit <- lmFit(v, design)
    })
  } else if (apply_logit) {
    # ロジット変換
    eps <- 1e-6
    cts_matrix <- pmax(pmin(cts_matrix, 1-eps), eps)
    cts_matrix <- log(cts_matrix/(1-cts_matrix))
    fit <- lmFit(cts_matrix, design)
  } else {
    # 通常データ
    fit <- lmFit(cts_matrix, design)
  }
  
  # eBayes適用（trend, robustパラメータを使用）
  fit <- eBayes(fit, trend=limma_trend, robust=limma_robust)
  
  # 結果テーブル初期化
  all_results <- data.frame(row.names = rownames(cts))
  
  # コントラスト行列の作成と比較の実行
  groups <- levels(coldata$condition)
  cat("Creating contrasts for groups:", paste(groups, collapse=", "), "\n")
  
  # リファレンスグループ以外の各グループとの比較
  for (group in groups) {
    if (group != ref_group) {
      # コントラスト作成
      contrast_name <- paste0(group, "-", ref_group)
      cat("Creating contrast:", contrast_name, "\n")
      
      tryCatch({
        contrast.matrix <- makeContrasts(contrasts=contrast_name, levels=design)
        cat("Contrast matrix created successfully\n")
        
        # コントラスト適用
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2, trend=limma_trend, robust=limma_robust)
        
        # 結果抽出
        top_table <- topTable(fit2, number=Inf)
        
        # 比較名の形式設定
        if (deseq2) {
          comp_name <- paste(group, ref_group, sep="_vs_")
        } else {
          comp_name <- paste(ref_group, group, sep="_vs_")
        }
        
        # DESeq2のスタイルに合わせて結果を保存
  all_results[paste(comp_name, 'log2FC', sep='.')] <- top_table$logFC[match(rownames(all_results), rownames(top_table))]
  all_results[paste(comp_name, 'pvalue', sep='.')] <- top_table$P.Value[match(rownames(all_results), rownames(top_table))]
  all_results[paste(comp_name, 'adj.pvalue', sep='.')] <- top_table$adj.P.Val[match(rownames(all_results), rownames(top_table))]
        
        # MA-plot作成
        png(paste0(res_dir, "/MA_", comp_name, ".png"))
        limma::plotMA(fit2, main=comp_name)
        dev.off()
        
      }, error = function(e) {
        cat("Error processing contrast", contrast_name, ":", e$message, "\n")
      })
    }
  }
  
  # 結果保存
  write.table(all_results, paste0(res_dir, "/limma_res.tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
  
  cat("limma eBayes analysis completed\n")
  sink()
}

#
calc_gam <- function() {
  suppressMessages(library(mgcv))
  suppressMessages(library(parallel))
  
  sink(paste0(res_dir, "/gam_output.txt"))
  cat("Running GAM analysis with", dist_short, "distribution\n")
  
  # グローバル環境から変数を確実に取得
  cts <- get("cts", envir = .GlobalEnv)
  cts_matrix <- as.matrix(cts)
  cat("Dimensions of cts_matrix:", dim(cts_matrix), "\n")
  
  # coldataとref_groupの確認
  cat("Checking coldata structure:\n")
  print(head(coldata))
  cat("Reference group:", ref_group, "\n")
  
  # 分布ファミリーに応じたデータ前処理
  if (dist_short == "beta") {
    # 0-1の境界調整
    cts_matrix <- pmax(pmin(cts_matrix, 1-epsilon), epsilon)
  } else if (dist_short == "poisson" || dist_short == "nb") {
    # カウントデータの整数化
    cts_matrix <- round(cts_matrix)
  }
  
  # 分布ファミリーの設定
  get_family <- function() {
    if (dist_short == "beta") {
      return(betar())
    } else if (dist_short == "gaussian") {
      return(gaussian())
    } else if (dist_short == "poisson") {
      return(poisson())
    } else if (dist_short == "nb") {
      return(negbin(theta = nb_theta))
    }
  }
  
  # 並列クラスターのセットアップ
  cl <- makeCluster(n_cores)
  
  # 各ワーカーノードで必要なパッケージを読み込む
  clusterEvalQ(cl, {
    library(mgcv)
    NULL  # 明示的な戻り値
  })
  
  # 関数内で使用する全ての変数を確実にエクスポート
  clusterExport(cl, c("cts_matrix", "coldata", "dist_short", "use_batch", "get_family", "ref_group", "deseq2"), 
                envir = environment())
  
  if (dist_short == "beta") {
    clusterExport(cl, "epsilon", envir = environment())
  } else if (dist_short == "nb") {
    clusterExport(cl, "nb_theta", envir = environment())
  }
  
  # coldataの因子型状態を保持するため、追加設定
  clusterEvalQ(cl, {
    if (exists("coldata")) {
      coldata$condition <- factor(coldata$condition, levels = levels(coldata$condition))
      if ("batch" %in% names(coldata)) {
        coldata$batch <- factor(coldata$batch, levels = levels(coldata$batch))
      }
    }
    NULL
  })
  
  # 遺伝子ごとの処理関数
  process_gene <- function(i) {
    # データフレームの準備
    gene_data <- cts_matrix[i,]
    df <- data.frame(y=gene_data, condition=coldata$condition)
    if (use_batch && "batch" %in% colnames(coldata)) {
      df$batch <- coldata$batch
    }
    
    # モデル式の構築
    if (use_batch && "batch" %in% colnames(df)) {
      formula_str <- "y ~ condition + batch"
    } else {
      formula_str <- "y ~ condition"
    }
    
    # 分布ファミリーの取得
    family_to_use <- get_family()
    
    # モデル適合を試行
    tryCatch({
      # GAMモデルのフィット
      model <- gam(as.formula(formula_str), family=family_to_use, data=df, method="REML")
      
      # 係数の抽出
      coef_summary <- summary(model)$p.table
      
      # 各条件の比較結果を抽出
      results <- c()
      
      # リファレンス以外の条件に対する結果を抽出
      for (cond in levels(df$condition)) {
        if (cond != ref_group) {
          coef_name <- paste0("condition", cond)
          if (coef_name %in% rownames(coef_summary)) {
            # ログフォールドチェンジ、p値などを抽出
            log2FC <- coef_summary[coef_name, "Estimate"]
            pvalue <- coef_summary[coef_name, "Pr(>|t|)"]
            
            # 結果を追加
            results <- c(results, log2FC, pvalue)
          }
        }
      }
      
      return(results)
      
    }, error=function(e) {
      # エラーが発生した場合はNAを返す
      return(rep(NA, 2*(length(levels(df$condition))-1)))
    })
  }
  
  # 並列処理を実行
  cat("Processing", nrow(cts_matrix), "genes with", n_cores, "cores\n")
  results_list <- parLapply(cl, 1:nrow(cts_matrix), process_gene)
  
  # クラスターを終了
  stopCluster(cl)
  
  # 結果を整形
  all_results <- data.frame(row.names = rownames(cts))
  
  # 全ての条件の組み合わせについて結果を格納
  comps <- c()
  for (other_cond in setdiff(levels(coldata$condition), ref_group)) {
    # DESeq2スタイルかHomerスタイルかによって名前を調整
    if (deseq2) {
      comp_name <- paste(other_cond, ref_group, sep="_vs_")
    } else {
      comp_name <- paste(ref_group, other_cond, sep="_vs_")
    }
    comps <- c(comps, comp_name)
    
    # インデックスの計算
    cond_idx <- which(levels(coldata$condition) == other_cond) - 1
    
    # 結果の抽出と保存
    log2FC_values <- sapply(results_list, function(x) {
      if (length(x) >= 2*cond_idx) return(x[2*cond_idx-1]) else return(NA)
    })
    pvalue_values <- sapply(results_list, function(x) {
      if (length(x) >= 2*cond_idx) return(x[2*cond_idx]) else return(NA)
    })
    
    all_results[paste(comp_name, 'log2FC', sep='.')] <- log2FC_values
    all_results[paste(comp_name, 'pvalue', sep='.')] <- pvalue_values
    
    # 多重検定補正
    all_results[paste(comp_name, 'adj.pvalue', sep='.')] <- p.adjust(pvalue_values, method="BH")
  }
  
  cat("Processed comparisons:", paste(comps, collapse=", "), "\n")
  
  # 結果を保存
  write.table(all_results, paste0(res_dir, "/glm_", dist_short, "_res.tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
  
  cat("GLM analysis completed\n")
  sink()
  
  return(all_results)
}

# Beta Regressionの実装
calc_betareg <- function() {
  suppressMessages(library(betareg))
  suppressMessages(library(lmtest))
  suppressMessages(library(parallel))
  
  sink(paste0(res_dir, "/betareg_output.txt"))
  cat("Running Beta Regression analysis\n")

    # 行列に変換
  cts_matrix <- as.matrix(cts)
  cat("Dimensions of cts_matrix:", dim(cts_matrix), "\n")
  
  # 0-1の境界調整
  cts_matrix <- pmax(pmin(cts_matrix, 1-epsilon), epsilon)
  
# クラスターセットアップと変数エクスポート部分を修正
cl <- makeCluster(n_cores)

# 各ワーカーノードで必要なパッケージを読み込む
clusterEvalQ(cl, {
  library(betareg)
  library(lmtest)
  NULL  # 明示的な戻り値
})

# 関数内で使用する全ての変数を確実にエクスポート
# coldata内の変数も正しく処理されるようにする
clusterExport(cl, c("cts_matrix", "coldata", "epsilon", "use_batch", "ref_group", "deseq2"), 
              envir = environment())

# coldataの因子型状態を保持するため、追加設定
clusterEvalQ(cl, {
  if (exists("coldata")) {
    coldata$condition <- factor(coldata$condition, levels = levels(coldata$condition))
    if ("batch" %in% names(coldata)) {
      coldata$batch <- factor(coldata$batch, levels = levels(coldata$batch))
    }
  }
  NULL
})
  # 全遺伝子を並列処理する関数
  process_gene <- function(i) {
    # データフレームの準備
    gene_data <- cts_matrix[i,]
    df <- data.frame(y=gene_data, coldata)
    
    # モデル式の構築
    if (use_batch && "batch" %in% colnames(df)) {
      formula_str <- "y ~ condition + batch"
    } else {
      formula_str <- "y ~ condition"
    }
    
    # モデル適合を試行
    tryCatch({
      model <- betareg(as.formula(formula_str), data=df)
      
      # 係数の抽出
      coef_summary <- summary(model)$coefficients$mean
      
      # 各条件の比較結果を抽出
      results <- c()
      
      # リファレンス以外の条件に対する結果を抽出
      for (cond in levels(df$condition)) {
        if (cond != ref_group) {
          coef_name <- paste0("condition", cond)
          if (coef_name %in% rownames(coef_summary)) {
            # ログフォールドチェンジ、p値、標準誤差などを抽出
            log2FC <- coef_summary[coef_name, "Estimate"]
            pvalue <- coef_summary[coef_name, "Pr(>|z|)"]
            
            # 結果を追加
            results <- c(results, 
                        log2FC=log2FC, 
                        pvalue=pvalue)
          }
        }
      }
      
      return(results)
      
    }, error=function(e) {
      # エラーが発生した場合はNAを返す
      return(rep(NA, 2*(length(levels(df$condition))-1)))
    })
  }
  
  # 並列処理を実行
  cat("Processing", nrow(cts_matrix), "genes with", n_cores, "cores\n")
  results_list <- parLapply(cl, 1:nrow(cts_matrix), process_gene)
  
  # クラスターを終了
  stopCluster(cl)
  
  # 結果を整形
  all_results <- data.frame(row.names = rownames(cts))
  
  # 全ての条件の組み合わせについて結果を格納
  for (ref_cond in ref_group) {
    for (other_cond in setdiff(levels(coldata$condition), ref_cond)) {
      # DESeq2スタイルかHomerスタイルかによって名前を調整
      if (deseq2) {
        comp_name <- paste(other_cond, ref_cond, sep="_vs_")
      } else {
        comp_name <- paste(ref_cond, other_cond, sep="_vs_")
      }
      
      # インデックスの計算
      cond_idx <- which(levels(coldata$condition) == other_cond)
      if (cond_idx > which(levels(coldata$condition) == ref_cond)) {
        cond_idx <- cond_idx - 1
      }
      
      # 結果の抽出と保存
      idx_offset <- 2 * (cond_idx - 1)
      log2FC_values <- sapply(results_list, function(x) x[idx_offset + 1])
      pvalue_values <- sapply(results_list, function(x) x[idx_offset + 2])
      
      all_results[paste(comp_name, 'log2FC', sep='.')] <- log2FC_values
      all_results[paste(comp_name, 'pvalue', sep='.')] <- pvalue_values
      
      # 多重検定補正
      all_results[paste(comp_name, 'adj.pvalue', sep='.')] <- p.adjust(pvalue_values, method="BH")
    }
  }
  
  # 結果を保存
  write.table(all_results, paste0(res_dir, "/betareg_res.tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
  
  cat("Beta Regression analysis completed\n")
  sink()
  
  return(all_results)
}
