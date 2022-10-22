##################################################################
#################### Quantitative DE analysis ####################

library("BiocStyle")
library("DEP")
library("dplyr")

my_data <- as.data.frame(read.table("intensities.txt", header = TRUE, sep = '\t'))

my_data$Gene_names %>% duplicated() %>% any()
my_data_unique <- make_unique(my_data, "Gene_names", "ID", delim = ";") 
my_data$name %>% duplicated() %>% any()

Intensity_columns <- grep("iBAQ_", colnames(my_data_unique))
my_experimental_design <- as.data.frame(read.table("sample_info.txt", header = TRUE, sep = '\t'))
data_se <- make_se(my_data_unique, Intensity_columns, my_experimental_design)

Intensity_columns <- grep("iBAQ_", colnames(my_data_unique))
data_se_parsed <- make_se_parse(my_data_unique, Intensity_columns)

pdf(file="frequency.pdf", height = 3, width = 5)
plot_frequency(data_se)
dev.off()

data_filt <- filter_missval(data_se, thr = 0)

pdf(file="numbers.pdf", height = 3, width = 5)
plot_numbers(data_filt)
dev.off()

pdf(file="coverage.pdf", height = 5, width = 3)
plot_coverage(data_filt)
dev.off()

data_norm <- normalize_vsn(data_filt)

pdf(file="normalization.pdf", height = 5, width = 7)
plot_normalization(data_filt, data_norm)
dev.off().

pdf(file="missval.pdf", height = 5, width = 5)
plot_missval(data_filt)
dev.off()

pdf(file="detect.pdf", height = 5, width = 4)
plot_detect(data_filt)
dev.off()

data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

pdf(file="imputation.pdf", height = 5, width = 4)
plot_imputation(data_norm, data_imp_man)
dev.off()

data_diff <- test_diff(data_imp_man, type = "control", control = "control")
data_diff_all_contrasts <- test_diff(data_imp_man, type = "all")

dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

pdf(file="pca.pdf", height = 7, width = 7)
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4) 
dev.off()

pdf(file="cor.pdf", height = 7, width = 7)
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
dev.off()

pdf(file="heatmap.pdf", height = 7, width = 7)
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
dev.off()

data_results <- get_results(dep)
data_results %>% filter(significant) %>% nrow()
colnames(data_results)
write.table(data_results,"quantitative_DE_results.txt", row.names = FALSE, sep = "\t")

df_wide <- get_df_wide(dep)
df_long <- get_df_long(dep)

save(data_se, data_norm, data_imp_man, data_diff, dep, file = "data.RData")
load("data.RData")

############################################################
#################### Binary DE analysis ####################

library(metagenomeSeq)

my_intensities <- as.data.frame(read.table("intensities.txt", header = TRUE, row.names = 1, sep = '\t'))
my_intensities <- newMRexperiment(my_intensities)

my_experimental_design <- as.data.frame(read.table("sample_info.txt", header = TRUE, sep = '\t'))

FET <- fitPA(my_intensities, my_experimental_design$condition)
write.table(FET,"FET.txt", sep = "\t")


######################################################
#################### Volcano plot ####################

library(ggplot2)

de <- read.table("quantitative_DE_results.txt", header = TRUE, row.names = 1, check.names=FALSE, sep = "\t")
ggplot(data=de, aes(x=case_vs_control_ratio, y=case_vs_control_p.val)) + geom_point()

p <- ggplot(data=de, aes(x=case_vs_control_ratio, y=-log10(case_vs_control_p.val))) + geom_point()

p <- ggplot(data=de, aes(x=case_vs_control_ratio, y=-log10(case_vs_control_p.val))) + geom_point() + theme_minimal()

p2 <- p + geom_vline(xintercept=c(-6, 6), col="red") +
    geom_hline(yintercept=-log10(1.0E-07), col="red") 
	
de$correlated <- "No"

de$correlated[de$case_vs_control_ratio > 6 & de$case_vs_control_p.val < 1.0E-07] <- "Positive"

de$correlated[de$case_vs_control_ratio < -6 & de$case_vs_control_p.val < 1.0E-07] <- "Negative"

p <- ggplot(data=de, aes(x=case_vs_control_ratio, y=-log10(case_vs_control_p.val), col=correlated)) + geom_point() + theme_minimal()

p2 <- p + geom_vline(xintercept=c(-6, 6), col="red") +
        geom_hline(yintercept=-log10(1.0E-07), col="red")
		
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("Negative", "Positive", "No")
p3 <- p2 + scale_colour_manual(values = mycolors)

de$delabel <- NA
de$delabel[de$correlated != "No"] <- de$Gene_name_single[de$correlated != "No"]

ggplot(data=de, aes(x=case_vs_control_ratio, y=-log10(case_vs_control_p.val), col=correlated, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text()
	
library(ggrepel)

pdf(file="volcano.pdf", height = 10, width = 10)

ggplot(data=de, aes(x=case_vs_control_ratio, y=-log10(case_vs_control_p.val), col=correlated, label=delabel)) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-6, 6), col="red") +
        geom_hline(yintercept=-log10(1.0E-07), col="red") +
		theme_bw() +
		theme(legend.position="none")
dev.off()