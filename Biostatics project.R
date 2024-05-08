Replication of results from the paper -
  Installing packages and libraries
install.packages("reshape")
install.packages("BiocManager")
install.packages("ggrepel")
install.packages("corrplot")
install.packages("pheatmap")
BiocManager::install("ComplexHeatmap")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db") # for converting gene accession numbers to gene names
BiocManager::install('PCAtools')
install.packages("rlang")
install.packages("gplots")
install.packages("factoextra")
install.packages("corrr")
library('corrr')
install.packages("ggcorrplot")
install.packages(“readr”)
library(ggcorrplot)
install.packages("FactoMineR")
library("FactoMineR")
library(PCAtools)
library(reshape)
library("org.Hs.eg.db")
library(tibble)
library(dplyr)
library(ggplot2)
library(gplots)
library(ggrepel)
library(pheatmap)
library(dplyr, warn.conflicts = FALSE)
install.packages("rlang")
library(rlang)
library(factoextra)
library(vctrs)
library(readr)
Reading infect.txt and plotting: -
  # ------------Code for Ratio of the infectivity vs different cell types plot--------
setwd("C:/Users/ABHIS/Downloads/biostats/datasets")
infect_data <- read.delim(file = "infect.txt", header = FALSE)
infect_data <- infect_data[order(infect_data$V2),]
infect_ratio <- infect_data$V2
names(infect_ratio) <- infect_data$V3
color_given <- function(value){
  if (value < 5) return("red")
  else if ( value < 25) return("green")
  else return("blue")
}
par(mar = c(6, 8, 4, 4))
color_vector <- unlist(lapply(infect_data$V2, color_given))
infect_data$V4 <- color_vector
barplot(infect_ratio, cex.names=0.6, horiz = TRUE,las=1, xlab = "Nef+/Nef- infectivity ratio",
        col=color_vector,
        main = "Ratio of the infectivity vs different cell types", border=0)\
Reading htseq files and merging to expression_df: -
  #Reading RNASeq Files and compile it into one using for loop
  for (i in 1:15) {
    file_name <- paste0("SRR21666", i + 23, ".htseq")
    RNASeq_Files[[paste0("C", i)]] <- read.table(file_name, header = FALSE)
  }
View(RNASeq_Files$C1)
View(Expression)
#Expression Table. Each row represents different Cell line. C1 represent sum of all, C2 =
no.of.reads of SERINC5, C3 = RPM
Expression[1,1] <- sum(RNASeq_Files$C3$V2)
Expression[1,2] <- RNASeq_Files$C3$V2[RNASeq_Files$C3$V1 == "ENSG00000164300"]
Expression[1,3] <- Expression[1,2]/(Expression[1,1]/1000000)
Expression[2,1] <- sum(RNASeq_Files$C11$V2)
Expression[2,2] <- RNASeq_Files$C11$V2[RNASeq_Files$C11$V1 == "ENSG00000164300"]
Expression[2,3] <- Expression[2,2]/(Expression[2,1]/1000000)
Expression[3,1] <- sum(RNASeq_Files$C9$V2)
Expression[3,2] <- RNASeq_Files$C9$V2[RNASeq_Files$C9$V1 == "ENSG00000164300"]
Expression[3,3] <- Expression[3,2]/(Expression[3,1]/1000000)
Expression[4,1] <- sum(RNASeq_Files$C1$V2)
Expression[4,2] <- RNASeq_Files$- Expression[1,2]/(Expression[1,1]/1000000)
C1$V2[RNASeq_Files$C1$V1 == "ENSG00000164300"]
Expression[4,3] <- Expression[4,2]/(Expression[4,1]/1000000)
Expression[5,1] <- sum(RNASeq_Files$C6$V2)
Expression[5,2] <- RNASeq_Files$C6$V2[RNASeq_Files$C6$V1 == "ENSG00000164300"]
Expression[5,3] <- Expression[5,2]/(Expression[5,1]/1000000)
Expression[6,1] <- sum(RNASeq_Files$C12$V2)
Expression[6,2] <- RNASeq_Files$C12$V2[RNASeq_Files$C12$V1 == "ENSG00000164300"]
Expression[6,3] <- Expression[6,2]/(Expression[6,1]/1000000)
Expression[7,1] <- sum(RNASeq_Files$C7$V2)
Expression[7,2] <- RNASeq_Files$C7$V2[RNASeq_Files$C7$V1 == "ENSG00000164300"]
Expression[7,3] <- Expression[7,2]/(Expression[7,1]/1000000)
Expression[8,1] <- sum(RNASeq_Files$C4$V2)
Expression[8,2] <- RNASeq_Files$C4$V2[RNASeq_Files$C4$V1 == "ENSG00000164300"]
Expression[8,3] <- Expression[8,2]/(Expression[8,1]/1000000)
Expression[9,1] <- sum(RNASeq_Files$C2$V2)
Expression[9,2] <- RNASeq_Files$C2$V2[RNASeq_Files$C2$V1 == "ENSG00000164300"]
Expression[9,3] <- Expression[9,2]/(Expression[9,1]/1000000)
Expression[10,1] <- sum(RNASeq_Files$C5$V2)
Expression[10,2] <- RNASeq_Files$C5$V2[RNASeq_Files$C5$V1 == "ENSG00000164300"]
Expression[10,3] <- Expression[10,2]/(Expression[10,1]/1000000)
Expression[11,1] <- sum(RNASeq_Files$C10$V2)
Expression[11,2] <- RNASeq_Files$C10$V2[RNASeq_Files$C10$V1 == "ENSG00000164300"]
Expression[11,3] <- Expression[11,2]/(Expression[11,1]/1000000)
Expression[12,1] <- sum(RNASeq_Files$C8$V2)
Expression[12,2] <- RNASeq_Files$C8$V2[RNASeq_Files$C8$V1 == "ENSG00000164300"]
Expression[12,3] <- Expression[12,2]/(Expression[12,1]/1000000)
Expression[13,1] <- sum(RNASeq_Files$C14$V2)
Expression[13,2] <- RNASeq_Files$C14$V2[RNASeq_Files$C14$V1 == "ENSG00000164300"]
Expression[13,3] <- Expression[13,2]/(Expression[13,1]/1000000)
Expression[14,1] <- sum(RNASeq_Files$C13$V2)
Expression[14,2] <- RNASeq_Files$C13$V2[RNASeq_Files$C13$V1 == "ENSG00000164300"]
Expression[14,3] <- Expression[14,2]/(Expression[14,1]/1000000)
Expression[15,1] <- sum(RNASeq_Files$C15$V2)
Expression[15,2] <- RNASeq_Files$C15$V2[RNASeq_Files$C15$V1 == "ENSG00000164300"]
Expression[15,3] <- Expression[15,2]/(Expression[15,1]/1000000)
Plotting the expression of SERINC5: -
  plot(infect_data$V2, Expression$V3,
       xlim = c(0, 40),
       xlab = "Nef+/Nef– infectivity ratio",
       ylab = "SERINC5 Expression (RPM)",
       cex = 2,
       col = c(rep("orangered", 8), rep("purple", 5), rep("darkblue", 2)),
       pch = 16)
# Add points
points(infect_data$V2, Expression$V3, col = "black", pch = 1, cex = 2)
# Fit linear regression model
fit <- lm(Expression$V3 ~ infect_data$V2)
# Get intercept and slope
intercept <- coef(fit)[1]
slope <- coef(fit)[2]
# Add trendline
abline(a = intercept, b = slope, col = "blue")
# Linear regression and Trendline
lm(infect_data$V2 ~ Expression$V3)
# Displaying Results
cor(infect_data$V2, Expression$V3)
summary(fit)$r.square
# r, R-square and p-value
text(x = 2.5, y = 75, paste("r:", round(sqrt(summary(fit)$r.square), digits = 4)), pos = 4)
text(x = 2.5, y = 70, paste("R^2:", round(summary(fit)$r.square, digits = 5)), pos = 4)
text(x = 2.5, y = 65, paste("p-value:", round(summary(fit)$coefficients[2, 4], digits = 7)), pos = 4)
HYPOTHESIS 1: - 
# Define protein names
protein_names <-
  c("ITGA1","ITGA2","ITGA3","ITGA4","ITGA5","ITGA6","ITGA7","ITGA8","ITGA9","ITGA10","ITG
A11","ITGA2B","ITGAD","ITGAE","ITGAL","ITGAM","ITGAV","ITGX")
# Define a function to calculate RPM
calculate_RPM <- function(gene_id, RNASeq_File) {
  total_reads <- sum(RNASeq_File$V2) # Total number of reads for all genes
  gene_reads <- RNASeq_File$V2[RNASeq_File$V1 == gene_id] # Reads for the gene
  rpm <- (gene_reads / total_reads) * 1e6 # Calculate RPM
  return(rpm)
}
# Calculate RPM for all ITAGs across all cell lines
genes <-
  c("ENSG00000213949","ENSG00000164171","ENSG00000005844","ENSG00000115232","EN
SG00000161638","ENSG00000091409","ENSG00000135424","ENSG00000077943","ENSG00
000144668","ENSG00000143127","ENSG00000137809","ENSG00000005961","ENSG0000015
6886","ENSG00000083457","ENSG00000005844","ENSG00000169896","ENSG00000138448"
    ,"ENSG00000140678")
# Create a list to store RPM values for each gene across all cell lines
rpm_values <- list()
for (i in 1:15) {
  gene_rpm <- sapply(genes, function(gene_id) {
    calculate_RPM(gene_id, RNASeq_Files[[paste0("C", i)]])
  })
  rpm_values[[i]] <- gene_rpm
}
# Create a data frame to store RPM values for each protein and cell line
rpm_df <- data.frame(Cell_Line = rep(cell_lines, each = length(protein_names)),
                     Protein = rep(protein_names, times = length(cell_lines)),
                     RPM = unlist(rpm_values))
# Save the data frame to a CSV file
write.csv(rpm_df, "RPM_values_all_cell_lines.csv", row.names = FALSE)
Plotting of graph for 18 ITGAs:-
  # Define colors for each graph
  colors <- c("lightblue", "coral",
              "lightpink","lightgreen","darkblue","darkgreen","yellow","seagreen","skyblue","brown","gold","grey
","maroon","orange","purple","cyan","ivory","khaki")
# Set custom margin
par(mar = c(5, 4, 4, 2)) # margin: bottom, left, top, right
# Increase the plot region size
par(plt = c(0.1, 0.9, 0.1, 0.9)) # plot region: xmin, xmax, ymin, ymax
# Plot fewer plots in each row and column to increase the size of individual plots
par(mfrow = c(9, 2)) # Arrange plots in a 9x2 grid
# Find the maximum RPM value across all plots to set the same y-axis limit
max_rpm <- max(sapply(rpm_values, max))
for (j in 1:length(protein_names)) {
  protein <- protein_names[j]
  gene_rpm_values <- sapply(rpm_values, `[`, j)
  # Plot each barplot with the same y-axis limits and title on y-axis
  barplot(gene_rpm_values, names.arg = cell_lines, las = 2, col = colors[j],
          ylab = paste(protein, "RPM"), xlab = "", ylim = c(0, max_rpm))
}
Preliminary tests
Normalizing data -> getting RPM values
For ITGA3
Shapiro test
# Running Shapiro-Wilk test
shapiro_test_result <- shapiro.test(infect_data$V2)
# Printing the result
print(shapiro_test_result)
Correlation Plot -
  # Load necessary libraries
  library(ggplot2)
# Create a data frame with RPM values for ITGA3 and SERINC5
rpm_data <- data.frame(ITGA3 = itga3_rpm, SERINC5 = serinc5_rpm)
# Calculate correlation coefficient
correlation <- cor(itga3_rpm, serinc5_rpm, method = "spearman")
# Plot correlation plot
ggplot(rpm_data, aes(x = ITGA3, y = SERINC5)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  RAMOS 5.62010816 0.02312802 11.1245762
WI38 0.67850697 15.4925757 0.28271124
labs(title = paste("Correlation Plot (Spearman Correlation =", round(correlation, 2), ")"),
     x = "ITGA3 RPM", y = "SERINC5 RPM")
#------
# Load necessary libraries
library(corrplot)
# Create a data frame with RPM values for ITGA3 and SERINC5
rpm_data <- data.frame(ITGA3 = itga3_rpm, SERINC5 = serinc5_rpm)
# Calculate the correlation matrix
correlation_matrix <- cor(rpm_data, method = "spearman")
# Plot correlation matrix heatmap
corrplot(correlation_matrix, method = "color", type = "upper",
         addCoef.col = "black", tl.col = "black", tl.srt = 45)
HYPOTHESIS 3 
# Define protein names
protein_names <- c("CD9")
# Define a function to calculate RPM
calculate_RPM <- function(gene_id, RNASeq_File) {
  total_reads <- sum(RNASeq_File$V2) # Total number of reads for all genes
  gene_reads <- RNASeq_File$V2[RNASeq_File$V1 == gene_id] # Reads for the gene
  rpm <- (gene_reads / total_reads) * 1e6 # Calculate RPM
  return(rpm)
}
# Calculate RPM for all ITAGs across all cell lines
genes <- c("ENSG00000010278")
# Create a list to store RPM values for each gene across all cell lines
rpm_values <- list()
for (i in 1:15) {
  gene_rpm <- sapply(genes, function(gene_id) {
    calculate_RPM(gene_id, RNASeq_Files[[paste0("C", i)]])
  })
  rpm_values[[i]] <- gene_rpm
}
# Create a data frame to store RPM values for each protein and cell line
rpm_df <- data.frame(Cell_Line = rep(cell_lines, each = length(protein_names)),
                     Protein = rep(protein_names, times = length(cell_lines)),
                     RPM = unlist(rpm_values))
# Save the data frame to a CSV file
write.csv(rpm_df, "RPM_values_all_cell_lines.csv", row.names = FALSE)
# Define colors for each graph
colors <- c("lightpink")
# Set custom margin
par(mar = c(5, 4)) # margin: bottom, left, top, right
# Increase the plot region size
par(plt = c(0.1, 0.9)) # plot region: xmin, xmax, ymin, ymax
# Find the maximum RPM value across all plots to set the same y-axis limit
max_rpm <- max(sapply(rpm_values, max))
for (j in 1:length(protein_names)) {
  protein <- protein_names[j]
  gene_rpm_values <- sapply(rpm_values, '[', j)
  # Plot each barplot with the same y-axis limits and title on y-axis
  barplot(gene_rpm_values, names.arg = cell_lines, las = 2, col = colors[j],
          ylab = paste(protein, "RPM"), xlab = "", ylim = c(0,max_rpm))
}
Correlation Plot:-
  # Calculate Spearman correlation between SERINC5 and CD9 RPM values
  spearman_result <- cor.test(serinc5_rpm, cd9_rpm, method = "spearman")
# Print the result
print(spearman_result)
library(ggplot2)
# Create a data frame with SERINC5 and CD9 RPM values
rpm_data <- data.frame(serinc5 = serinc5_rpm, cd9 = cd9_rpm)
# Calculate Spearman correlation coefficient
spearman_cor <- cor(serinc5_rpm, cd9_rpm, method = "spearman")
# Create the scatter plot
ggplot(rpm_data, aes(x = serinc5, y = cd9)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = paste("Spearman Correlation Plot (Correlation =", round(spearman_cor, 2), ")"),
       x = "SERINC5 RPM", y = "CD9 RPM")
HYPOTHESIS 4

# Define protein names
protein_names <- c("CD9","CD4")
# Define a function to calculate RPM
calculate_RPM <- function(gene_id, RNASeq_File) {
  total_reads <- sum(RNASeq_File$V2) # Total number of reads for all genes
  gene_reads <- RNASeq_File$V2[RNASeq_File$V1 == gene_id] # Reads for the gene
  rpm <- (gene_reads / total_reads) * 1e6 # Calculate RPM
  return(rpm)
}
# Calculate RPM for all ITAGs across all cell lines
genes <- c("ENSG00000010278","ENSG00000010610")
# Create a list to store RPM values for each gene across all cell lines
rpm_values <- list()
for (i in 1:15) {
  gene_rpm <- sapply(genes, function(gene_id) {
    calculate_RPM(gene_id, RNASeq_Files[[paste0("C", i)]])
  })
  rpm_values[[i]] <- gene_rpm
}
# Define colors for each graph
colors <- c("lightpink", "red")
# Set custom margin
par(mar = c(5, 4, 4, 2)) # margin: bottom, left, top, right
# Increase the plot region size
par(plt = c(0.1, 0.9, 0.1, 0.9)) # plot region: xmin, xmax, ymin, ymax
# Plot fewer plots in each row and column to increase the size of individual plots
par(mfrow = c(9, 2)) # Arrange plots in a 9x2 grid
# Find the maximum RPM value across all plots to set the same y-axis limit
max_rpm <- max(sapply(rpm_values, max))
for (j in 1:length(protein_names)) {
  protein <- protein_names[j]
  gene_rpm_values <- sapply(rpm_values, '[', j)
  # Plot each barplot with the same y-axis limits and title on y-axis
  barplot(gene_rpm_values, names.arg = cell_lines, las = 2, col = colors[j],
          ylab = paste(protein, "RPM"), xlab = "", ylim = c(0,max_rpm))
}
# Combine RPM values for different proteins into a single data frame
merged_rpm <- data.frame(Cell_Line = ordered_cell_lines)
for (j in 1:length(protein_names)) {
  protein <- protein_names[j]
  gene_rpm_values <- sapply(rpm_values, '[', j)
  merged_rpm[[protein]] <- gene_rpm_values
}
# Perform Shapiro-Wilk test on each protein's RPM values
shapiro_results <- lapply(protein_names, function(protein) {
  shapiro_test <- shapiro.test(merged_rpm[[protein]])
  return(c(Protein = protein, W = shapiro_test$statistic, p_value = shapiro_test$p.value))
})
# Convert Shapiro-Wilk results to a data frame
shapiro_df <- do.call(rbind, shapiro_results)
shapiro_df
# Save Shapiro-Wilk test results to a CSV file
write.csv(shapiro_df, "shapiro_results.csv", row.names = FALSE)
# Save merged RPM data to a CSV file
write.csv(merged_rpm, "merged_rpm_data.csv", row.names = FALSE)
Correlation plot for cd4 and cd9: -
  # Load necessary libraries
  library(ggplot2)
# Assuming 'merged_rpm' contains the RPM values for CD4 and CD9 across all cell lines
# Extract RPM values for CD4 and CD9
cd4_rpm <- merged_rpm$CD4
cd9_rpm <- merged_rpm$CD9
# Calculate correlation coefficient
correlation <- cor(cd4_rpm, cd9_rpm)
# Create a data frame with CD4 and CD9 RPM values
rpm_data <- data.frame(cd4 = cd4_rpm, cd9 = cd9_rpm)
# Plot correlation plot
ggplot(rpm_data, aes(x = cd4, y = cd9)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = paste("Correlation Plot (Correlation =", round(correlation, 2), ")"),
       x = "CD4 RPM", y = "CD9 RPM")

HYPOTHESIS 5 

Plotting
plot(infect_data$V2,Expression$V3, xlim=c(0,40), xlab = "Nef+/Nef– infectivity ratio", ylab =
       "BST2 Expression (RPM)", cex=2, col = c(rep("orangered",8),rep("purple",5),rep("darkblue",2)),
     pch = 16)
points(infect_data$V2,Expression$V3, col = "black", pch = 1, cex =2)
# Perform linear regression
lm_model <- lm(infect_data$V2 ~ Expression$V3)
# Extract R-squared value and p-value
r_square <- summary(lm_model)$r.squared
p_value <- summary(lm_model)$coefficients[2,4]
# Plot the data
pdf("Figure1d.pdf", width = 10, height=8)
par(mar = c(10, 10, 10, 10)) #margins
plot(infect_data$V2, Expression$V3, xlim=c(0,40), xlab = "Nef+/Nef– infectivity ratio",
     ylab = "BST2 Expression (RPM)", cex=2, col =
       c(rep("orangered",8),rep("purple",5),rep("navy",2)), pch = 16)
points(infect_data$V2, Expression$V3, col = "black", pch = 1, cex =2)
# Add the trendline
abline(lm_model)
# Add text with r-square and p-value
text(x = 2.5, y = 75, paste("R^2:", round(r_square, digits = 5)), pos = 4)
text(x = 2.5, y = 70, paste("p-value:", format(p_value, scientific = TRUE, digits = 2)), pos = 4)
dev.off()
Plotting the correlation plot:-
  # Calculate correlation coefficient
  correlation <- cor(bst2_rpm, serinc5_rpm)
# Create a data frame with CD4 and CD9 RPM values
rpm_data <- data.frame(bst2 = bst2_rpm, serinc5 = serinc5_rpm)
# Plot correlation plot
ggplot(rpm_data, aes(x = bst2, y = serinc5)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = paste("Correlation Plot (Correlation =", round(correlation, 2), ")"),
       x = "BST2 RPM", y = "SERINC5 RPM")