# TAXONOMY - GENERA

library(pheatmap)
library(ggplot2)
library(RColorBrewer)



# Read the CSV file
data <- read.csv("amr-expression-taxa.csv", row.names = 1, check.names = FALSE)

# Prepare data for the heatmap
mat_data <- data[, 3:76]  # Subset the expression values
Location <- data$Location  # Extract location information
Location_df <- as.data.frame(Location)
row.names(Location_df) <- rownames(data)



annotation_colors <- list(Location = c("Cooler" = "blue", "Hotbox" = "red"))

my_palette <- brewer.pal(9, "Set3")  # Using the BuPu palette from RColorBrewer

mat_data <- t(mat_data)

p <- pheatmap(mat_data,
         density.info="none",
         trace="none",
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         margins =c(12,15),
         treeheight_row=0,
         treeheight_col=0,
         Colv="NA",
         col = my_palette,
         cexRow=1,
         sepwidth=c(0.01,0.01),
         border_color="black",
         lhei = c(2,1),
         lwid = c(2,2),
         cellwidth = 40,
         cellheight = 20,
         colsep=1:ncol(mat_data),
         display_numbers = TRUE,
         number_format = "%.1f",
         number_color = "gray0",
         rowsep=1:nrow(mat_data),
         main = "Antimicrobial Resistance",
		 annotation_col = Location_df,
		 annotation_colors = annotation_colors,
		 show_rownames = TRUE,
		 show_colnames = TRUE
         )
		 
# Save the plot as a PNG file
ggsave("Heatmap_Taxa_AMR_Annotations.png", plot = p, width = 30, height = 30, dpi = 300)