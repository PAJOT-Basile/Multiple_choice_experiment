# Libraries
library(tidyverse)

# Import the number of the library from the arguments
args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]

# Import the table of all the barcodes for all the libraries
all_barcodes_library <- read.table(paste0(wd, "/barcodes/barcodes_all.tsv"), sep = "\t",
                           header = TRUE)

for (library_id in c(1, 2, 3)){
    for (index in all_barcodes_library %>% pull(code_index) %>% unique){
        for (lane in all_barcodes_library %>% pull(Lane) %>% unique){
                        temp_file <- all_barcodes_library %>%
                filter(code_index == index,
                       ID_Banque == paste0("Jaera_", library_id),
                       Lane == lane) %>%
                select(barcode_lettre, ID_DNA_RAD) 
                
            if (nrow(temp_file) == 0) next
            else{
                temp_file %>%
                    write.table(paste0(wd, "/barcodes/barcodes_", library_id, "_", index, "-", lane, ".tsv"),
                                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

            }
        }
    }
}
