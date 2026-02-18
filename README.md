# RudmanLab_GLM_Parallelization
GLM Parallelization in R

## 1. Structure your data
<img width="385" height="121" alt="image" src="https://github.com/user-attachments/assets/f5977af5-50b1-465f-a990-7c8832597c40" />

## 2. Long pivot your data
<img width="286" height="240" alt="image" src="https://github.com/user-attachments/assets/dad620ed-59d2-4f08-9bc2-3634c9dcb8db" />

## 3. Filter data with AF thresholds
````
DGRP_Spino3_long_filtered <- DGRP_Spino3_long %>%
  group_by(chrom, pos) %>%
  mutate(overall_mean_AF = mean(AF, na.rm = TRUE)) %>%
  filter(overall_mean_AF > 0.05) %>%
  ungroup()
````

## 4. Load in packages & remove CHR4
````
library(dplyr)
library(foreach)
library(doSNOW)

# Convert to regular data frame and ensure standard column types
DGRP_Spino3_long_filtered <- DGRP_Spino3_long_filtered %>%
  as.data.frame() %>%
  mutate(chrom = as.character(chrom),
         pos = as.numeric(pos)) %>%
  filter(!chrom %in% c("4"))

# Check whether you have all the chromosomes of interest
unique_chroms <- unique(DGRP_Spino3_long_filtered$chrom)

````

## 5. Run the parallelization pipeline
### Key notes: 1. Data autosaves each chromosome GLM | 2. Outputs & RAM consumption is cleared out <<< Speeds up processing immensely and makes it more resource efficient
````
# Set up parallel cluster once (outside the chromosome loop)
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
registerDoSNOW(my.cluster)

# Loop through each chromosome
for (current_chrom in unique_chroms) {
  
  print(paste("Processing chromosome:", current_chrom))
  
  # Subset data for current chromosome
  chrom_data <- DGRP_Spino3_long_filtered %>%
    filter(chrom == current_chrom)
  
  # Get unique positions for this chromosome
  unique_combination <- chrom_data %>%
    distinct(chrom, pos)
  
  print(paste("Number of positions:", nrow(unique_combination)))
  
  # Set up progress bar
  pb <- txtProgressBar(max = nrow(unique_combination), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Run parallel GLM for this chromosome
  results_list <- foreach(i = 1:nrow(unique_combination), 
                          .packages = c("dplyr", "stats"),
                          .errorhandling = 'stop', 
                          .combine = 'rbind',
                          .options.snow = opts) %dopar% {
                            
                            current_comb <- unique_combination[i, ]
                            
                            subset_data <- chrom_data %>%
                              filter(pos == current_comb$pos)
                            
                            glm_model <- glm(AF ~ Treatment,
                                             data = subset_data,
                                             family = quasibinomial, 
                                             na.action = na.exclude)
                            
                            model_summary <- summary(glm_model)
                            p_values <- coef(model_summary)[, 4]
                            treatment_p_value <- p_values["TreatmentT"] # <<<< Change this p_value["VARIABLE"] to your Column name, and Variable combo
                            
                            data.frame(chrom = current_comb$chrom,
                                       pos = current_comb$pos,
                                       p_value = treatment_p_value,
                                       stringsAsFactors = FALSE)
                          }
  
  close(pb)
  
  # Save results for this chromosome
  output_file <- paste0("DGRP_Spino3_GLM_results_", current_chrom, ".rds")
  saveRDS(results_list, output_file)
  print(paste("Saved results to:", output_file))
  print(paste("Completed chromosome:", current_chrom))
  print("---")
  
  # Clear memory
  rm(chrom_data, unique_combination, results_list)
  gc()
}

# Stop cluster after all chromosomes are done
parallel::stopCluster(cl = my.cluster)
closeAllConnections()
````

## 6. Load in data:
````
practice2L <- readRDS("DGRP_Spino3_GLM_results_2L.rds")
````
<img width="242" height="109" alt="image" src="https://github.com/user-attachments/assets/20df70bc-fdd8-4f30-b47e-23aa4de3666f" />



### Parallelization Time Comparison:
## What typically happens when you run GLM in R (1 core):
<img width="824" height="300" alt="image" src="https://github.com/user-attachments/assets/a06bfd08-22c9-426c-ab52-adb86e22c2ad" />

## What happens when you parallelize with multiple cores (3 cores):
<img width="828" height="313" alt="image" src="https://github.com/user-attachments/assets/5961803e-b5c3-4f4d-93e2-f2a49d64a8ca" />

### 1.1 minutes vs 0.3 minutes, a 3.666667 increase in processing speed (WITH ONLY 3x, imagine 60 cores!!!)
