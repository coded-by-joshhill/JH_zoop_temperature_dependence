# Estimating carbon weights from body lengths
# Josh Hill
# 06/01/2025



# Libraries and helpers ----
library(tidyverse)


# Read in the data ----
dat <- read_csv("https://www.dropbox.com/scl/fi/b1pl3iys1kqtuiwbfd99h/IngClear_dat.csv?rlkey=4ce0im2s1latych74hos8uors&st=8c3k443n&dl=1", 
                skip = 1) %>%
  mutate(ref_no = paste0("Hill_", row_number()),
         taxa = str_squish(taxa)) %>% # create a unique identifier (e.g., Hill_row#)
  relocate(ref_no, .before = everything()) # move it before all columns
glimpse(dat)



# Estimate carbon body mass on the basis of body length
larvDat <- dat %>% 
  select(ref_no, primRef, body_length_mm, BM_C, weight_unit) %>% 
  filter(primRef %in% c("Lombard2009", "Broms2003", "Aguirre2006", "DadonPilosof2023")) %>% 
  arrange(ref_no) %>% 
  mutate(BM_C = NA, # remove from normal script
         weight_unit = NA) %>% # I assigned BMC to NA here but will need to remove values from CSVs
  mutate(BM_C = calc_BMC(body_length_mm),
         weight_unit = "mg") %>% 
  select(-body_length_mm) 


dat2 <- dat %>% 
  left_join(larvDat, by = c("ref_no", "primRef", "BM_C", "weight_unit")) %>% 
  arrange(as.numeric(str_extract(ref_no, "\\d+")))

glimpse(dat2)

         
         




