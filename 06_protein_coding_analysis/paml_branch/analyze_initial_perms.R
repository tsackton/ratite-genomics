#analyze initial RER permutations

library(tidyverse)

data_path <- "../init_perms/" # path to the data
files <- dir(data_path, pattern = "*.out")

data <- data_frame(filename = files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_tsv(file.path(data_path, .))) # a new data column
  ) %>% unnest

pvals <- data %>% group_by(filename) %>% summarize(p_count = sum(!is.na(P) & P < 0.01)/sum(!is.na(P)))

pvals_ext <- pvals %>% filter(grepl("extended", filename)) %>% separate(filename, into=c("test", "target_sp", "tree", "drop", "drop2")) %>% select(-drop, -drop2, -tree)

pvals_ext %>% filter(test == "random") %>% ggplot(aes(p_count)) + geom_density()

pvals_ext %>% group_by(test) %>% summarize(median_count = mean(p_count))
