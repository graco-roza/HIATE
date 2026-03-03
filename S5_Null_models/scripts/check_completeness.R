library(dplyr)

# Folder path
folder <- "~/Downloads/null_output"

# List all files matching pattern ou_x.txt (x=1..165)
files <- list.files(folder, pattern = "^out_\\d+\\.txt$", full.names = TRUE)

# Function to extract last iteration number from last line of a file
get_last_iteration <- function(file) {
  lines <- readLines(file, warn = FALSE)
  if(length(lines) == 0) return(NA_integer_)
  last_line <- tail(lines, 1)
  iter <- sub(".*Processing iteration (\\d+) of 999.*", "\\1", last_line)
  as.integer(iter)
}

# Create dataframe
df <- data.frame(
  Id = as.integer(gsub("^out_(\\d+)\\.txt$", "\\1", basename(files))),
  iteration = unlist(sapply(files, get_last_iteration))
)

df <- arrange(df, Id)

df |> as_tibble() |>  filter(iteration <  999) |>  arrange(iteration)


