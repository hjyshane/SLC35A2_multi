# Basic Anchors
"^"      # Start of string
"$"      # End of string
"^file$" # Exact match for "file"

# Quantifiers
"*"      # Zero or more (0+)
"+"      # One or more (1+)
"?"      # Zero or one (0 or 1)
"{n}"    # Exactly n times
"{n,}"   # n or more times
"{n,m}"  # Between n and m times

# Character Classes
"[abc]"  # Any single character from a, b, or c
"[^abc]" # Any single character except a, b, or c
"[A-Z]"  # Any single character from A to Z
"[a-z]"  # Any single character from a to z
"[0-9]"  # Any single digit
"\\w"    # Word character [A-Za-z0-9_]
"\\d"    # Digit [0-9]
"\\s"    # Whitespace
"\\."    # . 
"csv$"   # matches "csv" at the end ($) of the string give me any other guides that can be used for pattern

# Examples for file matching: 
".*\\.csv$"                    # Any file ending in .csv
"^data_.*\\.csv$"             # Files starting with "data_" and ending in .csv
"[0-9]{4}\\.txt$"             # Files like "1234.txt"
"[A-Za-z]+_[0-9]+\\.csv$"     # Files like "Test_123.csv"
"(PTZ|SAL).*\\.csv$"          # Files containing PTZ or(|) SAL and ending in .csv

# Real-world examples for your data:
"GO_results_.*_(PTZ|SAL)_\\d+hr\\.csv$"    # GO results for PTZ/SAL with hour suffix
".*_(24|1)hr_.*\\.csv$"                    # Files containing 24hr or 1hr
"^GO_results_[A-Za-z]+_.*\\.csv$"          # GO results starting with letters

# Complex Examples
# Match specific cell types
"GO_results_(Neurons|Astrocytes)_.*\\.csv$"

# Match specific timepoints
".*_(24hr|1hr)_(PTZ|SAL)\\.csv$"

# Match files with numeric ranges
".*_\\d{1,2}hr_.*\\.csv$"     # 1-99 hours

# Combine multiple patterns
pattern <- paste0(
  "GO_results_",            # Prefix
  "(Neurons|Astrocytes)_",  # Cell types
  "(PTZ|SAL)_",            # Conditions
  "\\d+hr",                # Time
  "\\.csv$"                # Extension
)

# Match files with date patterns
".*_\\d{4}-\\d{2}-\\d{2}\\.csv$"    # YYYY-MM-DD
".*_\\d{8}\\.csv$"                   # YYYYMMDD

# Match specific experimental conditions
".*_(control|treated)_.*\\.csv$"

# Match replicate numbers
".*_rep[1-3]\\.csv$"                 # rep1, rep2, rep3

# Combine multiple conditions
pattern <- paste(
  "GO_results_",
  "([A-Za-z]+)_",                  # Cell type
  "(PTZ|SAL)_",                    # Treatment
  "(24|1)hr",                      # Time
  "\\.csv$",
  sep=""
)

# Test your patterns
test_files <- c(
  "GO_results_Neurons_PTZ_24hr.csv",
  "GO_results_Astrocytes_SAL_1hr.csv"
)
grep(pattern, test_files, value=TRUE)

# Use grep(pattern, files, value=TRUE) to test patterns
matched_files <- grep(pattern, files, value=TRUE)

# Use basename() to remove path information before matching

# Use list.files() with full.names=TRUE to get complete paths
files <- list.files(path=".", pattern=".csv$")

print(matched_files)
# Test your pattern
pattern <- "GO_results_.*_(PTZ|SAL)_\\d+hr\\.csv$"


# Use print(pattern) to verify your constructed pattern
# Example 1: Single file path
full_path <- "/igm/home/hxy008/PTZ_ATAC_scRNA_072024/WIP/File/GO_results_Astrocytes_PTZ_24hr.csv"
filename <- basename(full_path)
print(filename)
# Output:
# [1] "GO_results_Astrocytes_PTZ_24hr.csv"

# Example 2: Multiple file paths
files <- c(
  "/igm/home/hxy008/PTZ_ATAC_scRNA_072024/WIP/File/GO_results_Astrocytes_PTZ_24hr.csv",
  "/igm/home/hxy008/PTZ_ATAC_scRNA_072024/WIP/File/GO_results_Neurons_SAL_1hr.csv",
  "/igm/home/hxy008/PTZ_ATAC_scRNA_072024/WIP/File/GO_results_OPCs_24hrvs1hr_PTZ.csv"
)
filenames <- basename(files)
print(filenames)
# Output:
# [1] "GO_results_Astrocytes_PTZ_24hr.csv"
# [2] "GO_results_Neurons_SAL_1hr.csv"
# [3] "GO_results_OPCs_24hrvs1hr_PTZ.csv"

# Example 3: Using list.files()
files <- list.files(path = "~/PTZ_ATAC_scRNA_072024/WIP/File/", 
                    pattern = "GO_results", 
                    full.names = TRUE)
filenames <- basename(files)
print(filenames)
# Output might look like:
# [1] "GO_results_Astrocytes_PTZ_24hr.csv"
# [2] "GO_results_Astrocytes_SAL_1hr.csv"
# [3] "GO_results_Neurons_PTZ_24hr.csv"
# [4] "GO_results_Neurons_SAL_1hr.csv"
# [5] "GO_results_OPCs_24hrvs1hr_PTZ.csv"
# [6] "GO_results_OPCs_24hrvs1hr_SAL.csv"

# Example 4: Print with numbering
for(i in seq_along(filenames)) {
  print(paste(i, filenames[i]))
}
# Output:
# [1] "1 GO_results_Astrocytes_PTZ_24hr.csv"
# [1] "2 GO_results_Astrocytes_SAL_1hr.csv"
# [1] "3 GO_results_Neurons_PTZ_24hr.csv"
# [1] "4 GO_results_Neurons_SAL_1hr.csv"
# [1] "5 GO_results_OPCs_24hrvs1hr_PTZ.csv"
# [1] "6 GO_results_OPCs_24hrvs1hr_SAL.csv"

