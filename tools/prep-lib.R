# Download pre-built Rust library for spopt
# Called by configure/configure.win scripts

# Get package version from DESCRIPTION
version <- read.dcf("DESCRIPTION", fields = "Version")[[1]]

# Detect target platform
detect_target <- function() {
  os <- Sys.info()[["sysname"]]
  arch <- Sys.info()[["machine"]]

  if (os == "Windows") {
    return("x86_64-pc-windows-gnu")
  } else if (os == "Darwin") {
    if (arch == "arm64") {
      return("aarch64-apple-darwin")
    } else {
      return("x86_64-apple-darwin")
    }
  } else if (os == "Linux") {
    if (arch == "aarch64") {
      return("aarch64-unknown-linux-gnu")
    } else {
      return("x86_64-unknown-linux-gnu")
    }
  }

  stop("Unsupported platform: ", os, "-", arch)
}

target <- detect_target()
message("Detected target: ", target)

# Read checksums file
sums_file <- "tools/lib-sums.tsv"
if (!file.exists(sums_file)) {
  stop("Checksums file not found: ", sums_file)
}

sums <- read.delim(sums_file, header = TRUE, stringsAsFactors = FALSE)

# Find entry for our target
entry <- sums[sums$target == target, ]
if (nrow(entry) == 0) {
  stop("No pre-built binary available for target: ", target)
}

url <- entry$url[1]
expected_sha256 <- entry$sha256[1]

message("Downloading from: ", url)

# Download to tools directory
dest_file <- "tools/libspopt.a"
download_result <- tryCatch({
  download.file(url, dest_file, mode = "wb", quiet = FALSE)
  TRUE
}, error = function(e) {
  message("Download failed: ", e$message)
  FALSE
})

if (!download_result) {
  stop("Failed to download pre-built binary")
}

# Verify checksum
if (!is.na(expected_sha256) && nchar(expected_sha256) > 0) {
  actual_sha256 <- tools::md5sum(dest_file)  # R doesn't have sha256sum built-in

  # Use system sha256sum if available
  sha256_result <- tryCatch({
    if (Sys.info()[["sysname"]] == "Windows") {
      # Windows: use certutil
      output <- system2("certutil", c("-hashfile", dest_file, "SHA256"),
                        stdout = TRUE, stderr = TRUE)
      # Extract hash from second line
      hash <- gsub(" ", "", output[2])
    } else {
      # Unix: use shasum or sha256sum
      if (nzchar(Sys.which("sha256sum"))) {
        output <- system2("sha256sum", dest_file, stdout = TRUE)
      } else {
        output <- system2("shasum", c("-a", "256", dest_file), stdout = TRUE)
      }
      hash <- strsplit(output, " ")[[1]][1]
    }
    hash
  }, error = function(e) {
    message("Warning: Could not verify checksum: ", e$message)
    NA
  })

  if (!is.na(sha256_result) && sha256_result != expected_sha256) {
    unlink(dest_file)
    stop("Checksum mismatch!\n  Expected: ", expected_sha256, "\n  Got: ", sha256_result)
  }

  message("Checksum verified: ", sha256_result)
}

message("Successfully downloaded pre-built library to: ", dest_file)
