#' @title create_tmp_directory
#' @description non-exportable function. creates a temporary directory on the
#' user's Home directory in order to store netMHCpan file and do a netMHCpan call,
#' the directory is removed once the function is completed.
#' @return creates a temporary directory on the user's home directory, is part
#' of the checknetMHCpan() function and automatically deletes itself once the
#' function is completed
#' @details In order to do a netMHCpan call, this function creates a temporary
#' directory which would store the netMHCpan file, and then the results, before
#' it gets deleted once the function is complete.
#' @noRd
#' @keywords internal

create_tmp_directory <- function(){
  user_home_dir <- Sys.getenv("HOME")
  unique_identifier<- grep(pattern="\\d+", Sys.time(), value=TRUE)
  unique_identifier<- gsub(pattern="\\:", "", unique_identifier)
  unique_identifier<- gsub(pattern=" ", "--", unique_identifier)
  unique_temp_name<- paste0(unique_identifier, "_RHF_tmp_files")
  path_for_tmp_files <- file.path(user_home_dir, unique_temp_name)
  dir.create(path_for_tmp_files)
  return (path_for_tmp_files)
}
