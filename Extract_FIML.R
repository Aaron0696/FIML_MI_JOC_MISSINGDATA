# Function to extract relevant data from output files
extract <- function(segment)
{
  df <- data.frame()
  # loop through each core
  for (core in 1:numCore)
  {
    if (length(segment[[core]]) == 0)
    {
      print(paste("Error in Condition ",
                  condition,
                  ", Core ",
                  core,
                  sep = ""))
    } else
    {
      # loop through the lavaan objects within each core
      for (a in seq(1,length(segment[[core]]), by = 2))
      {
        df <- rbind(df,segment[[core]][[a]])
        names(df) <- names(segment[[core]][[a]])
        
        paste("Core ",
              core,
              " Number ",
              a,
              " Done",
              sep = "")
      }
    }
  }
  return(df)
}
