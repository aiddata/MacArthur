#Calculate the raw distance to every chinese project from every grid cell.
#Return the updated dataframe with new columns.
RDist <- function(cells, aid_locs)
{
  Dmatrix <- distm(coordinates(aid_locs), coordinates(cells), fun=distHaversine)
  return(Dmatrix)
}