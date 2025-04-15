##Get dat weight matrix 
weight_mat <- function(ssn.object,
         addfunccol = NULL,
         response.col # response variable
         )
{
  data <- ssn.object@obspoints@SSNPoints[[1]]@point.data
  data <- cbind(data, ssn.object@obspoints@SSNPoints[[1]]@point.coords)
  xcol <- "coords.x1"
  ycol <- "coords.x2"
  
  nIDs <- sort(as.integer(as.character(unique(data[,"netID"])))) ## Just gies all the unique netIDs
  dist.junc <- matrix(0, nrow = length(data[,1]), ncol = length(data[,1])) ## Prepares a matrix similar to the distance one that gives u the weights
  net.zero <-  matrix(0, nrow = length(data[,1]), ncol = length(data[,1])) ## Not sure-
  nsofar <- 0 ## For the for loop
  distord <- order(data[,"netID"],data[,"pid"])
  names(distord) <- rownames(data)[distord]
  for(i in nIDs) { ## This for loop puts all the distance matrices into on ebig one like I drew out
    workspace.name <- paste("dist.net", i, ".RData", sep = "")
    path <- file.path(ssn.object@path, "distance", "obs",
                      workspace.name)
    if(!file.exists(path)) {
      stop("Unable to locate required distance matrix")
    }
    file_handle <- file(path, open="rb")
    distmat <- unserialize(file_handle) ## Gets the pre-computed dist mat
    
    ordpi <- order(as.numeric(rownames(distmat))) ##
    close(file_handle)
    ni <- length(distmat[1,])
    dist.junc[(nsofar + 1):(nsofar + ni),(nsofar + 1):
                (nsofar + ni)] <- distmat[ordpi, ordpi, drop = F] ## Put distance matrix in debido lugar
    ## Make dist inf if not on same network
    dist.junc[(nsofar + 1):(nsofar + ni),(nsofar + 1):
                (nsofar + ni)]
    
    net.zero[(nsofar + 1):(nsofar + ni),(nsofar + 1):(nsofar + ni)] <- 1
    nsofar <- nsofar + ni
  }
  
  n.all <- length(data[,1])
  a.mat <- NULL
  b.mat <- NULL
  a.mat.data <- NULL
  b.mat.data <- NULL
  dist.hydro <- NULL
  dist.hydro.data <- NULL
  w.matrix.data <- NULL
  # ind <- dataXY.out$indvecs$ind.allxy # Equivelent below
  n.all <- length(data[,1])
  # create a vector of all TRUE values
  ind.all <- rep(TRUE, times = n.all)
  
  #### REMEMBER: if there is a bug, it may be because these include covariates with NA-see dataXY function for how to get 
  # rid of those
  ind.allcov <- ind.all
  ind.allxy <- ind <- ind.allcov & !is.na(data[,response.col])
  
  ## create any necessary matrices from distance and flow matrices
  if(!is.null(dist.junc) ) {
    ## maximum distance to common junction between two sites
    a.mat <- pmax(dist.junc,t(dist.junc))
    a.mat.data <- a.mat[ind,ind]
    ## minimum distance to common junction between two sites
    b.mat <- pmin(dist.junc,t(dist.junc))
    b.mat.data <- b.mat[ind,ind]
    ## hydrological distance
    dist.hydro <- as.matrix(dist.junc + t(dist.junc))
    ## subset stream distance to observed locations only
    dist.hydro.data <- dist.hydro[ind, ind]
  }

  if(missing(addfunccol) || is.null(addfunccol) ||
     length(addfunccol) == 0 || !(addfunccol %in% colnames(data))){
      stop("The specified value for addfunccol was invalid")
  }
  flow.con.mat <- 1 - (b.mat > 0)*1
  
  ################################
  ################################
  ## for n.all, get data from the first lines on this function, and in dataXY function, theres a line that
  ## Defines this as 
  # total number of observations
  # addfunccol is the column in data you want to base the weights on, that is h20landwhatever

  # addfunccol <- "computed.afv"
  # distord is above, as well as net.zero
  
  w.matrix <- sqrt(pmin(outer(data[distord,addfunccol], ## HERE IT IS~~THE WEIGHT MATRIX!!
                              rep(1, times = n.all)),
                        t(outer(data[distord,addfunccol],rep(1, times = n.all)))) /
                     pmax(outer(data[distord,addfunccol],rep(1, times = n.all)),
                          t(outer(data[distord,addfunccol],rep(1, times = n.all)))))*
    flow.con.mat*net.zero
  w.matrix.data <- w.matrix[ind, ind]

  net.zero.data <- net.zero[ind,ind]
  return(list('weight' = w.matrix.data, 'dist_from_junct' = dist.junc, 'flow_con' = flow.con.mat, 'net.zero' = net.zero.data,
  'a.mat' = a.mat.data, 'b.mat' = b.mat.data))
}
