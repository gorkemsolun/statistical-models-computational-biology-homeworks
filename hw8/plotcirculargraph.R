# input symmetric matrix with row and column names.
# quick circular plot for our derived network
plot_circular_graph = function(M){

    # 1) Extract protein names
    protein.names <- colnames(M)
    n <- length(protein.names)
    
    # 2) Compute coordinates on a circle
    #    We'll distribute the nodes evenly around the circle.
    angles <- seq(0, 2*pi, length.out = n+1)[1:n]
    xcoords <- cos(angles)
    ycoords <- sin(angles)
    
    # 3) Create an empty plot with equal aspect ratio
    plot(
      xcoords, ycoords,
      type = "n",            # don't plot the points yet
      xlim = c(-1.2, 1.2),   # extra space for labels
      ylim = c(-1.2, 1.2),
      asp  = 1,              # ensure x/y have same scale
      axes = FALSE,          # no axes
      xlab = "", ylab = ""
    )
    
    # 4) Draw edges for all non-zero entries in M
    #    We assume M is symmetric; if not, you can decide how to handle that.
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        if(M[i, j] != 0) {
          segments(xcoords[i], ycoords[i], xcoords[j], ycoords[j], col = "blue")
        }
      }
    }
    
    # 5) Draw points (nodes)
    points(xcoords, ycoords, pch = 19, col = "red", cex = 1.5)
    
    # 6) Add labels slightly offset from the points
    text(
      xcoords * 1.15, 
      ycoords * 1.15, 
      labels = protein.names, 
      cex = 1
    )
}    
