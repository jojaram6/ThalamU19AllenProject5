#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
#     THALAMUS IN THE MIDDLE - THEORETICAL ANATOMY - WASSERSTEIN DISTANCE      #
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#


#––––––––––––––––––––––––––––––––––––––––––––––––––––#
#                      Packages                      #
#––––––––––––––––––––––––––––––––––––––––––––––––––––#
library ( lattice ) 
library ( latticeExtra ) 
library ( PrettyCols ) 
library ( colors3d )
library ( viridis )
library ( grid ) 
library ( plotly ) 
library ( htmlwidgets )
library ( transport )
library ( parallel )
library ( combinat )
library ( vegan )



#––––––––––––––––––––––––––––––––––––––––––––––––––––#
#                      Database                      #
#––––––––––––––––––––––––––––––––––––––––––––––––––––#

CSVpath <- '~/_Julia_/ThalamusInTheMiddle/data/exp_raw/CCF_CellCoordinates_260109_3/'
CSVs <- list.files ( CSVpath ) 

#––––––––––––––––––––––––––––––––––––––––––––––––––––#
#                    Extract Data                    #
#––––––––––––––––––––––––––––––––––––––––––––––––––––#

### Single dataframe ––––––––––––––––––––––––––––––––#

CoordDF <- do.call ( 
  rbind , 
  lapply ( CSVs , \( x ) {
    tb <- read.csv ( paste ( CSVpath , x , sep = '/' ) ) 
    if ( nrow ( tb ) != 0 ) {
      df <- as.data.frame ( tb ) 
      colnames ( df ) <- c ( 'AP' , 'DV' , 'ML' ) 
      nm <- unlist ( strsplit ( x , '_' ) ) 
      df$mouse <- nm[1]
      df$nuclei <- nm[2]
      df$channel <- nm[3]
      df
    }
  } ) 
) 


# head ( droplevels ( subset ( CoordDF1 , mouse == '691386' & nuclei == 'MD' & ML > 230 & channel == '488') ) )
# head ( droplevels ( subset ( CoordDF2 , mouse == '691386' & nuclei == 'MD' & ML > 230 & channel == '488') ) )
# head ( droplevels ( subset ( CoordDF2 , mouse == '691386' & nuclei == 'MD' & ML > 230 & channel == '488' & AP == 250 & DV == 141 ) ) )

nucAll <- c ( 'MD' , 'CL' , 'PCN' , 'PT' , # MD ans shell
              'LH' , 'MH' , # Habenula (dorsal to MD)
              'PVT' , 'IMD' , 'CM' , 'IAM' , 'RH' , 'RE' , 'Xi' , # Midline group
              'VM' , 'VAL' , 'VPM' , 'VPL' , # Ventral group
              'AD' , 'IAD' , 'PO' ,  'AV' , 'AM' , # Anterior group
              'LD' , 'LP' , # Dorsal group
              'SMT' , 'PR' , # Ventral medial group
              'PF' , # Postorior group
              'RT' # Reticulate
)


str ( CoordDF ) 
CoordDF <- CoordDF[ , c ( 'AP' , 'ML' , 'DV' , 'mouse' , 'nuclei' , 'channel' ) ] # Reorder columns so that x = ~AP , y = ~ML , z = ~DV 
CoordDF$mouse <- factor ( CoordDF$mouse , unique ( CoordDF$mouse ) ) 
CoordDF$nuclei <- factor ( CoordDF$nuclei , levels = nucAll ) 
CoordDF$channel <- factor ( CoordDF$channel , levels = c ( '488' , '445' , '561' ) ) # Green > Blue > red injection site from medial to lateral for mouse 691386 694513
CoordDF$uninj <- with ( CoordDF , paste ( mouse , channel , sep = '_' ) )
CoordDF$uninj <- as.factor ( CoordDF$uninj )

#––––––––––––––––––––––––––––––––––––––––––––––––––––#



### Mirror flip the left injection cells ––––––––––––#
# !!! midline between ML227 and ML228 !!! Everything below 227 included needs to be flipped

FlipInj <- c ( '696669_445' , '696669_561' , '757190_561' ) # After visual check of labelling in the left MD

CoordDF[ which ( CoordDF$uninj %in% FlipInj ), 'ML'] <- 227 + ( 227 - CoordDF[ which ( CoordDF$uninj %in% FlipInj ), 'ML'] )

# plot_ly ( subset ( CoordDF , nuclei == 'MD' ) , type = 'scatter3d' , x = ~AP , y = ~ML , z = ~DV , # plot check
#           color = ~uninj , mode = 'markers' , marker = list ( size = 5 , symbol = 'circle' ) )

#––––––––––––––––––––––––––––––––––––––––––––––––––––#




### Proper injection site distance matrix –––––––––––#
InjSiteCoordDF <-  as.data.frame ( read.csv ( '~/_Julia_/ThalamusInTheMiddle/data/exp_raw/InjectionSiteEstimates_FlatmapProjected_260108.csv' ) )
InjSiteCoordDF <- InjSiteCoordDF[,-c ( 1:2 ) ]
colnames ( InjSiteCoordDF )[1:2] <- c ( 'mouse' , 'channel' )
InjSiteCoordDF$uninj <- with ( InjSiteCoordDF , paste ( mouse , channel , sep = '_' ) )
InjSiteCoordDF$uninj <- as.factor ( InjSiteCoordDF$uninj )

InjDistMat3D <- as.matrix ( dist ( InjSiteCoordDF[,c ( 'AP' , 'DV' , 'ML' ) ] ) )  
dimnames ( InjDistMat3D ) <- list ( InjSiteCoordDF$uninj , InjSiteCoordDF$uninj )

InjDistMatFlat <- as.matrix ( dist ( InjSiteCoordDF[,c ( 'flat_u' , 'flat_v' ) ] ) )  
dimnames ( InjDistMatFlat ) <- list ( InjSiteCoordDF$uninj , InjSiteCoordDF$uninj )

# Keep only the AAVrg injections
AAvrg <- droplevels ( subset ( InjSiteCoordDF , Virus == 'AAVrg-XFP' )$uninj )


#––––––––––––––––––––––––––––––––––––––––––––––––––––#



# Low N injections ––––––––––––––––––––––––––––––––––#

TooSmallInj <- sapply ( split ( df , ~uninj ) , nrow )
TooSmallInj <- TooSmallInj[ which ( TooSmallInj < 100 ) ]


#––––––––––––––––––––––––––––––––––––––––––––––––––––#




#––––––––––––––––––––––––––––––––––––––––––––––––––––#
#                Wasserstein distance                #
#––––––––––––––––––––––––––––––––––––––––––––––––––––#


### Choose nuclei system ––––––––––––––––––––––––––––#
nucList <- list()
# MD and surroundings
nucList[[1]] <- c ( 'MD' , 'LH' , 'CL' , 'PCN' , 'CM' , 'IMD' , 'PVT' , 'PT' , 'AD' , 'IAD' , 'PF' )

# VM and surroundings
nucList[[2]] <- c ( 'VM' , 'VAL' , 'VPM' , 'VPL' )

# RT
nucList[[3]] <- c ( 'RT' )

#––––––––––––––––––––––––––––––––––––––––––––––––––––#



### Wasserstein distance - MD and surroundings ––––––#

for ( n in 1:length ( nucList ) ) {
  
  nuc <- nucList[[n]]
  
  ##### All nuclei #####
  
  # Wasserstein distance matrix
  
  # df <- droplevels ( subset ( CoordDF , nuclei %in% nuc & uninj %in% InjSiteCoordDF$uninj & ML > 230 ) )
  df <- droplevels ( subset ( CoordDF , nuclei %in% nuc & uninj %in% AAvrg & ML > 230 & ! uninj %in% TooSmallInj ) )
  
  df$nuclei <- factor ( df$nuclei , levels = nuc )
  
  plot_ly (df , type = 'scatter3d' , x = ~AP , y = ~-ML , z = ~-DV , color = ~uninj , colors = mako ( nlevels ( df$uninj ) ) , mode = 'markers' , marker = list ( size = 5 , symbol = 'circle' ) )
  
  
  wppList <- 
    lapply ( split ( df , ~uninj ) , \(s){
      s <- s[,c ( 'AP' , 'ML' , 'DV' )]
      wpp ( s , rep ( 1/nrow ( s ) , nrow ( s ) ) )
    })
  
  cbn <- combn2 ( names ( wppList ) )
  
  system.time(
    distVec <- unlist (
      
      mclapply ( seq_len ( nrow ( cbn ) ) , function ( r ){
        # r <- 2
        i <- cbn[r, 1]
        j <- cbn[r, 2]
        d <- wasserstein ( wppList[[i]] , wppList[[j]] , p = 2 )
        d
      } , mc.cores = detectCores() - 1 )
      
    )
  )
  
  wdMat <- matrix ( 0 , length ( wppList ) , length ( wppList ) )
  wdMat[lower.tri ( wdMat )] <- distVec
  wdMat[upper.tri ( wdMat )] <- t ( wdMat )[upper.tri ( wdMat )]
  
  dimnames ( wdMat ) <- list ( names ( wppList ) , names ( wppList ) )
  
  
  # Centroid distance matrix
  ctrDF <- 
    do.call ( rbind , lapply ( split ( df , ~ uninj ) , \(x){
      # x <- split ( df , ~ uninj )[[1]]
      d <- x[1,c ( 'AP' , 'ML' , 'DV' , 'uninj')]
      d[1, c ( 'AP' , 'ML' , 'DV' )] <- colMeans ( x[,c ( 'AP' , 'ML' , 'DV' )])
      d
    } 
    )
    )
  
  CtrDistMat <-  as.matrix ( dist ( ctrDF[,c ( 'AP' , 'DV' , 'ML' ) ] ) )  
  dimnames ( CtrDistMat ) <- list ( ctrDF$uninj , ctrDF$uninj )
  
  # Within/between animal matrix
  withinMat <- expand.grid ( sapply ( strsplit ( rownames ( wdMat ) , '_' ) , \(x) x[1] ) , 
                             sapply ( strsplit ( colnames ( wdMat ) , '_' ) , \(x) x[1] ) )
  
  withinMat <- matrix ( withinMat$Var1 == withinMat$Var2 , nrow (wdMat) , ncol ( wdMat ) )
  
  
  # All distances dataframe
  
  InjDistMat3D <- InjDistMat3D[ rownames ( wdMat ) , colnames ( wdMat ) ]
  InjDistMatFlat <- InjDistMatFlat[ rownames ( wdMat ) , colnames ( wdMat ) ]
  CtrDistMat <- CtrDistMat[ rownames ( wdMat ) , colnames ( wdMat ) ]
  
  
  # Principal Coordinate Analysis (PCoA)
  
  wdPCOA <- cmdscale ( wdMat , 2 , T )
  cdPCOA <- cmdscale ( CtrDistMat , 2 , T )
  injPCOA <- cmdscale ( InjDistMatFlat , 2 , T )
  if ( !all ( rownames ( wdPCOA$points ) == rownames ( injPCOA$points ) ) ) stop ( 'Row names must match!')  # quick check
  wdPTall <- protest ( wdPCOA$points , injPCOA$points )
  cdPTall <- protest ( cdPCOA$points , injPCOA$points )
  
  wdDFall <- data.frame ( 
    'Wdist' = wdMat[lower.tri ( wdMat , diag = F )] ,
    'Cdist' = CtrDistMat[lower.tri ( CtrDistMat , diag = F )] ,
    'Idist3D' = InjDistMat3D[lower.tri ( InjDistMat3D , diag = F )] ,
    'IdistFlat' = InjDistMatFlat[lower.tri ( InjDistMatFlat , diag = F )] , 
    'Within' = withinMat[lower.tri ( withinMat , diag = F )] ,
    'Nuclei' = factor ( 'All' , levels = c ( nuc , 'All' ) ) 
  )
  
  # xypall <- 
  #   xyplot ( Wdist ~ IdistFlat , wdDF , 
  #            xlab = list ( 'Injection distance (µm)' , cex = cex.axs1 ) , 
  #            ylab = list ( 'Wasserstein distance (a.u.)' , cex = cex.axs1 ) , 
  #            xlim = c ( NA , 650 ) , ylim = c ( NA , 45 ) ,
  #            scales = list ( tck = - 0.8 , cex = cex.axs1 ) , 
  #            panel = function ( x , y , ... ){
  #              panel.xyplot ( x , y , ... )
  #              m <- lm ( Wdist ~ IdistFlat , wdDF )
  #              r <- summary ( m )$adj.r.squared
  #              panel.abline ( coef ( m ) )
  #              panel.loess ( x , y , col = 'red' , ... )
  #              panel.text ( 550 , 10 , bquote ( r^2 == .( round ( r , 2 ) ) ) )
  #              panel.text ( 10 , 40 , 'All' )
  #            })
  
  
  
  ##### Individual nuclei #####
  
  system.time(
    wdDFList <- lapply ( split ( df , ~ nuclei ) , \(sn){ 
      # sn <- split ( df , ~ nuclei )[[3]]
      wppList <- 
        lapply ( split ( sn , ~uninj , drop = T ) , \(s){
          # s <- split ( sn , ~uninj )[[3]]
          s <- s[,c ( 'AP' , 'ML' , 'DV' )]
          wpp ( s , rep ( 1/nrow ( s ) , nrow ( s ) ) )
        })
      
      cbn <- combn2 ( names ( wppList ) )
      
      distVec <- unlist (
        
        mclapply ( seq_len ( nrow ( cbn ) ) , function ( r ){
          # r <- 2
          i <- cbn[r, 1]
          j <- cbn[r, 2]
          d <- wasserstein ( wppList[[i]] , wppList[[j]] , p = 2 )
          d
        } , mc.cores = detectCores() - 1 )
        
      )
      
      wdMat <- matrix ( 0 , length ( wppList ) , length ( wppList ) )
      wdMat[lower.tri ( wdMat )] <- distVec
      wdMat[upper.tri ( wdMat )] <- t ( wdMat )[upper.tri ( wdMat )]
      
      dimnames ( wdMat ) <- list ( names ( wppList ) , names ( wppList ) )
      
      # Centroid distance matrix
      ctrDF <- 
        do.call ( rbind , lapply ( split ( sn , ~ uninj ) , \(x){
          # x <- split ( sn , ~ uninj )[[1]]
          d <- x[1,c ( 'AP' , 'ML' , 'DV' , 'uninj')]
          d[1, c ( 'AP' , 'ML' , 'DV' )] <- colMeans ( x[,c ( 'AP' , 'ML' , 'DV' )])
          d
        } 
        )
        )
      
      CtrDistMat <-  as.matrix ( dist ( ctrDF[,c ( 'AP' , 'DV' , 'ML' ) ] ) )  
      dimnames ( CtrDistMat ) <- list ( ctrDF$uninj , ctrDF$uninj )
      
      # Within/between animal matrix
      withinMat <- expand.grid ( sapply ( strsplit ( rownames ( wdMat ) , '_' ) , \(x) x[1] ) , 
                                 sapply ( strsplit ( colnames ( wdMat ) , '_' ) , \(x) x[1] ) )
      
      withinMat <- matrix ( withinMat$Var1 == withinMat$Var2 , nrow (wdMat) , ncol ( wdMat ) )
      
      
      # All distances dataframe    
      InjDistMat3D <- InjDistMat3D[ rownames ( wdMat ) , colnames ( wdMat ) ]
      InjDistMatFlat <- InjDistMatFlat[ rownames ( wdMat ) , colnames ( wdMat ) ]
      CtrDistMat <- CtrDistMat[ rownames ( wdMat ) , colnames ( wdMat ) ]
      
      # Principal Coordinate Analysis (PCoA)
      
      wdPCOA <- cmdscale ( wdMat , 2 , T )
      cdPCOA <- cmdscale ( CtrDistMat , 2 , T )
      injPCOA <- cmdscale ( InjDistMatFlat , 2 , T )
      if ( !all ( rownames ( wdPCOA$points ) == rownames ( injPCOA$points ) ) ) stop ( 'Row names must match!')  # quick check
      wdPT <- protest ( wdPCOA$points , injPCOA$points )
      cdPT <- protest ( cdPCOA$points , injPCOA$points )
      
      wdDF <- data.frame ( 
        'Wdist' = wdMat[lower.tri ( wdMat , diag = F )] ,
        'Cdist' = CtrDistMat[lower.tri ( CtrDistMat , diag = F )] ,
        'Idist3D' = InjDistMat3D[lower.tri ( InjDistMat3D , diag = F )] ,
        'IdistFlat' = InjDistMatFlat[lower.tri ( InjDistMatFlat , diag = F )] , 
        'Within' = withinMat[lower.tri ( withinMat , diag = F )] ,
        'Nuclei' = factor ( unique ( sn$nuclei ) , levels = c ( nuc , 'All' ) )
      )
      
      list ( wdDF , wdPT , cdPT )
      
    }
    )
    
  )
  
  # str(wdDFList)
  
  # Relist for Wasserstein Distance Procrustes test
  wdPTList <- lapply ( wdDFList , \(x) x[[2]] )
  wdPTList$All <- wdPTall
  
  # Relist for Centroid Distance Procrustes test
  cdPTList <- lapply ( wdDFList , \(x) x[[3]] )
  cdPTList$All <- cdPTall
  
  
  # Relist wdDFList
  wdDFList <- do.call ( rbind , lapply ( wdDFList , \(x) x[[1]] ) )
  wdDFList <- rbind ( wdDFList , wdDFall )
  
  
  
  #––––––––––––––––––––––––––––––––––––––––––––––––––––#
  
  
  #––––––––––––––––––––––––––––––––––––––––––––––––––––#
  #                      Figures                       #
  #––––––––––––––––––––––––––––––––––––––––––––––––––––#
  
  ### Graphical parameters ––––––––––––––––––––––––––––#
  
  cex <- 1
  cex.title <- 2
  cex.lab <- 1.5
  cex.strip <- 0.8
  cex.axs1 <- 1.1
  cex.axs2 <- 0.5
  cex.txt <- 0.7
  cex.leg <- 0.8
  asp <- 0.8
  
  colblind <- trellis.par.get ( 'superpose.symbol' ) $col
  # '#0072B2' '#E69F00' '#009E73' '#D55E00' '#56B4E9' '#F0E442' '#CC79A7'
  
  #––––––––––––––––––––––––––––––––––––––––––––––––––––#
  
  
  
  # Plot ––––––––––––––––––––––––––––––––––––––––––––––#
  
  if ( 'MD' %in% nuc ){
    
    ylimList <- lapply ( 1:nlevels ( wdDFList$Nuclei ) , \(x) c ( -5 , 40 ) ) ; ylimList[[7]][2] <- 80
    
    xyp1 <- 
      xyplot ( Wdist ~ IdistFlat | Nuclei , wdDFList , groups = Within ,
               xlab = list ( 'Injection distance (µm)' , cex = cex.lab ) , 
               ylab = list ( 'Wasserstein distance (a.u.)' , cex = cex.lab ) , 
               xlim = c ( -25 , 650 ) , ylim = ylimList , pch = 20 ,
               scales = list ( tck = - 0.8 , cex = cex.axs1 , relation = 'free' ) , 
               as.table = T , layout = c ( 4 , 3 ) , aspect = asp , 
               key = list ( 'points' = list ( col = c ( '#0072B2' , '#E69F00' ) , pch = 20 , cex = 2 ) , 
                            'text' = list ( c ( 'Between Animal' , 'Within animal' ) , cex = cex.leg ) , 
                            corner = c ( 0.4 , 1.04 ) , columns = 2 , padding.text = 2 ) ,
               legend = list ( 'inside' = list ( fun = draw.key ( key = list (
                 'lines' = list ( col = c ( 'black' , 'red' ) , lwd = 2 , size = 3 ) ,
                 'text' = list ( c ( 'Linear fit' , 'LOESS fit' ) , cex = cex.leg ) , columns = 2 , padding.text = 2 ) ) , 
                 corner = c ( 1.01 , 1.04 ) ) ) ,
               panel = function ( x , y , ... ){
                 panel.xyplot ( x , y , ... )
                 l <- lapply ( split ( wdDFList , ~Nuclei ) , \(s){
                   m<- lm ( Wdist ~ IdistFlat , s )
                   list ( coef = coef ( m ) , 
                          r = summary ( m )$adj.r.squared )
                 } )
                 panel.loess ( x , y , col = 'red' , lwd = 2 , ... )
                 panel.abline ( l[[panel.number()]]$coef , lwd = 2 )
                 # procruste stats
                 rpt <- sqrt ( 1 - wdPTList[[panel.number()]]$ss )
                 pval <- wdPTList[[panel.number()]]$signif
                 # text
                 if ( panel.number() == 7 ){
                   panel.text ( 690 ,  14, bquote ( r == .( format ( round ( ( sqrt ( l[[panel.number()]]$r ) ) , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , 4 , bquote ( r[proc] == .( format ( round ( rpt , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , -5 , bquote ( p[proc] == .( format ( round ( pval , 3 ) , nsmall = 3 ) ) ) , adj = c ( 1 , NA ) )
                 } else {
                   panel.text ( 690 , 6 , bquote ( r == .( format ( round ( sqrt ( l[[panel.number()]]$r ) , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , 0 , bquote ( r[proc] == .( format ( round ( rpt , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , -5 , bquote ( p[proc] == .( format ( round ( pval , 3 ) , nsmall = 3 ) ) ) , adj = c ( 1 , NA ) )
                 }
                 
               })
    
    pdf ( '~/_Julia_/ThalamusInTheMiddle/plots/Theoretical Anatomy/WassersteinInjectionSiteDistance_MDandSurrounding_v1.2_260130.pdf' , 10 , 7 )
    print ( xyp1 )
    dev.off()
    
    xyp2 <- 
      xyplot ( Cdist ~ IdistFlat | Nuclei , wdDFList , groups = Within , 
               xlab = list ( 'Injection distance (µm)' , cex = cex.lab ) , 
               ylab = list ( 'Centroid distance (a.u.)' , cex = cex.lab ) , 
               xlim = c ( -25 , 650 ) , ylim = ylimList , pch = 20 ,
               scales = list ( tck = - 0.8 , cex = cex.axs1 , relation = 'free' ) , 
               as.table = T , layout = c ( 4 , 3 ) , aspect = asp , 
               key = list ( 'points' = list ( col = c ( '#0072B2' , '#E69F00' ) , pch = 20 , cex = 2 ) , 
                            'text' = list ( c ( 'Between Animal' , 'Within animal' ) , cex = cex.leg ) , 
                            corner = c ( 0.4 , 1.04 ) , columns = 2 , padding.text = 2 ) ,
               legend = list ( 'inside' = list ( fun = draw.key ( key = list (
                 'lines' = list ( col = c ( 'black' , 'red' ) , lwd = 2 , size = 3 ) ,
                 'text' = list ( c ( 'Linear fit' , 'LOESS fit' ) , cex = cex.leg ) , columns = 2 , padding.text = 2 ) ) , 
                 corner = c ( 1.01 , 1.04 ) ) ) ,
               panel = function ( x , y , ... ){
                 panel.xyplot ( x , y , ... )
                 l <- lapply ( split ( wdDFList , ~Nuclei ) , \(s){
                   m<- lm ( Cdist ~ IdistFlat , s )
                   list ( coef = coef ( m ) , 
                          r = summary ( m )$adj.r.squared )
                 } )
                 panel.loess ( x , y , col = 'red' , lwd = 2 , ... )
                 panel.abline ( l[[panel.number()]]$coef , lwd = 2 )
                 # procruste stats
                 rpt <- sqrt ( 1 - cdPTList[[panel.number()]]$ss )
                 pval <- cdPTList[[panel.number()]]$signif
                 # text
                 if ( panel.number() == 7 ){
                   panel.text ( 690 ,  14, bquote ( r == .( format ( round ( ( sqrt ( l[[panel.number()]]$r ) ) , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , 4 , bquote ( r[proc] == .( format ( round ( rpt , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , -5 , bquote ( p[proc] == .( format ( round ( pval , 3 ) , nsmall = 3 ) ) ) , adj = c ( 1 , NA ) )
                 } else {
                   panel.text ( 690 , 6 , bquote ( r == .( format ( round ( sqrt ( l[[panel.number()]]$r ) , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , 0 , bquote ( r[proc] == .( format ( round ( rpt , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , -5 , bquote ( p[proc] == .( format ( round ( pval , 3 ) , nsmall = 3 ) ) ) , adj = c ( 1 , NA ) )
                 }
                 
               })
    
    pdf ( '~/_Julia_/ThalamusInTheMiddle/plots/Theoretical Anatomy/CentroidInjectionSiteDistance_MDandSurrounding_v1.2_260130.pdf' , 10 , 7 )
    print ( xyp2 )
    dev.off()
    
  }
  
  
  if ( 'VM' %in% nuc ){
    
    ylimList <- lapply ( 1:nlevels ( wdDFList$Nuclei ) , \(x) c ( -5 , 60 ) ) ; ylimList[[1]][2] <- 32 ; ylimList[[2]][2] <- 32
    
    xyp1 <- 
      xyplot ( Wdist ~ IdistFlat | Nuclei , wdDFList , groups = Within ,
               xlab = list ( 'Injection distance (µm)' , cex = cex.lab ) , 
               ylab = list ( 'Wasserstein distance (a.u.)' , cex = cex.lab ) , 
               xlim = c ( -25 , 650 ) , ylim = ylimList , pch = 20 ,
               scales = list ( tck = - 0.8 , cex = cex.axs1 , relation = 'free' ) , 
               as.table = T , layout = c ( 3 , 2 ) , aspect = asp , 
               key = list ( 'points' = list ( col = c ( '#0072B2' , '#E69F00' ) , pch = 20 , cex = 1.5 ) , 
                            'text' = list ( c ( 'Between Animal' , 'Within animal' ) , cex = cex.leg ) , 
                            corner = c ( 0. , 1.06 ) , columns = 2 , padding.text = 2 ) ,
               legend = list ( 'inside' = list ( fun = draw.key ( key = list (
                 'lines' = list ( col = c ( 'black' , 'red' ) , lwd = 2 , size = 3 ) ,
                 'text' = list ( c ( 'Linear fit' , 'LOESS fit' ) , cex = cex.leg ) , columns = 2 , padding.text = 2 ) ) , 
                 corner = c ( 1.02 , 1.06 ) ) ) ,
               panel = function ( x , y , ... ){
                 panel.xyplot ( x , y , ... )
                 l <- lapply ( split ( wdDFList , ~Nuclei ) , \(s){
                   m<- lm ( Wdist ~ IdistFlat , s )
                   list ( coef = coef ( m ) , 
                          r = summary ( m )$adj.r.squared )
                 } )
                 panel.loess ( x , y , col = 'red' , lwd = 2 , ... )
                 panel.abline ( l[[panel.number()]]$coef , lwd = 2 )
                 # procruste stats
                 rpt <- sqrt ( 1 - wdPTList[[panel.number()]]$ss )
                 pval <- wdPTList[[panel.number()]]$signif
                 # text
                 if ( panel.number() %in% 1:2 ){
                   panel.text ( 690 ,  4, bquote ( r == .( format ( round ( ( sqrt ( l[[panel.number()]]$r ) ) , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , -1 , bquote ( r[proc] == .( format ( round ( rpt , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , -5 , bquote ( p[proc] == .( format ( round ( pval , 3 ) , nsmall = 3 ) ) ) , adj = c ( 1 , NA ) )
                 } else {
                   panel.text ( 690 , 11 , bquote ( r == .( format ( round ( sqrt ( l[[panel.number()]]$r ) , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , 2 , bquote ( r[proc] == .( format ( round ( rpt , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , -5 , bquote ( p[proc] == .( format ( round ( pval , 3 ) , nsmall = 3 ) ) ) , adj = c ( 1 , NA ) )
                 }
                 
               })
    
    pdf ( '~/_Julia_/ThalamusInTheMiddle/plots/Theoretical Anatomy/WassersteinInjectionSiteDistance_VMandSurrounding_v1.2_260130.pdf' , 7.5 , 4.66 )
    print ( xyp1 )
    dev.off()
    
    xyp2 <- 
      xyplot ( Cdist ~ IdistFlat | Nuclei , wdDFList , groups = Within , 
               xlab = list ( 'Injection distance (µm)' , cex = cex.lab ) , 
               ylab = list ( 'Centroid distance (a.u.)' , cex = cex.lab ) , 
               xlim = c ( -25 , 650 ) , ylim = ylimList , pch = 20 ,
               scales = list ( tck = - 0.8 , cex = cex.axs1 , relation = 'free' ) , 
               as.table = T , layout = c ( 3 , 2 ) , aspect = asp , 
               key = list ( 'points' = list ( col = c ( '#0072B2' , '#E69F00' ) , pch = 20 , cex = 1.5 ) , 
                            'text' = list ( c ( 'Between Animal' , 'Within animal' ) , cex = cex.leg ) , 
                            corner = c ( 0. , 1.06 ) , columns = 2 , padding.text = 2 ) ,
               legend = list ( 'inside' = list ( fun = draw.key ( key = list (
                 'lines' = list ( col = c ( 'black' , 'red' ) , lwd = 2 , size = 3 ) ,
                 'text' = list ( c ( 'Linear fit' , 'LOESS fit' ) , cex = cex.leg ) , columns = 2 , padding.text = 2 ) ) , 
                 corner = c ( 1.02 , 1.06 ) ) ) ,
               panel = function ( x , y , ... ){
                 panel.xyplot ( x , y , ... )
                 l <- lapply ( split ( wdDFList , ~Nuclei ) , \(s){
                   m<- lm ( Cdist ~ IdistFlat , s )
                   list ( coef = coef ( m ) , 
                          r = summary ( m )$adj.r.squared )
                 } )
                 panel.loess ( x , y , col = 'red' , lwd = 2 , ... )
                 panel.abline ( l[[panel.number()]]$coef , lwd = 2 )
                 # procruste stats
                 rpt <- sqrt ( 1 - cdPTList[[panel.number()]]$ss )
                 pval <- cdPTList[[panel.number()]]$signif
                 # text
                 if ( panel.number() %in% 1:2 ){
                   panel.text ( 690 ,  4, bquote ( r == .( format ( round ( ( sqrt ( l[[panel.number()]]$r ) ) , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , -1 , bquote ( r[proc] == .( format ( round ( rpt , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , -5 , bquote ( p[proc] == .( format ( round ( pval , 3 ) , nsmall = 3 ) ) ) , adj = c ( 1 , NA ) )
                 } else {
                   panel.text ( 690 , 11 , bquote ( r == .( format ( round ( sqrt ( l[[panel.number()]]$r ) , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , 2 , bquote ( r[proc] == .( format ( round ( rpt , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
                   panel.text ( 690 , -5 , bquote ( p[proc] == .( format ( round ( pval , 3 ) , nsmall = 3 ) ) ) , adj = c ( 1 , NA ) )
                 }
                 
               })
    
    pdf ( '~/_Julia_/ThalamusInTheMiddle/plots/Theoretical Anatomy/CentroidInjectionSiteDistance_VMandSurrounding_v1.2_260130.pdf' , 7.5 , 4.66 )
    print ( xyp2 )
    dev.off()
    
  }
  
} # end if nuclist loop

#––––––––––––––––––––––––––––––––––––––––––––––––––––#

cex.lab <- 1.2

xyplot ( Cdist ~ IdistFlat | Nuclei , wdDFList , groups = Within , subset = c ( Nuclei == 'MD') ,
         xlab = list ( 'Flatmap injection distance (µm)' , cex = cex.lab ) , 
         ylab = list ( 'Centroid distance (µm)' , cex = cex.lab ) , 
         xlim = c ( -25 , NA ) , ylim = c ( -5 , 40 ) , pch = 19 , #col = c ( col1 , col2 ) , 
         scales = list ( tck = - 0.8 , cex = cex.axs1 , relation = 'free' ) , 
         aspect = asp , 
         panel = function ( x , y , ... ){
           panel.xyplot ( x , y , ... )
           l <- lapply ( split ( wdDFList , ~Nuclei ) , \(s){
             m<- lm ( Cdist ~ IdistFlat , s )
             list ( coef = coef ( m ) , 
                    r = summary ( m )$adj.r.squared )
           } )
           panel.loess ( x , y , col = 'red' , lwd = 2 , ... )
           panel.abline ( l[[panel.number()]]$coef , lwd = 2 )
           # procruste stats
           rpt <- sqrt ( 1 - cdPTList[[panel.number()]]$ss )
           pval <- cdPTList[[panel.number()]]$signif
           # text
           tx <- 600 ; ty <- -5 ; dty <- 3
           panel.text ( tx , ty+2*dty , bquote ( r == .( format ( round ( sqrt ( l[[panel.number()]]$r ) , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
           panel.text ( tx , ty+dty , bquote ( r[proc] == .( format ( round ( rpt , 2 ) , nsmall = 2 ) ) ) , adj = c ( 1 , NA ) )
           panel.text ( tx , ty , bquote ( p[proc] == .( format ( round ( pval , 3 ) , nsmall = 3 ) ) ) , adj = c ( 1 , NA ) )
         })

require ( car )
m1 <- lm (  Wdist ~ IdistFlat , wdDFList , subset = c ( Nuclei == 'LH' ) )             
m2 <- lm (  IdistFlat ~ Within , wdDFList , subset = c ( Nuclei == 'LH' ) )             

summary ( m1 ) ; summary ( m2 ) ; 
qqPlot(m)
plot(m)
