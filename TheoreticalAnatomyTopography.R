#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
#          THALAMUS IN THE MIDDLE - THEORETICAL ANATOMY - 3D TOPOLOGY          #
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



### Order and colour by injection site position –––––#

InjSiteCoordDF <-  as.data.frame ( read.csv ( '~/_Julia_/ThalamusInTheMiddle/data/exp_raw/InjectionSiteEstimates_FlatmapProjected_260108.csv' ) )
InjSiteCoordDF <- InjSiteCoordDF[,-c ( 1:2 ) ]
colnames ( InjSiteCoordDF )[1:2] <- c ( 'mouse' , 'channel' )
InjSiteCoordDF$uninj <- with ( InjSiteCoordDF , paste ( mouse , channel , sep = '_' ) )
InjSiteCoordDF$uninj <- as.factor ( InjSiteCoordDF$uninj )

# Keep only the AAVrg injections
AAvrgInjCoordDF <- droplevels ( subset ( InjSiteCoordDF , Virus == 'AAVrg-XFP' & uninj %in% levels ( CoordDF$uninj ) ) )

# Available injections in CoordDF

AAvrgCoordDF <- droplevels ( subset ( CoordDF , uninj %in% levels ( AAvrgInjCoordDF$uninj ) ) )

# order injection sites on ML axis
str(InjSiteCoordDF)
MLord <- as.vector ( AAvrgInjCoordDF[ order ( AAvrgInjCoordDF$flat_u ) ,'uninj' ] )
APord <- as.vector ( AAvrgInjCoordDF[ order ( AAvrgInjCoordDF$flat_v ) ,'uninj' ] )

levels ( AAvrgCoordDF$uninj ) %in% MLord
MLord %in% levels ( AAvrgCoordDF$uninj )

#––––––––––––––––––––––––––––––––––––––––––––––––––––#



# Low N injections ––––––––––––––––––––––––––––––––––#

TooSmallInj <- sapply ( split ( df , ~uninj ) , nrow )
TooSmallInj <- TooSmallInj[ which ( TooSmallInj < 100 ) ]


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
cex.leg <- 1
asp <- 1.1

colblind <- trellis.par.get ( 'superpose.symbol' ) $col
# '#0072B2' '#E69F00' '#009E73' '#D55E00' '#56B4E9' '#F0E442' '#CC79A7'

#––––––––––––––––––––––––––––––––––––––––––––––––––––#




### 3D dynamic plot - ML Topography –––––––––––––––––#

nuc <- nucAll
# nuc <- c ( 'MD' )
# nuc <- c ( 'MD' , 'LH' , 'CL' , 'PCN' , 'CM' , 'IMD' , 'PVT' , 'PT' , 'AD' , 'IAD' , 'PF' )
# nuc <- c ( 'VM' )
# nuc <- c ( 'VM' , 'VAL' , 'VPM' , 'VPL' )
# nuc <- c ( 'MD' , 'LH' , 'CL' , 'PCN' , 'CM' , 'IMD' , 'PVT' , 'PT' , 'AD' , 'IAD' , 'PF' , 'VM' , 'VAL' , 'VPM' , 'VPL' )
# nuc <- c ( 'RT' )


# df <- droplevels ( subset ( AAvrgCoordDF , nuclei == 'MD' & ML > 230 ) )
df <- droplevels ( subset ( AAvrgCoordDF , nuclei %in% nuc & ML > 230 & ! uninj %in% TooSmallInj ) )
df$uninj <- factor ( df$uninj , levels = MLord )

# main data cloud
p3D <-
  plot_ly ( df , type = 'scatter3d' , x = ~AP , y = ~ML , z = ~DV ,
            color = ~uninj , colors = mako ( length ( MLord ) , direction = -1 ) ,
            scene = 'scene' , mode = 'markers' , marker = list ( size = 5 , symbol = 'circle' )
  )

p3D <-
  plot_ly ( df , type = 'scatter3d' , x = ~AP , y = ~ML , z = ~DV ,
          color = ~nuclei , colors = mako ( nlevels( df$nuclei ) , direction = -1 ) ,
          scene = 'scene' , mode = 'markers' , marker = list ( size = 5 , symbol = 'circle' )
)


cam <- list ( eye = list ( x = -0.8 , y = -1.2 , z = 1.2 ) )
# xax <- list ( title = 'Antero-posterior' , range = c ( 220 , 330 ) )
# yax <- list ( title = 'Medio-lateral' , range = c ( 200 , 280 ) , autorange = 'reversed' )
# zax <- list ( title = 'Dorso-ventral' , range = c ( 120 , 200 ) , autorange = 'reversed' )

xax <- list ( title = 'Antero-posterior' )
yax <- list ( title = 'Medio-lateral' , autorange = 'reversed' )
zax <- list ( title = 'Dorso-ventral' , autorange = 'reversed' )

aspm <- 'auto'
# asp <- list ( x = 0.7 , y = 0.7 , z = 0.7 )
# annot <- list ( list ( text = c ( 'MD' ) , x = c ( 300 ) , y = c ( 270 ) , z = c ( 140 ) , showarrow = F ) )

p3D <- layout ( p3D , scene = list ( xaxis = xax , yaxis = yax , zaxis = zax ,
                                     camera = cam , aspectmode = aspm #, #aspectratio = asp ,
                                     # annotations = annot 
                                     ) )

p3D


# saveWidget ( p3D, '/Users/loicmagrou/_Julia_/ThalamusInTheMiddle/plots/Theoretical Anatomy/ThalamusLabelledCells_Scatter3D_M691386_MD_v1.2_251023.html' )



#––––––––––––––––––––––––––––––––––––––––––––––––––––#




### 3D dynamic plot - AP Topography –––––––––––––––––#

df <- droplevels ( subset ( AAvrgCoordDF , nuclei %in% nuc & ML > 230 ) )
df$uninj <- factor ( df$uninj , levels = APord )

# main data cloud
p3D <-
  plot_ly ( df , type = 'scatter3d' , x = ~AP , y = ~ML , z = ~DV , 
            color = ~uninj , colors = magma ( length ( MLord ) , direction = -1 ) , 
            scene = 'scene' , mode = 'markers' , marker = list ( size = 5 , symbol = 'circle' ) 
  )


cam <- list ( eye = list ( x = -0.8 , y = -1.2 , z = 1.2 ) )
# xax <- list ( title = 'Antero-posterior' , range = c ( 220 , 330 ) )
# yax <- list ( title = 'Medio-lateral' , range = c ( 200 , 280 ) , autorange = 'reversed' )
# zax <- list ( title = 'Dorso-ventral' , range = c ( 120 , 200 ) , autorange = 'reversed' )

xax <- list ( title = 'Antero-posterior' )
yax <- list ( title = 'Medio-lateral' , autorange = 'reversed' )
zax <- list ( title = 'Dorso-ventral' , autorange = 'reversed' )

aspm <- 'auto'
# asp <- list ( x = 0.7 , y = 0.7 , z = 0.7 )
# annot <- list ( list ( text = c ( 'MD' ) , x = c ( 300 ) , y = c ( 270 ) , z = c ( 140 ) , showarrow = F ) )

p3D <- layout ( p3D , scene = list ( xaxis = xax , yaxis = yax , zaxis = zax ,
                                     camera = cam , aspectmode = aspm #, #aspectratio = asp ,
                                     # annotations = annot 
                                     ) )

p3D


# saveWidget ( p3D, '/Users/loicmagrou/_Julia_/ThalamusInTheMiddle/plots/Theoretical Anatomy/ThalamusLabelledCells_Scatter3D_M691386_MD_v1.2_251023.html' )



#––––––––––––––––––––––––––––––––––––––––––––––––––––#



