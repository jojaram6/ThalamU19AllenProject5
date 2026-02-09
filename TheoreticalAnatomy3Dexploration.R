#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
#        THALAMUS IN THE MIDDLE - THEORETICAL ANATOMY - SPATIAL ANALYSIS       #
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#


#––––––––––––––––––––––––––––––––––––––––––––––––––––#
#                      Packages                      #
#––––––––––––––––––––––––––––––––––––––––––––––––––––#
library ( lattice ) 
library ( latticeExtra ) 
library ( PrettyCols ) 
library ( grid ) 
library ( plotly ) 
library ( htmlwidgets )


#––––––––––––––––––––––––––––––––––––––––––––––––––––#
#                      Database                      #
#––––––––––––––––––––––––––––––––––––––––––––––––––––#

CSVpath <- '~/_Julia_/ThalamusInTheMiddle/data/exp_raw/CCF_CellCoordinates_250822'
CSVs <- list.files ( CSVpath ) 

#––––––––––––––––––––––––––––––––––––––––––––––––––––#
#                    Extract Data                    #
#––––––––––––––––––––––––––––––––––––––––––––––––––––#

### Single dataframe ––––––––––––––––––––––––––––––––#

CoordDF <- do.call ( 
  rbind , 
  lapply ( CSVs , \ ( x ) {
    tb <- read.csv ( paste ( CSVpath , x , sep = '/' ) ) 
    if ( nrow ( tb ) != 0 ) {
      df <- as.data.frame ( tb ) 
      colnames ( df ) <- c ( 'AP' , 'DV' , 'ML' ) 
      nm <- unlist ( strsplit ( x , "_" ) ) 
      df$mouse <- nm[1]
      df$nuclei <- nm[2]
      df$channel <- nm[3]
      df
    }
  } ) 
) 

str ( CoordDF ) 
CoordDF$mouse <- factor ( unique ( CoordDF$mouse ) ) 
CoordDF$nuclei <- factor ( CoordDF$nuclei , levels = c ( 'PVT' , 'MD' , 'IMD' , 'CL' , 'PCN' , 'CM' ) ) 
CoordDF$channel <- factor ( CoordDF$channel , levels = c ( '488' , '445' , '561' ) ) # Green > Blue > red injection site from medial to lateral for mouse 691386 694513
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
# "#0072B2" "#E69F00" "#009E73" "#D55E00" "#56B4E9" "#F0E442" "#CC79A7"

#––––––––––––––––––––––––––––––––––––––––––––––––––––#



### Multivariate plot –––––––––––––––––––––––––––––––#


splom ( ~ CoordDF[1:3] , CoordDF , subset = c ( mouse == '691386' & nuclei == 'MD' ) , groups = channel , col = colblind[ c ( 3 , 1 , 4 ) ] , pch = 20 )

cloud ( DV ~ AP * ML , CoordDF , subset = c ( mouse == '691386' & nuclei == 'MD' ) , groups = channel , 
        col = colblind[ c ( 3 , 1 , 4 ) ] , pch = 20 ,
        scales = list ( arrows = F )
)





#––––––––––––––––––––––––––––––––––––––––––––––––––––#






##### MOUSE 691386 ––––––––––––––––––––––––––––––#####




### 3D dynamic plot –––––––––––––––––––––––––––––––––#

ThalNuc <- levels ( CoordDF$nuclei ) [c ( 5 , 6 , 3 , 4 , 2 , 1 ) ]

plotList <-
  lapply ( 1:length ( ThalNuc ) , \ ( x ) {
    p <- plot_ly ( subset ( CoordDF , mouse == '691386' & nuclei == ThalNuc[x] & ML > 230 ) ,
                   type = 'scatter3d' , x = ~AP , y = ~ML , z = ~DV , color = ~channel ,
                   colors = colblind[ c ( 3 , 1 , 4 ) ] , scene = if ( x == 1 ) "scene" else paste0( "scene", x ) ,
                   mode = "markers" , marker = list ( size = 5 , symbol = "circle" )
                   
    )
    p
  } )


plotList[[ length ( plotList ) + 1 ]] <-
  plot_ly ( subset ( CoordDF , mouse == '691386' & nuclei %in% ThalNuc & ML > 230 ) ,
            type = 'scatter3d' , x = ~AP , y = ~ML , z = ~DV , color = ~channel ,
            colors = colblind[ c ( 3 , 1 , 4 ) ] , scene = 'scene7' ,
            mode = "markers" , marker = list ( size = 5 , symbol = "circle" )
            
  )


cam <- list ( eye = list ( x = 2 , y = 0.1 , z = 0.2 ) )
xax <- list ( title =  'AnteroPosterior' , range = c ( 220 , 300 ) )
yax <- list ( title = 'MedioLateral' , range = c ( 230 , 280 ) , autorange = 'reversed' )
zax <- list ( title = 'Dorsoventral' , range = c ( 120 , 200 ) , autorange = 'reversed' )
aspm <- 'manual'
asp <- list ( x = 0.8 , y = 0.5 , z = 0.8 )
annot <- function ( txt ){
  list ( list ( text = txt ,
                x = 230 , y = 270 , z = 195 ,
                xref = 'AP' ,
                yref = 'ML' ,
                zref = 'DV' ,
                showarrow = F
  ) )}

fig <- subplot ( plotList[1:7] , nrows = 2 ) %>%
  layout (
    # First row
    scene = list ( xaxis = xax , yaxis = yax , zaxis = zax ,
                   domain = list ( x = c ( 0 , 0.25 ) , y = c ( 0 , 0.5 ) ) ,
                   camera = cam , aspectmode = aspm , aspectratio = asp ,
                   annotations = annot ('PCM') ) , # [1 , 1]
    scene2 = list ( xaxis = xax , yaxis = yax , zaxis = zax ,
                    domain = list ( x = c ( 0.25 , 0.5 ) , y = c ( 0 , 0.5 ) ) ,
                    camera = cam , aspectmode = aspm , aspectratio = asp ,
                    annotations = annot ('CM') ) , # [1 , 2]
    scene3 = list ( xaxis = xax , yaxis = yax , zaxis = zax ,
                    domain = list ( x = c ( 0.5 , 0.75 ) , y = c ( 0 , 0.5 ) ) ,
                    camera = cam , aspectmode = aspm , aspectratio = asp ,
                    annotations = annot ('IMD') ) , # [1 , 3]
    # [1 , 4] is intentionally skipped ( no scene4 )
    
    # Second row
    scene4 = list ( xaxis = xax , yaxis = yax , zaxis = zax ,
                    domain = list ( x = c ( 0 , 0.25 ) , y = c ( 0.5 , 1 ) ) ,
                    camera = cam , aspectmode = aspm , aspectratio = asp ,
                    annotations = annot ('CL') ) , # [2 , 1]
    scene5 = list ( xaxis = xax , yaxis = yax , zaxis = zax ,
                    domain = list ( x = c ( 0.25 , 0.5 ) , y = c ( 0.5 , 1 ) ) ,
                    camera = cam , aspectmode = aspm , aspectratio = asp ,
                    annotations = annot ('MD') ) , # [2 , 2]
    scene6 = list ( xaxis = xax , yaxis = yax , zaxis = zax ,
                    domain = list ( x = c ( 0.5 , 0.75 ) , y = c ( 0.5 , 1 ) ) ,
                    camera = cam , aspectmode = aspm , aspectratio = asp ,
                    annotations = annot ('PVT') ) , # [2 , 3]
    scene7 = list ( xaxis = xax , yaxis = yax , zaxis = zax ,
                    domain = list ( x = c ( 0.75 , 1 ) , y = c ( 0.25 , 0.75 ) ) ,
                    camera = cam , aspectmode = aspm , aspectratio = asp ,
                    annotations = annot ('All') ) # [2 , 4]
  )


fig <- fig %>%
  htmlwidgets::onRender (
    "function ( x , el ) {
 x.on ( 'plotly_relayout' , function ( d ) {
 const camera = Object.keys ( d ) .filter ( ( key ) => /\\.camera$/.test ( key ) ) ;
 if ( camera.length ) {
 const scenes = Object.keys ( x.layout ) .filter ( ( key ) => /^scene\\d*/.test ( key ) ) ;
 const new_layout = {};
 scenes.forEach ( key => {
 new_layout[key] = {...x.layout[key] , camera: {...d[camera]}};
 } ) ;
 Plotly.relayout ( x , new_layout ) ;
 }
 } ) ;
 }" )

fig


saveWidget ( fig, '/Users/loicmagrou/_Julia_/ThalamusInTheMiddle/plots/Theoretical Anatomy/ThalamusLabelledCells_Scatter3D_M691386_v1.0_250825.html' )

#––––––––––––––––––––––––––––––––––––––––––––––––––––#





##### MOUSE 694513 ––––––––––––––––––––––––––––––#####




### 3D dynamic plot –––––––––––––––––––––––––––––––––#

# ThalNuc <- levels ( CoordDF$nuclei ) [c ( 1 , 2 , 4 , 3 , 6 , 5 ) ]
# 
# plotList <- 
#   lapply ( 1:length ( ThalNuc ) , \ ( x ) {
#     p <- plot_ly ( subset ( CoordDF , mouse == '694513' & nuclei == ThalNuc[x] & ML > 230 ) , 
#                    type = 'scatter3d' , x = ~AP , y = ~ML , z = ~DV , color = ~channel , 
#                    colors = colblind[ c ( 3 , 1 , 4 ) ] , scene = if ( x == 1 ) "scene" else paste0( "scene", x ) , 
#                    mode = "markers" , marker = list ( size = 5 , symbol = "circle" ) 
#                    
#     ) 
#     p
#   } ) 
# 
# 
# plotList[[ length ( plotList ) + 1 ]] <- 
#   plot_ly ( subset ( CoordDF , mouse == '694513' & nuclei %in% ThalNuc & ML > 230 ) , 
#             type = 'scatter3d' , x = ~AP , y = ~ML , z = ~DV , color = ~channel , 
#             colors = colblind[ c ( 3 , 1 , 4 ) ] , scene = 'scene7' , 
#             mode = "markers" , marker = list ( size = 5 , symbol = "circle" ) 
#             
#   ) 
# 
# 
# cam <- list ( eye = list ( x = 2 , y = 0.1 , z = 0.2 ) )
# xax <- list ( title =  'AnteroPosterior' , range = c ( 220 , 300 ) )
# yax <- list ( title = 'MedioLateral' , range = c ( 230 , 280 ) )
# zax <- list ( title = 'Dorsoventral' , range = c ( 120 , 200 ) ) 
# aspm <- 'manual'
# asp <- list ( x = 0.8 , y = 0.5 , z = 0.8 )
# annot <- function ( txt ){
#   list ( list ( text = txt ,
#                 x = 230 , y = 270 , z = 195 ,
#                 xref = 'AP' ,
#                 yref = 'ML' ,
#                 zref = 'DV' ,
#                 showarrow = F
#   ) )}
# 
# fig <- subplot ( plotList[1:7] , nrows = 2 ) %>%
#   layout ( 
#     # First row
#     scene = list ( xaxis = xax , yaxis = yax , zaxis = zax , 
#                    domain = list ( x = c ( 0 , 0.25 ) , y = c ( 0 , 0.5 ) ) , 
#                    camera = cam , aspectmode = aspm , aspectratio = asp , 
#                    annotations = annot ('PVT') ) , # [1 , 1]
#     scene2 = list ( xaxis = xax , yaxis = yax , zaxis = zax , 
#                     domain = list ( x = c ( 0.25 , 0.5 ) , y = c ( 0 , 0.5 ) ) ,
#                     camera = cam , aspectmode = aspm , aspectratio = asp ,
#                     annotations = annot ('MD') ) , # [1 , 2]
#     scene3 = list ( xaxis = xax , yaxis = yax , zaxis = zax , 
#                     domain = list ( x = c ( 0.5 , 0.75 ) , y = c ( 0 , 0.5 ) ) ,
#                     camera = cam , aspectmode = aspm , aspectratio = asp ,
#                     annotations = annot ('CL') ) , # [1 , 3]
#     # [1 , 4] is intentionally skipped ( no scene4 ) 
#     
#     # Second row
#     scene4 = list ( xaxis = xax , yaxis = yax , zaxis = zax , 
#                     domain = list ( x = c ( 0 , 0.25 ) , y = c ( 0.5 , 1 ) ) ,
#                     camera = cam , aspectmode = aspm , aspectratio = asp ,
#                     annotations = annot ('IMD') ) , # [2 , 1]
#     scene5 = list ( xaxis = xax , yaxis = yax , zaxis = zax , 
#                     domain = list ( x = c ( 0.25 , 0.5 ) , y = c ( 0.5 , 1 ) ) ,
#                     camera = cam , aspectmode = aspm , aspectratio = asp ,
#                     annotations = annot ('CM') ) , # [2 , 2]
#     scene6 = list ( xaxis = xax , yaxis = yax , zaxis = zax , 
#                     domain = list ( x = c ( 0.5 , 0.75 ) , y = c ( 0.5 , 1 ) ) ,
#                     camera = cam , aspectmode = aspm , aspectratio = asp ,
#                     annotations = annot ('PCN') ) , # [2 , 3]
#     scene7 = list ( xaxis = xax , yaxis = yax , zaxis = zax , 
#                     domain = list ( x = c ( 0.75 , 1 ) , y = c ( 0.25 , 0.75 ) ) ,
#                     camera = cam , aspectmode = aspm , aspectratio = asp ,
#                     annotations = annot ('All') ) # [2 , 4]
#   ) 
# 
# 
# fig <- fig %>% 
#   htmlwidgets::onRender ( 
#     "function ( x , el ) {
#  x.on ( 'plotly_relayout' , function ( d ) {
#  const camera = Object.keys ( d ) .filter ( ( key ) => /\\.camera$/.test ( key ) ) ;
#  if ( camera.length ) {
#  const scenes = Object.keys ( x.layout ) .filter ( ( key ) => /^scene\\d*/.test ( key ) ) ;
#  const new_layout = {};
#  scenes.forEach ( key => {
#  new_layout[key] = {...x.layout[key] , camera: {...d[camera]}};
#  } ) ;
#  Plotly.relayout ( x , new_layout ) ;
#  }
#  } ) ;
#  }" ) 
# 
# fig
# 
# 
# saveWidget ( fig, '/Users/loicmagrou/_Julia_/ThalamusInTheMiddle/plots/Theoretical Anatomy/ThalamusLabelledCells_Scatter3D_M694513_v1.0_250918.html' )

#––––––––––––––––––––––––––––––––––––––––––––––––––––#
