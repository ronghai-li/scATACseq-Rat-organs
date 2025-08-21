#---------------------------------------------------------------------------------------#
#########------------------------Define functions-------------------------------#########
#---------------------------------------------------------------------------------------#
##Part of the function references：
#https://github.com/winstonbecker/scCRC_continuum?tab=readme-ov-file
#https://github.com/GreenleafLab/scScalpChromatin/blob/main/

{
  sample_name <- "all_samples"
  sample_list <- c("Heart_AN49", "Heart_AN50", 
                     "Kidney_AN10", "Kidney_AN7", "Kidney_AN8", 
                     "Liver_AN1", "Liver_AN2", "Liver_AN3", 
                     "Lung_AN19", "Lung_AN20", "Lung_AN21", 
                     "Ovary_AN64", "Ovary_AN65", 
                     "Pancreas_AN37", "Pancreas_AN38", "Pancreas_AN39",
                     "Spleen_AN13", "Spleen_AN14", "Spleen_AN15", 
                     "Thymus_AN34", "Thymus_AN35", "Thymus_AN36", 
                     "Thyroid_AN43", "Thyroid_AN44", "Thyroid_AN45")
  
  sample_list <- c(
                   "Thyroid_AN43"="#a61e4d",
                   "Thyroid_AN44"="#a61e4d", 
                   "Thyroid_AN45"="#a61e4d",
                   "Thymus_AN34"="#E64640",
                   "Thymus_AN35"="#E64640",
                   "Thymus_AN36"="#E64640",
                   "Heart_AN49"="#FB9649",
                   "Heart_AN50"="#FB9649",
                   "Lung_AN19"="#FFE680",
                   "Lung_AN20"="#FFE680", 
                   "Lung_AN21"="#FFE680",
                   "Liver_AN1"="#53C292",
                   "Liver_AN2"="#53C292",
                   "Liver_AN3"="#53C292",
                   "Spleen_AN13"="#0b7285",
                   "Spleen_AN14"="#0b7285",
                   "Spleen_AN15"="#0b7285",
                   "Kidney_AN10"="#1864ab",
                   "Kidney_AN7"="#1864ab",
                   "Kidney_AN8"="#1864ab",
                   "Pancreas_AN37"="#605CB8",
                   "Pancreas_AN38"="#605CB8",
                   "Pancreas_AN39"="#605CB8",
                   "Ovary_AN64"="#862e9c",
                   "Ovary_AN65"="#862e9c"
                   )
  
  Organs_list <- c("Thyroid"="#a61e4d",
                 "Thymus"="#E64640",
                 "Heart"="#FB9649",
                 "Lung"="#FFE680",
                 "Liver"="#53C292",
                 "Spleen"="#0b7285",
                 "Kidney"="#1864ab",
                 "Pancreas"="#605CB8",
                 "Ovary"="#862e9c")
  
  Major_list <- c( "Endocrine"="#F66216",
                   "Endothelial"="#CA5D7C",
                   "Epithelial"="#FCD06A",
                   "Immune"="#615985",
                   "Muscle"="#95DDDF",
                   "Unknown"="#02117E",
                   "Stromal"="#68A3A2")
  
  Celltypes  <- c(
                 "Skeletal_Myh2"= "#b13a63", #Thyroi
                 "Skeletal_Myh4"= "#bc567a", #Thyroid
                 "Thyroid_Endothelial" = "#c77290",#Thyroid
                 "Neuroendocrine" = "#d38fa6",#Thyroid
                 "Skeletal_Myh7"= "#deabbc",#Thyroid
                 "Thyroid_Stromal"= "#e9c7d3", #Thyroid
                 "Follicular"= "#f4e3e9",#Thyroid
                
                 "Immature_T_cell"= "#e95d58",#Thymus
                 "DN4_thymocyte_1"= "#ec7470", #Thymus
                 "Thymic_progenitors"= "#ef8b88",#Thymus
                 "DN4_thymocyte_2"= "#f3a3a0", #Thymus
                 "Thymic_Stromal"= "#f6bab7",#Thymus
                 "Thymocyte"= "#f9d1cf", #Thymus
                 "Thymic_T_Cell"= "#fce8e7",#Thymus

                 "Cardiac_Endothelial"= "#fca360", #Heart
                 "Cardiac_Stromal"= "#fcb077", #Heart
                 "Cardiomyocyte"= "#fdbd8d", #Heart
                 "Cardiac_Pericyte" = "#fdcba4",#Heart
                 "Cardiac_Macrophage"= "#fed8bb", #Heart
                 "Cardiac_Dendritic_cell"= "#fee5d2",#Heart
                 "Cardiac_Unknown" = "#fff2e8",#Heart
                 
                 "AT2"= "#fcd476",#Lung
                 "Pulmonary_Endothelial"= "#fdd883",#Lung
                 "AT1"= "#fddc8f", #Lung
                 "Pulmonary_Ciliated"= "#fde09c",#Lung
                 "Pulmonary_Myofibroblast"= "#fde4a8", #Lung
                 "Pulmonary_NKT"= "#fee8b5", #Lung
                 "Pulmonary_Macrophage"= "#feebc1",#Lung
                 "Goblet"= "#feefcd", #Lung
                 "Pulmonary_B_cell"= "#fef3da",#Lung
                 "Pulmonary_Epithelial"= "#fff7e6", #Lung
                 "Pulmonary_Smooth_muscle_cell"= "#fffbf3", #Lung
                 
                 "Portal_Hep"= "#66c99e",#Liver
                 "Central_Hep"= "#79d0aa", #Liver
                 "Liver_Endothelial"= "#8cd6b6", #Liver
                 "Kupffer_cell"= "#9fddc2",#Liver
                 "Hepatic_Stellate"= "#b3e4cf", #Liver
                 "Midzonal_Hep"= "#c6ebdb", #Liver
                 "Liver_NK_cell"= "#d9f1e7",#Liver
                 "Liver_B_cell"= "#ecf8f3",#Liver
                 
                 "Splenic_B_cell"= "#238091",#Spleen
                 "Splenic_T_cell"= "#3c8e9d", #Spleen
                 "Splenic_Macrophage"= "#549caa",#Spleen
                 "Splenic_Pre_B_cell"= "#6daab6", #Spleen
                 "Splenic_Pre_T_cell"= "#85b9c2", #Spleen
                 "Splenic_Stromal"= "#9dc7ce",#Spleen
                 "Splenic_Endothelial"= "#b6d5da",#Spleen
                 "Plasma_BC"= "#cee3e7",#Spleen
                 "Cycling"= "#e7f1f3",#Spleen
                 
                 "Proximal_Tubule_S2"= "#276eb1",#Kidney
                 "Proximal_Tubule_S1"= "#3779b6", #Kidney
                 "Proximal_Tubule_S3"= "#4683bc", #Kidney
                 "Ascending_LOH"= "#568dc1",#Kidney
                 "Renal_Endothelial"= "#6598c7", #Kidney
                 "Connecting_Tubule"= "#74a2cd", #Kidney
                 "Distal_Convoluted_Tubule"= "#84acd2",#Kidney
                 "Intercalated"= "#93b7d8",#Kidney
                 "Mesangial"= "#a3c1dd",#Kidney
                 "Principal"= "#b2cbe3", #Kidney
                 "Podocyte"= "#c1d6e9",#Kidney
                 "Renal_Macrophage"= "#d1e0ee",#Kidney
                 "Renal_T_cell"= "#e0eaf4",#Kidney
                 "Renal_Epithelial"= "#f0f5f9", #Kidney
                 
                 "Acinar"= "#807dc6",#Pancreas
                 "Delta"= "#a09dd4", #Pancreas
                 "Pancreatic_Stellate"= "#bfbee3",#Pancreas
                 "Pancreatic_Macrophage"= "#dfdef1",#Pancreas
                 
                 "Theca"= "#9141a5",#Ovary
                 "Luteal_cell"= "#9c54ae", #Ovary
                 "Ovarian_Endothelial"= "#a767b7",#Ovary
                 "Ovarian_Stromal"= "#b27ac0", #Ovary
                 "Granulosa"= "#bd8dc9", #Ovary
                 "Pre_Luteal_cell"= "#c8a0d2",#Ovary
                 "Ovarian_Monocyte"= "#d3b3db",#Ovary
                 "Ovarian_Dendritic_cell"= "#dec6e4",#Ovary
                 "Surface_epithelial"= "#e9d9ed",#Ovary
                 "Ovarian_Macrophage"= "#f4ecf6" #Ovary
  )
  
  
  qc_plots <- function(proj, samplename, groupBy, name){
    pdf(paste0("./violin_", name, "_", samplename, "_.pdf"), height = 3, width = 5, onefile=F)
    p <- plotGroups(
      ArchRProj = proj,
      groupBy = groupBy,
      colorBy = "cellColData",
      name = name,
      plotAs = "violin",
      alpha = 1,
      addBoxPlot = F,
      pal = sample_list
    ) +
      geom_violin(aes(fill = groupBy), color = "black", alpha = 0.3, scale = "area") +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 1, position = position_dodge(width = 0.9), 
                   color = "black", fill = "black") +
      geom_point(stat = "summary", fun = "median", color = "white", size = 0.3) +
      theme_ArchR() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    print(p)
    dev.off()
  }
  

    ArchR_standard_clustering <- function(proj, iterations_times, LSI_resolution,dimsToUse,cluster,Cluster_resolution){
      proj <- addIterativeLSI(
        ArchRProj = proj,
        useMatrix = "TileMatrix",
        name = "IterativeLSI",
        iterations = iterations_times,
        clusterParams = list(
          resolution = LSI_resolution,#cluster_resolution is a vector
          sampleCells = 10000,
          n.start = 10 ),
        varFeatures = 25000,
        dimsToUse = dimsToUse,
        force = TRUE
      )
      
      if (cluster){
        proj <- addClusters(
          input = proj,
          reducedDims = "IterativeLSI",
          method = "Seurat",
          name = "Clusters",
          resolution = Cluster_resolution,
          dimsToUse = dimsToUse,
          force = TRUE
        )
      }
      
      proj <- addUMAP(
        ArchRProj = proj,
        reducedDims = "IterativeLSI",
        name = "UMAP",
        nNeighbors = 60,
        minDist = 0.6,
        dimsToUse = dimsToUse,
        metric = "cosine",
        force = TRUE
      )
      return(proj)
    }

    cmaps_BOR <- list(
      # Many of these adapted from ArchR ColorPalettes.R by Jeff Granja or colors.R from BuenColors
      # https://github.com/GreenleafLab/ArchR/blob/master/R/ColorPalettes.R
      # https://github.com/caleblareau/BuenColors/blob/master/R/colors.R
      
      ## Sequential colormaps:
      solarExtra = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', 
                     '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D'), #buencolors
      sunrise = c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                  "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"),
      horizon = c('#000075', '#2E00FF', '#9408F7', '#C729D6', '#FA4AB5', 
                  '#FF6A95', '#FF8B74', '#FFAC53', '#FFCD32', '#FFFF60'),
      horizonExtra =c("#000436", "#021EA9", "#1632FB", "#6E34FC", "#C732D5",
                      "#FD619D", "#FF9965", "#FFD32B", "#FFFC5A"),
      blueYellow = c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                     "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"),
      sambaNight = c('#1873CC','#1798E5','#00BFFF','#4AC596','#00CC00',
                     '#A2E700','#FFFF00','#FFD200','#FFA500'), #buencolors
      wolfgang_basic = c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", 
                         "#1D91C0", "#225EA8", "#253494", "#081D58"), #buencolors
      wolfgang_extra = c("#FFFFFF", "#FCFED3", "#E3F4B1", "#ABDEB6", "#60C1BF", 
                         "#2A9EC1", "#206AAD", "#243996", "#081D58"), #buencolors
      whitePurple = c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6',
                      '#8c6bb1','#88419d','#810f7c','#4d004b'),
      whiteBlue = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf',
                    '#3690c0','#0570b0','#045a8d','#023858'),
      whiteViolet = c('#FFF7F3', '#FDE0DD', '#FCC5C0', '#FA9FB5', '#F768A1', 
                      '#DD3497', '#AE017E', '#7A0177', '#49006A'),
      comet = c("#E6E7E8","#3A97FF","#8816A7","black"),
      rocket = c("#E6E7E8","#F4855EFF","#BD1656FF","black"),
      rocket1 = c("#E6E7E8","#F4855EFF","#BD1656FF","#67000d"),
      
      flame_flame = c('#000033', '#0000A5', '#1E00FB', '#6F00FD', '#C628D6', 
                      '#FE629D', '#FF9B64', '#FFD52C', '#FFFF5F'), # buencolors
      
      flame_short = c('#000033', '#0000A5', '#1E00FB', '#6F00FD', '#C628D6', 
                      '#FE629D', '#FF9B64', '#FFD52C'), # Stop short of yellow (better for tracks, etc.)
      
      #7-colors
      greenBlue = c('#e0f3db','#ccebc5','#a8ddb5','#4eb3d3','#2b8cbe',
                    '#0868ac','#084081'),
      
      #6-colors
      beach = c("#87D2DB","#5BB1CB","#4F66AF","#F15F30","#F7962E","#FCEE2B"),
      
      #5-colors
      fireworks = c("white","#2488F0","#7F3F98","#E22929","#FCB31A"),
      greyMagma = c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF"),
      fireworks2 = c("black", "#2488F0","#7F3F98","#E22929","#FCB31A"),
      purpleOrange = c("#581845", "#900C3F", "#C70039", "#FF5744", "#FFC30F"),
      beach = c("#87D2DB","#5BB1CB","#4F66AF","#F15F30","#F7962E","#FCEE2B"),
      zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"), #wesanderson
      darjeeling = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"), #wesanderson
      rushmore = c("#E1BD6D", "#EABE94", "#0B775E","#35274A" , "#F2300F"), #wesanderson
      FantasticFox1 = c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20"), #wesanderson
      BottleRocket2 = c("#FAD510", "#CB2314", "#273046", "#354823", "#1E1E1E"), #wesanderson
      Moonrise3 = c("#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B"), #wesanderson
      fireworks = c("white","#2488F0","#7F3F98","#E22929","#FCB31A"),
      
      # Divergent sequential:
      coolwarm = c("#4858A7", "#788FC8", "#D6DAE1", "#F49B7C", "#B51F29"),
      brewer_yes = c("#053061", "#2971B1", "#6AACD0","#C1DDEB", "#F7F7F7", 
                     "#FACDB5", "#E58267", "#BB2933", "#67001F"), #buencolors
      brewer_celsius = c("#313695", "#5083BB", "#8FC3DD", "#D2ECF4", "#FFFFBF", 
                         "#FDD384", "#F88D51", "#DE3F2E", "#A50026"), #buencolors
      flame_blind = c("#0DB2AA", "#0AD7D3", "#00FFFF", "#B1FFFE", "#FFFFFF", 
                      "#FFA3EC", "#FF00D8", "#BD00EC", "#5F00FF"), #buencolors
      solar_flare = c('#3361A5', '#2884E7', '#1BA7FF', '#76CEFF', '#FFFFFF', 
                      '#FFE060', '#FA8E24', '#DA2828', '#A31D1D'), #buencolors
      brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', 
                     '#FACDB5', '#E58267', '#BB2933', '#67001F'), #buencolors
      
      ## Qualitative colormaps:
      
      # see: https://carto.com/carto-colors/
      cartoPrism = c('#7F3C8D', '#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310', 
                     '#008695', '#CF1C90', '#F97B72', '#4B4B8F'),
      cartoSafe = c('#88CCEE', '#CC6677', '#DDCC77', '#117733', '#332288', '#AA4499', '#44AA99',
                    '#999933', '#882255', '#661100', '#6699CC'),
      cartoBold = c('#7F3C8D' ,'#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310',
                    '#008695', '#CF1C90', '#f97b72', '#4b4b8f'),
      cartoAntique = c('#855C75', '#D9AF6B', '#AF6458', '#736F4C', '#526A83', '#625377', '#68855C',
                       '#9C9C5E', '#A06177', '#8C785D', '#467378'),
      cartoPastel = c('#66C5CC', '#F6CF71', '#F89C74', '#DCB0F2', '#87C55F', '#9EB9F3', '#FE88B1',
                      '#C9DB74', '#8BE0A4', '#B497E7', '#D3B484'),
      cartoVivid = c('#E58606', '#5D69B1', '#52BCA3', '#99C945', '#CC61B0', '#24796C', '#DAA51B',
                     '#2F8AC4', '#764E9F', '#ED645A', '#CC3A8E'),
      # 15 color
      circus = c("#D52126", "#88CCEE", "#FEE52C", "#117733", "#CC61B0", "#99C945", "#2F8AC4", "#332288", 
                 "#E68316", "#661101", "#F97B72", "#DDCC77", "#11A579", "#89288F", "#E73F74"),
      iron_man = c('#371377','#7700FF','#9E0142','#FF0080', '#DC494C',"#F88D51","#FAD510","#FFFF5F",'#88CFA4',
                   '#238B45',"#02401B","#0AD7D3","#046C9A", "#A2A475", 'grey35'),
      # The following 3 were designed by Ryan Corces.
      stallion = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B", "#FEE500","#8A9FD1","#C06CAB", "#D8A767",
                   "#90D5E4", "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8", "#6E4B9E","#0C727C", "#7E1416", "#E6C2DC"),
      calm = c("#7DD06F", "#844081", "#688EC1", "#C17E73", "#484125", "#6CD3A7", "#597873","#7B6FD0", "#CF4A31", "#D0CD47",
               "#722A2D", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D", "#403552", "#76D73C", "#96CED5", "#CE54D1", "#C48736"),
      kelly = c("#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A",
                "#FF7A5C", "#53377A", "#FF8E00","#B32851", "#F4C800", "#7F180D", "#93AA00", "#593315", "#F13A13")
    )
    
    theme_BOR <- function(base_size=14, base_family="Helvetica", border = TRUE) {
      library(grid)
      library(ggthemes)
      # Should plots have a bounding border?
      if(border){
        panel.border <- element_rect(fill = NA, color = "black", size = 0.7)
        axis.line <- element_blank()
      }else{
        panel.border <- element_blank()
        axis.line <- element_line(color = "black", size = 0.5)
      }
      
      (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = panel.border,
                axis.title = element_text(size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line = axis.line,
                axis.ticks = element_line(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                legend.direction = "vertical",
                legend.key.size= unit(0.5, "cm"),
                legend.spacing = unit(0, "cm"),
                legend.title = element_text(),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text()
        ))
    }

    fractionXbyY <- function(x, y, add_total=FALSE, xname="x", yname="y", ylab="proportion"){
      # Returns a melted dataframe with the proportion of x by y
      # e.g. if x is cluster and y is samples, will return the proportion of each sample 
      # making up the total amount of each cluster (i.e. the proportional contribution of each y to each x)
      # x, y: paired vectors
      # add_total: bool indicating whether you would like to add a 'total' entry (i.e. proportion of y 
      # across all x)
      # xname, yname: char with x and y labels
      # ylab: char with proportion label
      XbyYdf <- data.frame("xgroup" = x, "ygroup" = y) %>% 
        group_by(xgroup, ygroup) %>% # Group by both x and y
        summarize(n = n()) %>% ungroup() %>% # summarize (i.e. tally for each group(s))
        pivot_wider(names_from=xgroup, values_from=n, values_fill= list(n=0)) %>% # Expand into multiple columns
        as.data.frame()
      
      #sort lable sequence (Modify as needed)
      clusters <- unique(proj$Clusters)
      colOrder <- clusters[order(as.numeric(sub("C", "", clusters)))]
      #XbyYdf$ygroup <- factor(XbyYdf$ygroup, levels = colOrder)
     # XbyYdf <- XbyYdf[order(sapply(XbyYdf$ygroup, function(x) which(x == colOrder))), ]
      #XbyYdf<-XbyYdf[,c("ygroup",seq(0,8,by = 1))]
      XbyYdf<-XbyYdf[,c("ygroup",colOrder)]
      
      #colOrder <- names(Location)
      #XbyYdf$ygroup <- factor(XbyYdf$ygroup, levels = colOrder)
      #XbyYdf <- XbyYdf[order(sapply(XbyYdf$ygroup, function(x) which(x == colOrder))), ]
      #XbyYdf<-XbyYdf[, c("ygroup","S","I","C")]
      
      rownames(XbyYdf) <- XbyYdf[,1]
      XbyYmat <- as.matrix(XbyYdf[,-1]) %>% t()
      
      if(add_total){
        rnms <- rownames(XbyYmat)
        XbyYmat <- rbind(XbyYmat, colSums(XbyYmat))
        rownames(XbyYmat) <- c(rnms, "total")
      }
      XbyYmat <- (XbyYmat/rowSums(XbyYmat)) %>% reshape2::melt()
      colnames(XbyYmat) <- c(xname, yname, ylab) 
      # Force cluster to be qualitative
      XbyYmat[,xname] <- as.factor(XbyYmat[,xname])
      return(XbyYmat)
    }
    
    pairwiseColorInterpolations <- function(cols, colorspace = "Lab"){
      # Get all pairwise interpolations between a vector of colors
      rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
      interpolate <- function(c1, c2, colorspace){
        rgb2hex(colorRamp(c(c1, c2), space = colorspace)(0.5))
      }
      paired <- sapply(cols, function(x) sapply(cols, function(y) interpolate(x, y, colorspace)))
      unique(as.vector(paired))
    }
    
    getColorMap <- function(cmap, n, type='qualitative'){
      stopifnot(n >= 1)
      # Return a character vector of n colors based on
      # the provided colormap. If n > length(cmap), do
      # some smart interpolation to get enough colors
      names(cmap) <- NULL # Having names on colors causes problems for some plotting routines
      
      if(type == 'qualitative'){
        # If qualitative colormap, do 'mostDifferent' interpolation
        if(length(cmap) < n){
          cmap <- mostDifferentColors(
            pairwiseColorInterpolations(cmap), 
            colorspace = "Apple RGB", n = n, startingCols = cmap
          )
        }
        
      }else{
        # Otherwise, return sequential colors based on provided palette
        colfunc <- colorRampPalette(cmap)
        cmap <- colfunc(n)
      }
      cmap[1:n]
    }
    
    mostDifferentColors <- function(cols, n=20, colorspace="Lab", startingCols=NULL){
      stopifnot(length(cols) > n)
      rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue=255)
      
      # Convert sRGB to another colorspace (more 'perceptually uniform' colorspace, e.g. "Lab")
      rgbCols <- t(col2rgb(cols))
      conv <- grDevices::convertColor(rgbCols, from="sRGB", to=colorspace, scale.in=255)
      
      # Now select n 'furthest neighbors' colors
      # This performs an iterative procedure for picking colors that maximize
      # 'distance' to already selected colors. The first color is picked randomly.
      # If starting cols provided, add these to the list of picked cols
      if(!is.null(startingCols)){
        stConv <- grDevices::convertColor(t(col2rgb(startingCols)), from="sRGB", to=colorspace, scale.in=255)
        pickedColors <- list()
        for(i in seq_len(nrow(stConv))){
          pickedColors[[i]] <- stConv[i,]
        }
        remainingColors <- conv
      }else{
        idx <- sample(1:nrow(conv), 1)
        pickedColors <- list(conv[idx,])
        remainingColors <- conv[-idx,]
      }
      pickedLen <- length(pickedColors)
      
      # Iteratively add the furthest color from the selected colors
      for(i in seq(pickedLen, n - 1)){
        distList <- list()
        for(j in seq_along(pickedColors)){
          colJ <- pickedColors[[j]]
          distMat <- dist(rbind(colJ, remainingColors), method="euclidean") %>% as.matrix
          distList[[j]] <- distMat[2:nrow(distMat),1]
        }
        # Maximize the minimum distance between each color
        distMat <- do.call(cbind, distList)
        distMins <- apply(distMat, 1, FUN = min)
        idx <- which(max(distMins) == distMins)
        pickedColors[[i + 1]] <- remainingColors[idx,]
        remainingColors <- remainingColors[-idx,]
      }
      pickedLab <- do.call(rbind, pickedColors)
      pickedRgb <- round(grDevices::convertColor(pickedLab, from = colorspace, to = "sRGB", scale.out = 255),0)
      hex <- apply(pickedRgb, 1, rgb2hex)
      hex
    }
    
    stackedBarPlot <- function(df, xlab=NULL,ylab=NULL,xcol = 1, fillcol = 2, ycol = 3, cmap = NULL, border_color="black", covarLabel = "", namedColors=FALSE, barwidth=0.5){
      # Plot a stacked bar plot
      # Expects a 'melted' dataframe/matrix as input
      nsamp <- length(unique((df[,xcol])))
      # Assume that we want to show all xaxis labels
      xID <- unique((df[,xcol]))
      
      # Fix colormap if provided
      if(!namedColors){
        if(!is.null(cmap)){
          cmap <- getColorMap(cmap, n = length(unique((df[,fillcol]))))
        }else{
          cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique((df[,fillcol]))))
        }
      }
      
      p <- (
        ggplot(df, aes(x=df[,xcol], y=df[,ycol], fill=df[,fillcol]))
        + geom_bar(stat = "identity", position="fill", width=barwidth, color=border_color)
        + xlab(xlab)
        + ylab(ylab)
        + theme_BOR(border=FALSE)
        + theme(panel.grid.major=element_blank(), 
                panel.grid.minor= element_blank(), 
                plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
                aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
                axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) 
        + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
      )
      
      # If colormap provided, update colors
      if(namedColors){
        # Named colormap corresponding to discrete values in third column
        p <- p + scale_color_manual(values = cmap, limits = names(cmap), name = covarLabel)
        p <- p + scale_fill_manual(values = cmap, limits = names(cmap), name = covarLabel)
      }else{
        p <- p + scale_color_manual(values = cmap, name = covarLabel)
        p <- p + scale_fill_manual(values = cmap, name = covarLabel)
      }
      p 
    }
    
    qcBarPlot <- function(df, cmap = "#bfbea4", border_color="black", barwidth=0.5){
      # Plot a bar plot (df is a 2+ column dataframe with column 1 = x and column 2 = y)
      p <- (
        ggplot(df, aes(x=df[,1], y=df[,2]))
        + geom_bar(stat = "identity", fill = cmap, width=barwidth, color=border_color)
        + scale_fill_manual(values = cmap)
        + xlab(colnames(df)[1])
        + ylab(colnames(df)[2])
        + theme_BOR(border=FALSE)
        + theme(panel.grid.major=element_blank(), 
                panel.grid.minor= element_blank(), 
                plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
                #aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
                axis.text.x = element_text(angle = 90, hjust = 1)) 
        + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
      )
      p
    }
    
    ArchR_standard_plot <- function(proj,sample_name,clusterName="Clusters",embedding="UMAP"){
      
      message("plotEmbedding...")
      p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = embedding,verbose = FALSE,size = 0.03)
      p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = embedding)
      p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = embedding,plotAs = 'points')
      p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "nFrags", embedding = embedding,plotAs = 'points')
      p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Organs",pal = Organs_list,embedding =embedding,plotAs = 'points')
      p7 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DoubletEnrichment", embedding = embedding,plotAs = 'points')
      p8 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DoubletScore", embedding = embedding,plotAs = 'points')
      plotPDF(p1,p2,p3,p4,p5,p7,p8,name =paste0("./umap_plot_", sample_name ,".pdf"), ArchRProj = proj, addDOC = FALSE, width = 12, height = 12)
      
      message("Calculate which samples reside in which clusters...")
      cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
      cM <- cM / Matrix::rowSums(cM)
      p9 <- pheatmap::pheatmap(
        mat = as.matrix(cM),
        color = paletteContinuous("whiteBlue"),
        border_color = "black"
      )
      plotPDF(p9,name = paste0("./sample_with_cluster_", sample_name ,".pdf"), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
      
      message("Bar plot cluster counts...")
      clustVec <- getCellColData(proj)[[clusterName]] %>% gsub("[^[:digit:].]", "", .) %>% as.numeric()
      tabDF <- base::table(clustVec) %>% as.data.frame
      colnames(tabDF) <- c("Clusters", "count")
      
      pdf(paste0("cluster_Bar_Plot_",sample_name,".pdf"))
      print(qcBarPlot(tabDF))
      dev.off()
      
      #message("Stacked bar plot fraction samples in clusters...")
      # clustBySamp <- fractionXbyY(clustVec, proj$Sample, add_total=TRUE, xname="Cluster", yname="Sample")
      # 
      #pdf(paste0("Clust_By_Sample_BarPlot_",sample_name,".pdf"),width = 15, height = 10)
      # print(stackedBarPlot(clustBySamp,xlab="Cluster",ylab="Sample",barwidth=0.9))
      # dev.off()
      
      #message("Stacked bar plot fraction group in clusters...")
      # clustBySamp <- fractionXbyY(clustVec, proj$condition, add_total=TRUE, xname="Cluster", yname="condition")
      
      # pdf(paste0("Clust_By_condition_BarPlot_",sample_name,".pdf"),width = 15, height = 10)
      # print(stackedBarPlot(clustBySamp,xlab="Cluster",ylab="condition",barwidth=0.9))
      # dev.off()
      
      return(proj)
    }
    
    # Cluster visualization helpers
    relabelClusters <- function(proj, clusterName="Clusters"){
      # Relabel clusters to be ordered by cluster size
      
      ogClusters <- getCellColData(proj)[[clusterName]]
      tabDF <- base::table(ogClusters) %>% as.data.frame
      colnames(tabDF) <- c("Clusters", "count")
      tabDF["NewClusters"] <- rank(-tabDF$count)
      swapVec <- paste0("C", tabDF$NewClusters)
      names(swapVec) <- tabDF$Clusters
      
      # Now replace cluster names
      newClust <- sapply(ogClusters, function(x) swapVec[x]) %>% unname()
      proj <- addCellColData(proj, data=newClust, name=clusterName, cells=getCellNames(proj), force=TRUE)
      return(proj)
    }
    
    ArchR_feature_plot <- function(proj, sample_name, embedding, cell_type, markers,colorBy = "GeneScoreMatrix"){
      p <- plotEmbedding(ArchRProj = proj,
                         colorBy = colorBy, 
                         name = markers, 
                         embedding = embedding,
                         imputeWeights = getImputeWeights(proj),
                         log2Norm =  T,
                         plotAs = 'points',
                         pal=ArchRPalettes$greyMagma
      )
      
      plotPDF(plotList = p, 
              name = paste0("feature_plot_",sample_name,"_",cell_type,".pdf"),
              ArchRProj = proj, 
              addDOC = FALSE, 
              width = 5, 
              height = 5,
      )
    }
    
    jaccardIndex <- function(mat, i, j){
      # Calculate Jaccard Index between row i and column j in matrix mat
      # (Matrix is an intersection matrix of categories in rows i and columns j)
      AiB <- mat[i,j]
      AuB <- sum(mat[i,]) + sum(mat[,j]) - AiB
      AiB/AuB
    }
    
    
    prettyOrderMat <- function(mat, scale=TRUE, cutOff=1, lmat=NULL, clusterCols=TRUE){
      # Reorder mat in a prettier way for plotting
      # Adapted from Jeff's ArchR .binarySort
      # mat = matrix (like) object to sort
      # scale = should mat be scaled before building logical mat
      # cutOff = cutoff for lmat
      # lmat = logical matrix for ordering rows (binary sorting)
      # clusterCols = should columns be clustered?
      mat <- as.matrix(mat)
      
      if(is.null(lmat)){
        # Compute row Z-scores
        if(scale){
          lmat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
        }else{
          lmat <- mat
        }
        # Logical matrix of values passing cutoff 
        lmat <- lmat >= cutOff
      }
      
      # Transpose:
      mat <- t(mat)
      lmat <- t(lmat)
      
      # Identify column ordering:
      if(clusterCols){
        hc <- hclust(dist(mat))
        colIdx <- hc$order
        mat <- t(mat[colIdx,])
        lmat <- t(lmat[colIdx,])
      }else{
        mat <- t(mat)
        lmat <- t(lmat)
        hc <- NULL
      }
      
      # Identify row ordering:
      rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
      mat <- mat[rowIdx,]
      
      return(list(mat=mat, hclust=hc))
    }
    
    # Heatmap wrapper:
    BORHeatmap <- function(
    mat, # Data to plot (matrix or dataframe)
    limits = NULL, # Enforced limits for colormap (2 dimensional array)
    clusterCols = TRUE, # Should columns be clustered
    clusterRows = TRUE, # Should rows be clustered
    labelCols = FALSE, # Should columns be labeled
    labelRows = FALSE, # Should rows be labeled
    dataColors = NULL, # Colormap for plotting data
    dataColorMidPoint = NULL, # The data value to be the middle of the color map
    customRowLabel = NULL,
    customRowLabelIDs = NULL,
    customColLabel = NULL,
    customColLabelIDs = NULL,
    customLabelWidth = 0.15,
    useRaster = TRUE, # Should heatmap be rasterized
    rasterDevice = "CairoPNG",
    rasterQuality = 5, # Raster quality. Higher is {better?}
    fontSize = 6, # Font size for labels
    showColDendrogram = FALSE, # Should the column dendrogram be shown
    showRowDendrogram = FALSE, # Should the row dendrogram be shown
    borderColor = NA, # Color for lines between cells
    mapname = " ", # 'Name' to give heatmap
    legendTitle = " ", # Name of legend
    ...
    ){
      
      #Packages
      suppressPackageStartupMessages(require(ComplexHeatmap))
      suppressPackageStartupMessages(require(circlize))
      
      # Make sure mat is actually a matrix
      if(!is.matrix(mat)){
        message("'mat' needs to be a matrix. Converting...")
        mat <- as.matrix(mat)
      }
      
      # Prepare color function
      if(!is.null(limits)){
        ll <- limits[1]
        ul <- limits[2]
      }else{
        ll <- min(mat, na.rm=TRUE)
        ul <- max(mat, na.rm=TRUE)
      }
      # If no colormap provided, use solarExtra
      if(is.null(dataColors)){
        dataColors <- c("1"='#3361A5', "2"='#248AF3', "3"='#14B3FF', 
                        "4"='#88CEEF', "5"='#C1D5DC', "6"='#EAD397', 
                        "7"='#FDB31A', "8"= '#E42A2A', "9"='#A31D1D')
      }
      dataColFun <- makeColFun(ll, ul, dataColors, midpoint = dataColorMidPoint)
      
      message("Preparing Heatmap...")
      hm <- Heatmap(
        # Main components:
        matrix = mat,
        name = mapname,
        col = dataColFun,
        
        # Legend options:
        heatmap_legend_param = list(
          color_bar = "continuous",
          legend_direction = "vertical",
          legend_width = unit(1, "cm"),
          title = legendTitle
        ),
        rect_gp = gpar(col = borderColor), 
        
        # Column options:
        show_column_names = labelCols,
        cluster_columns = clusterCols,
        show_column_dend = showColDendrogram,
        clustering_method_columns = "ward.D2",
        #column_names_gp = gpar(fontsize = fontSize), 
        
        # Row options:
        show_row_names = labelRows,
        cluster_rows = clusterRows,
        show_row_dend = showRowDendrogram,
        clustering_method_rows = "ward.D2",
        #row_names_gp = gpar(fontsize = fontSize), 
        
        # Raster info:
        use_raster = useRaster,
        raster_device = rasterDevice,
        raster_quality = rasterQuality,
        
        # Other
        ...
      )
      
      # Add row labels if provided:
      if(!is.null(customRowLabel)){
        if(is.null(customRowLabelIDs)){
          customRowLabelIDs <- rownames(mat)[customRowLabel]
        }
        hm <- hm + rowAnnotation(
          link = anno_mark(at = customRowLabel, labels = customRowLabelIDs, labels_gp = gpar(fontsize = fontSize)),
          width = unit(customLabelWidth, "cm") + max_text_width(customRowLabelIDs)
        )
      }
      
      return(hm)
    }
    
    # This is used primarily for making colormaps for ComplexHeatmap
    makeColFun <- function(start, end, cmap, midpoint = NULL){
      # Make a color ramp function from provided start and end breaks,
      # and optionally a midpoint
      cmapLen <- length(cmap)
      if(!is.null(midpoint)){
        interpolate <- function(c1, c2, colorspace = "Lab"){
          rgb2hex(colorRamp(c(c1, c2), space = colorspace)(0.5))
        }
        if(length(cmap) %% 2 == 0){
          # Interpolate middle colors if necessary to get midpoint
          preMidIdx <- floor(cmapLen / 2)
          midCol <- interpolate(cmap[preMidIdx], cmap[preMidIdx + 1])
          cmap <- c(cmap[1:preMidIdx], midCol, cmap[(preMidIdx + 1):cmapLen])
          cmapLen <- length(cmap)
        }
        midIdx <- ceiling(cmapLen / 2)
        breaks <- c(seq(start, midpoint, length.out = midIdx), seq(midpoint, end, length.out = midIdx)[2:midIdx])
      } else {
        breaks <- seq(start, end, length.out = cmapLen)
      }
      colorRamp2(breaks, cmap)
    }
    
    grLims <- function(gr){
      # Get the minimum and maximum range from a GR
      if(length(gr) == 0){
        return(NA)
      }
      starts <- start(gr)
      ends <- end(gr)
      c(min(starts, ends), max(starts, ends))
    }
    
    plotClusterQC <- function(obj,pointSize=1.0, barwidth=0.9, sampleCmap=NULL, diseaseCmap=NULL){
      
      # Plot basic clustering plots
      # Set colormap
      qualcmap <- cmaps_BOR$stallion
      quantcmap <- cmaps_BOR$medication
      namedSampCmap <- TRUE
      namedDiseaseCmap <- TRUE
      
      if(is.null(sampleCmap)){
        sampleCmap <- qualcmap
        namedSampCmap <- FALSE
      }
      if(is.null(diseaseCmap)){
        diseaseCmap <- qualcmap
        namedDiseaseCmap <- FALSE
      }
      
      ### Stacked bar plot fraction celltype in condition ###
      clustBySamp <- fractionXbyY(se.integrated$seurat_clusters, se.integrated$species, add_total=F, xname="clusters", yname="species")
      
      pdf(paste0("./species_clusters_fraction.pdf"),width = 10, height = 7,onefile=F)
      print(stackedBarPlot(clustBySamp, xlab="clusters",ylab="species",cmap=sampleCmap, namedColors=namedSampCmap,barwidth=0.9))
      dev.off()
      
      ### Stacked bar plot fraction location in clusters ###
      #clustByDisease <- fractionXbyY(obj$condition, obj$location, add_total=F, xname="group", yname="location")
      
      #pdf(paste0("./location_group_fraction.pdf"))
      #print(stackedBarPlot(clustByDisease,xlab="group",ylab="location",cmap=diseaseCmap, namedColors=namedDiseaseCmap, barwidth=0.9))
      #dev.off()
    }
    
    seurat_feature_plot <- function(MSCT, sample_name, reduction, cell_type, markers){
      # make seurat feature plots for multiple markers and arrange into a grid
      # colon: seurat_object
      # sample_name: included in the plot save name
      # reduction: reduction to plot e.g. UMAP
      # cell_type: the name of the cell type the markers correspond to, will be added to plot name
      # markers: list of arkers to plot
      p1 <- FeaturePlot(MSCT, features = markers, reduction = reduction, combine = FALSE, pt.size = 0.5)
      fix.sc <- scale_colour_gradientn(colours = ArchRPalettes$greyMagma)
      if (length(p1)==1){
        width <- 4
        height <- 4
      } else if (length(p1)==2){
        width <- 8
        height <- 4
      } else if (length(p1)<5){
        width <- 8
        height <- 8
      } else if (length(p1)<7){
        width <- 12
        height <- 8
      } else if (length(p1)<10){
        width <- 12
        height <- 12
      } else if (length(p1)<13){
        width <- 16
        height <- 12
      } else if (length(p1)<17){
        width <- 16
        height <- 16
      }
      pdf(paste0("./", reduction, "_feature_plot_", sample_name, "_", cell_type ,".pdf"), width = width, height = height)
      print(CombinePlots(lapply(p1, function (x) AugmentPlot(x + fix.sc))))
      dev.off()
    }
    
    calcTopGo <- function(
    allGenes, interestingGenes=NULL, pvals=NULL, geneSel=NULL,
    nodeSize=5, ontology="BP",
    alg="weight01", stat="fisher", topNodes=50
    ){
      # Calculate GO term enrichments using topGO on provided data
      # https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
      # allGenes: vector of genenames to be used in GO term search. Expects gene 'symbol' 
      # interestingGenes: predefined list of 'instersting' genes. Incompatible with supplying pvalues.
      # geneSel: function for selecting 'interesting' genes. Can only really be a p-value cutoff...
      # pvals: vector of pvalues corresponding to geneList. If not provided, will assign everything to 1
      # nodeSize: will prune terms that have less than nodeSize number of genes
      # ontology: which GO ontology to use (MF, BP, CC)
      # alg: algorithm to be used for testing GO terms (topGO default is 'weight01')
      # stat: test statistic to use for significant GO terms
      # topNodes: how many GO terms to return in result table
      
      # Prepare geneList as expected for topGO (i.e. value vector with names of genes)
      if(!is.null(interestingGenes)){
        message(sprintf("Running GO enrichments with %s genes in universe of %s...", 
                        length(interestingGenes), length(allGenes)))
        geneList <- factor(as.integer(allGenes %in% interestingGenes))
        names(geneList) <- allGenes
        # Create topGOdata object
        GOdata <- suppressMessages(new(
          "topGOdata",
          ontology = ontology,
          allGenes = geneList,
          annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol",
          nodeSize = nodeSize
        ))
      }else{
        geneList <- pvals
        names(geneList) <- allGenes
        message(sprintf("Running GO enrichments with %s genes in universe of %s...", 
                        sum(geneSel(geneList)), length(allGenes)))
        GOdata <- suppressMessages(new(
          "topGOdata",
          ontology = ontology,
          allGenes = geneList,
          geneSel = geneSel,
          annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol",
          nodeSize = nodeSize
        ))
      }
      
      # Test for enrichment using Fisher's Exact Test
      GOresult <- suppressMessages(runTest(GOdata, algorithm=alg, statistic=stat))
      GenTable(GOdata, pvalue=GOresult, topNodes=topNodes, numChar=1000)
    }
    
    topGObarPlot <- function(goRes, cmap = NULL, nterms=10, border_color="black", 
                             barwidth=0.85, title="", enrichLimits=c(0.0, 5.5), barLimits=NULL){
      # Plot GO results in bar plot form
      goRes$log2FoldEnrichment <- log2(goRes$Significant / goRes$Expected)
      goRes$log2FoldEnrichment <- ifelse(goRes$log2FoldEnrichment > enrichLimits[2], enrichLimits[2], goRes$log2FoldEnrichment)
      goRes$threshPval <- ifelse(goRes$pvalue == "< 1e-30", 1e-30, as.numeric(goRes$pvalue))
      goRes$log10pval <- -log10(goRes$threshPval)
      if(!is.null(barLimits)){
        goRes$log10pval <- ifelse(goRes$log10pval < barLimits[2], goRes$log10pval, barLimits[2])
      }
      
      # Only plot the top nterms (reverse order to plot most significant at top)
      goRes <- goRes[1:nterms,]
      goRes <- goRes[nrow(goRes):1,]
      
      if(is.null(cmap)){
        cmap <- cmaps_BOR$rocket
      }
      p <- (
        ggplot(goRes, aes(x=Term, y=log10pval, fill=log2FoldEnrichment))
        + geom_bar(stat="identity", width=barwidth, color=border_color)
        + scale_x_discrete(
          limits=goRes$Term, # Required to prevent alphabetical sorting of terms
          labels= function(x) str_wrap(x, width=40) # wrap long GO term labels
        ) 
        + scale_fill_gradientn(colors=cmap, limits=enrichLimits)
        + xlab("")
        + ylab("-log10 pvalue")
        + ggtitle(title)
        + theme_BOR(border=FALSE)
        + theme(panel.grid.major=element_blank(), 
                panel.grid.minor= element_blank(), 
                plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
                #aspect.ratio = 6/nterms, # What is the best aspect ratio for a bar chart?
                axis.text.x = element_text(angle = 90, hjust = 1)) 
        + coord_flip()
      )
      if(!is.null(barLimits)){
        p <- p + scale_y_continuous(limits=barLimits, expand=c(0,0))
      }else{
        p <- p + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
      }
      p
    }
    
    Volcano_plots <- function(geneCompared,log2FC,FDR,pCutoff,name=name,GO=F){
      keyvals <- ifelse(
        geneCompared$avg_log2FC < -log2FC & geneCompared$p_val_adj < FDR, '#fcd686',# c("#d34a27", "#fcd686", "#8A9FD1")
        ifelse(geneCompared$avg_log2FC > log2FC & geneCompared$p_val_adj < FDR, '#d34a27',
               'black'))
      keyvals[is.na(keyvals)] <- 'black'
      names(keyvals)[keyvals == '#d34a27'] <- 'highly-EXP in First group'
      names(keyvals)[keyvals == 'black'] <- 'Nodiff'
      names(keyvals)[keyvals == '#fcd686'] <- 'highly-EXP in Second group'
      
      pdf(paste0(paste("Volcano_plot", name, sep = "_"), ".pdf"), height = 6, width = 6)   
      p <- EnhancedVolcano(geneCompared,
                           lab = rownames(geneCompared),
                           x = 'avg_log2FC',
                           y = 'p_val_adj',
                           selectLab = rownames(geneCompared)[which(names(keyvals) %in% c('highly-EXP in First group', 'highly-EXP in Second group'))] %>% head(10),
                           xlab = bquote(~Log[2]~ 'fold change'),
                           ylab = bquote(~-Log[10] ~ "FDR"),
                           title = name,
                           pCutoff = pCutoff,
                           FCcutoff = log2FC,
                           pointSize = 1,
                           labSize = 3,
                           #               shape = c(6, 4, 2, 11),
                           colCustom = keyvals,
                           colAlpha = 1,
                           #              legendPosition = 'left',
                           legendLabSize = 8,
                           legendIconSize = 5.0,
                           drawConnectors = TRUE,
                           widthConnectors = 0.5,
                           colConnectors = 'black',
                           arrowheads = F,
                           gridlines.major = FALSE,
                           gridlines.minor = FALSE,
                           border = 'full',
                           borderWidth = 1.0,
                           borderColour = 'black',
      )
      print(p)
      dev.off()
      
      if (GO){
        # Get GO enrichments for highly-EXP genes in First group
        GOresults <- calcTopGo(rownames(geneCompared), interestingGenes=rownames(geneCompared)[which(names(keyvals) %in% c('highly-EXP in First group'))])
        
        # Plots of GO term enrichments:
        pdf(paste0( "Frist_group_highly_EXP_GO_",name,".pdf"), width=8, height=2.5)
        print(topGObarPlot(GOresults, cmap = cmaps_BOR$comet, 
                           nterms=4, border_color="black", 
                           barwidth=0.85, title=paste0(name,"_1st"),barLimits=c(0, 5)))
        dev.off()
        
        # Get GO enrichments for highly-EXP genes in Second group
        GOresults <- calcTopGo(rownames(geneCompared), interestingGenes=rownames(geneCompared)[which(names(keyvals) %in% c('highly-EXP in Second group'))])
        
        # Plots of GO term enrichments:
        pdf(paste0( "Second_group_highly_EXP_GO_",name,".pdf"), width=8, height=2.5)
        print(topGObarPlot(GOresults, cmap = cmaps_BOR$comet, 
                           nterms=4, border_color="black", 
                           barwidth=0.85, title=paste0(name,"_2nd"),barLimits=c(0, 5)))
        dev.off()
      }
    }
    
    dotPlot <- function(df, xcol, ycol, color_col, size_col, xorder=NULL, yorder=NULL, cmap=NULL, 
                        color_label=NULL, size_label=NULL, aspectRatio=NULL, sizeLims=NULL, colorLims=NULL){
      # Plot rectangular dot plot where color and size map to some values in df
      # (Assumes xcol, ycol, color_col and size_col are named columns)
      
      # If neither x or y col order is provided, make something up
      # Sort df:
      if(is.null(xorder)){
        xorder <- unique(df[,xcol]) %>% sort()
      }
      if(is.null(yorder)){
        yorder <- unique(df[,ycol]) %>% sort()
      }
      if(is.null(aspectRatio)){
        aspectRatio <- length(yorder)/length(xorder) # What is the best aspect ratio for this chart?
      }
      df[,xcol] <- factor(df[,xcol], levels=xorder)
      df[,ycol] <- factor(df[,ycol], levels=yorder)
      df <- df[order(df[,xcol], df[,ycol]),]
      
      # Make plot:
      p <- (
        ggplot(df, aes(x=df[,xcol], y=df[,ycol], color=df[,color_col], size=ifelse(df[,size_col] > 0, df[,size_col], NA)))
        + geom_point()
        + xlab(xcol)
        + ylab(ycol)
        + theme_BOR(border=TRUE)
        + theme(panel.grid.major=element_blank(), 
                panel.grid.minor= element_blank(), 
                plot.margin = unit(c(0.25,0,0.25,1), "cm"), 
                aspect.ratio = aspectRatio,
                axis.text.x = element_text(angle = 90, hjust = 1)) 
        + guides(
          fill = guide_legend(title=""), 
          colour = guide_colourbar(title=color_label, override.aes = list(size=5)),
          size = guide_legend(title=size_label)
        )
      )
      if(!is.null(cmap)){
        if(!is.null(colorLims)){
          p <- p + scale_color_gradientn(colors=cmap, limits=colorLims, oob=scales::squish, name = "")
        }else{
          p <- p + scale_color_gradientn(colors=cmap, name = "")
        }
      }
      if(!is.null(sizeLims)){
        p <- p + scale_size_continuous(limits=sizeLims)
      }
      p
    }
}
