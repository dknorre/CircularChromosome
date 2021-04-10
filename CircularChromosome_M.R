
library(tidyverse)
library(Rsamtools)
library(RColorBrewer)
library(Biostrings)
library(IRanges)
library(rtracklayer)
library(circlize)

  #Set working directory
  setwd("/home/dima/Desktop/BiocR")
  
  deletionCoord <- read.csv("rho_coordinates.csv", header = TRUE, sep = ";")  
  
  
   ############ #Analyze complete genome GCcontent v0.1
   ########### to fix: handling genome tail, overlapping window, ...
   
  
   yeastGenome <- readDNAStringSet("/home/dima/Desktop/NGS/References/S288C_with2micron.fsa")
   # window = 800 #set the window for GC content
   # step = 400 # 
   GenomeSize = length(yeastGenome[[18]])
   # steps <- GenomeSize %/% step
   # 
   # GCcontent =  data.frame(pos = as.numeric(), GC = as.numeric(), ATskew = as.numeric())
   # for (i in seq(1, steps))
   # {
   # position <- i*(step)
   # chunk <- yeastGenome[[18]][(position-(window/2)):(position+(window/2))]
   # GC <- letterFrequency(chunk, "GC", as.prob = TRUE)
   # Aa <- letterFrequency(chunk, "A")
   # Tt <- letterFrequency(chunk, "T")
   # ATskew <- (Aa-Tt)/(Aa+Tt) # AT skew = (A − T)/(A + T)
   # GCcontent <- rbind(GCcontent, data.frame(pos = position, GC = GC, ATskew = ATskew ))
   # }
   # 
   # GCcontent <-GCcontent %>% mutate(cumulativeATskew = cumsum(GCcontent$ATskew))
   # 
   #####import MtDNA genomic features
      # rtracklayer
   mtDNAfeatures <- readGFF("/home/dima/Desktop/NGS/References/mtdna.gff3")
   as.data.frame(mtDNAfeatures) %>% filter(!is.na(type)) %>% group_by(type) %>% summarise(n = n())
   mtDNAlength <- as.data.frame(mtDNAfeatures) %>% filter(type == "chromosome") %>% pull(end)

   # tidy mtDNAfeatures (Yeast-specific)
   mtDNAfeatures <- as.data.frame(mtDNAfeatures) %>% filter(!is.na(type)) %>% 
      select(type, start, end, strand, Name, X_Gene) %>% 
      mutate(coordinates = (start+end)/2) %>% filter(type != "chromosome") %>%
      mutate(X_Gene = ifelse(!is.na(X_Gene), X_Gene, Name))
   
   
   # MtGenes <- mtDNAfeatures %>% filter(!is.na(X_Gene)) %>% filter(type == "gene")
   # MtRNAa <- mtDNAfeatures  %>% filter(type %in% c("rRNA_gene", "ncRNA_gene", "tRNA_gene"))
   # MtORIs <- mtDNAfeatures %>% filter(type == "origin_of_replication")
   # 
   # length(levels(mtDNAfeatures$type))
   FeatureColor <- brewer.pal(length(levels(mtDNAfeatures$type)),"Set2") %>%
      set_names(levels(mtDNAfeatures$type))

   
 ################ Circular maps

   ### BioCircus
 
  circos.par("gap.degree" = 0, 
             "cell.padding" = c(0, 0, 0, 0), 
             "start.degree" = 90,
             "track.height" = 0.1)
  circos.initialize("a", xlim = c(0, mtDNAlength))
  circos.track(ylim = c(0,1),bg.border = "white")
    # circos.text(mtDNAfeatures$coordinates,0.7,mtDNAfeatures$X_Gene, niceFacing = TRUE, facing = "clockwise", cex = 0.5)

  BEDmtDNAfeatures <- mtDNAfeatures %>% mutate(chr = "a") %>% select(chr, start, end, X_Gene)
  circos.genomicLabels(BEDmtDNAfeatures, 
                       labels.column = 4, 
                       side = "outside", 
                       col = FeatureColor[mtDNAfeatures$type],
                       cex = 0.5)
  
  # circos.info()
  
  
  circos.track(ylim = c(0,1),bg.border = "white")
  for (row in 1:nrow(mtDNAfeatures)) {
                                circos.arrow(mtDNAfeatures$start[row],
                                             mtDNAfeatures$end[row], 
                                             y= 0.7,
                                             width = 0.2,
                                             arrow.head.width = 0.3,
                                             arrow.head.length =  100,
                                             arrow.position = ifelse(mtDNAfeatures$strand[row] == "+", "end", "start"),
                                             col = FeatureColor[mtDNAfeatures$type[row]])
  }


 # circos.track(ylim = c(0,0.4))
 # circos.lines(x = GCcontent$pos, y = GCcontent$GC, area = TRUE, type = 'h', col = 'darkblue', lwd = 2)
 circos.xaxis(labels.cex = 0.3)
 # circos.track(ylim = c(min(GCcontent$ATskew),max(GCcontent$ATskew)))
 # circos.lines(x = GCcontent$pos, y = GCcontent$ATskew, type = 'l', col = 'orange', lwd = 2)
  ## to plot Coverage from real data  
  # circos.track(ylim = c(0,4000))
  # circos.lines(x = Coverage$pos, y = Coverage$count + mCoverage$count, area = TRUE, type = 'h',col = 'orange')
   # circos.track(ylim = c(0,0.4))
  # circos.lines(x = Coverage$pos, y = Coverage$count, area = TRUE, type = 'h', col = 'darkblue')

  ##DELETION ARCS
  
  
  #auxilary calculations
  deletionCoord<-  deletionCoord %>% mutate(length = end-start) %>% 
     mutate(length = ifelse(length <0, length + mtDNAlength, length)) %>% arrange(length) 
  col_fun = colorRamp2(c(0, 50, 100), c("green", "black", "red"))  
  
  #plotting arcs
  circos.track(ylim = c(-10,0),bg.border = "white")
   for (row in seq(1:nrow(deletionCoord))) {
      if (deletionCoord$end[row]>deletionCoord$start[row]) {
         circos.segments(deletionCoord$start[row], 10-row, deletionCoord$end[row], 10-row,
                         col = col_fun(deletionCoord$suppressivity[row]))
         }
      
      else {
         circos.segments(deletionCoord$start[row], 10-row, mtDNAlength, 10-row,
                         col = col_fun(deletionCoord$suppressivity[row]))
         circos.segments(1, 10-row, deletionCoord$end[row], 10-row,  
                         col = col_fun(deletionCoord$suppressivity[row]))
      }
   }




  
