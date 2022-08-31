


PI3K <- as.data.frame(read.csv(file = "/Users/yeliz/Desktop/Circular/hsa04150.PI3K-Akt.signaling.pathway.Hs.txt", 
                                     stringsAsFactors = FALSE, header=FALSE, sep = "\t")) %>% 
  pull(V2) %>% 
  str_extract("[^;]*") %>% 
  unique() %>% 
  as.data.frame() %>% 
  dplyr::rename(HumanName='.') %>% 
  mutate(Akt=1)






mTOR <- as.data.frame(read.csv(file = "/Users/yeliz/Desktop/Circular/hsa04150.mTOR.signaling.pathway.Hs.txt", 
                               stringsAsFactors = FALSE, header=FALSE, sep = "\t")) %>% 
  pull(V2) %>% 
  str_extract("[^;]*") %>% 
  unique() %>% 
  as.data.frame() %>% 
  dplyr::rename(HumanName='.') %>% 
  mutate(mTOR=1)


RAS <- as.data.frame(read.csv(file = "/Users/yeliz/Desktop/Circular/hsa04014.Ras.signaling.pathway.Hs.txt", 
                               stringsAsFactors = FALSE, header=FALSE, sep = "\t")) %>% 
  pull(V2) %>% 
  str_extract("[^;]*") %>% 
  unique() %>% 
  as.data.frame() %>% 
  dplyr::rename(HumanName='.') %>% 
  mutate(RAS=1)


MAPK <- as.data.frame(read.csv(file = "/Users/yeliz/Desktop/Circular/hsa04010.MAPK.signaling.pathway.Hs.txt", 
                               stringsAsFactors = FALSE, header=FALSE, sep = "\t")) %>% 
  pull(V2) %>% 
  str_extract("[^;]*") %>% 
  unique() %>% 
  as.data.frame() %>% 
  dplyr::rename(HumanName='.') %>% 
  mutate(MAPK=1)


Glioma$HumanName %>% union(Cancer$HumanName) %>% 
  union(Apoptosis$HumanName) %>% 
  union(cycle$HumanName) %>% 
  union(PI3K$HumanName) %>% 
  union(VEGF$HumanName) %>% 
  union(mTOR$HumanName) %>% 
  union(RAS$HumanName) %>% 
  union(MAPK$HumanName) -> union.pathways


hm1 %>% filter (HumanName %in% union.pathways) -> x


x %>% left_join(Glioma, by="HumanName") %>% 
  left_join(Cancer, by="HumanName") %>%
  left_join(Apoptosis, by="HumanName") %>% 
  left_join(cycle, by="HumanName") %>% 
  left_join(PI3K, by="HumanName") %>% 
  left_join(VEGF, by="HumanName") %>% 
  left_join(mTOR, by="HumanName") %>% 
  left_join(RAS, by="HumanName") %>% 
  left_join(MAPK, by="HumanName") -> binary







  
binary[is.na(binary)] <- 0

lfcmax <- max(binary$Log2FC)
lfcmin <- min(binary$Log2FC)


binary1 <- binary

binary1 <- binary1[c(1,3,4,5,6,7,8,9,10,2)]


binary1 %>% filter(!(duplicated(HumanName)| duplicated(HumanName, fromLast=TRUE))) -> binary2

rownames(binary2) <- binary2$HumanName

binary2 <- binary2[,-1]

binary3 <- binary2


binary3 %>% filter(Glioma+Cancer+Apoptosis+cycle+Akt+VEGF+mTOR+RAS!=0) -> binary3

lfcmax <- max(binary3)
lfcmin <- min(binary3)


paletteLength <- 50

myColor <- c(colorRampPalette(c("paleturquoise4","slategray2"))(floor((paletteLength-1)/2)), 
             colorRampPalette("cornsilk")(1), 
             colorRampPalette(c("moccasin","darkorange3"))(floor((paletteLength-1)/2)))
MyBreaks <- c(seq(lfcmin, -.58, length.out = floor(paletteLength/2)), 
              seq(.58, lfcmax, length.out = ceiling(paletteLength/2)))

GOChord(binary3 %>% dplyr::rename(logFC=Log2FC), space = 0.010, gene.order = c('Log2FC') ,gene.space = 0.25, gene.size = 1, process.label=3, border.size = 0.001, 
        ribbon.col = c("#e46b69","#4984da","#66cec9","#c5d15b","#fed963","#f88a58","black","yellow"), lfc.max=lfcmax, lfc.min =lfcmin,
        lfc.col=c("red","white","blue")) + 
  scale_fill_gradientn(limits=c(lfcmin,lfcmax),colours=myColor, values=MyBreaks %>% scales::rescale(from=c(lfcmin,lfcmax)))






















binary_file <- read.csv("/Users/yeliz/Desktop/circular_plot/binary_file.csv")
X1_path <- read.csv("/Users/yeliz/Desktop/circular_plot/X1_Path.CSV")
X1_path <- X1_path[,1:2]
colnames(X1_path) <- c("query.term", "logFC") 

X1_binary <- merge(binary_file, X1_path, by="query.term")
X1_binary <- distinct(X1_binary)

rownames(X1_binary) <- X1_binary$query.term
chord_X1 <- X1_binary[,2:8]      



binary_file <- read.csv("/Users/yeliz/Desktop/circular_plot/binary_file.csv")
X3_path <- read.csv("/Users/yeliz/Desktop/circular_plot/X3_Path.CSV")
X3_path <- X3_path[,1:2]
colnames(X3_path) <- c("query.term", "logFC") 

X3_binary <- merge(binary_file, X3_path, by="query.term")
X3_binary <- distinct(X3_binary)

rownames(X3_binary) <- X3_binary$query.term
chord_X3 <- X3_binary[,2:8]  


lfcmax <- max(chord_X1$logFC,chord_X3$logFC)
lfcmin <- min(chord_X1$logFC,chord_X3$logFC)

paletteLength <- 50

myColor <- c(colorRampPalette(c("paleturquoise4","slategray2"))(floor((paletteLength-1)/2)), colorRampPalette("cornsilk")(1), colorRampPalette(c("moccasin","darkorange3"))(floor((paletteLength-1)/2)))
MyBreaks <- c(seq(lfcmin, -.58, length.out = floor(paletteLength/2)), seq(.58, lfcmax, length.out = ceiling(paletteLength/2)))



GOChord(chord_X1, space = 0.010, gene.order = c('logFC') ,gene.space = 0.25, gene.size = 2, process.label=7, border.size = 0.001, 
        ribbon.col = c("#e46b69","#4984da","#66cec9","#c5d15b","#fed963","#f88a58"), lfc.max=8.6, lfc.min = -2.9,
        lfc.col=c("red","white","blue")) + 
  scale_fill_gradientn(limits=c(lfcmin,lfcmax),colours=myColor, values=MyBreaks %>% scales::rescale(from=c(lfcmin,lfcmax)))

GOChord(chord_X3, space = 0.015, gene.order = 'logFC',gene.space = 0.25, gene.size = 6.5, process.label=7, border.size = 0.001, 
        ribbon.col =c("#e46b69","#4984da","#66cec9","#c5d15b","#fed963","#f88a58"), lfc.max=8.6, lfc.min = -2.9,
        lfc.col=c('red','white','blue')) + 
  scale_fill_gradientn(limits=c(lfcmin,lfcmax),colours=myColor, values=MyBreaks %>% scales::rescale(from=c(lfcmin,lfcmax)))



c("#E57373","#FFCC80","paleturquoise3","#9FA8DA","#81C784","lightpink")


