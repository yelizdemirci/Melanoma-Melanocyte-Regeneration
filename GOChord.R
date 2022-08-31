
library(ggdendro)
library(GOplot)
library(tidyverse)

# Read the gene annotation file 
gtf.Dr <- "Desktop/Submission_MM/Genome-GTF-files/Danio_rerio.GRCz11.93.gtf"
ann.Dr <- elementMetadata(import(gtf.Dr))

# Read the fold changes of time points from DESeq results
load("~/Desktop/Submission_MM/RData/RES_1dpa.RData")
load("~/Desktop/Submission_MM/RData/RES_7dpa.RData")
load("~/Desktop/Submission_MM/RData/RES_nevi.RData")
load("~/Desktop/Submission_MM/RData/RES_melanoma.RData")

res_1dpa <- res_1dpa[,c(2,7)] %>% as.data.frame() %>% rownames_to_column('IDs')
colnames(res_1dpa) <- c('IDs','log2FC-1dpa', 'gene_name')

res_7dpa <- res_7dpa[,c(2,7)] %>% as.data.frame() %>% rownames_to_column('IDs')
colnames(res_7dpa) <- c('IDs','log2FC-7dpa', 'gene_name')

res_nevi <- res_nevi[,c(2,7)] %>% as.data.frame() %>% rownames_to_column('IDs')
colnames(res_nevi) <- c('IDs','log2FC-nevi', 'gene_name')

res_melanoma <- res_melanoma[,c(2,7)] %>% as.data.frame() %>% rownames_to_column('IDs')
colnames(res_melanoma) <- c('IDs','log2FC-melanoma', 'gene_name')

# Generate a fold change table (all)
all <- full_join(res_1dpa, res_7dpa, by='IDs' ) %>% full_join(res_nevi, by='IDs') %>% full_join(res_melanoma, by='IDs')
all <- all[,c(1,2,4,6,8,9)]

# Read the genes in pathways
NC <- read.csv('~/Desktop/Submission_MM/AmiGO-GOs-for-comparisons/3-GO-0001755-neural crest cell migration.txt', 
               stringsAsFactors = FALSE, header=FALSE, sep = "\t") %>% 
  pull(V1) %>% 
  str_extract("[^;]*") %>% 
  unique() %>% 
  as.data.frame() %>% 
  dplyr::rename(gene_name='.') %>% 
  mutate(NC=1)

GSD <- read.csv('~/Desktop/Submission_MM/AmiGO-GOs-for-comparisons/10-GO-0010001-glial cell differentiation.txt', 
                stringsAsFactors = FALSE, header=FALSE, sep = "\t") %>% 
  pull(V1) %>% 
  str_extract("[^;]*") %>% 
  unique() %>% 
  as.data.frame() %>% 
  dplyr::rename(gene_name='.') %>% 
  mutate(GSD=1)

FR <- read.csv('~/Desktop/Submission_MM/AmiGO-GOs-for-comparisons/13-GO-0031101-fin regeneration.txt', 
               stringsAsFactors = FALSE, header=FALSE, sep = "\t") %>% 
  pull(V1) %>% 
  str_extract("[^;]*") %>% 
  unique() %>% 
  as.data.frame() %>% 
  dplyr::rename(gene_name='.') %>% 
  mutate(FR=1)

FM <- read.csv('~/Desktop/Submission_MM/AmiGO-GOs-for-comparisons/14-GO-0033334-fin-morphogenesis.txt', 
               stringsAsFactors = FALSE, header=FALSE, sep = "\t") %>% 
  pull(V1) %>% 
  str_extract("[^;]*") %>% 
  unique() %>% 
  as.data.frame() %>% 
  dplyr::rename(gene_name='.') %>% 
  mutate(FM=1)


NC$gene_name %>% union(GSD$gene_name) %>% 
  union(FR$gene_name) %>% 
  union(FM$gene_name) -> union.pathways


# Take the union genes in all lists and remove duplicated ones
union.pathways <- union.pathways[! union.pathways %in% c('adgrg1', 'cdh23', 'cdh30', 'celsr2', 'col6a3', 'dlg4b', 'evx1', 
                                                         'hmcn2', 'itgb1b.1', 'pcdh11', 'pcdh1gb2', 'pcdh20', 'pcdh2ab2', 
                                                         'plxnb2a', 'ptprsa', 'robo2', 'selp', 'si:ch211-182p11.1', 'si:ch211-215c18.3', 
                                                         'si:dkey-238d)18.7', 'si:dkey-238d18.7')]

union.pathways <- union.pathways %>% setdiff(str_subset(., "^[(si:)|(zmp:)]"))


# Subset pathway genes from fold change
all %>% filter (gene_name.y.y %in% union.pathways) -> x


# Change the column names in the table
colnames(x) <- c("ID","log2FC-1dpa", "log2FC-7dpa", "log2FC-nevi", "log2FC-melanoma", "gene_name")


# Generate binary table for Chord plot
x %>% left_join(NC, by="gene_name") %>%  
  left_join(GSD, by="gene_name") %>%
  left_join(FR, by="gene_name") %>% 
  left_join(FM, by="gene_name") -> binary

binary[is.na(binary)] <- 0

rownames(binary) <- binary$gene_name

binary <- binary[ -c(1,6) ]

binary1 <- binary
binary1 <- binary1[c(5,6,7,8,1,2,3,4)]

binary2 <- binary1
binary2 %>% filter(NC+GSD+FR+FM!=0) -> binary3 # remove zeros

# Set max and min for fold changes
lfcmax <- max(binary3)
lfcmin <- min(binary3)


# Set s colour palette
paletteLength <- 50

myColor <- c(colorRampPalette(c("darkblue","skyblue"))(floor((paletteLength-1)/2)), colorRampPalette("cornsilk")(1), colorRampPalette(c("lightpink","darkred"))(floor((paletteLength-1)/2)))
MyBreaks <- c(seq(min(lfcmin), -log2(1.2), length.out = floor(paletteLength/2)), seq(log2(1.2), max(lfcmax), length.out = ceiling(paletteLength/2)))

# Order and cluster binary file based on fold change
binary3 %>% dplyr::select(starts_with("log")) %>% 
  as.matrix() %>% 
  {rownames(.)[hclust(dist(.))$order]} -> order

binary3 <- binary3[order,]


# Change the column names of the binary file
colnames(binary3) <- c("Neural crest cell migration", 
                       "Glial cell differentiation",
                       "Fin regeneration", 
                       "Fin-morphogenesis",  
                       "log2FC-1dpa",   
                       "log2FC-7dpa",   
                       "log2FC-nevi", 
                       "log2FC-melanoma")


# Plot the data
GOChord(binary3 %>% dplyr::rename(logFC=`log2FC-1dpa`),
        space = 0.000000001, 
        gene.order = 'none',
        gene.space = 0.3, 
        gene.size = 2, 
        nlfc=4,
        process.label=10, 
        limit = c(0,0),
        border.size = 0.0000000001, 
        ribbon.col = c((alpha('#d9c989',0.9)), (alpha('#c8acc1',0.9)), (alpha("#5493a6",0.6)), (alpha("#6f5554",0.7))),     
        lfc.max=lfcmax, lfc.min =lfcmin,
        lfc.col=c("red","white","blue")) + 
  scale_fill_gradientn (limits=c(lfcmin,lfcmax),colours=myColor, values=MyBreaks %>% scales::rescale(from=c(lfcmin,lfcmax))) 

















