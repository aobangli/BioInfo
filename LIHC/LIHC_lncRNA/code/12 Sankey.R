rm(list = ls())
load('correlation_coefficient_matrix.Rdata')
load(file = './LIHC_lipid_limma_DEG.Rdata')

arachidonic_gene = c("PLA2G10","PLA2G2D","PLA2G2E","PLA2G3",
                     "PLA2G2F","PLA2G12A","PLA2G12B","PLA2G1B",
                     "PLA2G5","PLA2G2A","PLA2G2C","PLA2G4E","PLA2G4A",
                     "JMJD7-PLA2G4B","PLA2G4B","PLA2G4C","PLA2G4D",
                     "PLA2G4F","PLA2G6","PLB1","PTGS1","PLAAT3",
                     "PTGS2","PTGES","PTGES2","PTGES3","CBR1","CBR3",
                     "PRXL2B","TBXAS1","PTGDS","HPGDS","AKR1C3",
                     "PTGIS","ALOX5","LTA4H","CYP4F2","CYP4F3","LTC4S",
                     "GGT1","GGT5","GPX6","GPX7","GPX2","GPX3","GPX1",
                     "GPX5","GPX8","CYP2E1","CYP2J2","CYP2U1","CYP4A11",
                     "CYP2C19","CYP4F8","ALOX12","ALOX12B","ALOX15B","CYP2B6",
                     "CYP2C8","CYP2C9","EPHX2","ALOX15")

lncRNAs = row.names(need_DEG[need_DEG$change != 'NOT', ])

########## 计算筛选出来的lncRNA和花生四烯酸基因的相关性桑葚图 ##########
library(ggalluvial)
library(reshape2)
matrix = r_matrix[lncRNAs, ]
data = melt(matrix)
colnames(data) = c('lncRNA', "arachidonic", "cor")
data = data[data$cor > 0.5, ]

san = to_lodes_form(data[,c(1,2)], axes = 2:1, id = "Group")

mycol <-
  rep(c('#223D6C', '#D20A13', '#FFD121', '#088247', '#11AA4D', '#58CDD9', '#7A142C', '#5D90BA',
        '#029149', '#431A3D', '#91612D', '#6E568C', '#E0367A', '#D8D155', '#64495D', '#7CC767',
        '#223D6C', '#D20A13', '#FFD121', '#088247', '#11AA4D', '#58CDD9', '#7A142C', '#5D90BA',
        '#029149', '#431A3D', '#91612D', '#6E568C', '#E0367A', '#D8D155', '#64495D', '#7CC767',
        '#223D6C', '#D20A13', '#FFD121', '#088247', '#11AA4D', '#58CDD9', '#7A142C', '#5D90BA',
        '#029149', '#431A3D', '#91612D', '#223D6C', '#D20A13', '#FFD121', '#088247', '#11AA4D',
        '#58CDD9', '#7A142C', '#5D90BA', '#029149', '#431A3D', '#91612D', '#6E568C', '#E0367A',
        '#D8D155', '#64495D', '#7CC767', '#223D6C', '#D20A13', '#FFD121', '#088247', '#11AA4D',
        '#58CDD9', '#7A142C', '#5D90BA', '#029149', '#431A3D', '#91612D', '#6E568C', '#E0367A',
        '#D8D155', '#64495D', '#7CC767'), 3)
        

pdf(file="Sankey.pdf",width = 6,height = 13)
ggplot(san,   #定义流向后的数据集 
       aes(x = x, stratum = stratum, alluvium = Group, fill = stratum, label = stratum)) +  
  geom_flow(width = 1/3,aes.flow = "forward") +  
  #geom_flow()函数控制边的视觉通道映射设定，也就是线条的颜色；主要由alluvium和weight决定 
  #forward表示线条的颜色与前面的变量一致，而backward则表示与后者一致 
  geom_stratum(alpha = 0.4, width = 1/3, color = 'grey70') + 
  #geom_stratum()函数控制节点的视觉通道映射设定，主要由stratum和weight决定 
  geom_text(stat = "stratum", size = 3.5, color="black") + 
  #size = 3.5代表名字的大小 
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = mycol) +
  xlab("") + 
  ylab("") + 
  theme_bw() +  
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank()) + #去掉坐标轴 
  theme(panel.grid =element_blank()) +  
  theme(panel.border = element_blank()) +  
  theme(axis.line.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 20, color="black")) + #显示分组名字 
  ggtitle("") + 
  guides(fill = FALSE)
dev.off()