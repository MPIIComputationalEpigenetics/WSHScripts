library(ggplot2)
ct_hetero_file <- read.csv("output_of_density_coverage_statistic_CT_HETERO")
contamination_file <- read.csv("output_of_density_coverage_statistic_SAMPLE_PURITY")
asm_file <- read.csv("output_of_density_coverage_statistic_ASM")
erosion_file <- read.csv("output_of_density_coverage_statistic_EROSION")
msd_file <- read.csv("output_of_density_coverage_statistic_MSD")
plot.path <- getwd()

to.plot <- rbind(cbind(ct_hetero_file,Scenario=rep("Cell type heterogeneity",nrow(ct_hetero_file))),cbind(contamination_file,Scenario=rep("Cellular contamination",nrow(contamination_file))),cbind(asm_file,Scenario=rep("ASM",nrow(asm_file))),cbind(erosion_file,Scenario=rep("Methylation erosion",nrow(erosion_file))),cbind(msd_file,Scenario=rep("Methylation Switching Domains",nrow(msd_file))))
colnames(to.plot)[3] <- "Scenario"
mean(to.plot$AvgCovg[to.plot$Scenario %in% "Cell type heterogeneity"]) 
mean(to.plot$AvgCovg[to.plot$Scenario %in% "Cellular contamination"])  
mean(to.plot$AvgCovg[to.plot$Scenario %in% "ASM"]) 
mean(to.plot$AvgCovg[to.plot$Scenario %in% "Methylation erosion"]) 
mean(to.plot$AvgCovg[to.plot$Scenario %in% "Methylation Switching Domains"]) 
mean(to.plot$NSites[to.plot$Scenario %in% "Cell type heterogeneity"])/50 
mean(to.plot$NSites[to.plot$Scenario %in% "Cellular contamination"])/50  
mean(to.plot$NSites[to.plot$Scenario %in% "ASM"])/50 
mean(to.plot$NSites[to.plot$Scenario %in% "Methylation erosion"])/50 
mean(to.plot$NSites[to.plot$Scenario %in% "Methylation Switching Domains"])/50

plot <- ggplot(to.plot,aes(x=Scenario,y=NSites/50))+geom_boxplot()+ylab("Number of CpGs per kb")+theme_bw()+theme(panel.grid=element_blank())
ggsave(file.path(plot.path,"CpG_density_100bp.pdf"))
plot <- ggplot(to.plot,aes(x=Scenario,y=AvgCovg))+geom_boxplot()+ylab("Average coverage per dataset")+theme_bw()+theme(panel.grid=element_blank())
ggsave(file.path(plot.path,"average_coverage_100bp.pdf"))

