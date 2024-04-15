#This block calculates death time of single-cells in microscopic images.
for (j in 1:1) {
  #Load functions; change directories to desktop.
  source("/Users/apple/Desktop/Persistence & Spatial Paper/Single-Cell Microscopy/Nikon_death_functions.R")
  setwd("Desktop/")
  total_cluster_loc=read.delim("ES_suppl4_S_areas_total.txt")
  
  #To identify cell clusters in the first image.
  a=total_cluster_loc[total_cluster_loc$Label==1,]
  plot(density(a$Area))
  single_cluster_at_0h=total_cluster_loc[total_cluster_loc$Label==1 & total_cluster_loc$Area>25,] 
  write.table(single_cluster_at_0h,"single_cluster_at_0h_ES_Supple4_S_XX02202024.txt",sep='\t')
  
  #To measure death (disappearance) time of single cells.
  max_image=22
  single_cluster_area_over_time=cell_area_over_time(total_cluster_loc,max_image)
  single_cell_death_time=death_disappearance_time(single_cluster_area_over_time,max_image)  
  write.table(single_cell_death_time,"ES_Supple4_S_death_time_area_25cutoff_XX02202024.txt",sep='\t')
}    

#This block correlates E. coli death time with So cell biomass within given neighborhoods.
for (j in 1:1) {
  source("/Users/apple/Desktop/Single-Cell Microscopy/Nikon_death_functions.R")
  source("/Users/apple/Desktop/Single-Cell Microscopy/Nikon_death_functions.R")
  
  setwd(("/Users/apple/Desktop/"))
  Co_E_death_total=read.delim("ES_Supple4_E_death_time_area_25cutoff_curated_XX02202024.txt")
  Co_S_death_total=read.delim("ES_Supple4_S_death_time_area_25cutoff_curated_XX02202024.txt")
  
  sprm_rho_FracNei=c(); sprm_pval_FracNei=c()
  rad_size=c(50,100)
  rad_size=c(rad_size,seq(200,3600,200))
  for (i in rad_size) { #i denotes the radius of the neiborhood size
    #spearman_cor_by_ngbr=Corr_rho_ngbr_frac_E_death_alive_S(Co_S_death_total,Co_E_death_total,i)
    spearman_cor_by_ngbr=Corr_rho_ngbr_frac_E_death(Co_E_death_total,Co_S_death_total,i)
    #spearman_death_time_fracNei=cor.test(spearman_cor_by_ngbr$S_Fraction, spearman_cor_by_ngbr$Death_Time, method="spearman",exact=FALSE)
    spearman_death_time_fracNei=cor.test(spearman_cor_by_ngbr$S_Fraction, spearman_cor_by_ngbr$Death_Time, method="spearman",exact=FALSE)
    sprm_rho_FracNei=c(sprm_rho_FracNei,spearman_death_time_fracNei$estimate)
    sprm_pval_FracNei=c(sprm_pval_FracNei,spearman_death_time_fracNei$p.value)
  }
  lag_corr_S_frac=data.frame(Neighborhood_Size=rad_size,Spearman_rho_S_fraction=sprm_rho_FracNei,Spearman_p_val=sprm_pval_FracNei)
  write.table(lag_corr_S_frac, "predicting_S_death_w_neighborhood_info_Co2.txt", sep="\t")
}

ggplot(spearman_cor_by_ngbr,aes(x=S_Fraction,y=Death_Time))+geom_point()

a=read.delim("Co1_AMP_death_area_total.txt")
ggplot(a,aes(x=X,y=-Y,col=Species,label=Death_Time))+geom_point(shape=1)+geom_text(hjust=0, vjust=0,size=3)

lambda=2.17145
curve(lambda*exp(-lambda*x), from=1, to=10, n=300, xlab="x", ylab="PDF", 
        col="blue", lwd=2, main="Exponential Distribution"  )


#This block of code predicts single-cell death time for a focal species by calculating
#the relative distance to closest nearby cells in the partner species (DLAN mechanism).
for (variable in vector) {
  source("/Users/apple/Desktop/Persistence & Spatial Paper/Single-Cell Microscopy/Nikon_death_functions.R")
  setwd(("/Users/apple/Desktop/"))
  
  Co_E_death_total=read.delim("ES_Suppl5_E_death_time_curated_122023.txt")
  Co_S_death_total=read.delim("ES_Suppl5_S_death_time_curated_122023.txt")
  relative_distance_correlate_E_death_w_alive_closest_S=death_time_to_alive_closest_So_with_standardization(Co_E_death_total,Co_S_death_total)
  relative_distance_correlate_E_death_w_alive_closest_S=death_time_to_alive_closest_partner_without_standardization(Co_E_death_total,Co_S_death_total)
  
  cor_closest_alive_So=cor.test(relative_distance_correlate_E_death_w_alive_closest_S$Eo_Death_Time,relative_distance_correlate_E_death_w_alive_closest_S$Closest_dist,method="spearman",exact=FALSE)
  ggplot(relative_distance_correlate_E_death_w_alive_closest_S,aes(x=Closest_dist,y=Eo_Death_Time))+geom_point()
  cor_closest_alive_So$p.value
  cor_closest_alive_So$estimate
  write.table(relative_distance_correlate_E_death_w_alive_closest_S,"relative_distance_correlate_E_death_w_alive_closest_S_Mutual3.txt",sep='\t')
  
  relative_distance_correlate_S_death_w_alive_closest_E=death_time_to_alive_closest_So_with_standardization(Co_S_death_total,Co_E_death_total)
  relative_distance_correlate_S_death_w_alive_closest_E=death_time_to_alive_closest_partner_without_standardization(Co_S_death_total,Co_E_death_total)
  So_cor_by_closest_alive_Eo=cor.test(relative_distance_correlate_S_death_w_alive_closest_E$Eo_Death_Time,relative_distance_correlate_S_death_w_alive_closest_E$Closest_dist,method="spearman",exact=FALSE)
  
}

#This block correlates E. coli single-cell death time with initial proximity to nearest So.
for (j in 1:1) {
  source("/Users/apple/Desktop/Single-Cell Microscopy/Nikon_death_functions.R")
  setwd(("/Users/apple/Desktop/"))
  Co_E_death_total=read.delim("Co3_E_AMP_death_time_area_25cutoff_curated.txt")
  Co_S_death_total=read.delim("Co3_S_AMP_death_time_area_25cutoff_curated.txt")
  
  larger_than_3=Co_E_death_total[Co_E_death_total$Death_Time>=3,]
  smaller_than_3=Co_E_death_total[Co_E_death_total$Death_Time<3,]
  l3_dist_to_S=death_time_to_closest_So(larger_than_3,Co_S_death_total)
  mean(l3_dist_to_S$Closest_dist)
  sd(l3_dist_to_S$Closest_dist)
  s3_dist_to_S=death_time_to_closest_So(smaller_than_3,Co_S_death_total)
  mean(s3_dist_to_S$Closest_dist)
  sd(s3_dist_to_S$Closest_dist)
  l3_S_death=death_time_to_closest_So_death_time(larger_than_3,Co_S_death_total)
  mean(l3_S_death$Closest_So_Death_Time)
  sd(l3_S_death$Closest_So_Death_Time)
  s3_S_death=death_time_to_closest_So_death_time(smaller_than_3,Co_S_death_total)
  mean(s3_S_death$Closest_So_Death_Time)
  sd(s3_S_death$Closest_So_Death_Time)
  
  correlate_E_death_w_alive_closest_S=death_time_to_alive_closest_So(Co_S_death_total,Co_E_death_total)
  cor_closest_alive_So=cor.test(correlate_E_death_w_alive_closest_S$Eo_Death_Time,correlate_E_death_w_alive_closest_S$Closest_dist,method="spearman",exact=FALSE)
  ggplot(correlate_E_death_w_alive_closest_S,aes(x=Closest_dist,y=Eo_Death_Time))+geom_point()
  cor_closest_alive_So=cor.test(a$Eo_Death_Time,a$Closest_dist,method="spearman",exact=FALSE)
  
  correlate_E_death_w_closest_S=death_time_to_closest_So(Co_E_death_total,Co_S_death_total)
  cor_closest_So=cor.test(correlate_E_death_w_closest_S$Eo_Death_Time,correlate_E_death_w_closest_S$Closest_dist,method="spearman",exact=FALSE)
  write.table(correlate_E_death_w_closest_S,"correlate_E_death_w_closest_S_Co3.txt",sep='\t')
  ggplot(correlate_E_death_w_closest_S,aes(x=Closest_dist,y=Eo_Death_Time))+geom_point()
  
  correlate_E_death_w_closest_S_death=death_time_to_closest_So_death_time(Co_E_death_total,Co_S_death_total)
  ggplot(correlate_E_death_w_closest_S_death,aes(x=Closest_So_Death_Time,y=Eo_Death_Time))+geom_point()
  cor_closest_So_death_time=cor.test(correlate_E_death_w_closest_S_death$Eo_Death_Time,correlate_E_death_w_closest_S_death$Closest_So_Death_Time,method="spearman",exact=FALSE)
}
