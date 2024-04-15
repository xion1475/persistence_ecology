library("ggplot2")
library(dplyr)

#This block measures the lag time of all E. coli cell clusters, which it time it takes to gain 10% more of biomass.
for (j in 1:1) {
  #Load functions; change directories to desktop.
  source("/Users/apple/Desktop/Single-Cell Microscopy/Nikon_cell_cluster_growth_functions.R")
  setwd(("/Users/apple/Desktop/"))
  total_cluster_loc=read.delim("Co_C5_S_areas_total_frames.txt")
  
  #To identify cell clusters in the first image.
  single_cluster_at_0h=total_cluster_loc[total_cluster_loc$Label==1,] 
  write.table(single_cluster_at_0h,"single_S_cluster_at_0h_CoA1.txt",sep='\t')
  
  #To measure single-cell cluster area changes within a few images of time.
  max_image=10 #Set a total number of images we look at before clusters begin to fuse with each other.
  total_cluster_loc=total_cluster_loc[total_cluster_loc$Label<=max_image,]
  single_cluster_area_over_time=cluster_area_over_time(total_cluster_loc,max_image)
  
  #To remove cell clusters with NA as results.
  
  single_cluster_area_over_time_transposed=data.frame(t(single_cluster_area_over_time))
  single_cluster_area_over_time_transposed$X1[single_cluster_area_over_time_transposed$X1==1]=NA
  one_removed=t(single_cluster_area_over_time_transposed)
  write.table(one_removed,"Co1_single_cluster_area_w_na.txt",sep='\t')
  single_cluster_area_na_omitted=t(na.omit(single_cluster_area_over_time_transposed))
  write.table(single_cluster_area_na_omitted,"Co_C5_S_single_cluster_area_na_omitted.txt",sep='\t')
  
  single_cluster_area_na_omitted=read.delim("Co_C4_single_cluster_area_na_omitted.txt")
  #To find lag time by frame number where biomass just gained 10% of biomass.
  cluster_lag_time=find_cluster_lag_time(single_cluster_area_na_omitted,max_image)
  
  #single_cluster_area_na_omitted=read.delim("remainder_areas_calc.txt",header=FALSE)
  #To find lag time by fitting an exponential line to the growth phase and then calculate when clusters gain 10% of biomass.
  cluster_lag_time=lag_time_by_exp_fit(single_cluster_area_na_omitted,max_image,1/3)
  write.table(cluster_lag_time,"fit_result_Co_C5_S_single_cluster_area_lag_measurement_raw.txt",sep='\t')
}

#Distance between E. coli and the closest S. enterica
for (i in 1:1) {
  source("/Users/apple/Desktop/Single-Cell Microscopy/Nikon_cell_cluster_growth_functions.R")
  setwd(("/Users/apple/Desktop/"))
  Co_E_lag_total=read.delim("CoA1_single_cluster_area_lag_measurement_curated_graph.txt")
  S_location_total=read.delim("S_location_0h_Co_A1.txt")
  E_distance_to_closest_S=E_coli_distance_to_closest_S_enterica(Co_E_lag_total,S_location_total)
  write.table(E_distance_to_closest_S, "E_distance_to_closest_S_Co_A1.txt", sep="\t")
}
E_distance_to_closest_S=read.delim("E_distance_to_closest_S_Co_B5.txt")
dist_lag=cor.test(E_distance_to_closest_S$LagTime, E_distance_to_closest_S$Closest_dist, method="spearman",exact=FALSE)
ggplot(E_distance_to_closest_S,aes(x=Closest_dist,y=LagTime))+geom_point()


#Cluster_size_prediction
for (i in 1:1) {
  source("/Users/apple/Desktop/Single-Cell Microscopy/Nikon_cell_cluster_growth_functions.R")
  setwd(("/Users/apple/Desktop/"))
  E_location_total=read.delim("CoA1_single_cluster_area_lag_measurement_curated_graph.txt")
  S_location_total=read.delim("S_location_0h_Co_A1.txt")
  E_distance_to_closest_S=E_coli_distance_to_closest_S_enterica(E_location_total,S_location_total)
  write.table(E_distance_to_closest_S, "E_distance_to_closest_S_coA.txt", sep="\t")
}
a=cor.test(E_distance_to_closest_S$CellArea,E_distance_to_closest_S$Closest_dist,method="spearman",exact=FALSE)
ggplot(E_distance_to_closest_S,aes(x=Closest_dist,y=CellArea))+geom_point()

#Neighborhood size calculation for single-cluster lag time (Dal Co style plot)
for (j in 1:1) {
  source("/Users/apple/Desktop/Single-Cell Microscopy/Nikon_cell_cluster_growth_functions.R")
  setwd(("/Users/apple/Desktop/"))
  Co_E_lag_total=read.delim("CoA1_single_cluster_area_lag_measurement_curated_graph.txt")
  E_location_total=read.delim("single_cluster_at_0h_CoA1.txt")
  S_location_total=read.delim("S_location_0h_Co_A1.txt")
  
  sprm_rho_neiFrac=c(); sprm_p_val_neiFrac=c(); sprm_rho_NumNei=c()
  rad_size=55
  for (i in rad_size) { #i denotes the radius of the neiborhood size
    spearman_cor_by_ngbr=Corr_rho_ngbr_frac_E_lag_growth_10_percent(Co_E_lag_total,E_location_total,S_location_total,i)
    spearman_lag_time_fracNei=cor.test(spearman_cor_by_ngbr$S_Fraction, spearman_cor_by_ngbr$Lag_Time, method="spearman",exact=FALSE)
    sprm_rho_neiFrac=c(sprm_rho_neiFrac,spearman_lag_time_fracNei$estimate)
    sprm_p_val_neiFrac=c(sprm_p_val_neiFrac,spearman_lag_time_fracNei$p.value)
    
    spearman_lag_time_NumNei=cor.test(spearman_cor_by_ngbr$S_Area,spearman_cor_by_ngbr$Lag_Time, method="spearman",exact=FALSE)
    sprm_rho_NumNei=c(sprm_rho_NumNei,spearman_lag_time_NumNei$estimate)
  }
  lag_corr_S_frac=data.frame(Neighborhood_Size=rad_size,Spearman_rho_FracNei=sprm_rho_neiFrac,Spearman_p_val=sprm_p_val_neiFrac, Spearman_rho_NumNei=sprm_rho_NumNei)
  write.table(lag_corr_S_frac, "fit_DalCo_Co_A1_wo8_rho_S_area_fraction_edge_cell_removed_both.txt", sep="\t")
  write.table(spearman_cor_by_ngbr, "Elag_vs_S_fract_CoA1_55.txt", sep="\t")
  cell_never_exit_lag=Co_E_lag_total[Co_E_lag_total$Lag_Time==8,]
  ggplot(spearman_cor_by_ngbr,aes(x=S_Fraction,y=Lag_Time))+geom_point()
  a=cor.test(spearman_cor_by_ngbr$S_Fraction,spearman_cor_by_ngbr$Lag_Time,method="spearman",exact=FALSE)
}

#Growth time and growth rate tests for single clusters within first 3 hours.
for (k in 1:1) {
  setwd("/Users/apple/Desktop/")  
  source("/Users/apple/Desktop/Single-Cell Microscopy/Nikon_cell_cluster_growth_functions.R")

  area_per_cell_per_frame=read.delim("CoA1_single_cluster_area_na_omitted.txt")
  cell_number=ncol(area_per_cell_per_frame)
  gr_gt_data=matrix(0,cell_number,4)
  time_series=c(0,0.33,0.67,1,1.33,1.67,2,2.33)
  for (i in 1:cell_number) {
    gr_gt_data[i,1]=area_per_cell_per_frame[8,i]
    gr_gt_data[i,2]=area_per_cell_per_frame[9,i]
    
    #Find growth time of a cluster.
    individual_growth_time_point=1; current_area=area_per_cell_per_frame[1,i]; j=2
    while (j<=7) {
      if (area_per_cell_per_frame[j,i]>current_area) {
        current_area=area_per_cell_per_frame[j,i]
        j=j+1
        individual_growth_time_point=individual_growth_time_point+1
      } else {
        #If at a given time point no growth happens, then we consider the duration up to previous
        #time point the growth time.
        j=10000 
      }
    }
    gr_gt_data[i,3]=individual_growth_time_point
    growth_levels=area_per_cell_per_frame[1:individual_growth_time_point,i]
    growth_vs_time=data.frame(time_series[1:individual_growth_time_point],growth_levels)
    colnames(growth_vs_time)[1]="Time"  
    model=lm(growth_levels~Time,data=growth_vs_time)
    gr_gt_data[i,4]=model$coefficients[2]
  }
  write.table(gr_gt_data,"gr_gt_data_Co1.txt",sep='\t')
}
