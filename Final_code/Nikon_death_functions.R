library(dplyr)
library("ggplot2")

#This function tracks single-cell area over time.
cell_area_over_time<-function(total_cluster_loc,max_image) {
  
  #We only study cells that show up since frame 1.
  E_0h=total_cluster_loc[total_cluster_loc$Label==1 & total_cluster_loc$Area>25,]
  potential_cell_number=nrow(E_0h); total_area_over_time_data=matrix(0,max_image+2,potential_cell_number)
  total_area_over_time_data[1,]=E_0h[,3]
  
  for (cell_ID in 1:potential_cell_number) {
    x_coord=E_0h[cell_ID,4]; y_coord=E_0h[cell_ID,5]
    total_area_over_time_data[max_image+1,cell_ID]=x_coord; total_area_over_time_data[max_image+2,cell_ID]=y_coord; 
    for (i in 2:max_image) {
      E_later=total_cluster_loc[total_cluster_loc$Label==i,]
      neighbors=E_later[E_later$X>x_coord-20 & E_later$X<x_coord+20 & E_later$Y>y_coord-20 & E_later$Y<y_coord+20,]
      cell_num_in_nei=nrow(neighbors)
      if (cell_num_in_nei==0) {
        total_area_over_time_data[i,cell_ID]=NA
      } else if (cell_num_in_nei==1) {
        total_area_over_time_data[i,cell_ID]=neighbors[1,3]
      } else if (cell_num_in_nei>1) {
        high_dist=1000; correct_area=1000
        for (j in 1:cell_num_in_nei) {
          xj=neighbors[j,4]; yj=neighbors[j,5]
          current_dist=sqrt((x_coord-xj)^2+(y_coord-yj)^2)
          if (current_dist<high_dist) {
            high_dist=current_dist
            correct_area=neighbors[j,3]
          }
        }
        total_area_over_time_data[i,cell_ID]=correct_area
      }
    }
    
  } 
  return(total_area_over_time_data)
}

#This function correlates Eo single-cell death time with location in the Dal Co style plot.
Corr_rho_ngbr_frac_E_death <- function(Co_E_death_total,S_location_total,radius_ngbr){

  E_death_number=nrow(Co_E_death_total); death_time <- Co_E_death_total[,2]; 
  S_neighbor_number=c(); E_neighbor_number=c(); S_neighbor_frac=c()
  
  for (i in 1:E_death_number){
    E_lag_coor_x=Co_E_death_total[i,3]; E_lag_coor_y=Co_E_death_total[i,4]
    E_neighbor=Co_E_death_total[Co_E_death_total$X>E_lag_coor_x-radius_ngbr & Co_E_death_total$X<E_lag_coor_x+radius_ngbr & Co_E_death_total$Y>E_lag_coor_y-radius_ngbr & Co_E_death_total$Y<E_lag_coor_y+radius_ngbr,]
    E_in_neighbor=nrow(E_neighbor); c=0     #Count how many E cells there are around a chosen E cell.
    if (E_in_neighbor>0) {
      for (j in 1:E_in_neighbor) {
        other_E_x=E_neighbor[j,3]; other_E_y=E_neighbor[j,4]
        dist=sqrt((E_lag_coor_x-other_E_x)^2+(E_lag_coor_y-other_E_y)^2)
        if (dist<=radius_ngbr) 
          {c=c+E_neighbor[j,5]} #Calculate total Eo biomass within the neighborhood.
      } 
    }
    
    #Count how many S cells there are around a chosen E cell.
    S_neighbor=S_location_total[S_location_total$X>E_lag_coor_x-radius_ngbr & S_location_total$X<E_lag_coor_x+radius_ngbr & S_location_total$Y>E_lag_coor_y-radius_ngbr & S_location_total$Y<E_lag_coor_y+radius_ngbr,]
    S_in_neighbor=nrow(S_neighbor); k=0
    if (S_in_neighbor>0) {
      for (j in 1:S_in_neighbor) {
        S_cell_coor_x=S_neighbor[j,3]; S_cell_coor_y=S_neighbor[j,4]
        dist=sqrt((E_lag_coor_x-S_cell_coor_x)^2+(E_lag_coor_y-S_cell_coor_y)^2)
        if (dist<=radius_ngbr) {k=k+S_neighbor[j,5]} 
      }
    }
    
    S_neighbor_number=c(S_neighbor_number,k)
    E_neighbor_number=c(E_neighbor_number,c)
    S_neighbor_frac=c(S_neighbor_frac,k/(k+c))
  }
  
  E_cell_lag_time_corr_w_S_frac=data.frame(Death_Time=death_time, S_Fraction=S_neighbor_frac, S_Number=S_neighbor_number)
  return(E_cell_lag_time_corr_w_S_frac)
}  
#ggplot(E_cell_lag_time_corr_w_S_frac,aes(x=S_Fraction,y=Lag_Time))+geom_point()
#ggplot(E_cell_lag_time_corr_w_S_frac,aes(x=S_Fraction,y=log(Lag_Time)))+geom_point()
#cor.test(log(E_cell_lag_time_corr_w_S_frac$Lag_Time),E_cell_lag_time_corr_w_S_frac$S_Fraction,method="spearman",exact=FALSE)

#This function correlates Eo single-cell death time with fraction of alive So within a given neighborhood.
Corr_rho_ngbr_frac_E_death_alive_S <- function(Co_E_death_total,S_location_total,radius_ngbr){
  
  E_death_number=nrow(Co_E_death_total); death_time <- Co_E_death_total[,2]; 
  S_neighbor_number=c(); E_neighbor_number=c(); S_neighbor_frac=c()
  
  for (i in 1:E_death_number){
    E_lag_coor_x=Co_E_death_total[i,3]; E_lag_coor_y=Co_E_death_total[i,4]
    E_neighbor=Co_E_death_total[Co_E_death_total$X>E_lag_coor_x-radius_ngbr & Co_E_death_total$X<E_lag_coor_x+radius_ngbr & Co_E_death_total$Y>E_lag_coor_y-radius_ngbr & Co_E_death_total$Y<E_lag_coor_y+radius_ngbr,]
    E_in_neighbor=nrow(E_neighbor); c=0     #Count how many E cells there are around a chosen E cell.
    if (E_in_neighbor>0) {
      for (j in 1:E_in_neighbor) {
        other_E_x=E_neighbor[j,3]; other_E_y=E_neighbor[j,4]
        dist=sqrt((E_lag_coor_x-other_E_x)^2+(E_lag_coor_y-other_E_y)^2)
        if (dist<=radius_ngbr & E_neighbor[j,2]>=Co_E_death_total[i,2]) 
          {c=c+E_neighbor[j,5]} #Calculate total Eo biomass within the neighborhood.
      } 
    }
    
    #Count how many S cells there are around a chosen E cell.
    S_neighbor=S_location_total[S_location_total$X>E_lag_coor_x-radius_ngbr & S_location_total$X<E_lag_coor_x+radius_ngbr & S_location_total$Y>E_lag_coor_y-radius_ngbr & S_location_total$Y<E_lag_coor_y+radius_ngbr,]
    S_in_neighbor=nrow(S_neighbor); k=0
    if (S_in_neighbor>0) {
      for (j in 1:S_in_neighbor) {
        S_cell_coor_x=S_neighbor[j,3]; S_cell_coor_y=S_neighbor[j,4]
        dist=sqrt((E_lag_coor_x-S_cell_coor_x)^2+(E_lag_coor_y-S_cell_coor_y)^2)
        if (dist<=radius_ngbr & S_neighbor[j,2]>=Co_E_death_total[i,2]) 
          {k=k+S_neighbor[j,5]} 
      }
    }
    
    S_neighbor_number=c(S_neighbor_number,k)
    E_neighbor_number=c(E_neighbor_number,c)
    S_neighbor_frac=c(S_neighbor_frac,k/(k+c))
  }
  
  E_cell_lag_time_corr_w_S_frac=data.frame(Death_Time=death_time, S_Fraction=S_neighbor_frac, S_Number=S_neighbor_number)
  return(E_cell_lag_time_corr_w_S_frac)
}  

#This function finds lag time for individual cell area data analyzed previously.
death_disappearance_time <- function (single_cluster_area_over_time,max_image) {
  total_cluster_number=ncol(single_cluster_area_over_time)
  Cluster_ID=seq(1,total_cluster_number,1)
  X_location=single_cluster_area_over_time[max_image+1,]
  Y_location=single_cluster_area_over_time[max_image+2,]
  disappearance_time=c()
  
  for (i in 1:total_cluster_number) {
    image_count=1; yet_to_find_death_time=TRUE
    while (image_count<=max_image & yet_to_find_death_time) {
      image_count=image_count+1
      if (is.na(single_cluster_area_over_time[image_count,i])) {yet_to_find_death_time=FALSE}
    }
    disappearance_time=c(disappearance_time,image_count)    
  }
  total_death_time_location=data.frame(Cluster_ID=Cluster_ID,Death_Time=disappearance_time,X=X_location,Y=Y_location)
  return(total_death_time_location)
}

#This function correlates Eo death time with the distance to the closest, alive So cell.
death_time_to_alive_closest_So <- function(E_location_total,S_enterica_location){
  E_death_time <- E_location_total[,2]; Dist_closest_S=c()
  E_cell_number=nrow(E_location_total); S_cell_number=nrow(S_enterica_location)
  for (i in 1:E_cell_number) {
    min_dist=1000000
    alive_total_cell_num=sum(E_location_total$Death_Time>=E_location_total[i,2])+sum(S_enterica_location$Death_Time>=E_location_total[i,2])
    for (j in 1:S_cell_number) {
      current_dist=sqrt((E_location_total[i,3]-S_enterica_location[j,3])^2+(E_location_total[i,4]-S_enterica_location[j,4])^2)
      if (S_enterica_location[j,2]>=E_location_total[i,2]) {
        min_dist=min(current_dist,min_dist)
      }
    }
    Dist_closest_S=c(Dist_closest_S,min_dist)
  }
  table_cloest_dist=data.frame(Eo_Death_Time=E_death_time, Closest_dist=Dist_closest_S)
  return(table_cloest_dist)
}

#This function correlates Eo death time with the distance to the closest So cell.
death_time_to_closest_So <- function(E_location_total,S_enterica_location){
  E_death_time <- E_location_total[,2]; Dist_closest_S=c()
  E_cell_number=nrow(E_location_total); S_cell_number=nrow(S_enterica_location)
  for (i in 1:E_cell_number) {
    min_dist=1000000
    for (j in 1:S_cell_number) {
      current_dist=sqrt((E_location_total[i,3]-S_enterica_location[j,3])^2+(E_location_total[i,4]-S_enterica_location[j,4])^2)
      min_dist=min(current_dist,min_dist)
    }
    Dist_closest_S=c(Dist_closest_S,min_dist)
  }
  table_cloest_dist=data.frame(Eo_Death_Time=E_death_time, Closest_dist=Dist_closest_S)
  return(table_cloest_dist)
}

#This function calculates the death time of the single So cell that is closest to the focal Eo cell.
death_time_to_closest_So_death_time <- function(E_location_total,S_enterica_location){
  E_death_time <- E_location_total[,2]; Dist_closest_S_death_time=c()
  E_cell_number=nrow(E_location_total); S_cell_number=nrow(S_enterica_location)
  for (i in 1:E_cell_number) {
    min_dist=1000000
    for (j in 1:S_cell_number) {
      current_dist=sqrt((E_location_total[i,3]-S_enterica_location[j,3])^2+(E_location_total[i,4]-S_enterica_location[j,4])^2)
      if (current_dist<=min_dist) {
        min_dist=current_dist
        min_death_time_S=S_enterica_location[j,2]
      }
    }
    Dist_closest_S_death_time=c(Dist_closest_S_death_time,min_death_time_S)
  }
  table_cloest_dist=data.frame(Eo_Death_Time=E_death_time, Closest_So_Death_Time=Dist_closest_S_death_time)
  return(table_cloest_dist)
}


#This function uses average, max, or min So death time in a neighborhood of radius (radius_ngbr) to predict Eo death time within the same neighborhood.
corr_rho_So_death_for_Eo_death_at_arbitrary_radius <- function(Co_E_death_total,Co_S_death_total,radius_ngbr){
  
  #Co_E_death_total=Co_E_death_total[Co_E_death_total$X>200 & Co_E_death_total$X<3072-200 & Co_E_death_total$Y>200 & Co_E_death_total$Y<3072-200,]
  E_death_number=nrow(Co_E_death_total); Eo_death_time <- c(); min_S_death=c(); max_S_death=c(); mean_S_death=c()
  
  for (i in 1:E_death_number){
    
    E_lag_coor_x=Co_E_death_total[i,3]; E_lag_coor_y=Co_E_death_total[i,4]
    Co_S_death_total[,5]=0
    #Count how many S cells there are around a chosen E cell.
    S_neighbor=Co_S_death_total[Co_S_death_total$X>E_lag_coor_x-radius_ngbr & Co_S_death_total$X<E_lag_coor_x+radius_ngbr & Co_S_death_total$Y>E_lag_coor_y-radius_ngbr & Co_S_death_total$Y<E_lag_coor_y+radius_ngbr,]
    
    if (nrow(S_neighbor)>0) {
      for (j in 1:nrow(S_neighbor)) {
        S_cell_coor_x=S_neighbor[j,3]; S_cell_coor_y=S_neighbor[j,4]
        dist=sqrt((E_lag_coor_x-S_cell_coor_x)^2+(E_lag_coor_y-S_cell_coor_y)^2)
        if (dist<=radius_ngbr) {S_neighbor[j,5]=1} 
      }
      S_actual_neighbor=S_neighbor[S_neighbor[,5]==1,]; S_actual_neighbor_number=nrow(S_actual_neighbor)
      if (S_actual_neighbor_number>0) {
        Eo_death_time=c(Eo_death_time,Co_E_death_total[i,2])
        min_S_death=c(min_S_death,min(S_actual_neighbor[,2]))
        max_S_death=c(max_S_death,max(S_actual_neighbor[,2]))
        mean_S_death=c(mean_S_death,mean(S_actual_neighbor[,2]))
      }
    }
  }
  
  E_cell_lag_time_corr_w_S_frac=data.frame(Eo_Death_Time=Eo_death_time, Min_S_Death_Time=min_S_death, Max_S_Death_Time=max_S_death, Mean_S_Death_Time=mean_S_death)
  #ggplot(E_cell_lag_time_corr_w_S_frac,aes(x=Max_S_Death_Time,y=Eo_Death_Time))+geom_point()
  ggplot(E_cell_lag_time_corr_w_S_frac,aes(x=Mean_S_Death_Time,y=Eo_Death_Time))+geom_point()
  return(E_cell_lag_time_corr_w_S_frac)
}  

predict_Eo_death_time <- function (Co_E_death_total,Co_S_death_total) {
  sprm_rho_Min=c(); sprm_rho_Max=c(); sprm_rho_Mean=c()
  sprm_p_val_Min=c(); sprm_p_val_Max=c(); sprm_p_val_Mean=c()
  rad_size=seq(500,3000,500)
  rad_size=seq(10,400,10)
  rad_size=10
  for (i in rad_size) { #i denotes the radius of the neiborhood size
    #Call the correlation function.
    spearman_cor_by_ngbr=corr_rho_So_death_for_Eo_death_at_arbitrary_radius(Co_E_death_total,Co_S_death_total,i)
    
    #Record correlation between minimal So death time with focal Eo death time.
    spearman_death_time_Min=cor.test(spearman_cor_by_ngbr$Eo_Death_Time, spearman_cor_by_ngbr$Min_S_Death_Time, method="spearman",exact=FALSE)
    sprm_rho_Min=c(sprm_rho_Min,spearman_death_time_Min$estimate)
    sprm_p_val_Min=c(sprm_p_val_Min,spearman_death_time_Min$p.value)
    ggplot(spearman_cor_by_ngbr,aes(x=Min_S_Death_Time,y=Eo_Death_Time))+geom_point()
    
    #Record correlation between max So death time with focal Eo death time.
    spearman_death_time_Max=cor.test(spearman_cor_by_ngbr$Eo_Death_Time, spearman_cor_by_ngbr$Max_S_Death_Time, method="spearman",exact=FALSE)
    sprm_rho_Max=c(sprm_rho_Max,spearman_death_time_Max$estimate)
    sprm_p_val_Max=c(sprm_p_val_Max,spearman_death_time_Max$p.value)
    
    #Record correlation between mean So death time with focal Eo death time.
    spearman_death_time_Mean=cor.test(spearman_cor_by_ngbr$Eo_Death_Time, spearman_cor_by_ngbr$Mean_S_Death_Time, method="spearman",exact=FALSE)
    sprm_rho_Mean=c(sprm_rho_Mean,spearman_death_time_Mean$estimate)
    sprm_p_val_Mean=c(sprm_p_val_Mean,spearman_death_time_Mean$p.value)
    
  }
  lag_corr_S_frac=data.frame(Neighborhood_Size=rad_size,sprm_rho_Min,sprm_p_val_Min,sprm_rho_Max,sprm_p_val_Max,sprm_rho_Mean,sprm_p_val_Mean)
  return(lag_corr_S_frac)
}

#This function finds the weighted So death time by distance around a focal Eo.
Eo_death_time_weighted_by_So <- function (Co_E_death_total,Co_S_death_total) {
 
  E_death_time <- Co_E_death_total[,2]; Weighted_S_Death=c()
  E_cell_number=nrow(Co_E_death_total); S_cell_number=nrow(Co_S_death_total)
  for (i in 1:E_cell_number) {
    current_weighted_death_time=0
    for (j in 1:S_cell_number) {
      current_dist=sqrt((Co_E_death_total[i,3]-Co_S_death_total[j,3])^2+(Co_E_death_total[i,4]-Co_S_death_total[j,4])^2)
      current_weighted_death_time=current_weighted_death_time+Co_S_death_total[j,2]/current_dist
    }
    Weighted_S_Death=c(Weighted_S_Death,current_weighted_death_time)
  }
  table_cloest_dist=data.frame(Eo_Death_Time=E_death_time, Weighted_S_Death=Weighted_S_Death)
  return(table_cloest_dist)
}

Eo_death_time_predicted_by_x_nearest_So <- function (Co_E_death_total,Co_S_death_total,x) {
  
  E_death_time <- Co_E_death_total[,2]; Meaninx=c(); Mininx=c(); Maxinx=c(); Weighted=c()
  E_cell_number=nrow(Co_E_death_total); S_cell_number=nrow(Co_S_death_total)
  for (i in 1:E_cell_number) {
    current_weighted_death_time=0
    duplicate_Co_S_total=Co_S_death_total
    for (j in 1:S_cell_number) {
      current_dist=sqrt((Co_E_death_total[i,3]-Co_S_death_total[j,3])^2+(Co_E_death_total[i,4]-Co_S_death_total[j,4])^2)
      duplicate_Co_S_total[j,5]=current_dist
    }
    duplicate_Co_S_total=with(duplicate_Co_S_total, duplicate_Co_S_total[order(Hour0_Area),])
    first_x=duplicate_Co_S_total[1:x,]
    Meaninx=c(Meaninx,mean(first_x$Death_Time)); Mininx=c(Mininx, min(first_x$Death_Time)); Maxinx=c(Maxinx, max(first_x$Death_Time))
    
    current_weighted_death=0
    for (k in 1:x) {
      current_weighted_death=current_weighted_death+first_x[k,2]/first_x[k,5]
    }
    Weighted=c(Weighted,current_weighted_death)
    
  }
  table_cloest_dist=data.frame(Eo_Death_Time=E_death_time, Mean_in_x=Meaninx, Min_in_x=Mininx, Max_in_x=Maxinx, Weighted_by_distance=Weighted)
  return(table_cloest_dist)
}


#This function calculates all E. coli's distance to nearby So cell at the E. coli death time, and then standardizes
#that distance against the average cell-cell distance based on the number of total alive cells at that point.
death_time_to_alive_closest_So_with_standardization <- function(E_location_total,S_enterica_location){
  E_death_time <- E_location_total[,2]; Dist_closest_S=c()
  E_cell_number=nrow(E_location_total); S_cell_number=nrow(S_enterica_location)
  for (i in 1:E_cell_number) {
    min_dist=1000000
    alive_total_cell_num=nrow(E_location_total[E_location_total$Death_Time>=E_location_total[i,2],])+nrow(S_enterica_location[S_enterica_location$Death_Time>=E_location_total[i,2],])
    average_dist=sqrt(3052*3052/alive_total_cell_num)
    for (j in 1:S_cell_number) {
      current_dist=sqrt((E_location_total[i,3]-S_enterica_location[j,3])^2+(E_location_total[i,4]-S_enterica_location[j,4])^2)
      if (S_enterica_location[j,2]>=E_location_total[i,2]) {
        min_dist=min(current_dist,min_dist)
      }
    }
    Dist_closest_S=c(Dist_closest_S,min_dist/average_dist)
  }
  table_cloest_dist=data.frame(Eo_Death_Time=E_death_time, Closest_dist=Dist_closest_S)
  return(table_cloest_dist)
}

#This function calculates all cells distance of a focal species to nearby partner cells at the focus cell death time.
death_time_to_alive_closest_partner_without_standardization <- function(E_location_total,S_enterica_location){
  E_death_time <- E_location_total[,2]; Dist_closest_S=c()
  E_cell_number=nrow(E_location_total); S_cell_number=nrow(S_enterica_location)
  for (i in 1:E_cell_number) {
    min_dist=1000000
    for (j in 1:S_cell_number) {
      current_dist=sqrt((E_location_total[i,3]-S_enterica_location[j,3])^2+(E_location_total[i,4]-S_enterica_location[j,4])^2)
      if (S_enterica_location[j,2]>=E_location_total[i,2]) {
        min_dist=min(current_dist,min_dist)
      }
    }
    Dist_closest_S=c(Dist_closest_S,min_dist)
  }
  table_cloest_dist=data.frame(Eo_Death_Time=E_death_time, Closest_dist=Dist_closest_S)
  return(table_cloest_dist)
}

