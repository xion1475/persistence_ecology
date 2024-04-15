library(dplyr)
library("ggplot2")

E_coli_distance_to_closest_S_enterica <- function(Co_E_lag_total,S_enterica_location){
  lag_time=Co_E_lag_total[,2]; cell_number <- Co_E_lag_total[,1]; Dist_closest_S=c()
  E_cell_number=nrow(Co_E_lag_total); S_cell_number=nrow(S_enterica_location)
  for (i in 1:E_cell_number) {
    min_dist=1000000
    for (j in 1:S_cell_number) {
      current_dist=sqrt((Co_E_lag_total[i,3]-S_enterica_location[j,3])^2+(Co_E_lag_total[i,4]-S_enterica_location[j,4])^2)
      min_dist=min(current_dist,min_dist)
    }
    Dist_closest_S=c(Dist_closest_S,min_dist)
  }
  table_cloest_dist=data.frame(CellNumber=cell_number, LagTime=lag_time, Closest_dist=Dist_closest_S)
  return(table_cloest_dist)
}

cluster_area_over_time<-function(total_cluster_loc,max_image) {
  
  #We only study cells that show up since frame 1.
  E_0h=total_cluster_loc[total_cluster_loc$Label==1 & total_cluster_loc$X<970 & total_cluster_loc$Y>52,]
  potential_cell_number=nrow(E_0h); total_area_over_time_data=matrix(0,max_image+2,potential_cell_number)
  total_area_over_time_data[1,]=E_0h[,3]
  
  for (cell_ID in 1:potential_cell_number) {
    x_coord=E_0h[cell_ID,4]; y_coord=E_0h[cell_ID,5]
    total_area_over_time_data[max_image+1,cell_ID]=x_coord; total_area_over_time_data[max_image+2,cell_ID]=y_coord; 
    for (i in 2:max_image) {
      E_later=total_cluster_loc[total_cluster_loc$Label==i,]
      neighbors=E_later[E_later$X>x_coord-10 & E_later$X<x_coord+10 & E_later$Y>y_coord-10 & E_later$Y<y_coord+10,]
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

#This function finds lag time for individual cell area data, where lag time is defined as the minimal frame by which a cell cluster gains 10% biomass (area).
find_cluster_lag_time <- function (single_cluster_area_na_omitted,max_image) {
  total_cluster_number=ncol(single_cluster_area_na_omitted)
  Cluster_ID=seq(1,total_cluster_number,1)
  X_location=single_cluster_area_na_omitted[max_image+1,]
  Y_location=single_cluster_area_na_omitted[max_image+2,]
  lag_time=c()
  
  #To find time it takes for individual clusters to reach 110% of biomass, which we define as lag time.
  for (i in 1:total_cluster_number) {
    initial_biomass=single_cluster_area_na_omitted[1,i]
    image_count=1; yet_to_find_lag_time=TRUE
    while (image_count<=max_image & yet_to_find_lag_time) {
      image_count=image_count+1
      if (single_cluster_area_na_omitted[image_count,i]>=1.2*initial_biomass) {
        yet_to_find_lag_time=FALSE
      }
    }
    lag_time=c(lag_time,image_count)
  }
  
  #Assemble all info to form a data frame.
  total_lag_location=data.frame(Cluster_ID=Cluster_ID,Lag_Time=lag_time,X=X_location,Y=Y_location)
  return(total_lag_location)
}

#To calculate individual cell lag time. Here, lag time is calculated by fitting a exponential
#line to the phase of the growth curve that's bigger than the initial biomass, and then back-calculate
#the time it takes for the cluster to gain 10% of biomass.
lag_time_by_exp_fit <- function (single_cluster_area_na_omitted,max_image,image_interval) {
  total_cluster_number=ncol(single_cluster_area_na_omitted)
  Cluster_ID=seq(1,total_cluster_number,1)
  X_location=single_cluster_area_na_omitted[max_image+1,]
  Y_location=single_cluster_area_na_omitted[max_image+2,]
  lag_time=c()
  
  #To find out lag time for individual clusters.
  for (i in 1:total_cluster_number) {
    initial_biomass=single_cluster_area_na_omitted[1,i]
    current_cluster=single_cluster_area_na_omitted[1:max_image,i]
    #To convert the areas to exponential values.
    current_cluster=log10(current_cluster/initial_biomass)
    
    #To find out whether there exists a biomass that's bigger than 0.
    cluster_bigger_than_0.05=which(current_cluster>0.05)
    if (length(cluster_bigger_than_0.05)==0) {
      current_lag_time=(max_image)*image_interval
    } else {
      right_side_id=min(cluster_bigger_than_0.05)
      left_side_id=right_side_id-1
      growth_time=c((left_side_id-1)*image_interval,(right_side_id-1)*image_interval)
      growth_biomass=c(current_cluster[left_side_id],current_cluster[right_side_id])
      growth_simulation=data.frame(growth_time,growth_biomass)
      growth_fit=summary(lm(growth_biomass~growth_time,data=growth_simulation))
      k=growth_fit$coefficients[2,1]; b=growth_fit$coefficients[1,1]
      current_lag_time=(0.05-b)/k           
    }
    lag_time=c(lag_time,current_lag_time)
    #To fit a linear line to the log-transformed biomass data for the growing phase.
    #growing_phase=current_cluster[2:max_image]
    #growth=data.frame(growth_phase_time_points,growing_phase)
    #exponential_fit=summary(lm(growing_phase~growth_phase_time_points,data=growth))
    #k=exponential_fit$coefficients[2,1]; b=exponential_fit$coefficients[1,1]
    #current_lag=(0.04139-b)/k
    #growth_phase_time_points=seq(1/3,(max_image-1)/3,1/3)
    #if (current_lag>0) {
    #  lag_time=c(lag_time,current_lag)
    #} else {
    #  lag_time=c(lag_time,(max_image)/3)
    #}
  }
  #Assemble all info to form a data frame.
  total_lag_location=data.frame(Cluster_ID=Cluster_ID,Lag_Time=lag_time,X=X_location,Y=Y_location)
  return(total_lag_location)
}



Corr_rho_ngbr_frac_E_lag_growth_10_percent <- function(Co_E_lag_total,E_location_total,S_location_total,radius_ngbr){
  E_lag_number=nrow(Co_E_lag_total); lag_time <- Co_E_lag_total[,2]; 
  S_neighbor_number=c(); E_neighbor_number=c(); S_neighbor_frac=c()
  
  for (i in 1:E_lag_number){
    E_lag_coor_x=Co_E_lag_total[i,3]; E_lag_coor_y=Co_E_lag_total[i,4]
    E_neighbor=E_location_total[E_location_total$X>E_lag_coor_x-radius_ngbr & E_location_total$X<E_lag_coor_x+radius_ngbr & E_location_total$Y>=E_lag_coor_y-radius_ngbr & E_location_total$Y<=E_lag_coor_y+radius_ngbr,]
    E_in_neighbor=nrow(E_neighbor); c=0     #Count how many E cells there are around a chosen E cell.
    if (E_in_neighbor>0) {
      for (j in 1:E_in_neighbor) {
        other_E_x=E_neighbor[j,4]; other_E_y=E_neighbor[j,5]
        dist=sqrt((E_lag_coor_x-other_E_x)^2+(E_lag_coor_y-other_E_y)^2)
        if (dist<=radius_ngbr) {c=c+E_neighbor[j,3]} #Calculate total Eo biomass within the neighborhood.
      } 
    }
    
    #Count how many S cells there are around a chosen E cell.
    S_neighbor=S_location_total[S_location_total$X>E_lag_coor_x-radius_ngbr & S_location_total$X<E_lag_coor_x+radius_ngbr & S_location_total$Y>E_lag_coor_y-radius_ngbr & S_location_total$Y<E_lag_coor_y+radius_ngbr,]
    S_in_neighbor=nrow(S_neighbor); k=0
    if (S_in_neighbor>0) {
      for (j in 1:S_in_neighbor) {
        S_cell_coor_x=S_neighbor[j,3]; S_cell_coor_y=S_neighbor[j,4]
        dist=sqrt((E_lag_coor_x-S_cell_coor_x)^2+(E_lag_coor_y-S_cell_coor_y)^2)
        if (dist<=radius_ngbr) {k=k+S_neighbor[j,2]} 
      }
    }
    
    S_neighbor_number=c(S_neighbor_number,k)
    E_neighbor_number=c(E_neighbor_number,c)
    S_neighbor_frac=c(S_neighbor_frac,k/(k+c))
  }
  
  E_cell_lag_time_corr_w_S_frac=data.frame(Lag_Time=lag_time, S_Fraction=S_neighbor_frac, S_Area=S_neighbor_number)
  return(E_cell_lag_time_corr_w_S_frac)
}  
#ggplot(E_cell_lag_time_corr_w_S_frac,aes(x=S_Fraction,y=Lag_Time))+geom_point()
#ggplot(E_cell_lag_time_corr_w_S_frac,aes(x=S_Fraction,y=log(Lag_Time)))+geom_point()
#cor.test(log(E_cell_lag_time_corr_w_S_frac$Lag_Time),E_cell_lag_time_corr_w_S_frac$S_Fraction,method="spearman",exact=FALSE)

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


Corr_rho_E_lag_by_S_lag <- function(Co_E_lag_total,Co_S_lag_total,radius_ngbr){
  E_lag_number=nrow(Co_E_lag_total); E_lag_time <- Co_E_lag_total[,2]; S_neighbor_lag_weighted=c()

  for (i in 1:E_lag_number){
    E_lag_coor_x=Co_E_lag_total[i,3]; E_lag_coor_y=Co_E_lag_total[i,4]
    #Count how many S cells there are around a chosen E cell.
    S_neighbor=Co_S_lag_total[Co_S_lag_total$X>E_lag_coor_x-radius_ngbr & Co_S_lag_total$X<E_lag_coor_x+radius_ngbr & Co_S_lag_total$Y>E_lag_coor_y-radius_ngbr & Co_S_lag_total$Y<E_lag_coor_y+radius_ngbr,]
    S_in_neighbor=nrow(S_neighbor); weighted_S_lag=100
    if (S_in_neighbor>0) {
      for (j in 1:S_in_neighbor) {
        S_cell_coor_x=S_neighbor[j,3]; S_cell_coor_y=S_neighbor[j,4]
        dist=sqrt((E_lag_coor_x-S_cell_coor_x)^2+(E_lag_coor_y-S_cell_coor_y)^2)
        if (dist<=radius_ngbr) {
          min_S_lag=min(min_S_lag,S_neighbor[j,2])
        } 
      }
    }
    S_neighbor_lag_weighted=c(S_neighbor_lag_weighted,min_S_lag)
  }
  
  E_cell_lag_time_corr_w_S_frac=data.frame(Lag_Time=E_lag_time, S_Lag=S_neighbor_lag_weighted)
  return(E_cell_lag_time_corr_w_S_frac)
}  

#This function calculates the time it takes for a current cell to gain 10% biomass.
lag_time_calculator_10p_biomass <- function (current_lineage,current_tree_ID,imaging_interval){
  
  current_lineage=data.frame(current_lineage)
  lag_time_data=c(current_tree_ID)
  starting_cell_ID=current_lineage[1,2]
  track_length=nrow(current_lineage)
  last_cell_ID=current_lineage[track_length,2]
  cell_lag_time=10000 #Arbitrary inital value for lag time measurement.
  current_min_area=current_lineage[1,9]; pointer=2; yet_to_find_lag_time=TRUE
  
  while (pointer<=track_length && yet_to_find_lag_time) {
    current_area=current_lineage[pointer,9]
    if (current_area>current_min_area*1.1){ #To identify when cell gains that 10% biomass, if at all.
      
      right_side_time=(current_lineage[pointer,4]-1)*imaging_interval
      left_side_time=(current_lineage[pointer-1,4]-1)*imaging_interval
      right_side_relative_log_biomass=log10(current_lineage[pointer,9]/current_min_area)
      left_side_relative_log_biomass=log10(current_lineage[pointer-1,9]/current_min_area)
      if (left_side_time==right_side_time) {
        left_side_time=(current_lineage[pointer-2,4]-1)*imaging_interval
        left_side_relative_log_biomass=log10(current_lineage[pointer-2,9]/current_min_area)
      }
      
      growth_time=c(left_side_time,right_side_time)
      growth_biomass=c(left_side_relative_log_biomass,right_side_relative_log_biomass)
      growth_simulation=data.frame(growth_time,growth_biomass)
      growth_fit=summary(lm(growth_biomass~growth_time,data=growth_simulation))
      k=growth_fit$coefficients[2,1]; b=growth_fit$coefficients[1,1]
      cell_lag_time=(log10(1.1)-b)/k  
      
      yet_to_find_lag_time=FALSE #To confirm we found the lag time and exit while loop.
      
    } else if (current_area<current_min_area) {
      current_min_area=current_area
    }
    pointer=pointer+1
  }  
  
  lag_time_data=c(lag_time_data,cell_lag_time)
  lag_time_data=c(lag_time_data,current_lineage[1,6])
  lag_time_data=c(lag_time_data,current_lineage[1,7])
  return(lag_time_data)
}

lag_time_calculator_0p_biomass <- function (current_lineage,current_tree_ID,imaging_interval){
  
  current_lineage=data.frame(current_lineage)
  lag_time_data=c(current_tree_ID)
  starting_cell_ID=current_lineage[1,2]
  track_length=nrow(current_lineage)
  last_cell_ID=current_lineage[track_length,2]
  cell_lag_time=10000 #Arbitrary inital value for lag time measurement.
  current_min_area=current_lineage[1,9]; pointer=2; yet_to_find_lag_time=TRUE
  
  while (pointer<=track_length && yet_to_find_lag_time) {
    current_area=current_lineage[pointer,9]
    if (current_area>current_min_area){ #To identify when cell gains that 10% biomass, if at all.
      cell_lag_time=(current_lineage[pointer,4]-1)*imaging_interval
      yet_to_find_lag_time=FALSE #To confirm we found the lag time and exit while loop.
    } 
    pointer=pointer+1
  }  
  
  lag_time_data=c(lag_time_data,cell_lag_time)
  lag_time_data=c(lag_time_data,current_lineage[1,6])
  lag_time_data=c(lag_time_data,current_lineage[1,7])
  return(lag_time_data)
}


calculate_single_cell_growth_curve <- function (current_lineage,imaging_interval,max_frame) {
  time_data=seq(0,(max_frame-1)*imaging_interval,imaging_interval)
  biomass_data=c()
  for (i in 1:max_frame) {
    current_image=current_lineage[current_lineage$ND.T==i,]
    current_biomass=sum(current_image$Area..µm..)
    biomass_data=c(biomass_data,current_biomass)
  }
  c_total_biomass=data.frame(Time=time_data,Biomass=biomass_data)
  return(c_total_biomass)
}

lag_time_calculator_10p_biomass_with_frame_numbers <- function (current_lineage,current_tree_ID,imaging_interval){
  
  current_lineage=data.frame(current_lineage)
  lag_time_data=c(current_tree_ID)
  starting_cell_ID=current_lineage[1,2]
  track_length=nrow(current_lineage)
  last_cell_ID=current_lineage[track_length,2]
  cell_lag_time=10000 #Arbitrary inital value for lag time measurement.
  current_min_area=current_lineage[1,9]; pointer=2; yet_to_find_lag_time=TRUE
  
  max_frame=max(current_lineage$ND.T)
  growth_curve=calculate_single_cell_growth_curve(current_lineage,imaging_interval,max_frame)
  
  while (pointer<=track_length && yet_to_find_lag_time) {
    current_area=current_lineage[pointer,9]
    if (current_area>current_min_area*1.1){ #To identify when cell gains that 10% biomass, if at all.
      cell_lag_time=(current_lineage[pointer,4]-1)*imaging_interval
      yet_to_find_lag_time=FALSE #To confirm we found the lag time and exit while loop.
    }  else if (current_area<current_min_area) {
      current_min_area=current_area
    }
    pointer=pointer+1
  }  
  
  lag_time_data=c(lag_time_data,cell_lag_time)
  lag_time_data=c(lag_time_data,current_lineage[1,6])
  lag_time_data=c(lag_time_data,current_lineage[1,7])
  return(lag_time_data)
}

#This function calculates cell lag time defined by time it takes for it to gain 10% biomass.
single_cell_lag_time_biomass <- function (cell_data,imaging_interval) {
  
  lag_time <- c(); lineage_number <- c(); coo_x<- c();coo_y<- c(); 
  cells_present=TRUE;current_tree_ID=1
  while (cells_present) {
    current_lineage=cell_data[cell_data$TreeID==current_tree_ID,]
    track_start=current_lineage[1,4]
    
    #We only analyze cells that already appearred in the first frame. Cells that show up in later
    #frames are offsprings of other cells.
    if (track_start==1) {
      label_check=current_lineage[current_lineage$ND.T==1,]
      actual_cell_number=nrow(label_check)
      
      #Due to the tracking setting of the Nikon Elements software, two cells nearby that 
      #appear in frame 1 can be considered offsprings from the same lineage. Their division 
      #may have happened before the start of imaging, which we do not consider. Here, we still 
      #measure their division starting from imaging. So we will consider them separately.
      i=1; starting_row_number=0; end_row_number=0
      while (i<actual_cell_number) {
        starting_row_number=which(row.names(current_lineage)==row.names(label_check[i,]))
        end_row_number=which(row.names(current_lineage)==row.names(label_check[i+1,]))-1
        
        actual_current_cell=current_lineage[starting_row_number:end_row_number,]
        #lag_time_data=lag_time_calculator_10p_biomass(actual_current_cell,current_tree_ID+i*0.1,imaging_interval)
        lag_time_data=lag_time_calculator_10p_biomass_with_frame_numbers(actual_current_cell,current_tree_ID+i*0.1,imaging_interval)
        lineage_number <- c(lineage_number,lag_time_data[1])
        lag_time <- c(lag_time,lag_time_data[2])
        coo_x<- c(coo_x,lag_time_data[3])
        coo_y<- c(coo_y,lag_time_data[4])
        
        i=i+1
      }
      
      starting_row_number=end_row_number+1; end_row_number=nrow(current_lineage)
      lag_time_data=lag_time_calculator_10p_biomass_with_frame_numbers(current_lineage[starting_row_number:end_row_number,],current_tree_ID,imaging_interval)
      lineage_number <- c(lineage_number,lag_time_data[1])
      lag_time <- c(lag_time,lag_time_data[2])
      coo_x<- c(coo_x,lag_time_data[3])
      coo_y<- c(coo_y,lag_time_data[4])
      
      current_tree_ID=current_tree_ID+1
    } else {cells_present=FALSE}
  }
  c_res=data.frame(CellNumber=lineage_number, LagTime=lag_time, Coord_x_um=coo_x, Coord_y_um=coo_y)
  return(c_res)
}

#a=calculate_single_cell_growth_curve(current_lineage,1/3,max_frame)
#ggplot(a,aes(x=Time,y=Biomass))+geom_point()
