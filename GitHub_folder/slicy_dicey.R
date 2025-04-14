slice_sampler <- function(w, param_bound, cond_dens,start_val, params_for_cond_dens, movement_num = .1){

  yval <- log(runif(1)) + do.call(cond_dens, c(start_val, params_for_cond_dens))
  
  ## Set the width-because I know this is unimodel, not going to be too wide
  left <- start_val - runif(1)*w
  right <- start_val + runif(1)*w
  if(right >=param_bound[2]){
    right = param_bound[2]
  }
  if(left <= param_bound[1]){
    left = param_bound[1]
  }
  
  while(yval < do.call(cond_dens, c(left, params_for_cond_dens))){
    # print(left)
    # print(do.call(cond_dens, c(left, params_for_cond_dens)))
    left <- left - movement_num
    if(left <= param_bound[1]){
      left = param_bound[1]
      break
    }
    # }
    # print('in 1')
  }
  
  while(yval < do.call(cond_dens, c(right, params_for_cond_dens))){
    # print(cond_dens)
    # print(right)
    # print(do.call(cond_dens, c(right, params_for_cond_dens)))
    # print(yval)
    
    right <- right + movement_num
    # print(right)
    if(right >= param_bound[2]){
      right = param_bound[2]
      print('reached right bound')
      break
    }
    # print('in 2')
  }
  
  ## Now take a sample for next x from this width
  new_val <- runif(1,left,right)
  # print(new_val)
  
  # Keep if logpost > newval
  while(do.call(cond_dens, c(new_val,params_for_cond_dens)) < yval){
    if(new_val > start_val){
      right <- new_val
    }else{
      left <- new_val
    }
    new_val <- runif(1,left,right)
    # print('in 3')
  }
  
  start_val <- new_val
  return(start_val)
  # pis <- c(pis,start_val)
}
