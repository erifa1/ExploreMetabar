# return list of modalities or taxas of a phyloseq object
get_cat_colors <- function(type, list_fact, phy_object){
  if(type == "taxa"){
    colors_list <- c(lapply(1:length(list_fact), FUN = function(i){
      data.frame(unique(tax_table(phy_object)[, list_fact[i]]), row.names = NULL)
    }))
    return(colors_list)
  }else if(type == "fact"){
    colors_list <- c(lapply(1:length(list_fact), FUN = function(i){
      data.frame(unique(sample_data(phy_object)[, list_fact[i]]), row.names = NULL)
    }))
    return(colors_list)
  }
}

# return a named list of association modality-color or taxa-color
get_modality_colors <- function(type, list_fact, modality, palettes){
  available_pal <- palettes
  colors <- list()
  for(k in 1:length(list_fact)){
    m <- modality[[k]]
    nb_modality <- dim(m)[1]
    missing_val <- NULL
    if(anyNA(m)){
      m <- m[!is.na(m), ]
      nb_modality <- length(m)
      missing_val <- c("#000000")
      names(missing_val) <- c("NA")
    }
    if(nb_modality < 52){
      pal <- available_pal[available_pal$length >= nb_modality, ]
      pal <- pal[which(pal$length == min(pal$length)), ]
      npal <- ifelse(type == "fact", 1, length(pal$package))
      available_pal <- available_pal[available_pal$package != pal$package[npal] | available_pal$palette != pal$palette[npal], ]
      df_colors <- c(paletteer::paletteer_d(palette = paste(pal$package[npal], pal$palette[npal], sep = "::"), n = nb_modality))
      df_colors <- df_colors[1:nb_modality]
    }else{
      random_colors <- gplots::col2hex(grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)])
      random_colors <- random_colors[random_colors != "#FFFFFF"]
      df_colors <- sample(random_colors, nb_modality)
    }
    names(df_colors) <- sapply(1:nb_modality, FUN  = function(i){data.frame(m)[i, ]})
    colors[[k]] <- c(df_colors, missing_val)
  }
  names(colors) <- list_fact
  return(colors)
}

# return a named list of colors associated to variables or taxas read from a xlsx file
get_modality_file <- function(path_file, type, list_fact, num_fact = FALSE, num_palettes = NULL, phy_obj){
  sheet_names <- readxl::excel_sheets(path_file)
  valid_names <- sheet_names %in% list_fact
  sheet_names <- sheet_names[valid_names]
  if(length(sheet_names) > 0){
    list_colors <- lapply(1:length(sheet_names), FUN = function(k){
      df_fact_colors <- as.data.frame(readxl::read_excel(path = path_file, sheet = sheet_names[k], col_names = FALSE))
      if(type == "fact"){
        if(num_fact[[sheet_names[k]]]){
          df_colors <- df_fact_colors[1, 1]
          names(df_colors) <- NULL
          df_colors <- df_colors[df_colors %in% num_palettes]
        }else{
          m <- get_cat_colors(type = type, list_fact = sheet_names[k], phy_object = phy_obj)[[1]]
          df_colors <- t(df_fact_colors)[2, ]
          valid_colors <- plotfunctions::isColor(df_colors)
          df_colors <- df_colors[valid_colors]
          names(df_colors) <- t(df_fact_colors)[1, valid_colors]
          df_colors <- df_colors[names(df_colors) %in% m[, 1]]
        }
      }else{
        m <- get_cat_colors(type = type, list_fact = sheet_names[k], phy_object = phy_obj)[[1]]
        df_colors <- t(df_fact_colors)[2, ]
        valid_colors <- plotfunctions::isColor(df_colors)
        df_colors <- df_colors[valid_colors]
        names(df_colors) <- t(df_fact_colors)[1, valid_colors]
        df_colors <- df_colors[names(df_colors) %in% m[, 1]]
      }
      return(df_colors)
    })
    names(list_colors) <- sheet_names
  }else{
    list_colors <- NULL
  }
  list_colors <- list_colors[!sapply(list_colors, FUN = is.null)]
  return(list_colors)
}

# return a named list with association modality-color or taxa-color randomly completed if there are missing colors in xlsx file
get_complete_color_file <- function(mod_file, missing_colors){
  random_colors <- gplots::col2hex(grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)])
  random_colors <- random_colors[random_colors != "#FFFFFF"]
  new_colors <- lapply(1:length(mod_file), FUN = function(i){
    missing <- missing_colors[[i]]
    if(missing[1] == "No color is missing." | length(missing) == 0){
      color <- mod_file[[i]]
    }else{
      df_colors <- sample(random_colors, length(missing))
      names(df_colors) <- missing
      color <- c(mod_file[[i]], df_colors)
    }
    return(color)
  })
  names(new_colors) <- names(mod_file)
  return(new_colors)
}
