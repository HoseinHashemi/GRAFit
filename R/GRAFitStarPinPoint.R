# Author: Hosein Hashemi, ICRAR, Dec 2017.
GRAFitStarPinPoint <- function(image = NULL, InitSeg = NULL, blur_sigma = NULL, segPlot = segPlot, pixscale = pix_scale, 
                               magzero = magzero, header = header, tolerance = 3, sigma = 2, smooth = smooth , pixcut = pixcut, skycut = 1, ... ) {
  
  # source(paste(GRAFitlib,'/GRAFitMainFinder.R',sep=''))
  
  seg_diff = profoundProFound(profoundImDiff( image, sigma = blur_sigma ), plot = FALSE, ... )
  
  # diff_main = GRAFitMainFinder(src_list =  seg_diff$segstats, imDim = dim(image), main_pin = FALSE)
  # mergeIDs = seg_diff$segstats[,"segID"][seg_diff$segstats[,"segID"] != diff_main$segID]
  
   diff_seg_dilate = profoundMakeSegimDilate( image, seg_diff$segim, size = 3,  
                                              plot = FALSE, ... )
  
  finalSeg = profoundSegimMerge(image, InitSeg, diff_seg_dilate$segim)  
  # finalSegDilate = profoundSegimMerge(image, mainSegDilate, diff_seg_dilate$segim)
  
  if (segPlot) profoundSegimPlot(image, finalSeg)
  
  output = return(list(finalSeg = finalSeg))
}

# END
  