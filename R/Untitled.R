wrk_dir = '~/Desktop/GRAFit_test_Run/'
data_dir = '/Users/22111305/Desktop/PhD/HST/data/COSMOS'


inCat = read.csv('~/Desktop/GRAFit_test_Run/Input.csv')

CATAID = inCat$D10CATAID
RA = inCat$RA
DEC = inCat$DEC

frame_im_name = GRAFitFrameFinder_v2(GRAFitlib = GRAFitlib,
                                     data_dir = data_dir,
                                     target_loc = c(RA, DEC))

GRAFitMaster(wrk_dir = wrk_dir, data_dir = data_dir, threadMode = 0, ncores = 1, nComp= 1, optimMode = 'optim', object_list = inCat)
