#!/bin/bash

# DOES NOT DELETE ORIGINALS

# a = recursive + preserve special files and file metadata, v = verbose, h = human readable
# note with rsycn, the / at the end of the first path says copy all contents of the folder, not the folder itself
# I probably do not need the wildcard or the / at the end of the destination, but it works
# use -n or --dry-run to test

# see test runs at the end
	# MAIN LOOP:
	rootsrcDIR=/Users/coop558/mydata/icemodel/output/region
	rootdstDIR=/Volumes/GDRIVE/matlab/GREENLAND/runoff/icemodel/output/region

for i in {2011..2018}
do
	# enbal
	sourceDIR=$rootsrcDIR/$i/enbal
	destDIR=$rootdstDIR/$i/enbal
	rsync -avh $sourceDIR/* $destDIR/

	# ice 1
	sourceDIR=$rootsrcDIR/$i/ice1
	destDIR=$rootdstDIR/$i/ice1
	rsync -avh $sourceDIR/* $destDIR/

	# ice 2
	sourceDIR=$rootsrcDIR/$i/ice2
	destDIR=$rootdstDIR/$i/ice2
	rsync -avh $sourceDIR/* $destDIR/

	# opts
	sourceDIR=$rootsrcDIR/$i/opts
	destDIR=$rootdstDIR/$i/opts
	rsync -avh $sourceDIR/* $destDIR/

done
exit #
#
# # test run
# # sourceDIR=/Volumes/GDRIVE/matlab/GREENLAND/runoff/icemodel/model/experiment/v8/region/skinmodel/output/modis/t4/2009/enbal
# # destDIR=/Volumes/Samsung_T5/matlab/GREENLAND/runoff/icemodel/model/experiment/v8/region/skinmodel/output/modis/t4/2009/enbal
# # rsync -avh $sourceDIR/enbal_1.mat $destDIR/enbal_1.mat
#
# # test loop:
# # rootsrcDIR=/Volumes/GDRIVE/matlab/GREENLAND/runoff/icemodel/model/experiment/v8/region/skinmodel/output/modis/t4
# # rootdstDIR=/Volumes/Samsung_T5/matlab/GREENLAND/runoff/icemodel/model/experiment/v8/region/skinmodel/output/modis/t4
# # for i in {2009..2009}
# # do
# # 	# enbal
# # 	sourceDIR=$rootsrcDIR/$i/enbal
# # 	destDIR=$rootdstDIR/$i/enbal
# # 	rsync -n --remove-source-files -avh $sourceDIR/enbal_1.mat $destDIR/
# # done
# # exit #


# sourceDIR=/Users/coop558/mydata/icemodel/output/region
# destDIR=/Volumes/GDRIVE/matlab/GREENLAND/runoff/icemodel/output/region
#
# 	rsync -avh $sourceDIR/skinmodel_reference_mar_albedo_2009.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_mar_albedo_2010.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_mar_albedo_2011.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_mar_albedo_2012.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_mar_albedo_2013.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_mar_albedo_2014.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_mar_albedo_2015.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_mar_albedo_2016.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_mar_albedo_2017.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_mar_albedo_2018.mat $destDIR/
#
# 	rsync -avh $sourceDIR/icemodel_reference_mar_albedo_2009.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_mar_albedo_2010.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_mar_albedo_2011.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_mar_albedo_2012.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_mar_albedo_2013.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_mar_albedo_2014.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_mar_albedo_2015.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_mar_albedo_2016.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_mar_albedo_2017.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_mar_albedo_2018.mat $destDIR/
#
# 	rsync -avh $sourceDIR/skinmodel_reference_modis_albedo_2009.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_modis_albedo_2010.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_modis_albedo_2011.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_modis_albedo_2012.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_modis_albedo_2013.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_modis_albedo_2014.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_modis_albedo_2015.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_modis_albedo_2016.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_modis_albedo_2017.mat $destDIR/
# 	rsync -avh $sourceDIR/skinmodel_reference_modis_albedo_2018.mat $destDIR/
#
# 	rsync -avh $sourceDIR/icemodel_reference_modis_albedo_2009.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_modis_albedo_2010.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_modis_albedo_2011.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_modis_albedo_2012.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_modis_albedo_2013.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_modis_albedo_2014.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_modis_albedo_2015.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_modis_albedo_2016.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_modis_albedo_2017.mat $destDIR/
# 	rsync -avh $sourceDIR/icemodel_reference_modis_albedo_2018.mat $destDIR/