import argparse
import PoolSubjs
import flexiblehubs
import sys


# camcan data
#studyforrest=PoolSubjs.poolsubjs('studyforrest',
studyforrest=flexiblehubs.FlexibleHubs('studyforrest',
                                        #datapth='/camcan/cc700/mri/pipeline/release004/data_fMRI_Unsmooth_Craddock/data_fMRI_Unsmooth_Craddock/aamod_roi_extract_epi_00002/$subj/Movie/',
                                        datapth='/camcan/power_parc/',
                                        subjlist='/camcan/studyforrest/studyforrest/subjlistCC_mo.txt',
                                        filenametemplate='{subj}_shen_264ROI_movie.txt',
                                        fileseparator='\t',
                                        nsess=1,
                                        nroi=264,
                                        resultspth='/camcan/results_rep/new_res/',
                                        atlaspth='/camcan/studyforrest/studyforrest/power_networks.csv',
                                        drop_first_column=True,
                                        has_header_row=False,
                                        network_dict={
                                            'FPN': [174, 175, 176, 177, 178, 179, 180, 181, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202],
                                            'CON': [47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60],
                                            'SAN': [203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220],
                                            'VAN': [235, 236, 237, 238, 239, 240, 241, 242], 
                                            'DAN': [251, 252, 256, 257, 258, 259, 260, 261, 262, 263, 264],
                                            'DMN': [74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 137, 139, ],
                                            'Motor': [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46],
                                            'Aud.': [61, 62, 63, 64, 65, 66, 67, 68, 68, 70, 71, 72, 73],
                                            'Vis.': [143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 19, 170, 171, 172, 173],
                                            'Subc.': [222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234],
                                            'Whole Brain': [174, 175, 176, 177, 178, 179, 180, 181, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 235, 236, 237, 238, 239, 240, 241, 242, 251, 252, 256, 257, 258, 259, 260, 261, 262, 263, 264, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 137, 139, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 61, 62, 63, 64, 65, 66, 67, 68, 68, 70, 71, 72, 73, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 19, 170, 171, 172, 173, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234],
                                            }
                                        )




if len(sys.argv)<2:
   print('Provide chunklen as a parameter')
   raise
else:
    for i, arg in enumerate(sys.argv):
         print('Argument %d: %s'%(i,arg))
    chunklen=int(sys.argv[1])
    subj=(sys.argv[2])

# Run stage 1
#studyforrest.stage1_analysis(chunklen, subj, corr=0)

# Run stage 2
studyforrest.stage2_analysis(chunklen, subj, corr=0)


