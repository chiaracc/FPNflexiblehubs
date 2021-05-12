# Chiara Caldinelli & Rhodri Cusack
# Flexible Hubs Project
# caldinec@tcd.ie

import csv
import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import os.path
from os import path
import pandas as pd
import json
import glob
import math
import pickle
import ptitprince as pt
import seaborn as sns
from scipy.io import loadmat
from scipy.stats import pearsonr
from scipy import stats
import statsmodels.api as sm
from statsmodels.api import OLS
from statsmodels.formula.api import ols
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import statistics
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import statsmodels
from statsmodels.stats.multicomp import (pairwise_tukeyhsd,
                                         MultiComparison,
                                         ) 
from statsmodels.stats.anova import anova_lm

class poolsubjs:
    datapth=None
    atlaspth=None
    mainpth=None
    resultspth=None
    network_dict=None
    subjlist=None
    drop_first_column=None
    has_header_row=None

    def __init__(self,name,datapth=None,atlaspth=None,nsess=None,nroi=None,resultspth=None,network_dict=None,subjlist=None,filenametemplate=None,fileseparator='\t',drop_first_column=None,has_header_row=None):
        self.name=poolsubjs
        self.mainpath=resultspth
        self.datapth=datapth
        self.atlaspth=atlaspth
        self.subjlist=subjlist
        self.nroi=nroi
        self.resultspth=resultspth
        self.nsess=nsess
        self.network_dict=network_dict
        self.filenametemplate=filenametemplate
        self.fileseparator=fileseparator
        self.drop_first_column=drop_first_column
        self.has_header_row=has_header_row

    ############################################################################
    ####################### 3 Assign ROIs to networks ##########################
    ############################################################################

    def stage3_analysis(self,chunklen,corr):
        # Import parcellation atlas
        df_index = pd.read_csv(self.atlaspth) 
        # Remove colour code
        df_index.drop(df_index.columns[3], axis=1, inplace=True)
        # Create index  
        index = pd.MultiIndex.from_frame(df_index)
        index.get_level_values(0)

        # Create output dataframe for BVC and GVC
        BVCdfout=pd.DataFrame(columns=['subject','net1','net2','MeanBVC'])
        GVCdfout=pd.DataFrame(columns=['subject','net1','net2','MeanGVC'])
        print(BVCdfout)
        print(GVCdfout)
        

        for line in open(self.subjlist).readlines():
            line = line[:-1]
            
            
            # Write in a text file inside the results' folder what subject and what stage is starting
            now=datetime.datetime.now()
            file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
            file.write("\n Starting stage 3 ** Assign ROIs to networks ** for subject %s, chunklen %d"%(line,chunklen) + " at " + str(now))
            file.close()
            print ("Starting stage 3 ** Assign ROIs to networks ** for subject %s, chunklen %d"%(line,chunklen) + " at " + str(now))
            # Load results from the previous stage
            if corr==1:
                res=np.load(os.path.join(self.datapth,'stdevresults_corr_%d_%s'%(chunklen,line) + '.npy'), allow_pickle=True)    
                print(res)
            else:
                res=np.load(os.path.join(self.datapth,'stdevresults_reg_%d_%s'%(chunklen,line) + '.npy'), allow_pickle=True)
            # Create a pandas dataframe from the results from the previous stage
            dfres = pd.DataFrame(res)
            print(dfres)
            # Calculate BVC
            for net1, roilist1 in self.network_dict.items():
                for net2, roilist2 in self.network_dict.items():
                    if not net1==net2:
                        BVCvallist=[]
                        for roi1 in roilist1:
                            for roi2 in roilist2:
#                                print(roi1,roi2)
                                if not dfres.iloc[(roi1-1),(roi2-1)] is None:
                                    BVCvallist.append(dfres.iloc[(roi1-1),(roi2-1)])
                        # Calculate the average standard deviation and append it to the results' list
                        BVCdictout={'subject':line,'net1':net1,'net2':net2,'MeanBVC':np.nanmean(BVCvallist)}
                        BVCdfout=BVCdfout.append(BVCdictout,ignore_index=True)


            # Calculate GVC
            for net1, roilist1 in self.network_dict.items():
                for net2, roilist2 in self.network_dict.items():
                    GVCvallist=[]
                    for roi1 in roilist1:
                        for roi2 in roilist2:
                            # If ROI does not exist skip this stage
                            if not dfres.iloc[(roi1-1),(roi2-1)] is None:
                                GVCvallist.append(dfres.iloc[(roi1-1),(roi2-1)])
                    # Calculate the average standard deviation and append it to the results' list
                    GVCdictout={'subject':line,'net1':net1,'net2':net2,'MeanGVC':np.nanmean(GVCvallist)}
                    GVCdfout=GVCdfout.append(GVCdictout,ignore_index=True)
                    # Write in a text file inside the results' folder what subject and what stage is done
                now=datetime.datetime.now()
                file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
                file.write("\n Finished stage 3 ** Assign ROIs to networks ** for subject %s, chunklen %d"%(line,chunklen) + " at " + str(now))
                file.close()
        print(BVCdfout)
        print(GVCdfout)
        if corr==1:
            BVCdfout.to_csv(os.path.join(self.resultspth,'results-assignBVCcorr%d'%(chunklen) + '.csv'))
            GVCdfout.to_csv(os.path.join(self.resultspth,'results-assignGVCcorr%d'%(chunklen) + '.csv')) 
        else:
            BVCdfout.to_csv(os.path.join(self.resultspth,'results-assignBVCreg%d'%(chunklen) + '.csv'))
            GVCdfout.to_csv(os.path.join(self.resultspth,'results-assignGVCreg%d'%(chunklen) + '.csv'))


    ############################################################################
    ######################### 4 Calculate statistics ###########################
    ############################################################################			       

    def stage4_analysis(self,chunklen,corr):
        # Write in a text file inside the results' folder what stage is starting
        now=datetime.datetime.now()
        file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
        file.write("\n Starting stage 4 **Calculate statistics** at " + str(now))
        file.close()
        print ("Starting stage 4 **Calculate statistics** at " + str(now))

        # Import .csv file with BVC results
        if corr==1:
            df = pd.read_csv(os.path.join(self.resultspth,'results-assignBVCcorr%d'%(chunklen) + '.csv'))
            df.drop(df.columns[0],axis=1,inplace=True)
            print(df)  
        else:
            df = pd.read_csv(os.path.join(self.resultspth,'results-assignBVCreg%d'%(chunklen) + '.csv'))
            df.drop(df.columns[0],axis=1,inplace=True)
            print(df)
        
        # Calculate difference between each subject's mean and grand mean
        dfsubjmean_against_grandmean = df.groupby('subject').transform('mean') - df.mean()
        df2 = df - dfsubjmean_against_grandmean
        df['MeanBVC']=df2['MeanBVC']
        
        # Calculate average for each network
        dfmean=df.groupby(['subject','net2']).mean()
        dfmean = dfmean.reset_index(drop=False)
        print(dfmean)

        # Save the averaged values as .csv
        if corr==1:
            dfmean.to_csv(os.path.join(self.resultspth,'finalresultsBVCcorr%d'%(chunklen) + '.csv'))  
        else:
            dfmean.to_csv(os.path.join(self.resultspth,'finalresultsBVCreg%d'%(chunklen) + '.csv'))

        ###### ANOVA BVC ######
        results = ols('MeanBVC ~ C(net2)', data=dfmean).fit()
        table = sm.stats.anova_lm(results)
        print ("Anova results for BVC")
        print(table)
        if corr==1:
            np.savetxt((os.path.join(self.resultspth,'AnovaBVCcorr%d'%(chunklen) + '.txt')), table, fmt='%s')
        else:
            np.savetxt((os.path.join(self.resultspth,'AnovaBVCreg%d'%(chunklen) + '.txt')), table, fmt='%s')

        ###### TTEST BVC ######
        mod = MultiComparison(dfmean['MeanBVC'], dfmean['net2'])
        comp = mod.allpairtest(stats.ttest_rel, method='Holm')
        print ("T-test results for BVC")
        print(comp[0])
        # Save results
        if corr==1:
            np.savetxt((os.path.join(self.resultspth,'ResultsMultipleComparisonBVCcorr%d'%(chunklen) + '.txt')), comp, fmt='%s')
        else:
            np.savetxt((os.path.join(self.resultspth,'ResultsMultipleComparisonBVCreg%d'%(chunklen) + '.txt')), comp, fmt='%s')
        
        ###### GRAPHS BVC ######

        
        # Create a custom color palette for graphs
        colors = ["#F7EA23", "#B51DA3", "#000000", "#1DB526", "#5692BF", "#F31111", "#31A8F1", "#F131DF", "#1B56A6", "#950707"]
        # Set the custom color palette
        custom_p = sns.set_palette(sns.color_palette(colors))

        # Violin plot
        plt.figure()
        ax = sns.violinplot(x="net2", y="MeanBVC", data=dfmean, ci=68, palette=custom_p, order=["FPN", "CON", "SAN", "DAN", "VAN", "DMN", "Motor", "Aud.", "Vis.", "Subc."])
        ax.set_ylabel("Mean BVC", fontsize=16)
        plt.title('BVC')
        plt.xlabel('')
        plt.ylim(0, 0.7)
        ax.set_ylabel("Mean BVC", fontsize=16)
        plt.grid(False)
        if corr==1:
            outFileViolin = "ViolinPlotBVCcorr.png"
            plt.savefig((os.path.join(self.resultspth,'ViolinPlotBVCcorr%d'%(chunklen) + '.png')), dpi=200)
        else:
            outFileViolin = "ViolinPlotBVC.png"
            plt.savefig((os.path.join(self.resultspth,'ViolinPlotBVCreg%d'%(chunklen) + '.png')), dpi=200)
        print(("Figure saved as {0}".format(outFileViolin)))

        plt.figure()

        # Bar plot
        
        if corr==1:
            ax = sns.barplot(x="net2", y="MeanBVC", data=dfmean, ci=68, palette=custom_p, capsize=.2, order=["FPN", "CON", "SAN", "DAN", "VAN", "DMN", "Motor", "Aud.", "Vis.", "Subc."])
            plt.xlabel('')
            ax.set_ylabel("Mean BVC", fontsize=16)
            plt.title('BVC')
            plt.ylim(0, 0.7)
            plt.grid(False)
            outFile = "BarPlotBVCcorr.png"
            plt.savefig((os.path.join(self.resultspth,'BarPlotBVCcorr%d'%(chunklen) + '.png')), dpi=200)
        else:
            ax = sns.barplot(x="net2", y="MeanBVC", data=dfmean, ci=68, palette=custom_p, capsize=.2, order=["FPN", "CON", "SAN", "DAN", "VAN", "DMN", "Motor", "Aud.", "Vis.", "Subc."])
            plt.xlabel('')
            ax.set_ylabel("Mean BVC", fontsize=16)
            plt.title('BVC')
            plt.ylim(0, 0.7)
            plt.grid(False)
            outFile = "BarPlotBVC.png"
            plt.savefig((os.path.join(self.resultspth,'BarPlotBVCreg%d'%(chunklen) + '.png')), dpi=200)
        print(("Figure saved as {0}".format(outFile)))

        
        # Calculate average per network, across all subjects
        dfmeansubj=df.groupby(['net2']).mean()
        if corr==1:
            dfmeansubj.to_csv(os.path.join(self.resultspth,'finalresultpersubjBVCcorr%d'%(chunklen) + '.csv'))
        else:
            dfmeansubj.to_csv(os.path.join(self.resultspth,'finalresultpersubjBVCreg%d'%(chunklen) + '.csv'))
            
        # Calculate SD per network, across all subjects
        dfmeansubj=df.groupby(['net2']).std()
        if corr==1:
            dfmeansubj.to_csv(os.path.join(self.resultspth,'finalresultSDpersubjBVCcorr%d'%(chunklen) + '.csv'))
        else:
            dfmeansubj.to_csv(os.path.join(self.resultspth,'finalresultSDpersubjBVCreg%d'%(chunklen) + '.csv'))

############ GVC

        # Import .csv file with GVC results
        if corr==1:
            df = pd.read_csv(os.path.join(self.resultspth,'results-assignGVCcorr%d'%(chunklen) + '.csv'))
            df.drop(df.columns[0],axis=1,inplace=True)
            print(df)
        else:
            df = pd.read_csv(os.path.join(self.resultspth,'results-assignGVCreg%d'%(chunklen) + '.csv'))
            df.drop(df.columns[0],axis=1,inplace=True)
            print(df)
        
        
        # Calculate difference between each subject's mean and grand mean
        dfsubjmean_against_grandmean = df.groupby('subject').transform('mean') - df.mean()
        df2 = df - dfsubjmean_against_grandmean
        df['MeanGVC']=df2['MeanGVC']
        
        
        # Calculate average for each network
        dfmean=df.groupby(['subject','net2']).mean()
        dfmean = dfmean.reset_index(drop=False)
        print(dfmean)

        # Save the averaged values as .csv
        if corr==1:
            dfmean.to_csv(os.path.join(self.resultspth,'finalresultsGVCcorr%d'%(chunklen) + '.csv'))
        else:
            dfmean.to_csv(os.path.join(self.resultspth,'finalresultsGVCreg%d'%(chunklen) + '.csv'))

        ###### ANOVA GVC ######
        results = ols('MeanGVC ~ C(net2)', data=dfmean).fit()
        table = sm.stats.anova_lm(results)
        print ("Anova results for GVC")
        print(table)
        if corr==1:
            np.savetxt((os.path.join(self.resultspth,'AnovaGVCcorr%d'%(chunklen) + '.txt')), table, fmt='%s')
        else:
            np.savetxt((os.path.join(self.resultspth,'AnovaGVCreg%d'%(chunklen) + '.txt')), table, fmt='%s')

        ###### TTEST GVC ######
        mod = MultiComparison(dfmean['MeanGVC'], dfmean['net2'])
        comp = mod.allpairtest(stats.ttest_rel, method='Holm')
        print ("T-test results for GVC")
        print(comp[0])
        if corr==1:
            np.savetxt((os.path.join(self.resultspth,'ResultsMultipleComparisonGVCcorr%d'%(chunklen) + '.txt')), comp, fmt='%s')
        else:
            np.savetxt((os.path.join(self.resultspth,'ResultsMultipleComparisonGVCreg%d'%(chunklen) + '.txt')), comp, fmt='%s')

        
        ###### GRAPHS GVC ######
        # Violin and jitter
        ax = sns.violinplot(x="net2", y="MeanGVC", data=dfmean, ci=68, color=".8")
        ax = sns.stripplot(x="net2", y="MeanGVC", data=dfmean, edgecolor="white", size=2, jitter=1)
        plt.ylim(0, 0.65)
        plt.xlabel('')
        ax.set_ylabel("Mean GVC", fontsize=16)
        plt.title('GVC')
        if corr==1:
            outFileRainViolin = "RainViolinPlotGVC.png"
            plt.savefig((os.path.join(self.resultspth,'RainViolinPlotGVCcorr%d'%(chunklen) + '.png')), dpi=200)
        else:
            outFileRainViolin = "RainViolinPlotGVC.png"
            plt.savefig((os.path.join(self.resultspth,'RainViolinPlotGVCreg%d'%(chunklen) + '.png')), dpi=200)
        print(("Figure saved as {0}".format(outFileRaincloud)))
        
        # Create a custom color palette for graphs
        colors = ["#F7EA23", "#B51DA3", "#000000", "#1DB526", "#5692BF", "#F31111", "#31A8F1", "#F131DF", "#1B56A6", "#950707"]
        # Set the custom color palette
        custom_p = sns.set_palette(sns.color_palette(colors))

        # Violin plot
        plt.figure()
        ax = sns.violinplot(x="net2", y="MeanGVC", data=dfmean, ci=68, palette=custom_p, order=["FPN", "CON", "SAN", "DAN", "VAN", "DMN", "Motor", "Aud.", "Vis.", "Subc."])
        ax.set_ylabel("Mean GVC", fontsize=16)
        plt.title('GVC')
        plt.xlabel('')
        plt.grid(False)
        plt.ylim(0, 0.6)
        if corr==1:
            outFileViolin = 'ViolinPlotGVCcorr%d'%(chunklen) + '.png'
            plt.savefig((os.path.join(self.resultspth,'ViolinPlotGVCcorr%d'%(chunklen) + '.png')), dpi=200)
        else:
            outFileViolin = 'ViolinPlotGVC%d'%(chunklen) + '.png'
            plt.savefig((os.path.join(self.resultspth,'ViolinPlotGVCreg%d'%(chunklen) + '.png')), dpi=200)
        print("Figure saved as {0}".format(outFileViolin))

        # Bar plot   
        if corr==1:
            ax = sns.barplot(x="net2", y="MeanGVC", data=dfmean, ci=68, palette=custom_p, capsize=.2, order=["FPN", "CON", "SAN", "DAN", "VAN", "DMN", "Motor", "Aud.", "Vis.", "Subc."])
            plt.xlabel('')
            ax.set_ylabel("Mean GVC", fontsize=16)
            plt.title('GVC')
            plt.ylim(0, 0.6)
            plt.grid(False)
            outFile = 'BarPlotGVCcorr%d'%(chunklen) + '.png'
            plt.savefig((os.path.join(self.resultspth,'BarPlotGVCcorr%d'%(chunklen) + '.png')), dpi=200)
        else:
            #ax = sns.barplot(x="net2", y="MeanGVC", data=dfmean, ci=68, palette=custom_p, capsize=.2, order=["FPN", "CON", "SAN", "DAN", "VAN", "DMN", "Motor", "Aud.", "Vis.", "Subc."])
            ax=sns.stripplot(x="net2", y="MeanGVC", data=dfmean, palette=custom_p)
            plt.xlabel('')
            ax.set_ylabel("Mean GVC", fontsize=16)
            plt.title('GVC')
            plt.ylim(0, 0.6)
            plt.grid(False)
            outFile = 'BarPlotGVC%d'%(chunklen) + '.png'
            plt.savefig((os.path.join(self.resultspth,'BarPlotGVCreg%d'%(chunklen) + '.png')), dpi=200)
            print(("Figure saved as {0}".format(outFile)))


        
        # Calculate average per network, across all subjects
        dfmeansubj=df.groupby(['net2']).mean()
        if corr==1:
            dfmeansubj.to_csv(os.path.join(self.resultspth,'finalresultpersubjGVCcorr%d'%(chunklen) + '.csv'))
        else:
            dfmeansubj.to_csv(os.path.join(self.resultspth,'finalresultpersubjGVCreg%d'%(chunklen) + '.csv'))
        
        # Calculate SD per network, across all subjects
        dfmeansubj=df.groupby(['net2']).std()
        if corr==1:
            dfmeansubj.to_csv(os.path.join(self.resultspth,'finalresultSDpersubjGVCcorr%d'%(chunklen) + '.csv'))
        else:
            dfmeansubj.to_csv(os.path.join(self.resultspth,'finalresultSDpersubjGVCreg%d'%(chunklen) + '.csv'))

        # Write in a text file inside the results' folder what stage is done
        now=datetime.datetime.now()
        file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
        file.write("\n Finished stage 4 **Calculate statistics** at " + str(now))
        file.close()
