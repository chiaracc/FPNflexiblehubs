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

#chunklenght=15
#networks=11

class FlexibleHubs:
    datapth=None
    atlaspth=None
    mainpth=None
    resultspth=None
    network_dict=None
    subjlist=None
    drop_first_column=None
    has_header_row=None

    def __init__(self,name,datapth=None,atlaspth=None,nsess=None,nroi=None,resultspth=None,network_dict=None,subjlist=None, filenametemplate=None, fileseparator='\t', drop_first_column=None, has_header_row=None):
        self.name=FlexibleHubs
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
        #self.nsubj=nsubj
        


        # Make results path if it doesn't already exist
        if not os.path.exists(resultspth):
            os.makedirs(resultspth)


    
    ##########################################################################################
    ################ 1 Build and run the correlation / regression model  #####################
    ##########################################################################################

    def stage1_analysis(self,chunklen,subj,corr):
        # Run a correlation model
        if corr==1:
            # If the file that it is being computed during this stage is already there do not start
            if not os.path.exists(os.path.join(self.datapth,"allresults_corr_%d_%s"%(chunklen,subj) + ".npy")):
                # Write what stage and what subjects are starting in a text file in the results folder
                now=datetime.datetime.now()
                file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
                file.write("\n Starting stage 1 ** Build and run the correlation model ** for subject %s with chunklen %d"%(subj,chunklen) + " at " + str(now))
                file.close()
                # Set the right directory
                os.chdir(self.datapth)
                # Create a numpy array for the results of the correlation
                allresultscorr=np.ndarray((self.nsess,self.nroi,self.nroi),dtype=np.object)        
                # For every run
                for r in range(self.nsess):
                    # Print what stage and subject is being computed
                    now=datetime.datetime.now()
                    print("Starting stage 1 ** Build and run the regression or correlation model ** for subject %s with chunklen %d"%(subj,chunklen) + " session %d"%(r) + " at " + str(now))                    
                    # Get filename
                    filename=self.filenametemplate.replace('{subj}',subj).replace('{sess0}',str(r+1)).replace('{sess1}',str(r+1))
                    # Read parcellated data in text file for Power's ROIs
                    if self.has_header_row:
                        df=pd.read_csv(os.path.join(self.datapth,filename), header=0, sep=self.fileseparator, index_col=None)
                    else:
                        df=pd.read_csv(os.path.join(self.datapth,filename), header=None, sep=self.fileseparator, index_col=None)
                    if self.drop_first_column:
                        df = df.drop(columns=0)
                    print(df)
                    # Get number total number of columns of time series
                    self.nrois=len(df.columns)
                    # Calculated how many chunks there are 
                    nchunks = int(np.floor(len(df)/chunklen))
                    print ("Chunk number is %d"%(nchunks))
                    # Get one chunk per time, keeping all ROIs
                    resultscorrchunk=np.ndarray((nchunks,self.nrois,self.nrois),dtype=np.object)
                    # This is collapsing coefficients of each chunk together to store them in the results' file
                    for chunkind in range(nchunks): 
                        data=np.array(df.iloc[chunkind*chunklen : (chunkind+1)*chunklen,:])
                        resultscorrchunk[chunkind][:][:] = np.corrcoef(np.transpose(data))

                    # Assigning the right names to ROIs
                    for roi2 in range(self.nrois):
                    # If the ROI is not missing
                        if not pd.isna(df.iloc[0,roi2]): 
                            for roi1 in range(self.nrois):
                                listofcoef = []
                                # Append each chunk's result to a list to be stored
                                for chunkind in range(nchunks):  
                                    cor = (resultscorrchunk[chunkind][roi1][roi2])  
                                    listofcoef.append(cor)
                                        
                                    #Look up name of roi in dataframe
                                    #e.g., '1','2','3'... and convert to zero-based integer index for our numpy array (e.g., 0,1,2)
                                    allresultscorr[r][int(df.columns[roi1])-1][int(df.columns[roi2])-1] = np.array(listofcoef)
                print(allresultscorr)
                np.save("allresults_corr_%d_%s"%(chunklen,subj) + ".npy",allresultscorr)

            # Write what stage and what subjects are starting in a text file in the results folder
            now=datetime.datetime.now()
            file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
            file.write("\n Finished stage 1 ** Build and run the correlation model ** for subject %s  with chunklen %d"%(subj,chunklen) + " at " + str(now))
            file.close()
        
        # Run a regression model
        else:
            # If the file that it is being computed during this stage is already there do not start this subject
            if not os.path.exists(os.path.join(self.datapth,"allresults_reg_%d_%s"%(chunklen,subj) + ".npy")):
                # Write what stage and what subjects are starting in a text file in the results folder
                now=datetime.datetime.now()
                file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
                file.write("\n Starting stage 1 ** Build and run the regression model ** for subject %s  with chunklen %d"%(subj,chunklen) + " at " + str(now))
                file.close()
                # Set the right directory
                os.chdir(self.datapth)
                now=datetime.datetime.now()
                file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
                file.write("\n Starting stage 1 ** Still running for subj** for subject %s with chunklen %d"%(subj,chunklen) + " at " + str(now))
                file.close()
                # Create a numpy array for the results of the regression
                allresults=np.ndarray((self.nsess,self.nroi,self.nroi),dtype=np.object)
                # For every run
                for r in range(self.nsess):
                    now=datetime.datetime.now()
                    file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
                    file.write("\n Starting stage 1 ** Still running for subj** for subject %s with chunklen %d"%(subj,chunklen) + " at " + str(now))
                    file.close()
                    print ("Starting stage 1 ** Build and run the regression model ** for subject %s with chunklen %d"%(subj,chunklen) + " session %d"%(r) + " at " + str(now))
                    # Get filename
                    filename=self.filenametemplate.replace('{subj}',subj).replace('{sess0}',str(r+1)).replace('{sess1}',str(r+1))
                    # Read parcellated data in text file for Power's ROIs
                    if self.has_header_row:
                        df=pd.read_csv(os.path.join(self.datapth,filename), header=0, sep=self.fileseparator, index_col=None)
                    else:
                        df=pd.read_csv(os.path.join(self.datapth,filename), header=None, sep=self.fileseparator, index_col=None)
                    # Drop the first column as it is empty for timeseries from AFNI
                    if self.drop_first_column:
                        df = df.drop(columns=0)
                    print(df)
                    # Get number total number of columns of time series
                    self.nrois=len(df.columns)
                    # Calculated how many chunks
                    nchunks = int(np.floor(len(df)/chunklen))
                    print ("Chunk number is %d"%(nchunks))
                    for roi2 in range(self.nrois):
                        # If the ROI does not exist skip it
                        if not pd.isna(df.iloc[0,roi2]):
                            matrix = np.zeros((nchunks*chunklen,nchunks*2+1))
                            roi2index=int(roi2)
                            # For every chunk
                            for chunkind in range(nchunks):
                                # Create a matrix for the regression model
                                matrix[chunkind*chunklen : (chunkind+1)*chunklen, chunkind+1] = 1
                                matrix[chunkind*chunklen : (chunkind+1)*chunklen, nchunks+1+chunkind] = df.iloc[chunkind*chunklen : (chunkind+1)*chunklen,roi2index]
                                print(matrix[chunkind*chunklen : (chunkind+1)*chunklen, nchunks+1+chunkind])
                                # Build a regression model
                                print ("Starting run %d for ROI2 %d"%(r,roi2))
                                for roi1 in range(self.nrois):
                                    print(df.iloc[:,roi1])
                                    # Build a regression model  
                                    result=sm.OLS(df.iloc[:,roi1][:np.size(matrix,0)],matrix).fit()
                                    #Look up name of roi in dataframe
                                    #e.g., '1','2','3'... and convert to zero-based integer index for our numpy array (e.g., 0,1,2)
                                    allresults[r][int(df.columns[roi1])-1][int(df.columns[roi2])-1]=np.array(result.params)
                                    
                # Save results in the subject's folder
                print('this is allresults:')
                print(allresults)
                np.save("allresults_reg_%d_%s"%(chunklen,subj) + ".npy",allresults)

            # Write what stage and what subjects are done in a text file in the results folder
            now=datetime.datetime.now()
            file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
            file.write("\n Finished stage 1 ** Build and run the regression model ** for subject %s with chunklen %d"%(subj,chunklen) + " at " + str(now))
            file.close()

    ############################################################################
    ######################### 2 Summarize results  #############################
    ############################################################################

    def stage2_analysis(self,chunklen,subj,corr):
        # Summarize results from the correlation model
        if corr == 1:
            # If the file that it is being computed during this stage is already there do not start this subject
            if not os.path.exists(os.path.join(self.datapth,'stdevresults_corr_%d_%s'%(chunklen,subj) + '.npy')):

                # Write in a text file inside the results' folder what subject and what stage is starting
                now=datetime.datetime.now()
                file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
                file.write("\n Starting stage 2 ** Summarize results from the correlation ** for subject %s"%(subj) + " at " + str(now))
                file.close()
                print ("Starting stage 2 ** Summarize results from the correlation ** for subject %s"%(subj) + " at " + str(now))
                # Load the results from the previous stage
                res=np.load(os.path.join(self.datapth,"allresults_corr_%d_%s"%(chunklen,subj) + ".npy"), allow_pickle=True)
                print(res)
                # Create a numpy array for storing results of this stage
                stdevresultscorr=np.ndarray((self.nroi,self.nroi),dtype=np.object) 

                for roi1 in range(self.nroi):
                    for roi2 in range(self.nroi):
                        betaslist = []
                        for ses in range(self.nsess):
                            beta=(res[ses][roi1][roi2])
                            # If the current ROI is missing skip this stage
                            if beta is None:
                                file = open((os.path.join(self.datapth,'emptyROIs%s.txt'%(subj))),'a') 
                                file.write('Empty at %d, %d, %d \n' %(ses,roi1,roi2))
                                file.close()
                            else:
                                newlist=[i for i in beta if abs(1-i) > (10e-6)]
                                #  All betas across all sessions
                                betaslist.extend(newlist) 
                        if betaslist:   
                            print('betaslist')
                            print(betaslist)
                            # calculate the st dev across all betas and copy it into the results' matrix
                            stdevresultscorr[roi1,roi2] = np.std(betaslist) 
                        else:
                            print('Beta list empty for roi pair %d and %d for part %s'%(roi1,roi2,subj))
                # Save the results into the subject's folder
                print('this is stdevresultscorr:')
                print(stdevresultscorr)
                np.save(os.path.join(self.datapth,'stdevresults_corr_%d_%s'%(chunklen,subj) + '.npy'),stdevresultscorr)    
        
        # Summarize results from the regression model
        else:
            # If the file that it is being computed during this stage is already there do not start this subject
            if not os.path.exists(os.path.join(self.datapth,'stdevresults_reg_%d_%s'%(chunklen,subj) + '.npy')):
                                  
                # Write in a text file inside the results' folder what subject and what stage is starting
                now=datetime.datetime.now()
                file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
                file.write("\n Starting stage 2 ** Summarize results from the correlation ** for subject %s"%(subj) + " at " + str(now))
                file.close()
                print ("Starting stage 2 ** Summarize results from the correlation ** for subject %s"%(subj) + " at " + str(now))
                # Load the results from the previous stage
                res=np.load(os.path.join(self.datapth,"NEWallresults_reg_%d_%s"%(chunklen,subj) + ".npy"), allow_pickle=True)
                print(res)
                # Create a numpy array for storing results of this stage
                stdevresults=np.ndarray((self.nroi,self.nroi),dtype=np.object)

                for roi1 in range(self.nroi):
                    for roi2 in range(self.nroi):
                        betaslist = []
                        for ses in range(self.nsess):
                            beta=(res[ses][roi1][roi2])
                            # If ROI does not exist skip the following
                            if beta is None:
                                file = open((os.path.join(self.datapth,'emptyROIs%s.txt'%(subj))),'a') 
                                file.write('Empty at %d, %d, %d \n' %(ses,roi1,roi2))
                                file.close()
                            else:
                                # See how long the list of results is (how many betas we stored from the previous stage)
                                nbeta=len(beta)
                                # Select the right betas in the vector beta and append them to the list of betas across all sessions 
                                # Get only the betas we need - the last ones
                                betaPairwise=(int((nbeta-1)/2))
                                newlist=beta[betaPairwise+1:]
                                newlist=[i for i in beta if abs(1-i) > (10e-6)]
                                betaslist.extend(newlist)
                                
                        if betaslist:       
                            # Calculate standard deviation across all betas and copy it into the results' matrix
                            stdevresults[roi1,roi2] = np.std(betaslist) 
                        else:
                            print('Beta list empty for roi pair %d and %d for part %s'%(roi1,roi2,subj))
                # Save the results into the subject's folder

                np.save(os.path.join(self.datapth,'stdevresults_reg_%d_%s'%(chunklen,subj) + '.npy'),stdevresults)
                # Write in a text file inside the results' folder what subject and what stage is done
                now=datetime.datetime.now()
                file = open((os.path.join(self.resultspth,'Activity.txt')),'a') 
                file.write("\n Finished stage 2 ** Summarize results ** for subject %s"%(subj) + " at " + str(now))
                file.close()
