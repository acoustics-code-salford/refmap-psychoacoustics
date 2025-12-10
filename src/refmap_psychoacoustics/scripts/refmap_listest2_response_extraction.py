# -*- coding: utf-8 -*-

# script

# --------------
# %% Description
# --------------
"""
This script extracts response data from the individual test data files,
compiles and saves to CSV or XLSX files for further analysis. 
"""
# --------
# %% Setup
# --------

# import statements
import sys
import os
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QFileDialog, QApplication
from warnings import simplefilter

# suppress pandas performance warnings
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

# enable copy-on-write mode for Pandas (will be default from Pandas 3.0)
pd.options.mode.copy_on_write = True

if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance()

fileExt = "*.csv"
filelist = list(QFileDialog.getOpenFileNames(caption="Open individual response data files in '03 Experiment\Experiment 2\Test_files\Response_data\Individual'",
                                             filter=fileExt))[0]

# initialize empty dataframe
df = pd.DataFrame()
# loop through selected files
for file in filelist:
    print(f"Reading file: {os.path.basename(file)}")
    temp_df = pd.read_csv(file)
    df = pd.concat([df, temp_df], ignore_index=True)
# reset index
df.reset_index(drop=True, inplace=True)


# fill missing values
df['sourceStart'] = df['sourceStart'].fillna("Baseline")
df['sourceEvents'] = df['sourceEvents'].fillna(0)
df['sourceInterval'] = df['sourceInterval'].fillna(0)

# make sourceEvents and sourceInterval integer type
df['sourceEvents'] = df['sourceEvents'].astype(int)
df['sourceInterval'] = df['sourceInterval'].astype(int)

# create source fields by extracting the relevant substring from
# the source column, or assigning "Baseline" for baseline conditions
df['sourceType'] = df['source'].str.extract(r'(H520|T150)', expand=False)
df['sourceProximity'] = df['source'].str.extract(r'(Near|Far)', expand=False)
df['sourceMode'] = df['source'].str.extract(r'(Overflight|Delivery)', expand=False)
df.loc[df['sourceStart'] == "Baseline", 'sourceType'] = "Baseline"
df.loc[df['sourceStart'] == "Baseline", 'sourceProximity'] = "Baseline"
df.loc[df['sourceStart'] == "Baseline", 'sourceMode'] = "Baseline"

# simplify ambient condition names by removing everything after the first underscore
df['ambient'] = df['ambient'].str.split('_').str[0]

# create stimulus field
df['stimulus'] = (df['ambient'].str.split('_').str[0] + "_" +
                  df['source'].str.split("trajectory").str[1].str.replace("FadeCal", "") + "_" +
                  df['sourceStart'].astype(str) + "_" +
                  df['sourceEvents'].astype(str))
df.loc[df['sourceStart'] == "Baseline", 'stimulus'] = (df['ambient'].str.split('_').str[0] + "_baseline")

# drop unnecessary columns
df.drop(columns=['video', 'source'], inplace=True)

# rearrange the columns
cols = df.columns.tolist()
cols.insert(cols.index('trial') + 1, cols.pop(cols.index('stimulus')))
cols.insert(cols.index('ambient') + 1, cols.pop(cols.index('sourceProximity')))
cols.insert(cols.index('ambient') + 1, cols.pop(cols.index('sourceMode')))
cols.insert(cols.index('ambient') + 1, cols.pop(cols.index('sourceType')))
df = df.reindex(columns=cols)

# # assign a stimulus ID using unique stimulus identifiers
# df['stimID'] = df['stimulus'].astype('category').cat.codes + 1
# # rearrange the columns
# cols = df.columns.tolist()
# cols.insert(cols.index('participant') + 1, cols.pop(cols.index('stimulus')))
# cols.insert(cols.index('stimulus') + 1, cols.pop(cols.index('stimID')))
# df = df.reindex(columns=cols)

# create separate DataFrame for each response
# end responses
df_endAnnoy = df[df['response'] == "EndAnnoyance"]
df_endAnnoy = df_endAnnoy.loc[:, ~df_endAnnoy.columns.str.contains('^Unnamed')]
df_endAnnoy.rename(columns={"rating" : "annoyance"}, inplace=True)
df_endAnnoy.drop(columns=['response'], inplace=True)

df_endCircumplex = df[df['response'] == "EndCircumplex"]
df_endCircumplex = df_endCircumplex.loc[:, ~((df_endCircumplex.columns.str.contains('^Unnamed')) & (df_endCircumplex.isnull().all()))]
df_endCircumplex.rename(columns={"rating" : "pleasantness"}, inplace=True)
df_endCircumplex.rename(columns={df_endCircumplex.columns[df_endCircumplex.columns.str.contains('^Unnamed')][0] : "eventfulness"}, inplace=True)
df_endCircumplex.drop(columns=['response'], inplace=True)
# merge DataFrames
df_endResponse = pd.merge(df_endAnnoy, df_endCircumplex, how='inner')

# extract MomentAnnoyance data
df_momentAnnoy = df[df['response'] == "MomentAnnoyance"]
df_momentAnnoy.drop(columns=['response'], inplace=True)
df_momentAnnoy.reset_index(drop=True, inplace=True)
# rename rating column and all following unnamed columns as time "t_" with a
# number series indicating the milliseconds in x100 ms intervals
time_cols = {}
for ii, col in enumerate(df_momentAnnoy.columns[df_momentAnnoy.columns.get_loc('rating'):]):
    if "Unnamed" in col:
        time_cols[col] = f"t_{(ii)}"
df_momentAnnoy.rename(columns={"rating" : "t_0"}, inplace=True)
df_momentAnnoy.rename(columns=time_cols, inplace=True)

# save compiled data to file

# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

outFilePath = QFileDialog.getExistingDirectory(caption="Choose output folder to save processed files in  in '03 Experiment\Experiment 2\Test_files\Response_data\Compiled'")

# save end response data
endOutFile = os.path.join(outFilePath, "refmap_listest2_endResponses.csv")
df_endResponse.to_csv(endOutFile, index=False)

# save momentary annoyance data
momentOutFile = os.path.join(outFilePath, "refmap_listest2_momentaryAnnoyance.csv")
df_momentAnnoy.to_csv(momentOutFile, index=False)


# questionnaire data

if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance()

fileExt = "*.xlsx"
filepath = list(QFileDialog.getOpenFileName(caption="Open questionnaire response data file in '03 Experiment\Experiment 2\Test_files\Questionnaire'",
                                            filter=fileExt))[0]

# read questionnaire data
df_questionnaire = pd.read_excel(filepath, index_col=0)
 # drop unnecessary columns
df_questionnaire.drop(columns=['Start time', 'Completion time', 'NoiseSensitive1', 'NoiseSensitive2',
                               'NoiseSensitive3', 'NoiseSensitive4',
                               'NoiseSensitive5', 'Comments'], inplace=True)
# fill na for AAMExperience (NA is "None"), for NativeLanguage (NA is "NotAnswered")
df_questionnaire['AAMExperience'] = df_questionnaire['AAMExperience'].fillna("None")
df_questionnaire['Nationality'] = df_questionnaire['Nationality'].fillna("NotAnswered")
df_questionnaire['NativeLanguage'] = df_questionnaire['NativeLanguage'].fillna("NotAnswered")
# remove substring following first space
df_questionnaire['AAMAttitude'] = df_questionnaire['AAMAttitude'].str.split(' ').str[0]
df_questionnaire['UKNational'] = df_questionnaire['UKNational'].str.replace("No answer", "NotAnswered")

# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

outFilePath = QFileDialog.getExistingDirectory(caption="Choose output folder to save processed files in  in '03 Experiment\Experiment 2\Test_files\Questionnaire'")

# save cleaned questionnaire data to file
questionOutFile = os.path.join(outFilePath, "refmap_listest2_questionnaireResponses.csv")

df_questionnaire.to_csv(questionOutFile, index=False)
