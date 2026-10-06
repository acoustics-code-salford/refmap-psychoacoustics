# -*- coding: utf-8 -*-

# script


# --------
# %% Setup
# --------

# import statements
import sys
import os
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QFileDialog, QApplication

# suppress pandas performance warnings
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

# enable copy-on-write mode for Pandas (will be default from Pandas 3.0)
pd.options.mode.copy_on_write = True

if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance()

# output path for saving data
outFilePath1 = QFileDialog.getExistingDirectory(caption="Choose output folder to save processed files in '03 Experiment\Experiment 1\Analysis\PostProcess'")
outFilePath2 = QFileDialog.getExistingDirectory(caption="Choose output folder to save processed files in '03 Experiment\Experiment 2\Analysis\PostProcess'")

# Load the data
dataBySubjTest1filepath = list(QFileDialog.getOpenFileName(caption="Open refmap_listest1_testdata_BySubj.csv file in '03 Experiment\Experiment 1\Analysis\PostProcess'",
                                            filter="refmap_listest1_testdata_BySubj.csv"))[0]
dataBySubjTest1 = pd.read_csv(dataBySubjTest1filepath, index_col=False)

dataByStimTest1filepath = list(QFileDialog.getOpenFileName(caption="Open refmap_listest1_testdata_ByStim.csv file in '03 Experiment\Experiment 1\Analysis\PostProcess'",
                                            filter="refmap_listest1_testdata_ByStim.csv"))[0]
dataByStimTest1 = pd.read_csv(dataByStimTest1filepath, index_col=0)

dataBySubjTest2filepath = list(QFileDialog.getOpenFileName(caption="Open refmap_listest2_testdata_BySubj.csv file in '03 Experiment\Experiment 2\Analysis\PostProcess'",
                                            filter="refmap_listest2_testdata_BySubj.csv"))[0]
dataBySubjTest2 = pd.read_csv(dataBySubjTest2filepath, index_col=False)

dataByStimTest2filepath = list(QFileDialog.getOpenFileName(caption="Open refmap_listest2_testdata_ByStim.csv file in '03 Experiment\Experiment 2\Analysis\PostProcess'",
                                            filter="refmap_listest2_testdata_ByStim.csv"))[0]
dataByStimTest2 = pd.read_csv(dataByStimTest2filepath, index_col=0)


# %%---------------------------
# subject-level data processing
# -----------------------------

# test 1
# ------
dataBySubjTest1Merge = dataBySubjTest1.copy()

test1_drop_cols = ["SessionPart", "Valence", "Arousal", "UAS_noticed",
                   "dValence", "dArousal", "CALBINRecFiles", "CALHEQRecFiles",
                   "UASLAeq", "AmbientLAeq", "SNRlevel", "Hearing_impaired",
                   "PANAS_positive", "PANAS_negative", "ResponseFactorLoud",
                   "ResponseFactorChar", "ResponseFactorDuration", "ResponseFactorQuantRep",
                   "ResponseFactorProximity", "ResponseFactorAmb", "ResponseFactorExper",
                   "ResponseFactorSens", "ResponseFactorComfort", "ResponseFactorPrivSafe",
                   "PartLoudMGSTPowAvg", "PartLoudMGST05Ex", "PartLoudMGLTPowAvg",
                   "PartLoudMGLT05Ex", "LoudQZ5323PowAvgMaxLR", "LoudQZ4182PowAvgMaxLR", "UASLoudQZ5323PowAvgMaxLR",
                   "UASLoudQZ4182PowAvgMaxLR", "AmbLoudQZ5323PowAvgMaxLR", "AmbLoudQZ4182PowAvgMaxLR"]
dataBySubjTest1Merge.drop(columns=test1_drop_cols, inplace=True)

dataBySubjTest1Merge['LAeqdiff'] = dataBySubjTest1Merge['UASLAeqMaxLR'] - dataBySubjTest1Merge['AmbLAeqMaxLR']

test1_rename_cols = {"StimFile": "Stimulus", "HighAnnoy": "HighlyAnnoyed",
                     "dHighAnnoy": "dHighlyAnnoyed",
                     "Home_Area": "HomeArea", "Area_soundscape": "SoundscapeArea",
                     "NSSTotal": "NoiseSensitivity", "AAM_attitude": "AAMAttitude"}
dataBySubjTest1Merge.rename(columns=test1_rename_cols, inplace=True)

dataBySubjTest1Merge["UKNational"] = dataBySubjTest1Merge["NationGeo"].apply(lambda x: "UK" if x == "UK" else "Non-UK")


# strip ".wav" substring from all rows in Stimulus column
dataBySubjTest1Merge['Stimulus'] = dataBySubjTest1Merge['Stimulus'].str.replace('.wav', '', regex=False)

# 1. Define specific renaming column ranges (start, end)
column_ranges = [
    ('LoudECMAPowAvgBin', 'ImpulsLoudWECMAPowAvgBin'),
    ('UASLoudECMAPowAvgBin', 'UASImpulsLoudWECMAPowAvgBin'),
    ('AmbLoudECMAPowAvgBin', 'AmbImpulsLoudWECMAPowAvgBin'),
    ('PartLoudSHMPowAvgBin', 'PartTonShpvBSHM05ExMaxLR'),
    ('dTonalECMAAvgMaxLR', 'dImpulsLoudWECMAPowAvgBin')
]

rename_dict = {}
cols = dataBySubjTest1Merge.columns.tolist()

for start, end in column_ranges:
    try:
        start_idx = cols.index(start)
        end_idx = cols.index(end)
        
        # Extract the block of names (inclusive)
        block = cols[start_idx : end_idx + 1]
        
        # Add the modified versions to our mapping dictionary
        for col in block:
            rename_dict[col] = col.replace('MaxLR', '').replace('Bin', '')

    except ValueError as e:
        print(f"Error occurred while processing range {start} to {end}: {e}")

# Apply the dictionary directly to rename
dataBySubjTest1Merge.rename(columns=rename_dict, inplace=True)

dataBySubjTest1Merge.to_csv(os.path.join(outFilePath1,
                                         "refmap_listest1_testdata_forValid_BySubj.csv"),
                            index=False)

# test 2
# ------
dataBySubjTest2Merge = dataBySubjTest2.copy()

test2_drop_cols = ["Pleasantness", "Eventfulness", "ProbHA20k", "ProbHA10k",
                   "dPleasantness", "dEventfulness", "dProbHA20k", "dProbHA10k",
                   "HATSRecFiles", "MA220MicRecFiles", "AmbientRef",
                   "UASProximity", "UASStart", "Nationality", "NativeLanguage"]
dataBySubjTest2Merge.drop(columns=test2_drop_cols, inplace=True)

test2_rename_cols = {"Trial": "TrialNumber"}
dataBySubjTest2Merge.rename(columns=test2_rename_cols, inplace=True)

dataBySubjTest2Merge['PartTrialNumber'] = dataBySubjTest2Merge['TrialNumber']
dataBySubjTest2Merge['StimDuration'] = 30

# find columns in dataBySubjTest1Merge not in dataBySubjTest2Merge and vice versa
cols1 = dataBySubjTest1Merge.columns.tolist()
cols2 = dataBySubjTest2Merge.columns.tolist()
print("Not in dataBySubjTest1Merge:" + str([col for col in cols2 if col not in cols1]))
print("Not in dataBySubjTest2Merge:" + str([col for col in cols1 if col not in cols2]))

# find any duplicates in cols1 or cols2
print("Duplicates in dataBySubjTest1Merge:" + str([col for col in cols1 if cols1.count(col) > 1]))
print("Duplicates in dataBySubjTest2Merge:" + str([col for col in cols2 if cols2.count(col) > 1]))

# %%------------------------------------------
# aligning participant and stimulus ID numbers
# --------------------------------------------

df1 = dataBySubjTest1Merge.copy()
df2 = dataBySubjTest2Merge.copy()

# --- Participant-level link tables (one row per linked participant) ---
links1 = df1[['ID', 'Exp2ID']].dropna().drop_duplicates().astype(int)
links2 = df2[['ID', 'Exp1ID']].dropna().drop_duplicates().astype(int)

# Each participant links to at most one participant in the other experiment
assert links1['ID'].is_unique and links1['Exp2ID'].is_unique, "Inconsistent Exp2ID links in df1"
assert links2['ID'].is_unique and links2['Exp1ID'].is_unique, "Inconsistent Exp1ID links in df2"

# Links agree in both directions
assert set(zip(links1['ID'], links1['Exp2ID'])) == set(zip(links2['Exp1ID'], links2['ID'])), \
    "Exp1ID / Exp2ID links are not reciprocal between datasets"

# --- Build mappings from the ORIGINAL IDs ---
offset = int(max(df1['ID'].max(), df2['ID'].max()))   # new df2 IDs all exceed every df1 ID

df2_to_df1 = dict(zip(links2['ID'], links2['Exp1ID']))          # original df2 ID -> original df1 ID
df1_to_df2 = dict(zip(links1['ID'], links1['Exp2ID'] + offset))  # original df1 ID -> NEW df2 ID

# --- (ii) Person-level ID: the experiment 1 ID is canonical ---
df1['PersonNum'] = df1['ID']
df2['PersonNum'] = df2['ID'].map(df2_to_df1).fillna(df2['ID'] + offset)

# --- Cross-reference columns, filled on every row of a linked participant ---
df1['Exp2ID'] = df1['ID'].map(df1_to_df2)
df2['Exp1ID'] = df2['ID'].map(df2_to_df1)

# --- (i) Make experiment-specific IDs unique across datasets ---
df2['ID'] = df2['ID'] + offset

# Nullable ints so NaNs don't force floats
for df, cols in [(df1, ['PersonNum', 'Exp2ID']), (df2, ['PersonNum', 'Exp1ID'])]:
    df[cols] = df[cols].astype('Int64')

# --- Checks ---
# ID values don't overlap across datasets
assert set(df1['ID']).isdisjoint(df2['ID'])

# Within each dataset, ID <-> PersonNum is one-to-one
for df in (df1, df2):
    assert df.groupby('ID')['PersonNum'].nunique().eq(1).all()
    assert df['PersonNum'].nunique() == df['ID'].nunique()

# Participants in both experiments share a PersonNum; everyone else's is unique to them
shared = set(df1['PersonNum']) & set(df2['PersonNum'])
assert len(shared) == len(links1)

# df1's Exp2ID values match df2's new ID values
linked_ids = df1['Exp2ID'].dropna().unique()
assert set(linked_ids) <= set(df2['ID'])

# drop the Exp1ID and Exp2ID columns, since we now have PersonNum as the canonical ID
df1.drop(columns=['Exp2ID'], inplace=True)
df2.drop(columns=['Exp1ID'], inplace=True)

# --- StimNum: unique stimulus number across both experiments ---
assert df1['StimID'].notna().all() and df2['StimID'].notna().all(), "Missing StimID values"

stims1 = sorted(df1['StimID'].unique())
stims2 = sorted(df2['StimID'].unique())

stim_map1 = {s: i for i, s in enumerate(stims1, start=1)}
stim_map2 = {s: i for i, s in enumerate(stims2, start=len(stims1) + 1)}

df1['StimNum'] = df1['StimID'].map(stim_map1)
df2['StimNum'] = df2['StimID'].map(stim_map2)

# --- Checks ---
assert set(df1['StimNum']).isdisjoint(df2['StimNum'])               # no overlap across experiments
assert df1['StimNum'].nunique() == df1['StimID'].nunique()          # one-to-one within each dataset
assert df2['StimNum'].nunique() == df2['StimID'].nunique()
assert df1['StimNum'].notna().all() and df2['StimNum'].notna().all()

dataBySubjTestMerge = pd.concat([df1.assign(Experiment=1), df2.assign(Experiment=2)], ignore_index=True)

# move the Experiment column to the front and the PersonNum column to the second position
# move StimNum before StimID
cols = dataBySubjTestMerge.columns.tolist()
cols.insert(0, cols.pop(cols.index('Experiment')))
cols.insert(1, cols.pop(cols.index('PersonNum')))
cols.insert(cols.index('StimID'), cols.pop(cols.index('StimNum')))
dataBySubjTestMerge = dataBySubjTestMerge[cols]

dataBySubjTestMerge.rename(columns={'ID': 'ParticipantNum'}, inplace=True)

# save merged data to CSV
dataBySubjTestMerge.to_csv(os.path.join(outFilePath2,
                                        "refmap_listest1_2_testdataMerge_BySubj.csv"),
                           index=False)


# %%----------------------------
# stimulus-level data processing
# ------------------------------

# test 1
# ------
dataByStimTest1Merge = dataByStimTest1.copy()

test1_drop_cols = ["SessionPart", "CALBINRecFiles", "CALHEQRecFiles",
                   "UASLAeq", "AmbientLAeq", "SNRlevel",
                   "PartLoudMGSTPowAvg", "PartLoudMGST05Ex", "PartLoudMGLTPowAvg",
                   "PartLoudMGLT05Ex", "LoudQZ5323PowAvgMaxLR", "LoudQZ4182PowAvgMaxLR",
                   "UASLoudQZ5323PowAvgMaxLR",
                   "UASLoudQZ4182PowAvgMaxLR", "AmbLoudQZ5323PowAvgMaxLR", "AmbLoudQZ4182PowAvgMaxLR"]

test1_drop_cols += list(dataByStimTest1Merge.columns[dataByStimTest1Merge.columns.get_loc('Arousal_1'):
                                                dataByStimTest1Merge.columns.get_loc('dHighAnnoy_19') + 1])
test1_drop_cols += list(dataByStimTest1Merge.columns[dataByStimTest1Merge.columns.get_loc('UAS_noticed_1'):
                                                dataByStimTest1Merge.columns.get_loc('NoticedPropCI_High') + 1])
test1_drop_cols.extend(["ArousalMedian", "ArousalMedianCI_Low", "ArousalMedianCI_High", "ArousalMean",
                        "ArousalMeanCI_Low", "ArousalMeanCI_High", "ValenceMedian", "ValenceMedianCI_Low",
                        "ValenceMedianCI_High", "ValenceMean", "ValenceMeanCI_Low", "ValenceMeanCI_High",
                        "dValenceMedian", "dValenceMedianCI_Low", "dValenceMedianCI_High", "dValenceMean",
                        "dValenceMeanCI_Low", "dValenceMeanCI_High", "dArousalMedian", "dArousalMedianCI_Low",
                        "dArousalMedianCI_High", "dArousalMean", "dArousalMeanCI_Low", "dArousalMeanCI_High"])

dataByStimTest1Merge.drop(columns=test1_drop_cols, inplace=True)

dataByStimTest1Merge['LAeqdiff'] = dataByStimTest1Merge['UASLAeqMaxLR'] - dataByStimTest1Merge['AmbLAeqMaxLR']

test1_rename_cols = {"StimFile": "Stimulus", "AnnoyMedian": "AnnoyanceMedian",
                     "AnnoyMedianCI_Low": "AnnoyanceMedianCI_Low",
                     "AnnoyMedianCI_High": "AnnoyanceMedianCI_High", "AnnoyMean": "AnnoyanceMean",
                     "AnnoyMeanCI_Low": "AnnoyanceMeanCI_Low", "AnnoyMeanCI_High": "AnnoyanceMeanCI_High",
                     "dAnnoyMedian": "dAnnoyanceMedian", "dAnnoyMedianCI_Low": "dAnnoyanceMedianCI_Low",
                     "dAnnoyMedianCI_High": "dAnnoyanceMedianCI_High", "dAnnoyMean": "dAnnoyanceMean",
                     "dAnnoyMeanCI_High": "dAnnoyanceMeanCI_High", "dAnnoyMeanCI_Low": "dAnnoyanceMeanCI_Low",
                     "HighAnnoyTotal": "HighlyAnnoyedTotal", "HighAnnoyProp": "HighlyAnnoyedProp",
                     "HighAnnoyPropCI_Low": "HighlyAnnoyedPropCI_Low",
                     "HighAnnoyPropCI_High": "HighlyAnnoyedPropCI_High",
                     "dHighAnnoyTotal": "dHighlyAnnoyedTotal", "dHighAnnoyProp": "dHighlyAnnoyedProp",
                     "dHighAnnoyPropCI_Low": "dHighlyAnnoyedPropCI_Low",
                     "dHighAnnoyPropCI_High": "dHighlyAnnoyedPropCI_High"}

dataByStimTest1Merge.rename(columns=test1_rename_cols, inplace=True)
dataByStimTest1Merge.rename(columns=rename_dict, inplace=True)

# test 2
# ------
dataByStimTest2Merge = dataByStimTest2.copy()

test2_drop_cols = ["HATSRecFiles", "MA220MicRecFiles", "AmbientRef",
                   "UASProximity", "UASStart"]
test2_drop_cols += list(dataByStimTest2Merge.columns[dataByStimTest2Merge.columns.get_loc('Annoyance_1'):
                                                dataByStimTest2Merge.columns.get_loc('dProbHA20k_9') + 1])
test2_drop_cols.extend(["EventfulnessMedian", "EventfulnessMedianCI_Low", "EventfulnessMedianCI_High", "EventfulnessMean",
                        "EventfulnessMeanCI_Low", "EventfulnessMeanCI_High", "PleasantnessMedian", "PleasantnessMedianCI_Low",
                        "PleasantnessMedianCI_High", "PleasantnessMean", "PleasantnessMeanCI_Low", "PleasantnessMeanCI_High",
                        "dPleasantnessMedian", "dPleasantnessMedianCI_Low", "dPleasantnessMedianCI_High", "dPleasantnessMean",
                        "dPleasantnessMeanCI_Low", "dPleasantnessMeanCI_High", "dEventfulnessMedian", "dEventfulnessMedianCI_Low",
                        "dEventfulnessMedianCI_High", "dEventfulnessMean", "dEventfulnessMeanCI_Low", "dEventfulnessMeanCI_High",
                        "ProbHA20kMean", "ProbHA20kMeanCI_Low", "ProbHA20kMeanCI_High", "ProbHA10kMean", "ProbHA10kMeanCI_Low",
                        "ProbHA10kMeanCI_High", "dProbHA20kMean", "dProbHA20kMeanCI_Low", "dProbHA20kMeanCI_High", "dProbHA10kMean",
                        "dProbHA10kMeanCI_Low", "dProbHA10kMeanCI_High", "AnnoyanceQuartile25", "AnnoyanceQuartile75", "EventfulnessQuartile25", "EventfulnessQuartile75", "PleasantnessQuartile25", "PleasantnessQuartile75",
                        "dAnnoyanceQuartile25", "dAnnoyanceQuartile75", "dEventfulnessQuartile25", "dEventfulnessQuartile75", "dPleasantnessQuartile25", "dPleasantnessQuartile75"])
dataByStimTest2Merge.drop(columns=test2_drop_cols, inplace=True)
dataByStimTest2Merge['StimDuration'] = 30


# find columns in dataByStimTest1Merge not in dataByStimTest2Merge and vice versa
cols1 = dataByStimTest1Merge.columns.tolist()
cols2 = dataByStimTest2Merge.columns.tolist()
print("Not in dataByStimTest1Merge:" + str([col for col in cols2 if col not in cols1]))
print("Not in dataByStimTest2Merge:" + str([col for col in cols1 if col not in cols2]))

# find any duplicates in cols1 or cols2
print("Duplicates in dataByStimTest1Merge:" + str([col for col in cols1 if cols1.count(col) > 1]))
print("Duplicates in dataByStimTest2Merge:" + str([col for col in cols2 if cols2.count(col) > 1]))


# %%--------------------------
# aligning stimulus ID numbers
# ----------------------------

df1 = dataByStimTest1Merge.copy()
df2 = dataByStimTest2Merge.copy()

# --- StimNum: unique stimulus number across both experiments ---
assert df1['StimID'].notna().all() and df2['StimID'].notna().all(), "Missing StimID values"

stims1 = sorted(df1['StimID'].unique())
stims2 = sorted(df2['StimID'].unique())

stim_map1 = {s: i for i, s in enumerate(stims1, start=1)}
stim_map2 = {s: i for i, s in enumerate(stims2, start=len(stims1) + 1)}

df1['StimNum'] = df1['StimID'].map(stim_map1)
df2['StimNum'] = df2['StimID'].map(stim_map2)

# --- Checks ---
assert set(df1['StimNum']).isdisjoint(df2['StimNum'])               # no overlap across experiments
assert df1['StimNum'].nunique() == df1['StimID'].nunique()          # one-to-one within each dataset
assert df2['StimNum'].nunique() == df2['StimID'].nunique()
assert df1['StimNum'].notna().all() and df2['StimNum'].notna().all()

dataByStimTestMerge = pd.concat([df1.assign(Experiment=1), df2.assign(Experiment=2)], ignore_index=True)# move StimNum before StimID

cols = dataByStimTestMerge.columns.tolist()
cols.insert(0, cols.pop(cols.index('Experiment')))
cols.insert(cols.index('StimID'), cols.pop(cols.index('StimNum')))
dataByStimTestMerge = dataByStimTestMerge[cols]

dataByStimTestMerge.rename(columns={'ID': 'ParticipantNum'}, inplace=True)

# save merged data to CSV
dataByStimTestMerge.to_csv(os.path.join(outFilePath2,
                                        "refmap_listest1_2_testdataMerge_ByStim.csv"),
                           index=False)