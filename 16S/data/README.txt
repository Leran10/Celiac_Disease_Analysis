{\rtf1\ansi\ansicpg1252\cocoartf2822
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fswiss\fcharset0 Helvetica-Bold;}
{\colortbl;\red255\green255\blue255;\red199\green200\blue201;\red26\green28\blue31;}
{\*\expandedcolortbl;;\cssrgb\c81961\c82353\c82745;\cssrgb\c13333\c14510\c16078;}
\margl1440\margr1440\vieww19700\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 metadata_20250717: an updated metadata file copied from ~/Handley Lab Dropbox/16S/Celiac/metadata/Updated_Metadata_with_Onset_Timelime.csv. Leran updated this file to include new timeline variables to consider the following: 
\fs30 \cf2 \cb3 \expnd0\expndtw0\kerning0
It has a new variable named \'93
\f1\b onset_timeline
\f0\b0 \'94, which treats the onset time points as \'93t0", and each time point prior to that will be \'93t0-6\'94 (means 6 months prior to onset), \'93t0-12" (means 12 months prior to onset), and so on. There are another two new variables names \'93onset_timeline_combined\'94 and \'93onset_timeline_numeric\'94\'a0 in the updated metadata. The onset_timeline_combined variable combines all the samples from over 42 months ahead of disease onset, because there are too few of samples at each of those time points (2 samples per time point in average). The \'93onset_timeline_numeric\'94 variable transformed the \'93onset_timeline\'94 variable into a numerical format, in case numerical is a better format for your model. The purpose of these variables is to help identify the potential viruses/bacterias that triggers the celiac disease, and from which time point ahead of disease onset.\
\
Also, comment on the phyloseq file.}