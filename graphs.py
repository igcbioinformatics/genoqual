#!/usr/bin/env python

import os
import sys
import commands
import getopt
import re
import sqlite3 as lite
import seaborn
import pandas
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import numpy as np
import datetime
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker

pandas.set_option('display.max_rows', 5000)
pandas.set_option('display.max_columns', 500)
pandas.set_option('display.width', 1000)


def mark_selection(row, myrun):

	if row['Id']==str(myrun):
		return 'Yes'

def millions(x, pos):
	'The two args are the value and tick position'
	return '%1.0fM' % (x*1e-6)

def print_help():
	print '''
The script reads in the runs database, and produces several graphs	that show run quality over time.
Arguments:
-r: current run number (will be highlighted on graphs)
-d: db file built with build_db.py
-o: output folder for images
	
	'''

if __name__ == '__main__':
	
	
	opts, args = getopt.getopt(sys.argv[1:],"hr:d:o:")
	output='./'
	myrun='0'
	meta_run = False
	nextseq_run = False
	for opt, arg in opts:
		if opt == '-h':
			print_help()
		if opt in ("-r"):
			myrun = arg
		if opt in ("-d"):
			db_path = arg
			results_path = os.path.dirname(arg) + '/results'
		if opt in ("-o"):
			output=arg
			if not os.path.exists(output):
				os.makedirs(output)
			
	con = None
	try:
		con = lite.connect(db_path)
		cur = con.cursor()
		cur.execute('SELECT SQLITE_VERSION()')
		data = cur.fetchone()
		print "SQLite version: %s"%data

	except lite.Error, e:

		print "Error %s:" % e.args[0]
		sys.exit(1)

	finally:
		query="SELECT Runs.Id, Runs.day, Runs.month, Runs.year, Runs.timestamp, Runs.length, Runs.meta, Runs.nextseq, Samples.SampleName, Samples.Direction, Samples.Lane, Samples.Q30_Percent, Samples.TotReads, Samples.TotBases_all, Samples.TotBases_Q30 FROM Runs INNER JOIN Samples on Runs.Id=Samples.Id"
		# cur.execute(query)
		# rows = cur.fetchall()
		fig, ax = plt.subplots()
		# for row in rows:
		# 	print row		
		#with open('%s/htmldummy'%output)
		df=pandas.read_sql(query, con)
		df['short_date']=(df['month'].map(str)+'/'+df['year'].map(str)).astype('category')
		df=df.sort_values(['year', 'month'],)
		df['selection'] = df.apply(mark_selection,axis=1,args=(myrun,))
		order= df.short_date.unique()
		# Q30% stripplot for R1
		fig, ax = plt.subplots()		
		df_R1=df[df.Direction=='R1']
		df_R1_selection=df_R1[df_R1.Id==myrun]
		#Flagging the run as meta if it's the case
		if df_R1_selection['meta'].iloc[0]:
			meta_run=True
		if df_R1_selection['nextseq'].iloc[0]:
			nextseq_run=True	
			print 'Nextseq!'
		# 
		df_R1_others=df_R1[df_R1.Id!=myrun]
		stripplotR1=seaborn.stripplot(x='short_date',y='Q30_Percent', data=df_R1_others, hue='length',jitter=True, order=order)
		stripplotR1_sel=seaborn.stripplot(x='short_date',y='Q30_Percent', marker='*', size=10, color='red', data=df_R1_selection, order=order, jitter=True)
		stripplotR1.set(ylim=(0, 100))
		for label in stripplotR1.get_xticklabels():
			label.set_rotation(45)
		stripplotR1.legend(loc='lower left', title='Length')
		stripplotR1.set_title('Per Sample %Q30 bases for R1\n', fontsize=16, fontweight='bold')
		plt.tight_layout()				
		plt.savefig('%s/fig1.png'%output)
		#seaborn.plt.show()
		
		
		# Q30% stripplot for R2
		fig, ax = plt.subplots()		
		df_R2=df[df.Direction=='R2']
		df_R2_selection=df_R2[df_R2.Id==myrun]
		df_R2_others=df_R2[df_R2.Id!=myrun]
		stripplotR2=seaborn.stripplot(x='short_date',y='Q30_Percent', data=df_R2_others, hue='length',jitter=True, order=order)
		stripplotR2_sel=seaborn.stripplot(x='short_date',y='Q30_Percent', marker='*', size=10, color='red', data=df_R2_selection, order=order, jitter=True)
		stripplotR2.set(ylim=(0, 100))
		for label in stripplotR2.get_xticklabels():
			label.set_rotation(45)
		stripplotR2.legend(loc='lower left', title='Length')
		stripplotR2.set_title('Per Sample %Q30 Bases for R2\n', fontsize=16, fontweight='bold')
		plt.tight_layout()				
		plt.savefig('%s/fig2.png'%output)
		#seaborn.plt.show()


		#Q30% boxplot
		fig, ax = plt.subplots()				
		df=df.sort_values(['Direction', 'length'],)
		df_selection=df[df.Id==myrun]
		df_others=df[(df.Id!=myrun) & (df.nextseq==nextseq_run)]
		if nextseq_run:
			boxplot=seaborn.boxplot(x="Lane", y="Q30_Percent", hue="length", data=df_others)			
		else:
			boxplot=seaborn.boxplot(x="Direction", y="Q30_Percent", hue="length", data=df_others)
		seaborn.despine(offset=20, trim=True)
		#boxplot.set(ylim=(0, 100))
		boxplot.legend(loc='lower left', title='Length')
		if df_selection.empty!=True:
			if nextseq_run:
				stripplot_sel=seaborn.stripplot(x='Lane',y='Q30_Percent', marker='*', size=10, color='red', data=df_selection, jitter=False)				
			else:
				stripplot_sel=seaborn.stripplot(x='Direction',y='Q30_Percent', marker='*', size=10, color='red', data=df_selection, jitter=False)
			stripplot_sel.legend(loc='lower left', title='Length')	
		boxplot.set_title('Per Sample %Q30 Bases\n', fontsize=16, fontweight='bold')		
		plt.tight_layout()				
		plt.savefig('%s/fig3.png'%output)		
		#seaborn.plt.show()


		#Total reads boxplot
		fig, ax = plt.subplots()	
		df_R1_others = df_R1_others[df.nextseq==nextseq_run]
		df_R1_others_grouped=df_R1_others.groupby(['timestamp', 'length'], as_index=False).sum()	
		print df_R1_others_grouped
		if df_selection.empty!=True:
			df_R1_selection_grouped=df_R1_selection.groupby(['timestamp', 'length'], as_index=False).sum()
			possible_lengths=sorted(df.length.unique())
			print df_R1_selection_grouped
			print possible_lengths
			for c in range(1,len(possible_lengths)+1): #adding fake rows for all possible lengths, so that the selection graph has the same categories as the 'others' graph overlayed on it
				df_R1_selection_grouped.loc[c]=[0,possible_lengths[c-1],0,0,0,0,0,0,900000000,0,0]	
			print df_R1_selection_grouped
			stripplot_sel=seaborn.stripplot(x='length',y='TotReads', marker='*', size=16, color='red', data=df_R1_selection_grouped, jitter=False, ax=ax)
			stripplot_sel.legend(loc='lower left', title='Length')
		boxplot=seaborn.boxplot(x="length", y="TotReads", data=df_R1_others_grouped)
		seaborn.despine(offset=20)			
		if nextseq_run:
			boxplot.set(ylim=(0, 400000000))	
		else:
			boxplot.set(ylim=(0, 30000000))			
		boxplot.set_title('Total Reads per Run\n', fontsize=16, fontweight='bold')
		boxplot.yaxis.set_major_formatter(ticker.FuncFormatter(millions))	
		boxplot.legend(loc='upper right', title='Length')
			
		plt.tight_layout()				
		plt.savefig('%s/fig4.png'%output)			
		#seaborn.plt.show()

		#Total reads line plot
		fig, ax = plt.subplots()		
		fig.set_size_inches(21.7, 7.27)
		df_R1=df_R1[df.nextseq==nextseq_run]
		df_R1['timestamp']=df_R1['timestamp'].astype(int)
		pointplot_bycolor=seaborn.pointplot(x=df_R1['timestamp'],y=df_R1['TotReads'], join=False, size=10, hue=df_R1["length"], ci=None, estimator=np.sum)
		pointplot_all=seaborn.pointplot(x=df_R1['timestamp'],y=df_R1['TotReads'], join=True, size=10, ci=None, estimator=np.sum)
		if df_selection.empty!=True:			
			index_sel=np.where(df_R1['timestamp'].unique()==df_R1_selection['timestamp'].unique())[0]		
			star=plt.plot(index_sel, [df_R1_selection['TotReads'].sum()], marker='*', color='r', markersize=21)
		labels= [item.get_text() for item in pointplot_all.get_xticklabels()]		
		pointplot_bycolor.legend(loc='upper right', title='Length')
		for i in range(0,len(labels)):
			labels[i]=datetime.datetime.fromtimestamp(int(labels[i])).strftime('%m/%y')
		pointplot_all.set_xticklabels(labels)
		for label in pointplot_bycolor.get_xticklabels():
			label.set_rotation(45)
			
		df_R1_grouped=df_R1.groupby(['timestamp', 'length','Id'], as_index=False).sum()
	
		for p in zip(pointplot_all.get_xticks(), df_R1_grouped['Id'], df_R1_grouped['TotReads']):
			pointplot_all.text(p[0], int(p[2])+1000000, p[1], color='gray', fontsize=10, zorder=103) 
		
		plt.setp(pointplot_all.lines, zorder=1)	
		plt.setp(pointplot_bycolor, zorder=100)	
		plt.setp(star, zorder=102)	
		
		pointplot_bycolor.set_title('Total Reads per Run\n', fontsize=16, fontweight='bold')		
		pointplot_bycolor.set_xlabel('\nDate', fontsize=12)
		pointplot_bycolor.set_ylabel('\nReads', fontsize=12)
		
		pointplot_bycolor.yaxis.set_major_formatter(ticker.FuncFormatter(millions))
		
		plt.tight_layout()		
		fig.savefig('%s/fig5.png'%output, figsize=(25,9))					
		#seaborn.plt.show()

			
		if meta_run==True:
			
			query="SELECT Runs.Id, Runs.day, Runs.month, Runs.year, Runs.timestamp, Runs.length, MetaSamples.SampleName, MetaSamples.User, MetaSamples.TotReads_assigned FROM Runs INNER JOIN MetaSamples on Runs.Id=MetaSamples.Id"
			# cur.execute(query)
			# rows = cur.fetchall()
			fig, ax = plt.subplots()
			# for row in rows:
			# 	print row		
			#with open('%s/htmldummy'%output)
			df_meta=pandas.read_sql(query, con)
			df_meta=df_meta.sort_values(['year', 'month'],)			
			df_meta['TotReads_assigned']=df_meta['TotReads_assigned'].astype(int)
			
			#Line plots
			fig, ax = plt.subplots()		
			fig.set_size_inches(21.7, 7.27)
			df_R1_meta=df_R1[df_R1.meta==1]
			pointplot_assigned_all=seaborn.pointplot(x=df_meta['timestamp'],y=df_meta['TotReads_assigned'], join=True, color='skyblue',size=10, zorder=1, ci=None, estimator=np.sum, label='Line 2')			
			
			df_R1_meta['timestamp']=df_R1_meta['timestamp'].astype(int)

			pointplot_all=seaborn.pointplot(x=df_R1_meta['timestamp'],y=df_R1_meta['TotReads'], join=True,size=10, zorder=3, ci=None, estimator=np.sum)
			pointplot_bycolor=seaborn.pointplot(x=df_R1_meta['timestamp'],y=df_R1_meta['TotReads'], join=False, size=10, hue=df_R1_meta["length"], zorder=10, ci=None, estimator=np.sum)
			pointplot_bycolor.set(ylim=(0,40000000))	
			if df_selection.empty!=True:			
				index_sel=np.where(df_R1_meta['timestamp'].unique()==df_R1_selection['timestamp'].unique())[0]		
				star=plt.plot(index_sel, [df_R1_selection['TotReads'].sum()], marker='*', color='r', markersize=21)
			labels= [item.get_text() for item in pointplot_all.get_xticklabels()]		

			for i in range(0,len(labels)):
				labels[i]=datetime.datetime.fromtimestamp(int(labels[i])).strftime('%m/%y')
			pointplot_all.set_xticklabels(labels)
			for label in pointplot_bycolor.get_xticklabels():
				label.set_rotation(45)
			df_R1_meta_grouped=df_R1_meta.groupby(['timestamp', 'length','Id'], as_index=False).sum()
			for p in zip(pointplot_all.get_xticks(), df_R1_meta_grouped['Id'], df_R1_meta_grouped['TotReads']):
				pointplot_all.text(p[0], int(p[2])+1000000, p[1], color='gray', fontsize=10, zorder=103) 
			pointplot_bycolor.set_title('Total Reads per Run (16s)\n', fontsize=16, fontweight='bold')
			
			handles, labels = pointplot_bycolor.get_legend_handles_labels()
			plt.legend(handles=handles+[plt.Line2D((0,1),(0,0), color='skyblue', marker='o', label='Assigned')])
			
			pointplot_bycolor.set_xlabel('\nDate', fontsize=12)
			pointplot_bycolor.set_ylabel('\nReads', fontsize=12)
			plt.setp(pointplot_bycolor.collections, zorder=100)		
			plt.setp(star, zorder=101)
			
			pointplot_bycolor.yaxis.set_major_formatter(ticker.FuncFormatter(millions))

			plt.tight_layout()		
			fig.savefig('%s/fig5meta.png'%output, figsize=(25,9))						
			
			
			#Total reads with assigned sample for GEU, boxplot meta
			fig, ax = plt.subplots()		
			fig.set_size_inches(23.7, 9.27)			
			df_GEU=df_meta[df_meta.User=='GEU']
			df_GEU['Type'] = np.where(np.logical_or(df_GEU.SampleName.str.contains("Pos"),df_GEU.SampleName.str.contains("pos")), 'Pos', 'Neg')
			
			df_GEU_grouped=df_GEU.groupby(['timestamp', 'Id', 'length', 'Type'], as_index=False).sum()	
			df_GEU_grouped_selection=df_GEU_grouped[(df_GEU_grouped.Id==myrun) & (df_GEU_grouped.Type=='Pos')]
			g = seaborn.pointplot(x="timestamp", y="TotReads_assigned", hue="Type", data=df_GEU_grouped)
			if df_GEU_grouped_selection.empty!=True:			
				index_sel=np.where(df_GEU['timestamp'].unique()==df_GEU_grouped_selection['timestamp'].unique())[0]				
				star=plt.plot(index_sel, [df_GEU_grouped_selection['TotReads_assigned'].sum()], marker='*', color='r', markersize=21)			
			labels= [item.get_text() for item in g.get_xticklabels()]		
			
			for i in range(0,len(labels)):
				labels[i]=datetime.datetime.fromtimestamp(int(labels[i])).strftime('%m/%y')
			g.set_xticklabels(labels)
			for label in g.get_xticklabels():
				label.set_rotation(45)	
					
			ids=sorted(set(list(df_GEU_grouped['Id'])), key=lambda x: list(df_GEU_grouped['Id']).index(x))
			for p in zip(g.get_xticks(), ids, df_GEU_grouped.iloc[::2, :]['TotReads_assigned']):
				g.text(p[0], int(p[2])+50000, p[1], color='gray', fontsize=10, zorder=103) 	
			g.set(ylim=(0,1000000))	
			g.set_title('Total Reads assigned to controls (16s)\n', fontsize=18, fontweight='bold')		
			g.set_xlabel('\nDate', fontsize=12)
			g.set_ylabel('\nReads', fontsize=12)			
			fig.savefig('%s/fig6meta.png'%output, figsize=(25,9))						
			plt.tight_layout()

			# Total assigned reads per user in the single run
			fig, ax = plt.subplots()	
			df_current=df_meta[df_meta.Id==myrun]
			df_current_grouped=df_current.groupby(['User'], as_index=False).sum()
			barplot=seaborn.barplot(x="User", y="TotReads_assigned", data=df_current_grouped)
			#barplot.yaxis.set_major_formatter(ticker.FuncFormatter(millions))
			barplot.set_ylabel('\nReads', fontsize=12)
			barplot.set_title('Total Reads assigned to each user (16s)\n', fontsize=14, fontweight='bold')					
			fig.savefig('%s/fig7meta.png'%output)
			
			# QIIME PCoA plot for controls
