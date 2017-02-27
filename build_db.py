#!/usr/bin/env python

import os
import sys
import commands
import getopt
import time
import datetime
import sqlite3 as lite
import glob

def atoi(text):
	return int(text) if text.isdigit() else text

def month_converter(month):
	months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
	return months.index(month) + 1

def print_help():
	print "Help"

if __name__ == '__main__':
	
	basepath_save=sys.argv[1] # Mauro: Temporary fix to test it, reads the official one but writes on my desktop
	basepath='/afs/igc.gulbenkian.pt/folders/gen-com/USERS/UBI/'
	print 'Basepath:', basepath
	print 'Basepath save', basepath_save
	#basepath='/home/mtruglio/Desktop/Test_UBI_folder'
	folders= [ name for name in os.listdir(basepath+'/results/') if (os.path.isdir(os.path.join(basepath+'/results/', name)) and name.startswith('Run')) ]


	con = None

	try:
		con = lite.connect('%s/runs_samples.db'%basepath_save) # Mauro: The temporary fix affects this too.

		cur = con.cursor()    
		cur.execute('SELECT SQLITE_VERSION()')

		data = cur.fetchone()

		print "SQLite version: %s" % data

	except lite.Error, e:

		print "Error %s:" % e.args[0]
		sys.exit(1)

	finally:
		cur.execute("CREATE TABLE IF NOT EXISTS Runs(Id TEXT PRIMARY KEY, day INT, month INT, year INT, timestamp INT, length INT, meta BOOLEAN, nextseq BOOLEAN, UNIQUE(Id))")
		cur.execute("CREATE TABLE IF NOT EXISTS Samples(Id TEXT, SampleName TEXT, Direction TEXT, Lane TEXT, Q30_Percent FLOAT, TotReads INT, TotBases_all INT, TotBases_Q30 INT, UNIQUE(Id, SampleName))")
		cur.execute("CREATE TABLE IF NOT EXISTS MetaSamples(Id TEXT, SampleName TEXT, User TEXT, TotReads_assigned INT, UNIQUE(Id, SampleName))")
		for f in folders:
			meta=0
			nextseq=0
			if '_' in f:
				run_n='_'.join(f.split('_')[1:])
			report_path=basepath+'/results/%s/report.txt'%f

			if os.path.exists(report_path) and os.path.exists(basepath+'/input/%s/SampleSheet.csv'%f):
				
				if os.path.exists(basepath+'/results/%s/Qiime'%f):
					meta=1
				if os.path.exists(basepath+'/results/%s/Merged_lanes'%f):
					nextseq=1		
	
				reads_descr=commands.getoutput("grep -a2 '\[Reads\]' %s/input/%s/SampleSheet.csv | tail -n2"%(basepath,f)).split()
				if len(reads_descr)==2:
					paired=1
				else:
					paired=0

				length=int(reads_descr[0].replace(',',''))
				if length==301: length=300
				if length==251: length=250
				if length==151: length=150				
				run_date=time.ctime(os.path.getmtime(report_path))
				month=month_converter(run_date.split()[1])
				year=run_date.split()[4]
				day=run_date.split()[2]
				timestamp=int(time.mktime(datetime.datetime.strptime('%s/%s/%s'%(day,month,year), "%d/%m/%Y").timetuple()))
				cur.execute("INSERT OR IGNORE INTO Runs VALUES('%s', %s, %s, %s, %s, %s, %s, %s)"%(run_n, day, month, year, timestamp, length, str(meta), str(nextseq)))
				oldreport=0
				stored_samplename=''
				for line in open(report_path).readlines():
					direction=''
					if line.startswith('Sample ID'):
						if not '%Q30' in line:
							oldreport=1
							#print 'NO Q30% ###########################'
						if not 'Filename' in line:
							break
						continue
					samplename=line.split('\t')[1]
					try:
						direction=samplename.split('_')[-2]
					except:
						if samplename==stored_samplename and line.split()[0]!='Total':
							direction='R2'
							stored_samplename=''
						elif line.split()[0]=='Total':
							continue
						else:
							direction='R1'
							stored_samplename=samplename
					else:
						if not direction in ['R1','R2']:
							#print('Bad report!')
							break
					
					
					try:
						lane=samplename.split('_')[-3]+'_'+direction
					except:
						break
					
					if oldreport==1:
						totbases_Q30=line.split('\t')[3].replace(',','')
						totbases_all=line.split('\t')[4].replace(',','')
						#print totbases_Q30
						#print totbases_all
						Q30_percent=100*float(totbases_Q30)/float(totbases_all)
						totreads=line.split('\t')[5].replace(',','')

					else:
						Q30_percent=line.split('\t')[4]
						totreads=line.split('\t')[6].replace(',','')
						totbases_all=line.split('\t')[5].replace(',','')
						totbases_Q30=line.split('\t')[3].replace(',','')

					#print run_n, samplename, direction, Q30_percent, totreads, totbases_all, totbases_Q30
					cur.execute("INSERT OR IGNORE INTO Samples VALUES('%s','%s','%s','%s', %s, %s, %s, %s)"%(run_n, samplename, direction, lane, Q30_percent, totreads, totbases_all, totbases_Q30))
				
				if meta == 1:
					users=[item for item in os.listdir(basepath+'/results/%s/Qiime'%f) if item not in ['params.txt', 'temp_slout', 'temp_slout_unfiltered']]
					for u in users:
						if u == 'Controls_analysis':
							continue
						if os.path.exists(basepath+'/results/%s/Qiime/%s/fastq/counts.txt'%(f, u)) and os.path.exists(basepath+'/results/%s/Qiime/%s/counts.txt'%(f, u)):
							counts=open(basepath+'/results/%s/Qiime/%s/fastq/counts.txt'%(f, u)).readlines()
							for line in counts:
								samplename=line.split()[0]
								totreads_assigned=int(line.split()[1])
								cur.execute("INSERT OR IGNORE INTO MetaSamples VALUES('%s','%s', '%s', %s)"%(run_n, samplename, u, totreads_assigned))
						else:
							print('Error for user %s, %s'%(u,f))
							continue	

		if con:
			con.commit()
			con.close()


