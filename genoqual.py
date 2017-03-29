# -*- coding: utf-8 -*-

from __future__ import print_function, division
import os
import sys
import logging
import argparse
import tempfile
import shutil
import gzip
import string
import urllib
import atexit
import urlparse
import csv
import datetime
from time import sleep
from Bio import SeqIO
from itertools import cycle, imap, izip
from contextlib import nested
from glob import glob
from subprocess import Popen, PIPE
from HTSeq import FastqReader
from numpy import std
from collections import OrderedDict


def hamming(str1, str2):
	"""Compute hamming distance between 2 strings

	>>> hamming('abc', 'abc')
	0
	>>> hamming('abc', 'abcd')
	1
	>>> hamming('abc', 'cbc')
	1
	>>> hamming('abc', 'bca')
	3
	"""
	return sum(imap(str.__ne__, str1, str2))


def get_by_distance_from_dict(d, key, size_matters=True, distance=hamming,
                              threshold=1):
	"""Return all values of given dictionary that have a smaller or equal
	distance to the given threshold.

	If size_matters=True keys must have the same size or they will be treated
	as above threshold.
	"""
	matches = []
	size = len(key)

	for k in d.iterkeys():
		# If key size is different, assume no possible match
		if size_matters and len(k) != size:
			continue

		if distance(k, key) <= threshold:
			matches.append(k)

	if len(matches) > 1:
		raise Error("Read with sequence {0} ambiguously matched more than "
		            "one tag {1}".format(key, matches))

	if matches:
		return matches[0]
	else:
		return None


class Error(Exception):
	pass


class ConfigParser(object):
	"""Reads the Config format (which is basically a CSV)
	"""
	def __init__(self, params):
		self.params = params
		self.configfile = self.params["config"]
		self.config = None

		# Ensure characters are valid to avoid problems afterwards
		self.transtable = string.maketrans(" .", "--")
		self.validchars = string.ascii_letters + string.digits + " .-_"

		# Columns expected to be present in the config file
		self._expected_columns = self.params["columns"]
		# Boolean columns (columns that should map to True/False
		self._booleans = self.params["booleans"]

		self.parse()

	def _booleanify(self, val):
		"""Try to convert given value to a boolean in a somewhat extensive way

		It will try to convert Yes/No 0/1 Y/N to a boolean. Other values will
		raise an exception.
		"""
		bools = {
		        "yes": True,
		        "no": False,
		        "y": True,
		        "n": False,
		        "1": True,
		        "0": False,
		        '': False,
		}

		try:
			return bools[val.lower()]
		except KeyError:
			raise Error("Value '{0}' is not a valid option for that column, "
			            "should be one of {1}".format(val, bools.keys()))

	def _sanitize(self, filename):
		fname = filename.translate(self.transtable)

		for char in fname:
			if char not in self.validchars:
				raise Error("Config file contains invalid characters in "
				            "column 'Sample ID'. Offending character is "
				            "{0!r} in '{1}'".format(char, filename))

		return fname

	def parse(self):
		"""Read all settings from given configfile
		"""
		if not os.path.isdir(self.params["input_path"]):
			raise Error("Folder '{0}' not found, did you specify the prefix "
			            "correctly?".format(self.params["input_path"]))

		if not os.path.isfile(self.configfile):
			raise Error("Config file not found at '{0}'"
			            .format(self.configfile))

		# Prevent parsing the file more than once if parse is called twice
		if self.config:
			log.info("Using cached config settings")
			return

		log.info("Parsing config file '%s'", os.path.basename(self.configfile))
		with open(self.configfile) as fh:
			for line in csv.reader(fh):
				# Ignore empty lines
				if not line:
					continue

				if self.config is None:
					self.config = []
					self._ids = line

					# Check that required columns are present in the file
					for col in self._expected_columns:
						if col not in self._ids:
							raise Error("Column '{0}' is missing from "
							            "'{0}'".format(col, self.configfile))

					continue

				data = {}

				for i, value in enumerate(zip(self._ids, line)):
					id, val = value
					data[id] = val

				# Keep an ID safe for use as filename prefix
				data["ID"] = self._sanitize(data["Sample_ID"])
				
				# Figure out what should be the filename from sample data
				data["prefix"] = "{0}_".format(data["ID"])

				# Create a variable that will control the execution of blastn for contaminants when necessary 
				data["skip_contaminants"] = False
				
				# Transform requested columns into a boolean
				for key, value in data.iteritems():
					if key in self._booleans:
						data[key] = self._booleanify(value)

				log.debug("Parsed config line as %s", data)
				self.config.append(data)

		if self.config is None:
			raise Error("Provided config file is empty")


class DataWrapper(object):
	"""Object used to create a temporary folder where data and temporary files
	will be created.

	Provided files are usually in gzip format. This wrapper also ensures they
	are decompressed in place.

	Needs a ConfigParser object containing a data entry with information about
	each sample used.
	"""
	def __init__(self, config):
		self.dir = tempfile.mkdtemp(prefix="genoqual_")
		atexit.register(lambda: self.cleanup)

		self.config = config
		self.data = []

	def cleanup(self):
		if os.path.isdir(self.dir):
			log.debug("Removing temporary folder and all its contents")
			shutil.rmtree(self.dir)

	def _get_files(self, prefix):
		"""Get all files that match the given pattern in the input folder
		"""
		files = glob(os.path.join(self.config.params["input_path"],
		                          "{0}*".format(prefix)))

		# Check that we have files with valid extensions
		valid = set([".gz", ".gzip", ".fastq"])
		valid_files = []
		for file in files:
			if os.path.splitext(file)[-1] in valid:
				valid_files.append(file)
			else:
				log.warn("Unexpected file extension in file '%s'. "
				         "Ignoring file.", file)
		return tuple(sorted(valid_files))

	def prepare(self):
		log.debug("Making sure final output folder is clean")

		final_outdir = self.config.params["output_path"]
		if os.path.isdir(final_outdir):
			if not glob(os.path.join(final_outdir, "*")):
				shutil.rmtree(final_outdir)
			else:
				raise Error("Output folder '{0}' exists and is not empty. "
				            "Remove it before re-launching the analysis."
				            .format(final_outdir))

		log.debug("Looking for files matching identifiers in the config file")

		for entry in self.config.config:
			prefix = entry["prefix"]
			files = self._get_files(prefix)
			log.debug("Files for prefix %s"%prefix)
			entry["original_files"] = files

			# Check that we have paired reads
			count = len(files)
			if not count:
				raise Error("No files found matching the prefix '{0}'"
				            .format(prefix))
			
			if count == 1:
				log.info("Assuming data as single-end on Sample ID '%s'.",
				         entry["Sample_ID"])

			elif count == 2 and self.config.params["meta"]==False:
				log.info("Assuming data as paired-end on Sample ID '%s'.",
				         entry["Sample_ID"])
			elif count == 4:
				log.info("Assuming data as NextSeq single-end on Sample ID '%s'.",
					 entry["Sample_ID"])	
				self.config.params["nextseq"] = True
			elif count == 8:
				log.info("Assuming data as NextSeq paired-end on Sample ID '%s'.",
						 entry["Sample_ID"])	
				self.config.params["nextseq"] = True
				
			elif self.config.params["meta"]==True:			
				if count == 3:
					log.info("All three files for metagenomics analysis were found for Sample ID '%s'.",
					         entry["Sample_ID"])
				elif count == 2 and '_I1_001.fastq' in ' '.join(files):
					log.info("Two files for metagenomics analysis were found for Sample ID '%s'. Assuming R1 and I1.",
						 entry["Sample_ID"])					
				else:
					raise Error("Unexpected number of files for metagenomics analysis matching the prefix "
					            "'{0}'. Matching files were: '{1}'"
					            .format(prefix, files))
			else:
				raise Error("Unexpected number of files matching the prefix "
				            "'{0}'. Matching files were: '{1}'"
				            .format(prefix, files))

			fasta_files = []
			for file in files:
				# Decompress each gzip file
				fname = os.path.basename(file)
				name, ext = os.path.splitext(fname)

				if ext in [".gz", ".gzip"]:
					out = os.path.join(self.dir, name)

				else:
					out = os.path.join(self.dir, fname)

				# Keep track of the files in the temporary folder
				fasta_files.append(out)

			entry["files"] = tuple(fasta_files)

			self.data.append(entry)


	def run(self):
		for entry in self.data:
			for file, out in zip(entry["original_files"], entry["files"]):
				# Decompress each gzip file
				fname = os.path.basename(file)
				name, ext = os.path.splitext(fname)

				if ext in [".gz", ".gzip"]:
					log.info("Decompressing '%s' into temporary folder", fname)
					os.popen('zcat {0} > {1}'.format(file, out))
					#with gzip.open(file) as infile:
						#with open(out, 'w') as outfile:
							#while True:
								## Read/Write 1 MB at a time
								#data = infile.read(1024 ** 2)
								#if not data:
									#break
								#outfile.write(data)
				else:
					log.info("Copying '%s' to temporary folder", fname)
					# This should only be fastq files
					shutil.copyfile(file, out)
			if self.config.params["nextseq"] == True:
				log.info("Joining lanes R1 for sample {0}".format(entry["Sample_ID"]))
				self.join_lanes(entry["Sample_ID"], "R1")
				if len(entry["files"])==8:
					log.info("Joining lanes R2 for sample {0}".format(entry["Sample_ID"]))
					self.join_lanes(entry["Sample_ID"], "R2")
					
			

	def join_lanes(self, sampleID, direction):
		if not os.path.exists(os.path.join(self.dir,'Merged_lanes')):
			os.makedirs(os.path.join(self.dir,'Merged_lanes'))
		os.system("cd {0}; cat {1}_*_{2}_* > {0}/Merged_lanes/{1}_{2}.fastq ".format(self.dir, sampleID, direction))
		
class DataCollector(object):
	"""Reaps all objects used in the pipeline and formats data in a friendly
	memory representation for text/html output
	"""
	def __init__(self, pipeline):
		self.pipe = pipeline
		self.output = []
		self.demux = {}
		self.links = set([])
		self.stats_all_lanes = {}
		self.collect()

	def collect(self):
		columns = (
		        "Sample ID",
		        "Filename",
		        "Q20",
		        "Q30",
		        "%Q30",
		        "Total bases",
		        "Total reads",
		        "FastQC",
		        "Primary organism",
		        "Contaminants",
		        "Total mapped bases",
		        "Mean coverage",
		        "Theoretical coverage Q20",
		        "Theoretical coverage Q30",
		        "Total mapped reads",
		        "Total reads (used)",
		        "Qualimap",
		        "Demux read mean",
		        "Demux read stdev",
		        "Demux reads failed",
		        "Demux",
		        "Merged reads",
		        "Merged reads failed",
		)

		self.output.append(columns)
		# These columns (index number) should be links in HTML output
		self.links.update([7, 9, 16, 20])#, 20])

		sum_q20 = 0
		sum_q30 = 0
		sum_bases = 0
		sum_reads = 0
		sum_mapped_bases = 0
		sum_mapped_reads = 0
		sum_total_reads_used = 0
		theoretical_cov_Q20_sum = 0
		theoretical_cov_Q30_sum = 0
		multiqc_file = self.pipe.data.config.params["output_url"]+'/multiqc_report.html'


		for entry in self.pipe.data.config.config:
			sample_id = entry["Sample_ID"]
			self.stats_all_lanes[sample_id]=OrderedDict()
			self.stats_all_lanes[sample_id]['q20'] = 0
			self.stats_all_lanes[sample_id]['q30'] = 0
			self.stats_all_lanes[sample_id]['percent_q30'] = 0
			self.stats_all_lanes[sample_id]['total_bases'] = 0
			self.stats_all_lanes[sample_id]['total_reads'] = 0
			self.stats_all_lanes[sample_id]['fastqc'] = '...'
			self.stats_all_lanes[sample_id]['primary'] = ''
			self.stats_all_lanes[sample_id]['contam'] = ''
			self.stats_all_lanes[sample_id]['mapped_bases'] = ''
			self.stats_all_lanes[sample_id]['mean_cov'] = ''
			self.stats_all_lanes[sample_id]['cov_q20'] = ''
			self.stats_all_lanes[sample_id]['cov_q30'] = ''
			self.stats_all_lanes[sample_id]['mapped_reads'] = ''					
			self.stats_all_lanes[sample_id]['used_reads'] = ''															
			self.stats_all_lanes[sample_id]['qualimap'] = ''
			self.stats_all_lanes[sample_id]['demux_mean'] = ''
			self.stats_all_lanes[sample_id]['demux_stdev'] = ''
			self.stats_all_lanes[sample_id]['demux_failed'] = ''
			self.stats_all_lanes[sample_id]['demux'] = ''
			self.stats_all_lanes[sample_id]['merged_reads'] = ''
			self.stats_all_lanes[sample_id]['merged_failed'] = ''
			
			
			for filename in entry["files"]:
				if filename.endswith("I1_001.fastq"):
					continue
				fname = os.path.basename(filename)

				quality = self.pipe.qual.data[filename]
				q20 = str(quality[20])
				q30 = str(quality[30])
				total_bases = str(quality["total_bases"])
				total_reads = str(quality["total_reads"])
				sum_q20+=int(q20)
				sum_q30+=int(q30)
				sum_bases+=int(total_bases)
				sum_reads+=int(total_reads)
				
				self.stats_all_lanes[sample_id]['q20'] += int(q20)
				self.stats_all_lanes[sample_id]['q30'] += int(q30)			
				self.stats_all_lanes[sample_id]['total_bases'] += int(total_bases)
				self.stats_all_lanes[sample_id]['total_reads'] += int(total_reads)
				
				percent_q30 = str(round((quality[30]/quality["total_bases"])*100, 2))

				fastqc_file = self.pipe.qc.data[filename].replace(
				        self.pipe.data.dir,
				        self.pipe.data.config.params["output_url"],
				)
				
				# jcostaDamage
				if filename in self.pipe.tax.data:
					try:
						cont_file, primary_org = self.pipe.tax.data[filename]
					except KeyError:
						log.debug("Filename '%s' has no contaminants "
						          "information, ignoring.", filename)
						contaminants_file = primary_org = ""
					else:
						if primary_org == "=":
							contaminants_file = primary_org
						else:
							contaminants_file = cont_file.replace(
							        self.pipe.data.dir,
							        self.pipe.data.config.params["output_url"],
							)
							self.stats_all_lanes[sample_id]['primary'] = primary_org
							if os.path.basename(contaminants_file)!='skipped':
								self.stats_all_lanes[sample_id]['contam'] = '<a href={0}>Report</a>'.format(contaminants_file)
							else:
								self.stats_all_lanes[sample_id]['contam'] = 'skipped'
							
							
				else:
					contaminants_file = '='
					primary_org = '='

				try:
					ref = entry["Reference"]
				except KeyError:
					ref = None

				if ref and filename in self.pipe.cov.data:
					cov = self.pipe.cov.data[filename]

					total_mapped_bases = cov["total_mapped_bases"]
					self.stats_all_lanes[sample_id]['mapped_bases'] = total_mapped_bases
					mean_coverage = cov["mean_coverage"]
					try:
						self.stats_all_lanes[sample_id]['mean_cov'] = "{0:.2f}".format(float(mean_coverage))
					except ValueError:
						self.stats_all_lanes[sample_id]['mean_cov'] = mean_coverage
					try:
						theoretical_cov_Q20 = "{0:.2f}".format(quality[20]/cov["ref_length"])
					except ValueError:
						theoretical_cov_Q20 = 'ERROR'	
					self.stats_all_lanes[sample_id]['cov_q20'] = theoretical_cov_Q20
					try:
						theoretical_cov_Q30 = "{0:.2f}".format(quality[30]/cov["ref_length"])
					except ValueError:
						theoretical_cov_Q30 = 'ERROR'	
					self.stats_all_lanes[sample_id]['cov_q30'] = theoretical_cov_Q30
					
					total_mapped_reads = cov["total_mapped_reads"]
					self.stats_all_lanes[sample_id]['mapped_reads'] = total_mapped_reads					
					total_reads_used = cov["total_reads_used"]
					self.stats_all_lanes[sample_id]['used_reads'] = total_reads_used															
					qualimap_file = cov["qualimap"]
					self.stats_all_lanes[sample_id]['qualimap'] = qualimap_file #get it as it is (tmp path)


					theoretical_cov_Q20_sum+=float(theoretical_cov_Q20)
					theoretical_cov_Q30_sum+=float(theoretical_cov_Q30)

					if not qualimap_file in ("=", "ERROR"):
						sum_mapped_bases+=int(total_mapped_bases)
						sum_mapped_reads+=int(total_mapped_reads.split()[0])
						sum_total_reads_used+=int(total_reads_used)

				else:
					total_mapped_bases = total_mapped_reads = ""
					total_reads_used = qualimap_file = mean_coverage = ""
					theoretical_cov_Q20 = ''
					theoretical_cov_Q30 = ''

				if qualimap_file != "=":
					qualimap_file = qualimap_file.replace(
					        self.pipe.data.dir,
					        self.pipe.data.config.params["output_url"],
					)
					if qualimap_file!='':
						self.stats_all_lanes[sample_id]['qualimap'] = "<a href=\"{0}\">Report</a>".format(qualimap_file) # convert it to proper address
				
				if filename in self.pipe.dmux.data:
					if self.pipe.dmux.data[filename] == "=":
						demux_mean = "="
						demux_stdev = "="
						demux_read_failed = "="
						demux_url = "="
					else:
						samples, demux = self.pipe.dmux.data[filename]

						dmux_read = []
						demux_mean = 0
						demux_stdev = 0
						demux_read_failed = 0

						total = len(demux)

						log.debug("Demux samples data is %s", samples)
						log.debug("Demux data is %s", demux)
						counts = demux["options"]["read_counts"]

						for key1 in counts:
							if key1 == "rejected":
								demux_read_failed += counts["rejected"]
								continue

							for key2 in counts[key1]:
								dmux_read.append(counts[key1][key2])

						demux_read_failed = str(demux_read_failed)
						demux_mean = "{0:.2f}".format(sum(dmux_read) / total)
						demux_stdev = "{0:.2f}".format(std(dmux_read, ddof=1))

						demux_file = os.path.join(
						        demux["options"]["samplepath"],
						        "demux.html"
						)

						self.collect_demux_merge(filename, demux_file)

						demux_url = demux_file.replace(
						        self.pipe.data.dir,
						        self.pipe.data.config.params["output_url"],
						)
				else:
					demux_mean = ""
					demux_stdev = ""
					demux_read_failed = ""
					demux_url = ""
				
				if filename in self.pipe.merge.data:
					if self.pipe.merge.data[filename] is "=":
						merged_reads = "="
						merged_reads_failed = "="
					else:
						merged_reads = 0
						merged_reads_failed = 0

						for val in self.pipe.merge.data[filename].itervalues():

							merged_reads += val["reads_combined"]
							merged_reads_failed += val["reads_notcombined"]
						
						merged_reads = str(merged_reads)
						merged_reads_failed = str(merged_reads_failed)
						self.stats_all_lanes[sample_id]['merged_reads'] = merged_reads
						self.stats_all_lanes[sample_id]['merged_failed'] = merged_reads_failed					
				else:
					merged_reads = ""
					merged_reads_failed = ""
					
				
				self.output.append([
				        sample_id,
				        fname,
				        q20,
				        q30,
				        percent_q30,
				        total_bases,
				        total_reads,
				        fastqc_file,
				        primary_org,
				        contaminants_file,
				        total_mapped_bases,
				        mean_coverage,
				        theoretical_cov_Q20,
				        theoretical_cov_Q30,
				        total_mapped_reads,
				        total_reads_used,
				        qualimap_file,
				        demux_mean,
				        demux_stdev,
				        demux_read_failed,
				        demux_url,
				        merged_reads,
				        merged_reads_failed,
				])




		self.output.append(['Total',
		                          '',
		                          str(sum_q20),
		                          str(sum_q30),
		                          "{0:.2f}".format(100*sum_q30/float(sum_bases)),
		                          str(sum_bases), str(sum_reads),
		                          multiqc_file, '', '',
		                          str(sum_mapped_bases),
		                          '','','',

		                          str(sum_mapped_reads),
		                          str(sum_total_reads_used),
		                          ''])
		

	def collect_demux_merge(self, filename, demux_file):
		columns = (
		        "Demux ID",
		        "Seq.F",
		        "Seq.R",
		        "Demuxed reads",
		        "Demux failed reads",
		        "Merged reads",
		        "Merge failed reads",
		)

		data = []

		data.append(columns)

		files, demux = self.pipe.dmux.data[filename]

		log.debug("Preparing demux data %s for file %s", demux, filename)
		samples = demux["options"]["samples"]
		stats = demux["options"]["read_counts"]

		for i, sample_id in enumerate(sorted(samples)):
			seqF, seqR = samples[sample_id]

			# Demuxed files
			demuxed_reads = str(stats[seqF][seqR])
			if i == 0:
				# First line of the table
				demux_failed_reads = str(stats["rejected"])
			else:
				demux_failed_reads = "="

			# Merged files
			if filename in self.pipe.merge.data:
				vals = self.pipe.merge.data[filename][sample_id]
				reads_combined = str(vals["reads_combined"])
				reads_notcombined = str(vals["reads_notcombined"])
			else:
				reads_combined = ""
				reads_notcombined = ""

			values = (sample_id, seqF, seqR, demuxed_reads, demux_failed_reads,
			          reads_combined, reads_notcombined)
			log.debug("Demux file values %s", values)
			data.append(values)

		self.demux[filename] = (demux_file, data)

class DataFormatter(object):
	"""Used to generate the HTML report that gets produced at the end of all
	steps
	"""
	def __init__(self, data, config):
		self.data = data.output
		self.config = config
		self.stats_all_lanes = data.stats_all_lanes
		self.links = data.links
		self.demux = data.demux

	def write_tabular(self, filename):
		"""Write a tab delimited version of the data
		"""
		with open(filename, 'w') as fh:
			for line in self.data:
				fh.write("\t".join(line) + "\n")

		for filename in self.demux:
			path, values = self.demux[filename]
			file = os.path.splitext(path)[0] + ".txt"

			with open(file, 'w') as fh:
				for line in values:
					fh.write("\t".join(line) + "\n")

	def encode_figures(self, output_path):
		encoded_figs=OrderedDict()
		for img in sorted(glob('%s/tmp_images/*.png'%output_path)):
			with open(img, "rb") as f:
				encoded_figs[os.path.basename(img)]=f.read().encode("base64")

		shutil.rmtree('%s/tmp_images/'%output_path)

		return encoded_figs

	def write_html(self, filename):
		"""Write an html version of the data
		"""
		self._html_file(filename, self.data, "Quality Control", self.links)

		for filename in self.demux:
			file, values = self.demux[filename]
			self._html_file(file, values, "Demultiplexing and merging", ())

	def _html_file(self, filename, values, title, links):

		encoded_figs=self.encode_figures(os.path.dirname(filename))
		row_color = cycle((
		        "style='background-color:white;'",
		        "style='background-color:#eefafc;'",
		))

		i = datetime.datetime.now()
		with open(filename, 'w') as fh:
			fh.write("<!DOCTYPE html>\n")
			fh.write("<html>\n<head>\n<title>{0}</title>\n".format(title))
			fh.write("<link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css' integrity='sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u' crossorigin='anonymous'>\n")
			fh.write("<script type='text/javascript' src='http://code.jquery.com/jquery-1.9.1.js'></script>\n")	
			fh.write("<script>\n")
			fh.write("$(document).ready(function(){\n")
			fh.write("\t$('tbody > tr:not(\".header\")').hide();\n")
			fh.write("\t$('.header').click(function(){\n")
			fh.write("\t\t$(this).toggleClass('expand').nextUntil('tr.header').slideToggle(100);\n")
			fh.write("\t});\n")
			fh.write("})\n")
			fh.write("</script>\n")
			fh.write("<style>\n")
			fh.write("body {margin: 0; padding: 20px;}\n")
			fh.write("table {width: 100%; border-collapse: collapse;}\n")
			fh.write("table, th, td {border: 1px solid #cecece; vertical-align:middle !important}\n")
			fh.write("th {text-align:center; font-size:13px; vertical-align:middle}\n")
			fh.write("td {vertical-align: middle; text-align: center; font-size:12px}\n")
			fh.write("tr.header.expand{cursor:pointer;}\n")
			fh.write(".header .sign:after{content: ''; background-image:url(http://icons.iconarchive.com/icons/icojam/blue-bits/16/math-minus-icon.png); background-size:10px 10px; display:inline-block; width:10px; height:10px; vertical-align:middle;}\n")
			fh.write(".header.expand .sign:after{content:''; background-image:url(http://icons.iconarchive.com/icons/icojam/blue-bits/16/math-add-icon.png); background-size:10px 10px; display:inline-block; width:10px; height:10px;vertical-align:middle; }\n")
			fh.write("</style>\n")
			fh.write("</head>\n<body>\n")
			fh.write("<h2><b>GenoQual Report </b><small>v.1.3</small></h2>\n")
			fh.write("<p>Report generated on %s/%s/%s, %s:%s for run %s</p><br>"%(i.day,i.month,i.year,i.hour,i.minute, args.run_folder))
			fh.write("<div class='table-responsive'>\n")
			fh.write("<table class='table table-striped'>\n")			
			for i, line in enumerate(values):
				if i == 0:
					fh.write("<tr class='header' style='background-color:#006d7a;color: white;'>") 
				elif i == len(values)-1:
					fh.write("<tr class='header' style='background-color:#fff7ea;'>")
						
				if self.config["nextseq"] == True:

					if '_L001_R1_' in line[1]: # Before L001_R1 of a sample, we write the summary line for all the sample's lanes
						sample=line[0].split('_')[0]
						fh.write("<tr class='header expand' {0}>".format(row_color.next()))
						fh.write("<td><span style='vertical-align:middle'>{0} </span><span class=\"sign\"></td>".format(sample))
						fh.write("<td></td>") # empty field for filename
						self.stats_all_lanes[sample]['percent_q30'] = round(100*float(self.stats_all_lanes[sample]['q30'])/self.stats_all_lanes[sample]['total_bases'],2)
						for e in self.stats_all_lanes[sample]:
							try:
								value = int(self.stats_all_lanes[sample][e])
							except ValueError:
								fh.write("<td>{0}</td>".format(self.stats_all_lanes[sample][e]))
							else:
								fh.write("<td>{0:,}</td>".format(value))								
						fh.write("</tr>\n")				
						#Clear cov and tax fields for L001_R1, as they are now in the summary line.
						line[8] = line[9] = '='
						for n in range(10,17):
							line[n] = ''
						for n in range(17, len(line)):
							if line[n] != '':
								line[n] = '='
								
					if i!=0 and i !=len(values)-1:
						fh.write("<tr {0}>".format(row_color.next()))
				else:
					if i!=0 and i !=len(values)-1:
						fh.write("<tr class='header' {0}>".format(row_color.next()))
				
				for j, elem in enumerate(line):
					if i == 0:
						fh.write("<th>")
					else:
						fh.write("<td>")
					if elem=="Total":
						fh.write('<b>')
					if elem and not elem in ["=","skipped"] and j in links and i != 0:
						name = os.path.basename(elem)
						if i==len(values)-1:
							fh.write("<a href='{0}'>MultiQC report</a>".format(elem))
						else:
							fh.write("<a href='{0}'>Report</a>".format(elem))

					else:
						try:
							elem = int(elem)
						except ValueError:
							fh.write("{0}".format(elem))
						else:
							# Use a comma as thousands delimiter
							fh.write("{0:,}".format(elem))
					if i == 0:
						fh.write("</th>")
					else:
						if elem=="Total":
							fh.write('</b>')
						fh.write("</td>")
				fh.write("</tr>\n")
			fh.write("</table>\n</div>\n")
			fh.write("<div class='table-responsive'>\n")
			fh.write("<table>\n")
			fh.write("<tr class='header'>")
			for c,img in enumerate(encoded_figs.keys()):
				if c<4 or c==7:
					if c%2==0 and c!=0:
						fh.write("</tr>\n<tr class='header'>")
					fh.write('<td><img src="data:image/png;base64,%s" class="img-rounded" width="80%%"></td>\n'%encoded_figs[img])
					
					if c==7:
						fh.write("<tr class='header'>")
						fh.write('<td><p style="font-size:18px; font-weight:bold;padding-top: 5px;margin-bottom: 15px;">Historical unweighted unifrac PCoA plot for controls</p>\n'
							 '<iframe src="%s/Qiime/Controls_analysis/bdiv/unweighted_unifrac_emperor_pcoa_plot/index.html" style="margin:0;width:100%%;height:450px">'
							 'Alternative text for browsers that do not understand IFrames.</iframe></td></tr>'%os.path.dirname(filename))			
				if c>=4 and c<7:
					fh.write("<tr class='header'>\n<td colspan=2>\n")					
					fh.write('<p style="text-align:center;"><img src="data:image/png;base64,%s" class="img-rounded" width="95%%"></p>\n'%encoded_figs[img])
					fh.write("</td></tr>\n")
					
			fh.write("</table>\n</div>\n")
			fh.write("</body>\n</html>\n")


class Process(object):
	def __init__(self, cmd):
		self.cmd = cmd

	def run(self, stdout=PIPE, stderr=PIPE, stdin=None, trust_exitcode=True,
	        env=None, piped=False):
		"""Launch or specified command in a subprocess

		:arg trust_exitcode: raises IOError if command exits with non-zero code
		:arg env: dictionary with shell environment to use
		:arg piped: don't wait for command to finish to allow piping output to
		            subsequent commands

		If piped is True trust_exitcode has no effect.
		"""
		log.debug("Launching subprocess %r", self.cmd)

		if isinstance(stdin, basestring):
			_in = PIPE
			input = stdin

			if piped:
				raise TypeError("Cannot mix stdin stdin='string' with "
				                "piped=True")
		else:
			_in = stdin
			input = None

		try:
			self.p = Popen(self.cmd, stdout=stdout, stderr=stderr, stdin=_in,
			               env=env)
		except OSError:
			log.error("Command %r failed to execute. Is the file in the path?",
			          self.cmd)
			raise

		# If stdin is a pipe from another process, close it to allow the first
		# process to receive SIGPIPE if the second process finishes before.
		try:
			name = stdin.name
		except AttributeError:
			pass
		else:
			if name == "<fdopen>":
				stdin.close()

		# Don't wait for process to finish if it's going to be piped with
		# another command.
		if piped:
			return

		out, err = self.p.communicate(input)
		self.finish()

		if trust_exitcode and self.p.returncode != 0:
			raise IOError("Command {0!r} failed with exit code {1}.\n"
			              "Stdout:\n{2}\n Stderr:\n{3}\n"
			              .format(self.cmd, self.p.returncode, out, err))

		log.debug("Subprocess finished with output:\n"
		          "Stdout:\n%s\nStderr:\n%s\n", out, err)

		return out, err

	def finish(self):
		# Wait for processes to finish to avoid zombies
		self.p.wait()


class BaseTool(object):
	"""Base skeleton for external tools or pipeline components.
	Implements a check phase to validate requirements prior to execution
	"""
	def __init__(self, datacfg):
		self.datacfg = datacfg
		self.data = {}

	def prepare(self):
		"""Perform pre-flight checks to ensure all necessary resources are
		available and settings are correctly defined
		"""
		raise NotImplementedError("Subclasses need to override this method")

	def run(self):
		raise NotImplementedError("Subclasses need to override this method")


class QualityStats(BaseTool):
	"""Perform calculation of q20 and q30 statistics
	"""
	def prepare(self):
		log.info("Preparing files for q20 and q30 analysis")
		for entry in self.datacfg.data:
			for file in entry["files"]:
				if file.endswith("I1_001.fastq"):
					continue
				self.data[file] = None

	def run(self):
		self._compute_q(20, 30)

	def _compute_q(self, *args):
		for infile in self.data:
			filename = os.path.basename(infile)

			log.info("Calculating quality scores for file '%s'", filename)
			result = {x: 0 for x in args}
			total_bases = 0
			total_reads = 0
			with open(infile) as fh:
				for read in FastqReader(fh):
					total_bases += read.qual.size
					total_reads += 1
					for val in result:
						result[val] += read.qual[read.qual >= val].size

			# Keep result of calculation
			result["total_bases"] = total_bases
			result["total_reads"] = total_reads
			log.debug("Quality scores for file '%s': %r", filename, result)

			self.data[infile] = result

class MultiQC(BaseTool):
	"""Run MultiQC and write a summary report in the main results folder
	"""
	def __init__(self, *args, **kwargs):
		super(MultiQC, self).__init__(*args, **kwargs)

		self.destdir = self.datacfg.dir

	def prepare(self):
		log.info("Performing pre-flight checks on MultiQC")
		cmd = ("multiqc", "--version")
		p = Process(cmd)
		out, err = p.run()

		if err or not out.startswith("multiqc, version"):
			raise IOError("Unexpected output from MultiQC wrapper.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))

		log.info("MultiQC pre-flight checks finished successfully")

	def run(self):
		log.info("Running MultiQC")
		try:
			cmd = [
			        "multiqc",
			        "--no-data-dir", "-o", self.destdir, self.destdir
			]
			p = Process(cmd)
			out, err = p.run()

		except:
			raise IOError("MultiQC error.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))

class Qiime(BaseTool):

	def __init__(self, *args, **kwargs):
		super(Qiime, self).__init__(*args, **kwargs)

		self.destdir = os.path.join(self.datacfg.dir, "Qiime")
		self.infile = None
		self.infile_rev = None
		self.infile_fwd = None
		self.indexfile = None
		self.metafile = None
		self.users = []
		self.inputdir = None        #Useful only in the merged R1+R2 case
		self.original_index = None  #

	def prepare(self):
		log.info("Performing pre-flight checks on Qiime")
		cmd = ("print_qiime_config.py")
		p = Process(cmd)
		out, err = p.run()

		if err or not out.startswith("\nSystem information"):
			raise IOError("Unexpected output from Qiime wrapper.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))

		cmd = ("biom")
		p = Process(cmd)
		out, err = p.run()

		if err or not out.startswith("Usage:"):
			raise IOError("Unexpected output from Qiime wrapper.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))


		for entry in self.datacfg.data:
			if not entry["Merge_reads"]:
				log.info('QIIME on R1')
				for file in entry["files"]:
					if file.endswith("I1_001.fastq"):
						self.indexfile = file
						self.original_index = file						
					elif file.endswith("R1_001.fastq"):
						self.infile = file
						self.infile_fwd = file
					elif file.endswith("R2_001.fastq"):
						self.infile_rev = file
			else:
				log.info('QIIME on merged')
				self.inputdir = os.path.join(self.datacfg.dir, "Merged_reads", entry["Sample_ID"])
				for file in entry["files"]:
					if file.endswith("R1_001.fastq"):
						self.infile = self.inputdir+'/'+os.path.basename(file+".extendedFrags.fastq")
						self.infile_fwd = file
						log.debug("Original R1 file should be:")
						log.debug(self.infile_fwd)						
						log.debug("Merged file should be:")
						log.debug(self.infile)
					elif file.endswith("I1_001.fastq"):
						self.original_index = file
						self.indexfile = self.inputdir+'/'+os.path.basename(file.replace(".fastq","_fixed.fastq"))
						log.debug("Index file should be:")
						log.debug(self.indexfile)
					elif file.endswith("R2_001.fastq"):						
						self.infile_rev = file
						log.debug("Original R2 file should be:")
						log.debug(self.infile_rev)						
						
		if not os.path.isfile(self.datacfg.config.params["metadata"]):
			raise Error("Metadata file not found at '{0}'"
			            .format(self.datacfg.config.params["metadata"]))
		else:
			self.metafile = self.datacfg.config.params["metadata"]
		
		out, err = self.check_metafile(self.metafile, os.path.join(self.datacfg.dir, "validate_metadata"))
		if err or not out.startswith("No errors"):
			log.error("Error while validating metadata file.")
			log_metadata = open(glob(self.datacfg.dir+"/validate_metadata/*.log")[0]).read()
			log.error(log_metadata)
			sys.exit(1)
		else:
			log.info("Metadata file validated")
			shutil.rmtree(self.datacfg.dir+"/validate_metadata/")
		log.info("Qiime pre-flight checks finished successfully")


	def run(self):
		log.info("Entering QIIME")
		if not os.path.isdir(self.destdir):
			os.mkdir(self.destdir)
			
		#Start the Qiime
		log.info("Writing parameters file")
		with open('%s/params.txt'%self.destdir, 'w') as params:
			params.write("pick_otus:enable_rev_strand_match\tTrue\n")
			params.write("biom-summarize-table:qualitative\tTrue\n")
			params.write("make_emperor:ignore_missing_samples\tTrue\n")
			params.write("beta_diversity_through_plots:ignore_missing_samples\tTrue\n")
			params.write("beta_diversity:metrics\tbray_curtis,euclidean,unweighted_unifrac,weighted_unifrac,binary_jaccard,abund_jaccard\n")
			params.write("alpha_diversity:metrics\tobserved_species,chao1,shannon,PD_whole_tree\n")
		paramfile = self.destdir+'/params.txt'

		if self.inputdir != None:
			log.info("Correcting barcodes in order to match the merged fastq file...")
			self.fix_header_pairs()

		if self.indexfile!=None and self.infile!=None and self.metafile!=None:
			log.debug(self.indexfile)
			log.debug(self.infile)
			log.debug(self.metafile)
			log.debug("Checking metafile for whitespaces...")
			self.metafile_correct(self.metafile)
			
			#Splitting libraries...
			try:
				out, err = self.split_libraries(self.infile, self.indexfile, self.metafile, self.destdir+"/temp_slout")
			except:
				log.info('Splitting failed. Trying to split without reversing barcodes...')
				shutil.rmtree(self.destdir+"/temp_slout")
				out, err = self.split_libraries(self.infile, self.indexfile, self.metafile, self.destdir+"/temp_slout", rev=0)				
			
			#Reads separation by user...
			log.info("Starting separation by user")
			cmd = [
				"split_sequence_file_on_sample_ids.py",
				"-o", self.destdir+"/temp_slout/out_fasta",
				"-i", self.destdir+"/temp_slout/seqs.fna",
				"--file_type", "fasta"]
			p = Process(cmd)
			out, err = p.run()
			
			fastapath = self.destdir+"/temp_slout/out_fasta/"
			fastafiles = glob(fastapath+"*.fasta")	
			for f in fastafiles:
				user = os.path.basename(f).split('.')[0]
				if user not in self.users:
					self.users.append(user)
			for user in self.users:
				user_fastapath = self.destdir+'/'+user+'/slout'
				user_fastqpath = self.destdir+'/'+user+'/fastq'				
				os.makedirs(user_fastapath)
				os.makedirs(user_fastqpath)				
				os.system('cat %s/%s*.fasta > %s/seqs.fna'%(fastapath,user,user_fastapath))
				log.info('QIIME analysis for user %s'%user)
				
				# OTUs picking...
				try:
					out,err = self.pick_otus(self.destdir+"/%s/slout/seqs.fna"%user, self.destdir+"/%s/otus"%user, paramfile)
				except:
					log.info('OTU picking failed for user %s, with error:'%user)
					log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))
					continue
				
				# Biom summarizing...
				try:
					out,err = self.biom_summarize(self.destdir+"/%s/otus/otu_table_mc2_w_tax_no_pynast_failures.biom"%user, 
					                    self.destdir+"/%s/counts.txt"%user)
				except:
					log.info('BIOM summarize failed for user %s, with error:'%user)
					log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))
					continue	
				
				# Finding minimum count for rarefaction...(min:1000)
				min_count=1000
				if os.path.exists(self.destdir+"/%s/counts.txt"%user):
					c = open(self.destdir+"/%s/counts.txt"%user).readlines()
					for line in c:
						if line.startswith('%s.'%user):
							min_count_tmp = int(float(line.split()[1]))
							if min_count_tmp > min_count:
								min_count = min_count_tmp
								break
						
				# Core analyses...
				try:
					out, err = self.core_analyses(self.destdir+"/%s/otus/otu_table_mc2_w_tax_no_pynast_failures.biom"%user,
					                   self.destdir+"/%s/cdout"%user, self.metafile,
					                   self.destdir+"/%s/otus/rep_set.tre"%user, str(min_count), paramfile)
				except:
					log.info('Core analyses failed for user %s, with error:'%user)
					log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))
					continue						

			#Splitting again, without any quality constraint, to give the user the fastq files to repeat
			#the analysis.
			
			try:
				out, err = self.split_libraries(self.infile_fwd, self.indexfile, self.metafile, self.destdir+"/temp_slout_unfiltered", n=999, q=0, fastq=1)			
			except:
				log.info('Splitting failed. Trying to split without reversing barcodes...')
				shutil.rmtree(self.destdir+"/temp_slout_unfiltered")
				out, err = self.split_libraries(self.infile_fwd, self.indexfile, self.metafile, self.destdir+"/temp_slout_unfiltered", n=999, q=0, rev=0, fastq=1)			

				
			if self.infile_rev != None:	#i.e. input is paired-end
				log.info('Pairing reads after splitting into fastq format...')
				self.fix_reads_pairs()
			log.info('Splitting fastq files for users...')			
			self.split_fastq()
			shutil.rmtree(self.destdir+"/temp_slout")
			shutil.copy(self.destdir+"/temp_slout_unfiltered/split_library_log.txt", self.destdir+"/demultiplex_log.txt")			
			shutil.rmtree(self.destdir+"/temp_slout_unfiltered")

		else:
			log.error("Could not find either the reads or the index file:"
			          "Reads:\n{0}\n Index:\n{1}\n".format(self.infile, self.indexfile))
			raise
		
	def split_libraries(self, infile, index, meta, dest, n=0,q=19,rev=1, fastq=0):
		cmd = [
		        "split_libraries_fastq.py",
		        "-o", dest,
		        "-i", infile,
		        "-b", index,
		        "-m", meta,
		        "-n", str(n),
		        "-q", str(q)]

		if rev==1:
			cmd = cmd+["--rev_comp_mapping_barcodes", "--rev_comp_barcode"]
		if fastq==1:
			cmd = cmd+["--store_demultiplexed_fastq"]
			
		p = Process(cmd)
		out, err = p.run()	
		return out, err
	
	def pick_otus(self, i, o, paramfile):
		log.info('Entered otu pick')
		cmd = [
                        "pick_open_reference_otus.py",
		        "-i", i,		        
                        "-o", o,
                        "-p", paramfile]	
		p = Process(cmd)
		out, err = p.run()	
		return out, err	

	def check_metafile(self, i, outdir):
		log.info('Checking metafile...')
		cmd = [
		        "validate_mapping_file.py",
		        "-m", i,		        
		        "-o", outdir]	
		p = Process(cmd)
		out, err = p.run()	
		return out, err	
	
	def biom_summarize(self, i, o):
		cmd = [
                        "biom", "summarize-table",
                        "-i", i,
                        "-o", o]		
		p = Process(cmd)
		out, err = p.run()	
		return out, err	
	
	def core_analyses(self, i, o, meta, tree, min_count, params):
		cmd = [
		        "core_diversity_analyses.py",
		        "-i", i,
		        "-o", o,
		        "-m", meta,
		        "-t", tree,
		        "-e", min_count,
		        "-p", params]
		p = Process(cmd)
		out, err = p.run()
		return out, err	
		
	def split_fastq(self):
		users = []
		records_f = SeqIO.parse(open(self.destdir+"/temp_slout_unfiltered/seqs.fastq", "rU"), "fastq")
		if os.path.exists(self.destdir+"/temp_slout_unfiltered/seqs_rev.fastq"):
			records_r = SeqIO.parse(open(self.destdir+"/temp_slout_unfiltered/seqs_rev.fastq", "rU"), "fastq")
			for (forward, reverse) in izip(records_f,records_r):
				forward.description=' '.join(forward.description.split()[1:3])
				user = forward.id.split('.')[0]
				if not user in users:
					users.append(user)
				sample = '_'.join(forward.id.split('_')[0].split('.')[1:])
				forward.id=forward.description
				forward.description=''
				assert forward.id.split()[0]==reverse.description.split()[0]
				user_fastqpath=self.destdir+'/%s/fastq'%user				
				handleout_f = open(user_fastqpath+"/%s_L001_R1_001.fastq"%sample, "ab")
				handleout_r = open(user_fastqpath+"/%s_L001_R2_001.fastq"%sample, "ab")
				SeqIO.write(forward, handleout_f, 'fastq')
				SeqIO.write(reverse, handleout_r, 'fastq')
				handleout_f.close()
				handleout_r.close()			
			
		else:
			log.debug("Found seqs.fastq only")			
			for forward in records_f:
				forward.description=' '.join(forward.description.split()[1:3])
				user = forward.id.split('.')[0]
				if not user in users:
					users.append(user)				
				sample = '_'.join(forward.id.split('_')[0].split('.')[1:])
				forward.id=forward.description
				forward.description=''				
				user_fastqpath=self.destdir+'/%s/fastq'%user				
				handleout_f = open(user_fastqpath+"/%s_L001_R1_001.fastq"%sample, "ab")
				SeqIO.write(forward, handleout_f, 'fastq')
				handleout_f.close()
				
		counts=open(self.destdir+"/temp_slout_unfiltered/split_library_log.txt").readlines()
		for l in counts:
			user = l.split('.')[0]
			if user in users:
				fileout=open(self.destdir+'/%s/fastq/counts.txt'%user,'a')
				fileout.write(l)
				fileout.close()		
		
		return
				

	def fix_header_pairs(self):
		handle = open(self.infile, "rU")
		headers=[]
		for record in SeqIO.parse(handle, "fastq"):
			headers.append(record.description)
		handle.close()
		headers=set(headers)
		handle = open(self.original_index, "rU")
		handleout = open(self.indexfile, "wb")
		for record in SeqIO.parse(handle, "fastq"):
			if record.description in headers:
				SeqIO.write(record, handleout, 'fastq')
		handle.close()
		handleout.close()
		return
	
	def fix_reads_pairs(self):
		handle = open(self.destdir+"/temp_slout_unfiltered/seqs.fastq", "rU")
		headers=[]
		for record in SeqIO.parse(handle, "fastq"):
			headers.append(record.description.split()[1])
		handle.close()
		headers=set(headers)
		handle = open(self.infile_rev, "rU")
		handleout = open(self.destdir+"/temp_slout_unfiltered/seqs_rev.fastq", "wb")
		for record in SeqIO.parse(handle, "fastq"):
			if record.description.split()[0] in headers:
				SeqIO.write(record, handleout, 'fastq')
		handle.close()
		handleout.close()
		return	
	
	def metafile_correct(self, filename):	
		#Replaces whitespace with dot
		filedata = None
		with open(filename, 'r') as file:
			filedata = file.read()
		# Replace the target string
		filedata = filedata.replace(' ', '.')
		filedata = filedata.replace('_', '.')		
		# Write the file out again
		with open(filename, 'w') as file:
			file.write(filedata)	
		return	
	
class QiimeOnControls(Qiime):
	def __init__(self, *args, **kwargs):
		super(QiimeOnControls, self).__init__(*args, **kwargs)
		self.results_path = self.datacfg.config.params["base_path"]+'results/'
		self.input_path = self.datacfg.config.params["base_path"]+'input/'
		self.current_run = self.datacfg.config.params["run_folder"]

	def run(self):
		log.info('Starting historical QIIME analysis of controls...')
		metaruns=[]
		# Create a big single metadata file for all controls. Duplicate names are avoided by adding 
		# the run ID to the samplename. Duplicate barcodes are not fixed.
		os.makedirs(self.destdir+'/Controls_analysis')
		control_metadata = open(self.destdir+'/Controls_analysis/GEU_metadata.csv', 'w')
		#control_metadata = open('/home/mtruglio/Desktop/GEU_metadata.csv', 'w')
		control_metadata.write('#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tControlType\tDescription\n')
		
		# Also, create a big single seqs.fna file from all GEU runs. The headers will be modified according to the
		# metadata, again in order to avoid duplicate names.
		#control_seqs = open('/home/mtruglio/Desktop/GEU_seqs.fna', 'wb')
		control_seqs = open(self.destdir+'/Controls_analysis/GEU_seqs.fna', 'wb')
		
		for item in os.listdir(self.results_path):
			if os.path.exists(os.path.join(self.results_path, item)+'/Qiime/GEU/slout/seqs.fna'):
				metaruns.append(item)
		metaruns += [self.current_run]
		for item in metaruns:
			metain = open(self.input_path+item+'/metadata.csv').readlines()
			
			for line in metain:
				if line.startswith('GEU'):
					samplename=line.split()[0]+'.'+item.split('_')[1] #make it unique
					if('Pos' in samplename) or ('pos' in samplename):
						controltype='Pos'							
					else:
						controltype='Neg'
						
					if item==self.current_run:
						controltype+='_current'					
					
					control_metadata.write(samplename+'\t'+'\t'.join(line.split()[1:3])+'\t'+controltype+'\t'+samplename+'\n')
			
			if item != self.current_run:
				records_fna = SeqIO.parse(open(os.path.join(self.results_path, item)+'/Qiime/GEU/slout/seqs.fna', "rU"), "fasta")
			else:
				records_fna = SeqIO.parse(open(self.destdir+'/GEU/slout/seqs.fna', "rU"), "fasta")
				
			for r in records_fna:
				entry_n=r.id.split('_')[-1]
				r.id=r.id.split('_')[0]+'.%s'%item.split('_')[1]+'_%s '%entry_n+' '.join(r.description.split()[1:])
				r.description=''
				SeqIO.write(r, control_seqs, 'fasta')
		control_metadata.close()
		control_seqs.close()
		
		try:
			out,err = self.pick_otus(self.destdir+'/Controls_analysis/GEU_seqs.fna', self.destdir+'/Controls_analysis/otus', self.destdir+'/params.txt')
		except:
			log.info('OTU picking failed, with error:')
			log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))
		
		
		try:
			out,err = self.biom_summarize(self.destdir+"/Controls_analysis/otus/otu_table_mc2_w_tax_no_pynast_failures.biom", 
				            self.destdir+"/Controls_analysis/counts.txt")
		except:
			log.info('BIOM summarize failed, with error:')
			log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))
		
		counts=open(self.destdir+"/Controls_analysis/counts.txt").readlines()
		current_run_counts=[]
		for line in counts:
			if line.startswith('GEU') and line.split()[0].split('.')[-1].strip(':') == self.current_run.split('_')[1]:
				c = int(float(line.split()[-1]))
				current_run_counts.append(c)
		min_count=min(current_run_counts)
		
		try:
			out,err = self.beta_diversity(self.destdir+"/Controls_analysis/otus/otu_table_mc2_w_tax_no_pynast_failures.biom",
			                    self.destdir+"/Controls_analysis/bdiv",
			                    self.destdir+"/Controls_analysis/otus/rep_set.tre", 
			                    self.destdir+'/Controls_analysis/GEU_metadata.csv', min_count, 
			                    self.destdir+'/params.txt')
		except:
			log.info('Beta diversity failed, with error:')
			log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))
			
		try:
			out,err = self.taxa_plots(self.destdir+"/Controls_analysis/otus/otu_table_mc2_w_tax_no_pynast_failures.biom",
		                            self.destdir+"/Controls_analysis/taxa_summary",
		                            self.destdir+'/Controls_analysis/GEU_metadata.csv', 
		                            self.destdir+'/params.txt')
		except:
			log.info('Taxa plotting failed, with error:')
			log.info("Stdout:\n{0}\n Stderr:\n{1}\n".format(out, err))			
			
	def beta_diversity(self, i, o, tree, meta, min_count, params):
		if min_count<80:
			min_count=80
		cmd = [
			"beta_diversity_through_plots.py",
			"-i", i,
			"-o", o,
			"-m", meta,
			"-t", tree,
			"-e", str(min_count),
			"-p", params]
		p = Process(cmd)
		out, err = p.run()		
		return out, err	
	
	def taxa_plots(self, i, o, meta, params):
		cmd = [
			"summarize_taxa_through_plots.py",
			"-i", i,
			"-o", o,
			"-m", meta,
			"-p", params]
		p = Process(cmd)
		out, err = p.run()		
		return out, err			
		
class FastQC(BaseTool):
	"""Run FastQC and store the results in folders with names matching the
	input files
	"""
	def __init__(self, *args, **kwargs):
		super(FastQC, self).__init__(*args, **kwargs)

		self.destdir = os.path.join(self.datacfg.dir, "FastQC")

	def prepare(self):
		log.info("Performing pre-flight checks on FastQC")
		cmd = ("fastqc", "-v")
		p = Process(cmd)
		out, err = p.run()

		if err or not out.startswith("FastQC v"):
			raise IOError("Unexpected output from FastQC wrapper.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))

		log.info("FastQC pre-flight checks finished successfully")

		log.info("Preparing files for FastQC analysis")
		for entry in self.datacfg.data:
			# Keep track of the directory where FastQC stored the output
			for file in entry["files"]:
				if file.endswith("I1_001.fastq"):
					continue
				fname = os.path.basename(file)
				report = os.path.join(self.destdir, os.path.splitext(fname)[0] + "_fastqc.html")
				self.data[file] = report

	def run(self):
		if not os.path.isdir(self.destdir):
			os.mkdir(self.destdir)

		for infile in self.data:
			# jcostaDamage
			try:
				fname = os.path.basename(infile)
				log.info("Running FastQC on file '%s'", fname)
				cmd = [
				        "fastqc",
				        "-t", self.datacfg.config.params["args"].threads,
				        "-o", self.destdir,
				        infile
				]
				p = Process(cmd)
				out, err = p.run()

				zip_file = os.path.dirname(self.data[infile]) + ".zip"
				log.debug("Removing FastQC zipped output '%s'", zip_file)
				os.remove(zip_file)
			except:
				continue


class TaxonomyClassifier(BaseTool):
	"""Use a Seqtk and Blast to estimate the taxonomical distribution of the
	sampled reads
	"""
	def __init__(self, *args, **kwargs):
		super(TaxonomyClassifier, self).__init__(*args, **kwargs)

		self.destdir = os.path.join(self.datacfg.dir, "Contaminants")

	def prepare(self):
		log.info("Performing pre-flight checks on Seqtk")
		cmd = ("seqtk",)
		p = Process(cmd)
		out, err = p.run(trust_exitcode=False)

		if "seqtk <command> <arguments>" not in err:
			raise IOError("Unexpected output from Seqtk.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))

		log.info("Seqtk pre-flight checks finished successfully")

		log.info("Performing pre-flight checks on blastn")
		cmd = ("blastn", "-version")
		p = Process(cmd)
		out, err = p.run()

		if err or not out.startswith("blastn: "):
			raise IOError("Unexpected output from blastn.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))

		log.info("blastn pre-flight checks finished successfully")



	def run(self):
		log.info("Checking if all samples have a file with forward reads")
		for entry in self.datacfg.data:
			if entry["Contaminants"] and entry["skip_contaminants"]==False:
				log.debug("Contaminants requested for '%s'", entry["ID"])
				for infile in entry["files"]:
					if not infile.endswith("L001_R1_001.fastq"):
						# R2 should have the same results as R1,other lanes as well
						# cont_file, primary_organism
						self.data[infile] = (None, "=")
						continue

					self.data[infile] = None
			else:
				if entry["skip_contaminants"]==True:
					log.info('Blast search will be skipped for sample {0} due to good alignment'.format(entry["Sample_ID"]))
					for infile in entry["files"]:
						if infile.endswith("L001_R1_001.fastq"):	
							self.data[infile] = ('skipped', os.path.splitext(entry["Reference"])[0])
							
				continue
					
		if not os.path.isdir(self.destdir):
			os.mkdir(self.destdir)

		for file in self.data:
			
			# jcostaDamage (except: at the end..daaaaa)
			try:
				# Skip R2 reads
				if self.data[file] is not None:
					continue

				filename = os.path.basename(file)
				name, ext = os.path.splitext(filename)

				log.info("Running '%s' for taxonomic distribution analysis",
				         filename)

				outsample = os.path.join(self.destdir, name + ".sample")
				# jcostaDamage
				# SAMPLE
				sample = "1000"
				log.debug("Sampling %s sequences from '%s' and writing to '%s'",
				          sample, filename, outsample)

				cmd = ("seqtk", "sample", file, sample)
				p = Process(cmd)

				with open(outsample, 'w') as fh:
					out, err = p.run(stdout=fh)

				outfasta = os.path.join(self.destdir, name + ".fasta")
				log.debug("Converting sample '%s' from FastQ to Fasta format and "
				          "writing to '%s'", outsample, outfasta)
				cmd = ("seqtk", "seq", "-A", outsample)
				p = Process(cmd)
				with open(outfasta, 'w') as fh:
					out, err = p.run(stdout=fh)

					outblast = os.path.join(self.destdir, name + ".txt")
					log.info("BLASTing sampled sequences against nt database. "
					         "This may take a while...")
					log.debug("Sampled file '%s' will be blasted against the nt "
					          "database located at '%s'",
					          outfasta, os.environ["BLASTDB"])

				cmd = (
				        "blastn",
				        "-task", "megablast",
				        "-outfmt", "6 sscinames sseqid",
				        "-ungapped",
				        "-max_target_seqs", "1",  # Keep only the best hit
				        "-max_hsps", "1",         # And only keep 1 HSP per hit
				        "-evalue", "1e-20",
				        "-num_threads", self.datacfg.config.params["args"].threads,
				        "-query", outfasta,
				        "-db", "nt",  # Look up at location set in BLASTDB env var
				)
				p = Process(cmd)

				with open(outblast, 'w') as fh:
					out, err = p.run(stdout=fh, env=os.environ)

				log.info("Binning BLAST results per taxonomic group.")

				taxs = {}
				with open(outblast) as fh:
					for line in fh:
						# First line contains the Scientific taxonomy name
						tax = line.rstrip().split("\t", 1)[0]
						try:
							taxs[tax] += 1
						except KeyError:
							taxs[tax] = 1

				outcounts = os.path.join(self.destdir, name + "_count.txt")
				log.debug("Writing binning result to '%s'", outcounts)

				with open(outcounts, 'w') as fh:
					# Sort by species count in descending order
					sortk = lambda x: x[1]
					total_hits = 0
					highest_hits = None
					for sp, count in sorted(taxs.items(), key=sortk, reverse=True):
						if highest_hits is None:
							highest_hits = sp

						fh.write("{0}\t{1}\n".format(count, sp))
						total_hits += count

					fh.write("{0}\t{1}{2}\n".format(
					        total_hits,
					        "Sequences that hit Blast nt database out of ",
					        sample))

					# If BLAST fails to get even one hit
					if highest_hits is None:
						highest_hits = "No hits"
			except:
				continue

			self.data[file] = (outcounts, highest_hits)


class CoverageAndQuality(BaseTool):
	"""Use Picard tools to count how many sequences align to a given reference
	"""
	def __init__(self, *args, **kwargs):
		super(CoverageAndQuality, self).__init__(*args, **kwargs)
		self.input = {}
		self.destdir = os.path.join(self.datacfg.dir, "Quality")
		self.tmpdir = os.path.join(self.datacfg.dir, "Tmp")
		self.indexdir = os.path.join(self.datacfg.config.params["ref_path"],
		                             "index")
		self.tot_ref_length = 0
		self.nextseq = False

	def prepare(self):
		self.nextseq = self.datacfg.config.params['nextseq']		
		log.info("Performing pre-flight checks on Bwa aligner")
		cmd = ("bwa")
		p = Process(cmd)
		out, err = p.run(trust_exitcode=False)

		if "bwa <command> [options]" not in err:
			raise IOError("Unexpected output from Bwa.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))

		log.info("Performing pre-flight checks on SamTools")

		cmd = ("samtools",)
		p = Process(cmd)
		out, err = p.run(trust_exitcode=False)

		if "Program: samtools" not in err:
			raise IOError("Unexpected output from Samtools.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))

		log.info("Performing pre-flight checks on qualimap")

		cmd = ("qualimap", "-h")
		p = Process(cmd)
		out, err = p.run()

		if "QualiMap v." not in out:
			raise IOError("Unexpected output from Qualimap.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))

		log.info("Preparing output files and confirming references")
		for entry in self.datacfg.data:
			params = self.datacfg.config.params

			log.debug("Checking sample '%s' for reference", entry["Sample_ID"])
			try:
				ref = entry["Reference"]
			except KeyError:
				ref = None

			if not ref:
				log.debug("Got no reference in sample '%s'",
				          entry["Sample_ID"])
				continue

			reference = os.path.join(params["ref_path"], ref)
			index = os.path.join(self.indexdir, ref)


			if not os.path.isfile(reference):
				raise Error("Reference file '{0}' doesn't exist in '{1}'"
				            .format(ref, reference))

			# Reference length, for theoretical coverage
			self.tot_ref_length = 0
			for seq_record in SeqIO.parse(str(reference), "fasta"):
				output_line = '%s\t%i' % (seq_record.id, len(seq_record))
				self.tot_ref_length+=len(seq_record)

			# Use the sample ID as folder location
			file = entry["ID"]
			fname = os.path.basename(file)
			folder = os.path.splitext(fname)[0]
			outputfolder = os.path.join(self.destdir, folder)
			samfile = os.path.join(self.tmpdir, entry["ID"] + ".sam")
			bamfile = os.path.join(self.tmpdir, entry["ID"] + ".bam")
			sampleid = entry["Sample_ID"]
			log.debug("Preparing to run %s against reference %s with index %s and length %s bp",
			          entry["files"], reference, index, self.tot_ref_length)
			log.debug("Using temporary files %s, %s and final files in %s",
			          samfile, bamfile, outputfolder)

			files=()

			for file in entry["files"]:
				if not (file.endswith("R1_001.fastq") or file.endswith("R2_001.fastq")):
					continue
				files = files + (file,)
			self.input[files] = (reference, index, samfile, bamfile,
			                     outputfolder, sampleid)


	def run(self):
		if not os.path.isdir(self.destdir):
			os.mkdir(self.destdir)

		if not os.path.isdir(self.indexdir):
			os.mkdir(self.indexdir)

		if not os.path.isdir(self.tmpdir):
			os.mkdir(self.tmpdir)

		for files in self.input:
			try:
				files_str = ", ".join(map(os.path.basename, files))
				if self.input[files] is None:
					log.info("Skipping QC of %s - missing a reference",
					         files_str)
					continue

				ref, index, sam, bam, folder,sampleid = self.input[files]
				ref_str = os.path.basename(ref)

				log.info("Running %s for quality control", files_str)

				index_file = index + ".bwt"
				if not os.path.isfile(index_file):
					log.info("Indexing reference '%s'", ref_str)

					cmd = ("bwa", "index", ref, "-p", index)
					p = Process(cmd)
					p.run()

				log.info("Aligning files %s against reference %s", files_str,
				         ref_str)
				
				if self.nextseq == True:
					nextseq_files = ()
					nextseq_files = nextseq_files + (os.path.join(self.datacfg.dir,'Merged_lanes/' + sampleid + "_R1.fastq"),)
					if len(files) == 8:
						nextseq_files = nextseq_files + (os.path.join(self.datacfg.dir,'Merged_lanes/' + sampleid + "_R2.fastq"),)
					cmd = (
						"bwa",
						"mem",
						"-t", self.datacfg.config.params["args"].threads,
						index) + nextseq_files
					p = Process(cmd)					
					
				else:
					cmd = (
						"bwa",
						"mem",
						"-t", self.datacfg.config.params["args"].threads,
						index) + files
					p = Process(cmd)

				with open(sam, 'w') as fh:
					p.run(stdout=fh)

				log.info("Converting formats SAM -> BAM and sorting")
				log.debug("SAM: '%s' -> BAM '%s'", sam, bam)
				os.popen("samtools view -@ {2} -Sb {0} > {1}".format(sam, bam.replace('.bam', '_unsorted.bam'), self.datacfg.config.params["args"].threads))
				os.popen("samtools sort -@ {2} -o {1} {0} ".format(bam.replace('.bam', '_unsorted.bam'), bam, self.datacfg.config.params["args"].threads))
				#cmd1 = ("samtools", "view", "-Sb", sam)
				#cmd2 = ("samtools", "sort", "-", os.path.splitext(bam)[0])
				#p1 = Process(cmd1)
				#p1.run()
				#p2 = Process(cmd2)
				#p2.run(stdin=p1.p.stdout)
				

				# Avoid zombies on the first process
				#p1.finish()

				log.debug("Removing SAM file")

				# Remove SAM file to reduce space
				os.remove(sam)
				os.remove(bam.replace('.bam', '_unsorted.bam'))

				log.info("Indexing BAM file")
				log.debug("Indexed file: '%s'", bam)

				cmd = ("samtools", "index", bam)
				p = Process(cmd)
				p.run()

				log.info("Running Qualimap on %s", files_str)
				cmd = (
				        "qualimap",
				        "bamqc",
				        "-nt", self.datacfg.config.params["args"].threads,
				        "-bam",
				        bam,
				        "-outdir", folder,
				)
				p = Process(cmd)
				# Before was: out, err = p.run()
				p.run()

				log.debug("Removing BAM files")

				# Remove BAM file and index
				for f in glob(bam + "*"):
					os.remove(f)

				log.info("Parsing Qualimap output for number of mapped bases")

				stats = {}

				with open(os.path.join(folder, "genome_results.txt")) as fh:
					for l in fh:
						l = l.strip()

						if l.startswith("number of reads ="):
							stats["total_reads"] = l.split()[-1].replace(",", "")
							continue

						if l.startswith("number of mapped reads ="):
							stats["mapped_reads"] = l.split()[-2].replace(",", "")
							continue

						if l.startswith("number of mapped bases ="):
							stats["mapped_bases"] = l.split()[-2].replace(",", "")
							continue

						if l.startswith("mean coverageData ="):
							stats["mean_coverage"] = l.split()[-1].replace("X", "")
							continue

				perc_mapped_reads = (
				        int(stats["mapped_reads"]) / int(stats["total_reads"]) * 100
				)
				log.info(perc_mapped_reads)
				if perc_mapped_reads >= 80:
					for entry in self.datacfg.data:
						if entry["Sample_ID"]==os.path.splitext(os.path.basename(sam))[0] and entry["Contaminants"]:
							log.info('More than 80% of the reads mapped to the given reference genome; contaminants analysis will be skipped.')							
							entry['skip_contaminants']=True
				self.data[files[0]] = {
				        "total_mapped_reads": "{0} ({1:.2f}%)".format(
				                stats["mapped_reads"], perc_mapped_reads),
				        "total_mapped_bases": stats["mapped_bases"],
				        "total_reads_used": stats["total_reads"],
				        "mean_coverage": stats["mean_coverage"],
				        "ref_length": self.tot_ref_length,
				        "qualimap": os.path.join(folder, "qualimapReport.html"),
				}

				if len(files) == 2:
					self.data[files[1]] = {
					        "total_mapped_reads": "=",
					        "total_mapped_bases": "=",
					        "total_reads_used": "=",
					        "mean_coverage": "=",
					        "ref_length": self.tot_ref_length,
					        "qualimap": "=",
					}
			except:
				self.data[files[0]] = {
				        "total_mapped_reads": "ERROR",
				        "total_mapped_bases": "ERROR",
				        "total_reads_used": "ERROR",
				        "mean_coverage": "ERROR",
				        "ref_length": self.tot_ref_length,
				        "qualimap": "ERROR",
				}

				if len(files) == 2:
					self.data[files[1]] = {
					        "total_mapped_reads": "=",
					        "total_mapped_bases": "=",
					        "total_reads_used": "=",
					        "mean_coverage": "=",
					        "ref_length": self.tot_ref_length,
					        "qualimap": "=",
					}
				continue

		os.rmdir(self.tmpdir)


class Demultiplexer(BaseTool):
	"""Use the first bases on each sequence and information contained in the
	demux_map file to split different sequences to each output.
	"""
	def __init__(self, *args, **kwargs):
		super(Demultiplexer, self).__init__(*args, **kwargs)

		self.destdir = os.path.join(self.datacfg.dir, "Demultiplex")

	def _parse_demux_file(self, filename, filepath):
		data = {}

		samplepath = os.path.join(self.destdir, os.path.splitext(filename)[0])
		os.makedirs(samplepath)

		# Add 2 files for reads that don't match any tag
		rejected1 = os.path.join(samplepath, "rejected_1.fastq")
		rejected2 = os.path.join(samplepath, "rejected_2.fastq")
		data["rejected"] = (rejected1, rejected2)
		log.debug("Rejected demux files %s", data["rejected"])

		# Extra information useful during and after demuxing
		data["options"] = {}
		data["options"]["samplepath"] = samplepath
		stats = data["options"]["read_counts"] = {}
		stats["rejected"] = 0
		samples = data["options"]["samples"] = {}

		# NOTE Keep track of files for following pipeline steps (mergereads)
		files = []

		with open(filepath) as fh:
			header = None

			prev_tagsize = None

			for line_no, line in enumerate(csv.reader(fh)):
				if header is None:
					header = line

					expected = 7
					cols = len(header)
					if cols != expected:
						raise Error("Demux file {0} has an unexpected {1} "
						            "number of columns. Expected {2}".format(
						                    filename, cols, expected))
					continue

				sample_id = line[0]

				if not sample_id:
					raise Error("Demux file {0} has invalid data. "
					            "Sample column cannot be empty - line {1}."
					            .format(filename, line_no + 1))

				sample1 = os.path.join(samplepath, sample_id + "_R1.fastq")
				sample2 = os.path.join(samplepath, sample_id + "_R2.fastq")
				log.debug("Demux files will be %s and %s", sample1, sample2)

				strandF, tagF, seqF, strandR, tagR, seqR = line[1:]

				# Check that indeed strands are in the correct order
				if not (strandF == "F" and strandR == "R"):
					raise Error("Demux file {0} should have Forward "
					            "information on columns 2-4 and Reverse "
					            "information on columns 5-7".format(filename))

				if len(seqF) == len(seqR):
					if prev_tagsize is None:
						pass

					elif len(seqF) != prev_tagsize:
						log.debug("Tag previous tag had size %s while current "
						          "has %s", prev_tagsize, len(seqF))
						raise Error("Tags on file {0} have different lengths. "
						            "current algorithm requires same size "
						            "tags".format(filename))

					prev_tagsize = len(seqF)

				else:
					log.debug("Tag forward and tag reverse are %s and %s in "
					          "size", len(seqF), len(seqR))
					raise Error("Forward and reverse tags on file {0} have "
					            "different lengths. Current algorithm "
					            "requires same size tags".format(filename))

				try:
					data[seqF][seqR] = (sample1, sample2)
					stats[seqF][seqR] = 0
				except KeyError:
					data[seqF] = {seqR: (sample1, sample2)}
					stats[seqF] = {seqR: 0}

				if sample_id in samples:
					raise Error("Duplicate sample_id '{0}' on demux file {1} "
					            ". Sample_id column must be unique and "
					            "non-empty.".format(sample_id, filename))

				samples[sample_id] = (seqF, seqR)

				files.append((sample1, sample2))

			if prev_tagsize is None:
				raise Error("No tag found on file {0}".format(filename))

		log.info("Tag size to use in algorithm is %s", prev_tagsize)
		data["options"]["tagsize"] = prev_tagsize

		return files, data

	def _open_filehandles(self, demux):
		"""Open files for demultiplexing. This method avoids having too many
		open files
		"""
		for k1 in demux.keys():
			if k1 == "rejected":
				demux[k1] = (open(demux[k1][0], 'w'),
				             open(demux[k1][1], 'w'))

			elif k1 == "options":
				pass

			else:
				for k2 in demux[k1].keys():
					demux[k1][k2] = (open(demux[k1][k2][0], 'w'),
					                 open(demux[k1][k2][1], 'w'))

	def _close_filehandles(self, demux):
		"""Close files for demultiplexing. This method avoids having too many
		open files
		"""
		for k1 in demux.keys():
			if k1 == "rejected":
				for fh in demux[k1]:
					fh.close()
				demux[k1] = map(lambda x: x.name, demux[k1])

			elif k1 == "options":
				pass

			else:
				for k2 in demux[k1].keys():
					for fh in demux[k1][k2]:
						fh.close()
					demux[k1][k2] = map(lambda x: x.name, demux[k1][k2])

	def _reject_seq(self, rejected, seq1, seq2):
		"""Used to write Fastq sequences to the rejected/unclassified files
		"""
		seq1.write_to_fastq_file(rejected[0])
		seq2.write_to_fastq_file(rejected[1])

	def prepare(self):
		for entry in self.datacfg.data:
			if not entry["Demux"]:
				continue

			demuxfile = entry["Demux"]
			demuxpath = os.path.join(self.datacfg.config.params["input_path"],
			                         demuxfile)
			log.debug("Using Demux at %s for sample %s", demuxpath,
			          entry["Sample_ID"])

			if not os.path.isfile(demuxpath):
				raise Error("Demux file {0} couldn't be found at {1}".format(
				        demuxfile, self.datacfg.config.params["input_path"]))

			log.info("Parsing demux file %s", demuxfile)

			for infile in entry["files"]:
				if not infile.endswith("R1_001.fastq"):
					self.data[infile] = "="
					continue

				demux_files, demux = self._parse_demux_file(demuxfile,
				                                            demuxpath)
				entry["demux_files"] = demux_files

				self.data[infile] = (entry["files"], demux)

	def run(self):
		for file in self.data:
			if self.data[file] == "=":
				continue

			files, demux = self.data[file]
			sample = os.path.basename(file).replace("_L001_R1_001.fastq", '')

			log.debug("Opening filehandles for demux of sample %s", sample)
			self._open_filehandles(demux)

			r1, r2 = files
			tagsize = demux["options"]["tagsize"]
			read_counts = demux["options"]["read_counts"]

			total_counts = {"classified": 0,
			                "unclassified": 0}

			with nested(open(r1), open(r2)) as (fh1, fh2):
				f1 = iter(FastqReader(fh1))
				f2 = iter(FastqReader(fh2))

				while True:
					try:
						seq1 = next(f1)
						seq2 = next(f2)
					except StopIteration:
						break

					key1 = seq1.seq[:tagsize]
					key2 = seq2.seq[:tagsize]

					match1 = get_by_distance_from_dict(demux, key1)

					if match1 is None:
						self._reject_seq(demux["rejected"], seq1, seq2)
						total_counts["unclassified"] += 1
						read_counts["rejected"] += 1
						continue

					match2 = get_by_distance_from_dict(demux[match1], key2)

					if match2 is None:
						self._reject_seq(demux["rejected"], seq1, seq2)
						total_counts["unclassified"] += 1
						read_counts["rejected"] += 1
						continue

					files = demux[match1][match2]

					# Remove the tag sequence from final output
					new_seq1 = seq1[tagsize:]
					new_seq2 = seq2[tagsize:]
					# Use the original name. [part] is appended when slicing
					new_seq1.name = seq1.name
					new_seq2.name = seq2.name

					new_seq1.write_to_fastq_file(files[0])
					new_seq2.write_to_fastq_file(files[1])

					total_counts["classified"] += 1
					read_counts[match1][match2] += 1

			log.debug("Results of demux are %s classified and %s unclassified",
			          total_counts["classified"], total_counts["unclassified"])
			demux["options"]["total_counts"] = total_counts

			log.debug("Closing filehandles for demux of sample %s", sample)
			self._close_filehandles(demux)


class MergeReads(BaseTool):
	"""Use FLASH to fuse the paired reads into a single read.

	This process only makes sense if the dataset is 16S or if the sequenced
	fragments are expected to be smaller than ~500bp (300+300-overlap)
	"""
	def __init__(self, *args, **kwargs):
		super(MergeReads, self).__init__(*args, **kwargs)
		self.nextseq = False		
		self.destdir = os.path.join(self.datacfg.dir, "Merged_reads")

	def prepare(self):
		log.info("Performing pre-flight checks on Flash")
		self.nextseq = self.datacfg.config.params['nextseq']
		cmd = ("flash", "--version")
		p = Process(cmd)
		out, err = p.run(trust_exitcode=False)

		if "FLASH v1.2.11" not in out:
			raise IOError("Unexpected output from Flash.\n"
			              "Stdout:\n{0}\n Stderr:\n{1}\n"
			              .format(out, err))

		log.info("Flash pre-flight checks finished successfully")
			
		for entry in self.datacfg.data:
			if not entry["Merge_reads"]:
				continue
			#if self.nextseq:
			#	pass
			#else:
			infile = None
			indexfile = None
			for file in entry["files"]:
				if file.endswith("L001_R1_001.fastq"):
						infile = file
				elif file.endswith("I1_001.fastq"):
					indexfile = file
					print(indexfile)
				else:
					self.data[file] = "="	
					
			if not entry["Demux"]:
				if self.nextseq:
					merged_lanes_dir = os.path.join(self.datacfg.dir,'Merged_lanes/')
					self.data[infile] = [tuple([merged_lanes_dir + entry["Sample_ID"] + "_R1.fastq", merged_lanes_dir + entry["Sample_ID"] + "_R2.fastq"])]
				else:
					self.data[infile] = [entry["files"]]
				if indexfile!=None:
					self.data[infile][0] = tuple(x for x in self.data[infile][0] if x!=indexfile)
			else:
				self.data[infile] = entry["demux_files"]
				
			print(infile, self.data[infile])

			final_files = []

			for f1, f2 in self.data[infile]:
				if not entry["Demux"] and not self.nextseq:
					prefix = os.path.basename(f1).replace("_R1_001.fastq", "")
				else:
					prefix = os.path.basename(f1).replace("_R1.fastq", "")						
				outputdir = os.path.join(self.destdir, entry["Sample_ID"])
				final_files.append((f1, f2, outputdir, prefix))

			self.data[infile] = final_files

	def run(self):
		for file in self.data:
			if self.data[file] == "=":
				continue

			results = {}

			for f1, f2, outdir, sample_id in self.data[file]:
				if not os.path.isdir(outdir):
					os.makedirs(outdir)

				logfile = os.path.join(outdir, sample_id + "_flash.log")

				log.info("Running Flash on sample files %s %s", f1, f2)
				cmd = (
				        "flash",
				        "-t", "1",
				        "-o", sample_id,
				        "-d", outdir,
				        "-M", "250",
				        f1, f2)
				p = Process(cmd)

				with open(logfile, 'w') as fh:
					p.run(stdout=fh)

				log.info("Collecting FLASH output files")
				merged_reads = os.path.join(outdir,
				                            sample_id + ".extendedFrags.fastq")
				failed_f1 = os.path.join(outdir,
				                         sample_id + ".notCombined_1.fastq")
				failed_f2 = os.path.join(outdir,
				                         sample_id + ".notCombined_2.fastq")
				histogram = os.path.join(outdir, sample_id + ".histogram")
				hist = os.path.join(outdir, sample_id + ".hist")

				# Parse flash's logfile to retrieve counts
				total_reads_combined = 0
				total_reads_notcombined = 0

				with open(logfile) as fh:
					for line in fh:
						if "Combined pairs:" in line:
							total_reads_combined = int(line.split()[-1])
						elif "Uncombined pairs:" in line:
							total_reads_notcombined = int(line.split()[-1])

				results[sample_id] = {
				        "merged_reads": merged_reads,
				        "not_combined_R1": failed_f1,
				        "not_combined_R2": failed_f2,
				        "hist": hist,
				        "histogram": histogram,
				        "logfile": logfile,
				        "reads_combined": total_reads_combined,
				        "reads_notcombined": total_reads_notcombined,
				}

			self.data[file] = results


class Pipeline(object):
	"""Sequentially execute the steps of the pipeline, performing checks first
	and preparing all steps and only afterwards executing each component
	"""
	def __init__(self, config):
		self.config = config

		self._steps = []
		self.step("data", DataWrapper(self.config))		
		self.step("qc", FastQC(self.data))
		self.step("qual", QualityStats(self.data))
		self.step("cov", CoverageAndQuality(self.data))		
		self.step("tax", TaxonomyClassifier(self.data))
		self.step("multiqc", MultiQC(self.data))
		self.step("dmux", Demultiplexer(self.data))
		self.step("merge", MergeReads(self.data))
		if self.config.params["args"].meta:
			self.step("meta", Qiime(self.data))
			self.step("metaControls", QiimeOnControls(self.data))

	def step(self, name, obj):
		setattr(self, name, obj)
		self._steps.append(obj)

	def prepare(self):
		try:
			for step in self._steps:
				step.prepare()
		except:
			log.fatal("Cleaning up after fatal error during prepare phase...")
			self.cleanup()
			raise

	def run(self):
		try:
			for step in self._steps:
				try:
					step.run()
				except:
					continue
		except:
			log.fatal("Cleaning up after fatal error during run phase...")
			self.cleanup()
			raise

	def cleanup(self):
		log.debug("Cleaning up Pipeline leftovers")
		self.data.cleanup()

	def finish(self):
		"""Export the results as HTML and move all result files to their final
		destination.
		"""
		log.debug("Removing .fastq files used in the analysis")
		for file in glob(os.path.join(self.data.dir, "*.fastq")):
			os.remove(file)

		log.debug("Preparing data for reporting")
		data = DataCollector(self)

		log.debug("Generating Tabular report")
		out = DataFormatter(data,self.config.params)
		out.write_tabular(os.path.join(self.data.dir,
		                               self.config.params["report_tab"]))

		#log.debug("Generating HTML report")
		#out.write_html(os.path.join(self.data.dir,
		#                            self.config.params["report_html"]))

		log.info("Moving all files to final destination")
		shutil.move(self.data.dir, self.config.params["output_path"])

		sleep(5)
		log.debug("Updating runs and samples db")
		os.system("%s/build_db.py %s"%(os.path.dirname(os.path.realpath(__file__)), os.path.join(self.config.params["base_path"], self.config.params["output_folder"])))

		log.debug("Creating images for quality history")
		os.system("%s/graphs.py -o %s/tmp_images -d %s/runs_samples.db -r %s"%(os.path.dirname(os.path.realpath(__file__)), self.config.params["output_path"],
		                                                                       os.path.join(self.config.params["base_path"], self.config.params["output_folder"]), args.run_folder))

		log.debug("Generating HTML report")
		log.info('Joining paths: %s %s'%(self.config.params["output_path"], self.config.params["report_html"] ))
		out.write_html(os.path.join(self.config.params["output_path"],
		                            self.config.params["report_html"]))


		log.debug("Writing report to requested output location")
		extra_output = self.config.params["args"].output
		if extra_output:
			shutil.copy(self.config.params["output_file_html"], extra_output)

		log.debug("Correcting AFS permissions for POSIX strict clients")
		for root, _, files in os.walk(self.config.params["output_path"]):
			# Folder permissions
			os.chmod(root, 0777)
			# File permissions
			for file in files:
				os.chmod(os.path.join(root, file), 0666)

		log.info("Pipeline finished successfully")
		log.info("Results can be found at '%s'",
		         self.config.params["output_path"])


def parse_args():
	parser = argparse.ArgumentParser(description="Perform Quality Control on"
	                                 " Illumina reads")
	parser.add_argument('run_folder',
	                    help="Folder suffix where files with reads are "
	                    "located. E.g. 16 if folder name is Run_16")
	parser.add_argument('--base_path',
	                    help="Location of input/reference/results folders",
	                    #default="/home/mtruglio/Desktop/Test_UBI_folder/") #Mauro: correggi qui
	                    default="/afs/igc.gulbenkian.pt/folders/gen-com/USERS/UBI/")
	                    
	parser.add_argument('--output', '-o',
	                    help="Also store a copy of the HTML report in given "
	                    "location")
	parser.add_argument('--meta', '-m',
	                    help="Perform basic QIIME run",
	                    action='store_true')
	parser.add_argument('--blastdb',
	                    help="Location of blast databases (to be set as the "
	                    "env variable BLASTDB",
	                    default="/home/opt/blast/ncbi-blast-2.4.0+/db") # Mauro correggi qui
	parser.add_argument('--url', '-u',
	                    help="Url to use in HTML output",
	                    default="minsk.igc.gulbenkian.pt:8080/static/genoqual_results/") # QUI
	parser.add_argument('--threads', '-t',
	                    help="Number of threads to use in subprocesses",
	                    default='10')
						#default=os.environ.get("GALAXY_SLOTS", "8"))
	parser.add_argument('--verbose', '-v', action="count",
	                    help="Verbosity level. -vvv is the highest level")
	return parser.parse_args()


def prepare_settings(args):
	d = {}

	d["base_path"] = args.base_path
	d["config_name"] = "QCconf.csv"
	d["columns"] = ("Sample_ID", "Description", "Contaminants", "Reference",
	                "Demux", "Merge_reads")
	d["booleans"] = {"Contaminants", "Merge_reads"}

	d["input_folder"] = "input"
	d["output_folder"] = "results"
	d["meta"] = args.meta
	d["nextseq"] = False
	d["ref_folder"] = "references"
	d["report_html"] = "report.html"
	d["report_tab"] = "report.txt"
	d["run_folder"] = "Run_{0}".format(args.run_folder)

	d["input_path"] = os.path.join(d["base_path"], d["input_folder"],
	                               d["run_folder"])
	d["output_path"] = os.path.join(d["base_path"], d["output_folder"],
	                                d["run_folder"])
	d["output_url"] = (
	        "http://" +
	        urllib.quote(urlparse.urljoin(args.url+d["run_folder"],''), ':/')
#        urllib.quote(urlparse.urljoin(args.url, d["run_folder"]))
	)
	d["ref_path"] = os.path.join(d["base_path"], d["ref_folder"])
	d["config"] = os.path.join(d["input_path"], d["config_name"])
	if args.meta:
		d["metadata"] = os.path.join(d["input_path"], "metadata.csv")
	d["output_file_html"] = os.path.join(d["output_path"], d["report_html"])
	d["output_file_tab"] = os.path.join(d["output_path"], d["report_tab"])
	d["args"] = args

	# Also prepare environment settings that are required by some tools
	os.environ["BLASTDB"] = args.blastdb

	log.debug("Prepared settings are %s", d)

	return d


def main(args):
	# Setting things up
	settings = prepare_settings(args)
	config = ConfigParser(settings)
	# Pre-flight check on all steps
	pp = Pipeline(config)
	pp.prepare()
	# Actually compute stuff
	pp.run()
	# Generate HTML report with links to files and copy files to final location
	pp.finish()


if __name__ == "__main__":
	args = parse_args()
	if args.verbose == 1:
		level = logging.WARN
	elif args.verbose == 2:
		level = logging.INFO
	elif args.verbose >= 3:
		level = logging.DEBUG
	else:
		level = logging.ERROR

	logging.basicConfig(
	        format="%(asctime)s - %(levelname)s - %(message)s",
	        level=level,
	)
	log = logging.getLogger(__name__)

	try:
		main(args)
	except Error as e:
		log.fatal(e)
		sys.exit(1)
	except Exception as e:
		log.exception(e)
		sys.exit(1)
	

# vim: ai sts=4 et sw=4
