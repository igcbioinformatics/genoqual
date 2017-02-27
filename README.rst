This pipeline is composed of the following steps:

* Accept 1 parameter (the number that identifies the RUN and matches the folder name)
* Read config file and perform pre-flight checks
    * Make sure all files exist and all necessary information is present
* Quality statistics on all read files
    * Compute q20 and q30 (number of bases above given quality value)
* Run FastQC on all fastq files.
    * Can we get q20 and q30 from FastQC directly? (one less go through all files)
* Sample each file (R1 reads only as R2 are the same sequence in opposite direction)
* Convert fastq files to fasta (strip quality)
* Blast (megablast) the sample against nt (include the taxonomy identifier in the output)
* Count the hits per taxonomy group or higher taxonomy level
* Export results in csv or tab delimited format as well as HTML format for readability

* If a reference is given, perform several measurements of quality and coverage
