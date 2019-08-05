:: This batch script runs CoNet on a WINDOWS node/server. Note that the 
:: Windows implementation does not allow multithreading. A corresponding
:: LINUX script is available in the "demo" folder of the CoNet download. 
:: 
:: Four steps are required to run the full CoNet procedure. However, each 
:: step generates a output file that can be loaded in order to repeat single
:: steps without running the fill pipeline. Set the task/step to false (or 
:: anything else but true) to skip it.  


@echo off

:::: -- set all necessary paths 

:: point to the CoNet jar file 
SET path_to_conet_jar=".\lib\CoNet.jar" 
:: point to the java distibution (I found that Java10 on windows no longer supports CoNet)
SET java18="C:\Program Files\Java\jdk1.8.0_162\bin\java.exe"

:: Input folder
set inputfolder=.\R\Data
:: Output folder - this will also be set as working directory
set outputfolder=.\R\Results\CoNet
:: input file name
set inputfile=MT-QIIME-RC-1000.txt
:: threshold file name 
set th_file=%outputfolder%\th_file-intersect-1000-R.C.txt

:: output file identifier (+step +merge_by +id)
set file=RC-1000
:: run ID
:: set /p id=What file ID? 
set id=001 




:::: -- activate tasks/steps: Set the task/step to false (or anything else but true) to skip it.

:: compute inital network 
set initial_network=true
:: compute permutation
set permutate=true
:: compute bootstraping
set bootstraps=true
:: compute final network 
set final_network=true



:::: -- computational settings:

:: general
set d_type=count

:: preprocessing
:: is done beforehead

:: Methods 
set methods=correl_pearson/correl_spearman/dist_bray/dist_kullbackleibler

::  manually chosen threshold will be used only if th_file path leads to a dead end
::  in that case a waring will be send.
set pearson=0.9
set spearman=0.9
set kw=0.01
set bc=0.01

:: Merge settings
set merge_by=intersection
set minsupport=3

:: permutation and bootstrapping
set edgethreshold=0.05
set iterations=100
set multicorr=none


:: set working directory
cd %outputfolder%
:: reset id
set set id=%merge_by%-%id%

:: if not yet there, create the config files
if not exist %inputfolder%/CoNet-Config.txt (
	copy nul > %inputfolder%/CoNet-Config.txt
	echo NOTIFICATION: did not find a config file. Created an empty one instead.
	pause
	)

:: if thresholds file exists, use thresholds file
set use_th_file=false
if exist %th_file% (
	set use_th_file=true
	echo NOTIFICATION: will use thresholds from file %th_file%
	pause
	)


:: FIRST STEP: calculate initial network 
if %initial_network% equ true (
	if %use_th_file% equ true ( 
		%java18% -Xmx%memory%m -cp %path_to_conet_jar% be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser ^
		-i %inputfolder%\%inputfile% -o %outputfolder%\inital-network-%file%-%id%.gdl --inputformat biomtable -t %d_type% -f gdl ^
		--filter row_minsum/rand --filterparameter 1.0 --stand col_norm --nantreatmentparam 1 --nantreatment none ^
		--method ensemble -E %methods% --ensembleparamfile %th_file% ^
		--multigraph --scoremergestrategy mean --networkmergestrategy %merge_by% --minsupport %minsupport% -Z %inputfolder%\CoNet-Config.txt
	) else ( 
		%java18% -Xmx%memory%m -cp %path_to_conet_jar% be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser ^
		-i %inputfolder%\%inputfile% -o %outputfolder%\inital-network-%file%-%id%.gdl --inputformat biomtable -t %d_type% -f gdl ^
		--filter row_minsum/rand --filterparameter 1.0 --stand col_norm --nantreatmentparam 1 --nantreatment none ^
		--method ensemble -E %methods% --ensembleparams sim_varlogratio~lowerThreshold=%logvar%\correl_pearson~upperThreshold=%pearson%\correl_pearson~lowerThreshold=-%pearson%/correl_spearman~upperThreshold=%spearman%/correl_spearman~lowerThreshold=-%spearman%/dist_bray~lowerThreshold=%bc%/dist_kullbackleibler~lowerThreshold=%kw% ^
		--multigraph --scoremergestrategy mean --networkmergestrategy %merge_by% --minsupport %minsupport% -Z %inputfolder%\CoNet-Config.txt
	)
)


:: SECOND STEP: estimate null distribution (permutation)
if %permutate% equ true (
	if %use_th_file% equ true (
		%java18% -Xmx%memory%m -cp %path_to_conet_jar% be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser ^
		-i %inputfolder%\%inputfile% -o %outputfolder%\permut-network-%file%-%id%.gdl --inputformat biomtable -t %d_type% -f gdl ^
		--filter row_minsum/rand --filterparameter 1.0 --stand col_norm --nantreatmentparam 1 --nantreatment none ^
		--method ensemble -E %methods% --ensembleparamfile %th_file% ^
		--multigraph --scoremergestrategy mean --networkmergestrategy %merge_by% --minsupport %minsupport% ^
		--pvaluemerge max --renorm --edgethreshold 0.05 --iterations %iterations% --resamplemethod shuffle_rows -I %outputfolder%\permut-%file%-%id%.txt -K edgeScores --renorm --scoreexport ^
		-Z %inputfolder%\CoNet-Config.txt
	) else (
		%java18% -Xmx%memory%m -cp %path_to_conet_jar% be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser ^
		-i %inputfolder%\%inputfile% -o %outputfolder%\permut-network-%file%-%id%.gdl --inputformat biomtable -t %d_type% -f gdl ^
		--filter row_minsum/rand --filterparameter 1.0 --stand col_norm --nantreatmentparam 1 --nantreatment none ^
		--method ensemble -E %methods% --ensembleparams sim_varlogratio~lowerThreshold=%logvar%/correl_pearson~upperThreshold=%pearson%/correl_pearson~lowerThreshold=-%pearson%/correl_spearman~upperThreshold=%spearman%/correl_spearman~lowerThreshold=-%spearman%/dist_bray~lowerThreshold=%bc%/dist_kullbackleibler~lowerThreshold=%kw% ^
		--multigraph --scoremergestrategy mean --networkmergestrategy %merge_by% --minsupport %minsupport% ^
		--pvaluemerge max --renorm --edgethreshold 0.05 --iterations %iterations% --resamplemethod shuffle_rows -I %outputfolder%\permut-%file%-%id%.txt -K edgeScores --renorm --scoreexport ^
		-Z %inputfolder%\CoNet-Config.txt
	)
)


:: THIRD STEP: estimate confidence (bootstraping) 
if %bootstraps% equ true (
	if %use_th_file% equ true (
		%java18% -Xmx%memory%m -cp %path_to_conet_jar% be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser ^
		-i %inputfolder%\%inputfile% -o %outputfolder%\boot-network-%file%-%id%.gdl --inputformat biomtable -t %d_type% -f gdl ^
		--filter row_minsum/rand --filterparameter 1.0 --stand col_norm --nantreatmentparam 1 --nantreatment none ^
		--method ensemble -E %methods% --ensembleparamfile %th_file% ^
		--multigraph --scoremergestrategy mean --networkmergestrategy %merge_by% --minsupport %minsupport% ^
		--pvaluemerge max --edgethreshold 0.05 --iterations %iterations% --resamplemethod bootstrap -I %outputfolder%\boot-%file%-%id%.txt -K edgeScores --scoreexport ^
		-Z %inputfolder%\CoNet-Config.txt
	) else (
		%java18% -Xmx%memory%m -cp %path_to_conet_jar% be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser ^
		-i %inputfolder%\%inputfile% -o %outputfolder%\boot-network-%file%-%id%.gdl --inputformat biomtable -t %d_type% -f gdl ^
		--filter row_minsum/rand --filterparameter 1.0 --stand col_norm --nantreatmentparam 1 --nantreatment none ^
		--method ensemble -E %methods% --ensembleparams sim_varlogratio~lowerThreshold=%logvar%/correl_pearson~upperThreshold=%pearson%/correl_pearson~lowerThreshold=-%pearson%/correl_spearman~upperThreshold=%spearman%/correl_spearman~lowerThreshold=-%spearman%/dist_bray~lowerThreshold=%bc%/dist_kullbackleibler~lowerThreshold=%kw% ^
		--multigraph --scoremergestrategy mean --networkmergestrategy %merge_by% --minsupport %minsupport% ^
		--pvaluemerge max --edgethreshold 0.05 --iterations %iterations% --resamplemethod bootstrap -I %outputfolder%\boot-%file%-%id%.txt -K edgeScores --scoreexport ^
		-Z %inputfolder%\CoNet-Config.txt
	)
)

:: FOURTH STEP: merge final network
if %final_network% equ true (
	if %use_th_file% equ true (
		%java18% -Xmx%memory%m -cp %path_to_conet_jar% be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser ^
		-i %inputfolder%\%inputfile% -o %outputfolder%\ReBoot-network-%file%-%id%.gdl --inputformat biomtable -t %d_type% -f gdl ^
		--filter row_minsum/rand/confidence_boot --filterparameter 1.0 --stand col_norm --nantreatmentparam 1 --nantreatment none ^
		--method ensemble -E %methods% --ensembleparamfile %th_file% ^
		--multigraph --scoremergestrategy mean --networkmergestrategy %merge_by% --minsupport %minsupport% ^
		--pvaluemerge max --edgethreshold 0.05 --multicorr %multicorr% --iterations %iterations% --resamplemethod bootstrap -I %outputfolder%/boot-%file%-%id%.txt -K edgeScores ^
		--nulldistribfile %outputfolder%\permut-%file%-%id%.txt -Z %inputfolder%\CoNet-Config.txt
	) else (
		%java18% -Xmx%memory%m -cp %path_to_conet_jar% be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser ^
		-i %inputfolder%\%inputfile% -o %outputfolder%\ReBoot-network-%file%-%id%.gdl --inputformat biomtable -t %d_type% -f gdl ^
		--filter row_minsum/rand/confidence_boot --filterparameter 1.0 --stand col_norm --nantreatmentparam 1 --nantreatment none ^
		--method ensemble -E %methods% --ensembleparams sim_varlogratio~lowerThreshold=%logvar%/correl_pearson~upperThreshold=%pearson%/correl_pearson~lowerThreshold=-%pearson%/correl_spearman~upperThreshold=%spearman%/correl_spearman~lowerThreshold=-%spearman%/dist_bray~lowerThreshold=%bc%/dist_kullbackleibler~lowerThreshold=%kw% ^
		--multigraph --scoremergestrategy mean --networkmergestrategy %merge_by% --minsupport %minsupport% ^
		--pvaluemerge max --edgethreshold 0.05 --multicorr %multicorr% --iterations %iterations% --resamplemethod bootstrap -I %outputfolder%/boot-%file%-%id%.txt -K edgeScores ^
		--nulldistribfile %outputfolder%\permut-%file%-%id%.txt -Z %inputfolder%\CoNet-Config.txt
	)
)

ECHO on 
echo ******** calculation is DONE **********
EXIT
