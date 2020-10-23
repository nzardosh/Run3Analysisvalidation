#!/bin/bash

bash clean.sh

CASE=5

DOCONVERT=1   # Convert AliESDs.root to AO2D.root.
DOQA=0        # Run the QA task with O2.
DORUN1=1      # Run the heavy-flavour tasks with AliPhysics.
DORUN3=1      # Run the heavy-flavour tasks with O2.
DOCOMPARE=1   # Compare AliPhysics and O2 output.

RUN5=0        # Use Run 5 input.
CONVSEP=1     # Convert ESD files separately.
PARALLELISE=0 # Parallelise O2 tasks.
DEBUG=0       # Print out more information.

# Default settings
JSON="$PWD/dpl-config_std.json"
JSONRUN5="$PWD/dpl-config_run5.json"
ISMC=0
TRIGGERSTRINGRUN2=""
TRIGGERBITRUN3=-1
NMAX=-1

if [ $CASE -eq 0 ]; then
  INPUTDIR="../twikiinput"
  STRING="AliESDs_ppK0starToyMC.root"
fi

if [ $CASE -eq 1 ]; then
  INPUTDIR="/mnt/temp/Run3data/data/LHC15o_246751/pass1"
  STRING="15000246751019.110/AliESDs.root"
  TRIGGERSTRINGRUN2="CV0L7-B-NOPF-CENT"
  TRIGGERBITRUN3=5 #FIXME
  NMAX=5
fi

if [ $CASE -eq 2 ]; then
  INPUTDIR="/data/Run3data/alice_sim_2015_LHC15k1a3_246391/246391"
  ISMC=1
  STRING="00*/AliESDs.root"
  NMAX=1
fi

if [ $CASE -eq 3 ]; then
  INPUTDIR="/data/Run3data/output"
  STRING="00*/AliESDs.root"
fi

if [ $CASE -eq 4 ]; then
  INPUTDIR="/data/Run3data/alice_sim_2018_LHC18a4a2_cent/282099"
  STRING="0*/AliESDs.root"
fi

if [ $CASE -eq 5 ]; then
  INPUTDIR="/mnt/temp/Run3data_Vit/LHC18a4a2_cent/282341"
  STRING="00*/AliESDs.root"
fi

#INPUTDIR="/data/Run3data/output" #K0* MC injected
#INPUTDIR="/data/Run3data/alice_sim_2018_LHC18a4a2_cent/282099" #D2H MC sample
#INPUTDIR="/data/Run3data/alice_sim_2015_LHC15k1a3_246391/246391" #HIJING MC PbPb

# Lists of input files
LISTFILESALI="list_ali.txt"
ls $INPUTDIR/$STRING > $LISTFILESALI
LISTFILESO2="listrun3.txt"
LISTFILESO2RUN5="listrun5.txt"

# Output files names
FILEOUTALI="JetRun1.root"
FILEOUTO2="AnalysisResults.root"
FILEOUTQA="AnalysisResultsQA.root"

# Steering commands
ENVALI="alienv setenv AliPhysics/latest -c"
ENVO2="alienv setenv O2/latest -c"
ENVALIO2="alienv setenv AliPhysics/latest,O2/latest -c"
CMDROOT="root -b -q -l"

# Adjust settings for Run5.
O2INPUT=$LISTFILESO2
if [ $RUN5 -eq 1 ]; then
  echo -e "\nUsing Run 5 settings and O2 input"
  O2INPUT=$LISTFILESO2RUN5
  JSON="$JSONRUN5"
fi


# Convert AliESDs.root to AO2D.root.
if [ $DOCONVERT -eq 1 ]; then
  [ -f "$LISTFILESALI" ] || { echo "Converting: Error: File $LISTFILESALI does not exist."; exit 1; }
  echo -e "\nConverting... ($(cat $LISTFILESALI | wc -l) files)"
  if [ $DOQA -eq 1 ]; then
    echo "Setting MC mode ON."
    ISMC=1
  fi
  if [ $CONVSEP -eq 1 ]; then
    echo "Converting files separately"
    $ENVALI bash convert_batch.sh $LISTFILESALI $LISTFILESO2 $ISMC $DEBUG || exit 1 # Run the batch script in the ALI environment.
  else
    LOGFILE="log_convert.log"
    rm -f $LOGFILE
    echo "logfile: $LOGFILE"
    $ENVALI $CMDROOT "convertAO2D.C(\"$LISTFILESALI\", $ISMC, $NMAX)" > $LOGFILE 2>&1 || { echo "Error"; exit 1; }
    echo "$PWD/AO2D.root" > $LISTFILESO2
    rm -f $FILEOUTO2
  fi
fi

# Perform simple QA studies with O2.
if [ $DOQA -eq 1 ]; then
  #LOGFILE="log_o2_qa.log"
  [ -f "$O2INPUT" ] || { echo "QA task: Error: File $O2INPUT does not exist."; exit 1; }
  echo -e "\nRunning the QA task with O2... ($(cat $O2INPUT | wc -l) files)"
  rm -f $FILEOUTO2 $FILEOUTQA
  O2ARGS="--shm-segment-size 16000000000 --configuration json://$JSON"
  if [ $PARALLELISE -eq 1 ]; then
    NPROC=3
    echo "Parallelisation ON ($NPROC)"
    O2ARGS="$O2ARGS --pipeline qa-tracking-kine:$NPROC,qa-tracking-resolution:$NPROC"
  fi
  O2EXEC="o2-analysis-qatask $O2ARGS -b"
  O2SCRIPT="script_o2_qa.sh"
  cat << EOF > $O2SCRIPT # Create a temporary script with the full O2 commands.
#!/bin/bash
$O2EXEC
EOF
  #$ENVO2 bash $O2SCRIPT > $LOGFILE 2>&1 || { echo "Error"; exit 1; } # Run the script in the O2 environment.
  #grep WARN $LOGFILE | sort -u
  $ENVO2 bash o2_batch.sh $O2INPUT $JSON $O2SCRIPT $DEBUG || exit 1 # Run the batch script in the O2 environment.
  rm -f $O2SCRIPT
  mv $FILEOUTO2 $FILEOUTQA
  rm -rf output_o2_qa
  mv output_o2 output_o2_qa
  mv log_o2.log log_o2_qa.log
fi

# Run the jet tasks with AliPhysics.
if [ $DORUN1 -eq 1 ]; then
  [ -f "$LISTFILESALI" ] || { echo "Jet tasks ALI: Error: File $LISTFILESALI does not exist."; exit 1; }
  echo -e "\nRunning the Jet tasks with AliPhysics... ($(cat $LISTFILESALI | wc -l) files)"
  #$ENVALI bash ali_batch.sh $LISTFILESALI $JSON $FILEOUTALI # Run the batch script in the ALI environment.
  $ENVALIO2 bash ali_batch.sh $LISTFILESALI $JSON $FILEOUTALI $DEBUG || exit 1 # Run the batch script in the ALI+O2 environment.
fi

# Run the heavy-flavour tasks with O2.
if [ $DORUN3 -eq 1 ]; then
  #LOGFILE="log_o2_jet.log"
  [ -f "$O2INPUT" ] || { echo "Jet tasks O2: Error: File $O2INPUT does not exist."; exit 1; }
  echo -e "\nRunning the Jet tasks with O2... ($(cat $O2INPUT | wc -l) files)"
  rm -f $FILEOUTO2
  # Option --configuration has priority over --aod-file.
#  O2ARGS="--shm-segment-size 16000000000 --configuration json://$PWD/dpl-config_std.json --aod-file $AOD3NAME"
  O2ARGS="--shm-segment-size 16000000000 --configuration json://$JSON"
  O2ARGS_JETFINDER="$O2ARGS"
  O2ARGS_JETSUBSTRUCTURE="$O2ARGS"
  O2ARGS_JETFINDERHADRONRECOIL="$O2ARGS"
  if [ $PARALLELISE -eq 1 ]; then
    NPROC=3
    echo "Parallelisation ON ($NPROC)"
    O2ARGS_JETFINDER="$O2ARGS_JETFINDER --pipeline jet-finder:$NPROC"
    O2ARGS_SUBSTRUCTURE="$O2ARGS_SUBSTRUCTURE --pipeline jet-substructure:$NPROC"
    O2ARGS_JETFINDERHADRONRECOIL="$O2ARGS_JETFINDERHADRONRECOIL --pipeline jet-finder-hadron-recoil:$NPROC"
  fi
  O2EXEC_JETFINDER="o2-analysis-jet-finder $O2ARGS_JETFINDER"
  O2EXEC_JETSUBSTRUCTURE="o2-analysis-jet-substructure $O2ARGS_JETSUBSTRUCTURE"
  O2EXEC_JETFINDERHADRONRECOIL="o2-analysis-jet-finder-hadron-recoil $O2ARGS_JETFINDERHADRONRECOIL"
  #O2EXEC="$O2EXEC_SKIM | $O2EXEC_PIDTPC | $O2EXEC_PIDTOF | $O2EXEC_CAND | $O2EXEC_SEL | $O2EXEC_TASK -b"
  O2EXEC="$O2EXEC_JETFINDER | $O2EXEC_JETSUBSTRUCTURE | $O2EXEC_JETFINDERHADRONRECOIL -b"
  #O2EXEC="$O2EXEC_JETFINDER | $O2EXEC_JETSUBSTRUCTURE  -b" 
  O2SCRIPT="script_o2_jet.sh"
  cat << EOF > $O2SCRIPT # Create a temporary script with the full O2 commands.
#!/bin/bash
$O2EXEC
EOF
  #$ENVO2 bash $O2SCRIPT > $LOGFILE 2>&1 || { echo "Error"; exit 1; } # Run the script in the O2 environment.
  #grep WARN $LOGFILE | sort -u
  $ENVO2 bash o2_batch.sh $O2INPUT $JSON $O2SCRIPT $DEBUG || exit 1 # Run the batch script in the O2 environment.
  rm -f $O2SCRIPT
  rm -rf output_o2_jet
  mv output_o2 output_o2_jet
  mv log_o2.log log_o2_jet.log
fi


# Compare AliPhysics and O2 output.
if [ $DOCOMPARE -eq 1 ]; then
  LOGFILE="log_compare.log"
  echo -e "\nComparing... (logfile: $LOGFILE)"
  ok=1
  for file in "$FILEOUTALI" "$FILEOUTO2"; do
    [ -f "$file" ] || { echo "Error: File $file does not exist."; ok=0; }
  done
  [ $ok -ne 1 ] && exit 1
  $ENVALI $CMDROOT "Compare.C(\"$FILEOUTO2\",\"$FILEOUTALI\")" > $LOGFILE 2>&1 || { echo "Error"; exit 1; }
fi

echo -e "\nDone"

rm $LISTFILESALI

exit 0
