#!/bin/bash
#
# solvation.sh: Reads the input file for the solvation program using
#               VMD selection style. Than runs VMD to build a temporary
#               pdb file (as namdenergy does) and writes a temporary
#               input file for the solvation program. Finally, runs 
#               the solvation program.
#
#  Run with: ./solvation.sh solvation.inp
#
#  L. Martinez, Institut Pasteur, Apr 22, 2008.
#
#  Version 17.224
#
# IMPORTANT:
# Path for solvation program:
#
solvation=PWD/bin/solvation
#
if [ ! -e $solvation ]; then
  echo " ERROR: The solvation executable is not in the specified "
  echo "        path. Modify the path to the executable in the "
  echo "        solvation.sh script. "
  exit
fi

# Read solvation input file and extract relevant data

file=$1
while read line; do
  keyword=`echo $line | awk '{print $1}'`
  case "$keyword" in
    psf) psf=`echo $line | cut -d' ' -f 2` ;;
    dcd) dcd=`echo $line | cut -d' ' -f 2` ;;
    output) output=`echo $line | cut -d' ' -f 2` ;;
    solute) sel1=`echo $line | cut -d' ' -f 2-` ;;
    solvent) sel2=`echo $line | cut -d' ' -f 2-` ;;
  esac
done < <(cat $file) 

# Write VMD input file

vmdfile=$output.vmdtemp
groups=$output.solvationtemp
echo "
set psf $psf
set dcd $dcd
set solute \"$sel1\"
set solvent \"$sel2\"
set groupfile $groups

mol new \$psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile \$dcd type dcd first 0 last 1 step 1 filebonds 1 autobonds 1 waitfor all
 
# Clear beta
set notsel [atomselect top all frame first]
\$notsel set beta 0
\$notsel set x 0
\$notsel set y 0
\$notsel set z 0
 
# Add 1.00 to beta field of the solute
set sel1 [ atomselect top \$solute ]
\$sel1 set beta 1
 
# Add 2.00 to beta field of the solvent
set sel2 [ atomselect top \$solvent ]
\$sel2 set beta 2
 
# Write temporary pdb file with beta values
\$notsel writepdb \$groupfile
 
exit " > $vmdfile  

# Run VMD to build a temporary pdb file with group definitions

echo " ####################################################"
echo " "
echo "   Running VMD to define selections ... " 
echo " "
echo " ####################################################"
vmd -dispdev text < $vmdfile > $vmdfile.log      

# Check if there are errors in the vmd log file

grep "ERROR" $vmdfile.log

# Write the temporary solvation input file pointing to the group definitions

solvationinput=$output.solvationinp
echo "#" > $solvationinput
while read line; do
  keyword=`echo $line | awk '{print $1}'`
  case $keyword in
    solute | reference ) echo "groups $groups" >> $solvationinput ;;
    solvent | groups ) ;;
    * ) echo $line >> $solvationinput ;;
  esac  
done < <(cat $file)   

# Run solvation program

$solvation $solvationinput

# Erase temporary files

rm -f $solvationinput
rm -f $vmdfile
rm -f $vmdfile.log
rm -f $groups

