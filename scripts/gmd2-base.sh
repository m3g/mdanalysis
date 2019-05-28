#!/bin/bash
#
# gmd2.sh: Reads the input file for the gmd2 program using
#         VMD selection style. Than runs VMD to build a temporary
#         pdb file (as namdenergy does) and writes a temporary
#         input file for the gmd program. Finally, runs 
#         the gmd program.
#
#  Run with: ./gmd2.sh gmd2.inp
#
#  L. Martinez, Institut Pasteur, Apr 22, 2008.
#
#  Version 19.148
#
# IMPORTANT:
# Path for gmd2 program:
#
gmd2=PWD/bin/gmd2
#
if [ ! -e $gmd2 ]; then
  echo " ERROR: The gmd2 executable is not in the specified "
  echo "        path. Modify the path to the executable in the "
  echo "        gmd2.sh script. "
  exit
fi

# Read gmd2 input file and extract relevant data

file=$1
while read line; do
  keyword=`echo $line | awk '{print $1}'`
  case "$keyword" in
    psf) psf=`echo $line | cut -d' ' -f 2` ;;
    dcd) dcd=`echo $line | cut -d' ' -f 2` ;;
    output) output=`echo $line | cut -d' ' -f 2` ;;
    solute) sel1=`echo $line | cut -d' ' -f 2-` ;;
    solvent) sel2=`echo $line | cut -d' ' -f 2-` ;;
    exclude) sel3=`echo $line | cut -d' ' -f 2-` ;;
  esac
done < <(cat $file) 

# If the "exclude" selection was not set, set it equal to solute

if [ -z "$sel3" ]; then 
  sel3=$sel1
fi

# Write VMD input file

vmdfile=$output.vmdtemp
groups=$output.gmd2temp
echo "
set psf $psf
set dcd $dcd
set solute \"$sel1\"
set solvent \"$sel2\"
set exclude \"$sel3\"
set groupfile $groups

mol new \$psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile \$dcd type dcd first 0 last 1 step 1 filebonds 1 autobonds 1 waitfor all
 
# Clear beta
set notsel [atomselect top all frame first]
\$notsel set beta 0
\$notsel set occupancy 0
\$notsel set x 0
\$notsel set y 0
\$notsel set z 0
 
# Add 1.00 to beta field of the solute
set sel1 [ atomselect top \$solute ]
\$sel1 set beta 1
 
# Add 2.00 to beta field of the solvent
set sel2 [ atomselect top \$solvent ]
\$sel2 set beta 2

# Add 1.00 to to occupancy field of the excluded zone
set sel3 [ atomselect top \$exclude ]
\$sel3 set occupancy 1

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

# Write the temporary gmd2 input file pointing to the group definitions

gmd2input=$output.gmd2inp
echo "#" > $gmd2input
while read line; do
  keyword=`echo $line | awk '{print $1}'`
  case $keyword in
    solute | reference ) echo "groups $groups" >> $gmd2input ;;
    solvent | groups ) ;;
    exclude | groups ) ;;
    * ) echo $line >> $gmd2input ;;
  esac  
done < <(cat $file)   

# Run gmd2 program

$gmd2 $gmd2input

# Erase temporary files

rm -f $gmd2input
rm -f $vmdfile
rm -f $vmdfile.log
rm -f $groups

