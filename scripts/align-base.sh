#!/bin/bash
#
# align.sh: Reads the input file for the align program using
#           VMD selection style. Than runs VMD to build a temporary
#           pdb file (as namdenergy does) and writes a temporary
#           input file for the align program. Finally, runs 
#           the align program to compute the alignment.
#
#  Run with: ./align.sh align.inp
#
#  L. Martinez, Institut Pasteur, Apr 02, 2008.
#
#  Version 19.323
#
# IMPORTANT:
# Path for align program: modify if not in the current directory
#
align=PWD/bin/align
#
if [ ! -e $align ]; then
  echo " ERROR: The align executable is not in the specified "
  echo "        path. Modify the path to the executable in the "
  echo "        align.sh script. "
  exit
fi

# Read aligninput file and extract relevant data

file=$1
while read line; do
  keyword=`echo $line | awk '{print $1}'`
  case "$keyword" in
    psf) psf=`echo $line | cut -d' ' -f 2` ;;
    dcd) dcd=`echo $line | cut -d' ' -f 2` ;;
    output) output=`echo $line | cut -d' ' -f 2` ;;
    reference) sel1=`echo $line | cut -d' ' -f 2-` ;;
    compute_rmsd_of) sel2=`echo $line | cut -d' ' -f 2-` ;;
  esac
done < <(cat $file) 

# Write VMD input file

vmdfile=$output.vmdtemp
groups=$output.groupstemp
echo "
set psf $psf
set dcd $dcd
set reference \"$sel1\"
set compute_reference_of \"$sel2\"
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
 
# Add 1.00 to beta field of first selection
set sel1 [ atomselect top \$reference ]
\$sel1 set occupancy 1
 
# Add 2.00 to beta field of second selection
set sel2 [ atomselect top \$compute_reference_of ]
\$sel2 set beta 2
 
# Write fake pdb file with beta values
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

# Write the temporary align input file pointing to the group definitions

aligninput=$output.aligninp
echo "#" > $aligninput
while read line; do
  keyword=`echo $line | awk '{print $1}'`
  case $keyword in
    reference ) echo "groups $groups" >> $aligninput ;;
    groups | compute_rmsd_of ) ;;
    * ) echo $line >> $aligninput ;;
  esac  
done < <(cat $file)   

# Run align program

$align $aligninput

# Erase temporary files

rm -f $aligninput
rm -f $vmdfile
rm -f $vmdfile.log
rm -f $groups

