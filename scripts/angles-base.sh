#!/bin/bash
#
# angles.sh: Reads the input file for the angles program using
#            VMD selection style. Than runs VMD to build a temporary
#            pdb file (as namdenergy does) and writes a temporary
#            input file for the angles program. Finally, runs 
#            the angles program.
#
#  Run with: ./angles.sh angles.inp
#
#  L. Martinez, Institut Pasteur, Apr 22, 2008.
#
#  Version 19.136
#
# IMPORTANT:
# Path for angles program:
#
angles=PWD/bin/angles
#
if [ ! -e $angles ]; then
  echo " ERROR: The angles executable is not in the specified "
  echo "        path. Modify the path to the executable in the "
  echo "        angles.sh script. "
  exit
fi

# Read angles input file and extract relevant data

file=$1
while read line; do
  keyword=`echo $line | awk '{print $1}'`
  case "$keyword" in
    psf) psf=`echo $line | cut -d' ' -f 2` ;;
    dcd) dcd=`echo $line | cut -d' ' -f 2` ;;
    output) output=`echo $line | cut -d' ' -f 2` ;;
    atom1) sel1=`echo $line | cut -d' ' -f 2-` ;;
    atom2) sel2=`echo $line | cut -d' ' -f 2-` ;;
    atom3) sel3=`echo $line | cut -d' ' -f 2-` ;;
    atom4) sel4=`echo $line | cut -d' ' -f 2-` ;;
  esac
done < <(cat $file) 

# Write VMD input file

vmdfile=$output.vmdtemp
groups=$output.anglestemp
echo "
set psf $psf
set dcd $dcd
set atom1 \"$sel1\"
set atom2 \"$sel2\"
set atom3 \"$sel3\"
set atom4 \"$sel4\"
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
 
# Add 1.00 to beta field of the first atom
set sel1 [ atomselect top \$atom1 ]
\$sel1 set occupancy 1
 
# Add 2.00 to beta field of the second atom
set sel2 [ atomselect top \$atom2 ]
\$sel2 set occupancy 2
 
# Add 3.00 to beta field of the third atom
set sel3 [ atomselect top \$atom3 ]
\$sel3 set beta 3
 
# Add 4.00 to beta field of the fourth atom
set sel4 [ atomselect top \$atom4 ]
\$sel4 set beta 4
 
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

# Write the temporary angles input file pointing to the group definitions

anglesinput=$output.anglesinp
echo "#" > $anglesinput
while read line; do
  keyword=`echo $line | awk '{print $1}'`
  case $keyword in
    atom1 ) echo "groups $groups" >> $anglesinput ;;
    atom2 ) ;;
    atom3 ) ;;
    atom4 ) ;;
    * ) echo $line >> $anglesinput ;;
  esac  
done < <(cat $file)   

# Run angles program

$angles $anglesinput

# Erase temporary files

rm -f $anglesinput
rm -f $vmdfile
rm -f $vmdfile.log
rm -f $groups

