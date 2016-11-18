#!/bin/bash
#
# Distance.sh: Reads the input file for the distance program using
#              VMD selection style. Than runs VMD to build a temporary
#              pdb file (as namdenergy does) and writes a temporary
#              input file for the distance program. Finally, runs 
#              the distance program to compute the distances.
#
#  Run with: ./distance.sh distance.inp
#
#  L. Martinez, Institut Pasteur, Apr 02, 2008.
#
#  Version 16.146
#
# IMPORTANT:
# Path for distance program: modify if not in the current directory
#
distance=PWD/bin/distance
#
if [ ! -e $distance ]; then
  echo " ERROR: The distance executable is not in the specified "
  echo "        path. Modify the path to the executable in the "
  echo "        distance.sh script. "
  exit
fi

# Read distance input file and extract relevant data

file=$1
while read line; do
  keyword=`echo $line | awk '{print $1}'`
  if [ "$keyword" = "psf" ]; then
    psf=`echo $line | cut -d' ' -f 2`
  fi
  if [ "$keyword" = "dcd" ]; then
    dcd=`echo $line | cut -d' ' -f 2`
  fi
  if [ "$keyword" = "output" ]; then
    output=`echo $line | cut -d' ' -f 2`
  fi
  if [ "$keyword" = "group1" ]; then
    sel1=`echo $line | cut -d' ' -f 2-`
  fi
  if [ "$keyword" = "group2" ]; then
    sel2=`echo $line | cut -d' ' -f 2-`
  fi
done < <(cat $file) 

# Write VMD input file

vmdfile=$output.vmdtemp
groups=$output.groupstemp
echo "
set psf $psf
set dcd $dcd
set group1 \"$sel1\"
set group2 \"$sel2\"
set groupfile $groups

mol new \$psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile \$dcd type dcd first 0 last 1 step 1 filebonds 1 autobonds 1 waitfor all
 
# Clear beta
set notsel [atomselect top all frame first]
\$notsel set beta 0
\$notsel set x 0
\$notsel set y 0
\$notsel set z 0
 
# Add 1.00 to beta field of first selection
set sel1 [ atomselect top \$group1 ]
\$sel1 set beta 1
 
# Add 2.00 to beta field of second selection
set sel2 [ atomselect top \$group2 ]
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

# Write the temporary distance input file pointing to the group definitions

distanceinput=$output.distanceinp
echo "#" > $distanceinput
while read line; do
  keyword=`echo $line | awk '{print $1}'`
  case $keyword in
    group1 ) echo "groups $groups" >> $distanceinput ;;
    group2 | groups ) ;;
    * ) echo $line >> $distanceinput ;;
  esac  
done < <(cat $file)   

# Run distance program

$distance $distanceinput

# Erase temporary files

rm -f $distanceinput
rm -f $vmdfile
rm -f $vmdfile.log
rm -f $groups

