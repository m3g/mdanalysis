#!/bin/bash
####################################################################################################

# Software name:

package=mdanalysis

# HTML file containing download list:

downloads=~/public_html/mdanalysis/download.shtml

# GIT URL:

giturl=https://github.com/leandromartinez98/mdanalysis

# Name of file containing version number

versionfile=./version.f90

####################################################################################################

year=`date +%y`
day=`date +%j`
version="${year:0:1}${year:1:1}.$day"

file="$package-$version.tar.gz" 
version_file=$version

i=2
while grep -q $file $downloads ; do
  file=$package-$version.$i.tar.gz
  version_file=$version.$i
  i=`expr $i + 1`
done
version=$version_file
file=$package-$version.tar.gz
echo "Will create file: $file"

#
# Update version number in source files
#

for file in `ls ./*`; do
  if [ $file == "./$package/src/common.f90" ] ; then
    sed -e "s/' Version 16.323
    \mv -f $file.temp $file
    sed -e "s/! Version 16.323
    \mv -f $file.temp $file
  elif [ $file == "./$package/src/algencan-pocket.f" ]; then
    sed -e "s/lm-Version 16.323
    \mv -f $file.temp $file
  else
    sed -e "s/Version 16.323
    \mv -f $file.temp $file
  fi
done
for file in `ls ../scripts/*`; do
  sed -e "s/Version 16.323
  \mv -f $file.temp $file
done
for file in `ls ../input/*`; do
  sed -e "s/Version 16.323
  \mv -f $file.temp $file
done
file=./Makefile
sed -e "s/Version 16.323
\mv -f $file.temp $file

git add -A .
git commit -m "Changed version file to $version"
git tag -a $version -m "Release $version"
git push origin master tag $version

newline="<tr><td width=190px valign=top><a href=$giturl/archive/$version.tar.gz> $file </a></td><td> Release $version </td></tr>"
htmlfile=$downloads

writeline=yes
while IFS= read -r line ; do
 
  if [ "$writeline" = "yes" ] ; then
    echo $line >> ./htmlfile_new_temp
  fi
  if [ $writeline = "no" ] ; then
    writeline="yes"
  fi

  if [[ $line == *NEW_VERSION_HERE* ]] ; then
    echo $newline >> ./htmlfile_new_temp
  fi

done < "$htmlfile"
rm $htmlfile
mv htmlfile_new_temp $htmlfile   

echo " Done. " 

