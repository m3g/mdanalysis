#!/bin/bash
####################################################################################################

# Software name:

package=mdanalysis

# HTML file containing download list:

downloads=~/public_html/mdanalysis/download.shtml

# GIT URL:

giturl=https://github.com/leandromartinez98/mdanalysis

# Name of file containing version number

versionfile=./common.f90

####################################################################################################

#git log --pretty=oneline 16.323...16.330 | awk '{$1=""; print "-"$0}'

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

src=$versionfile
sed -e "s/' Version.*/' Version $version '/" $src > $src.temp
\mv -f $src.temp $src
for src in `ls ../scripts/*`; do
  sed -e "s/Version.*/Version $version/" $src > $src.temp
  \mv -f $src.temp $src
done
for src in `ls ../input/*`; do
  sed -e "s/Version.*/Version $version/" $src > $src.temp
  \mv -f $src.temp $src
done
src=./Makefile
sed -e "s/Version.*/Version $version/" $src > $src.temp
\mv -f $src.temp $src

git add -A
git commit -m "Changed version file to $version"
git tag -a $version -m "Release $version"
git push origin master tag $version

today=`date +"%b %d, %Y"`
changelog="https://github.com/leandromartinez98/$package/releases/tag/$version"
newline="<tr><td width=190px valign=top><a href=$giturl/archive/$version.tar.gz> $file </a></td><td> Released on $today - <a target=newpage href=$changelog> [chan
ge log at github] </a></td></tr>"
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

echo "----------------------"
echo "CHANGE LOG:"
echo "----------------------"
range=`git tag | tail -n 2 | xargs | sed 's! !...!'`
git log --pretty=oneline $range | awk '{$1=""; print "-"$0}'
echo "----------------------"

echo " Done. " 

