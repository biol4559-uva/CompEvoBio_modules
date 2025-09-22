
ls -d *\$* | awk '{print $0}'


for i in $( ls -d *\$* );
do
    rname=$( sed 's/([A-Za-z0-9_]*)/$1/g' | sed 's/\$//g' )
    echo $i"\t'${rname}
done


cd /standard/BerglandTeach/mapping_output

ls -d * | while read line; do
  rname=$( echo ${line} | sed 's/([A-Za-z0-9_]*)/$1/g' | sed 's/\$//g' )
  nFiles=$( ls -lh $line | wc -l )
  op=$( echo $line $rname $nFiles)
  if [ "$nFiles" -lt 14 ]; then
    echo $line $nFiles
    mv $line ../badOutput_2
  fi
done


ls -ld ExpEvo* | while read line; do
  # line="drwxrws---+ 3 asn8nx  biol4020-aob2x 14 Sep 16 23:24 ExpEvo_PRJNA657615_L_pop2_gen59_1999-10-01"
  year=$( echo $line | cut -f8 -d' ' )

  if [[ "$year" -eq "2024" ]]; then
    echo $line $nFiles
    mv $line ../2024Output
  fi
done
