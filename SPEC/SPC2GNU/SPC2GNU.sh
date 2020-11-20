if [ $# -ne 2 ]
    then 
    echo "illegal number of parameters"
    echo 'USAGE: Script to remove header and prepare gnu file 
       CALL: SPC2GNU input.spc output.gnu'
    exit
fi

if [ -f "$1" ]; then
   echo
else 
   echo 'Input file do not exists'
   exit
fi


ne=$(awk  '/NE /{print $2}' $1)
nt=$(awk  '/NT /{print $2}' $1)
np=$(awk  '/NP /{print $2}' $1)
nohead=$(sed '1,/#######################/d' $1)

if [ $ne -eq 1 ] && [ $nt -gt 1  ] && [ $np -gt 1 ]
then
 echo 'Constant energy map recognised'
 sed '1,/#######################/d' $1  | awk -v n=$nt ' {print;} NR % n == 0 { print ""; }' > $2
elif [ $ne -gt 1 ] && [ $nt -gt 1 ] && [ $np -eq 1 ]
then
 echo 'E(K||) recognised'
 sed '1,/#######################/d' $1  | awk -v n=$nt ' {print;} NR % n == 0 { print ""; }' > $2
elif [ $ne -gt 1 ] && [ $nt -eq 1 ] && [ $np -eq 1 ]
 then
 echo 'I(E) recognised'
 sed '1,/#######################/d' $1 > $2
else
 echo 'FORMAT NOT RECOGNISED !!!!' 
 echo 'ONLY HEADER WILL BE REMOVED!!!'
 sed '1,/#Decription/d' $1 > $2
fi
 


