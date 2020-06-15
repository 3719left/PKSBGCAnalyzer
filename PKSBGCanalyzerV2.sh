#!/bin/bash
#Usage ./PKSBGCanalyzerV1.sh input_protein.faa(1) Clusterfinderresult(2) proteincompositionfile(3) outputfoldername(4) PfamA.hmm(5) mode(6)
#MODE: S strict =; R relaxed >=


##0 early stage checks
# check if required files exist
if [ ! -e $1 ]
then
    echo "Protein file is missing."
    exit 1
fi

if [ ! -e $2 ]
then
    echo "Clustrfinderesult is missing."
    exit 1
fi

if [ ! -e $3 ]
then
    echo "Proteincomposition file is missing."
    exit 1
fi

if [ ! -e $5 ]
then
    echo "PfamA file is missing."
    exit 1
fi

if [ -e $4 ]
then
    echo "Warning output file exists!."
    echo 'Do you want to overwrite it? (Y/N):'
    read choice
    if [ "$choice" = 'Y' -o "$choice" = 'y' ]
    then
      echo "Overwrite and Continue!"
    else
      exit 1
    fi
fi

#check if prameter number is right and set running mode to relaxed(R) if not offered.
if [ "$#" -eq 5 ]
then
    set $1 $2 $3 $4 $5 R
elif [ "$#" -eq 6 ]
then
    if [ $6 = S -o $6 = R ]
    then
        :
    else
        echo "Mode value should either be S or R!"
        exit 1
    fi
elif [ "$#" -lt 5 -o "$#" -gt 6 ]
then
    echo "There should be 6 parameters! and the last one is running mode (default:R--relaxed)"
    echo "You just entered $#"
    exit 1
fi

if [ ! -e $1.ssi ]
then
    esl-sfetch --index $1
fi

echo "Pass command check with $6 mode..."
#construct folder structure.
mkdir $4
mkdir $4/BGC $4/PKScore $4/Domains $4/Domains/KS
mkdir -p $4/intermediate/scripts
mkdir $4/intermediate/files

cat $3 | while read a
do
    if [ `echo $a | awk '{print $3}'` ]  && [ `echo $a | awk '{print $3}'` = E ]
    then
        d=`echo $a | awk '{print $1}'`
        mkdir $4/Domains/$d
    fi
done

#Selection between two mode
if [ $6 = S ] 
then
    SIG=-eq #strict mode have to be "="
else
    SIG=-ge #relaxed mode can be ">="
fi

## get BGC information
# get serial number of target BGC
echo "Finding target BGCs and Extracting BGC information..."
grep "\S" $2 > $4/intermediate/files/NA0
awk '{print $1}' $4/intermediate/files/NA0 | sort -u > $4/intermediate/files/TargetBGCnum0

#Reorganize Clusterfinderresult
echo "
import sys
inputfile1 = open(sys.argv[1])
inputfile2 = open(sys.argv[2])
outputfile = open(sys.argv[3], 'w')

lines1 = inputfile1.readlines()
lines2 = inputfile2.readlines()

inputfile1.close()
inputfile2.close()

#BGC information list10
list10=['']
Numlist=[]
for a in lines1:
    Numlist.append(int(a.strip().split()[0]))
    
Numlist2=Numlist[:]
Numlist2.reverse()
LEN=len(Numlist)
for x in lines2:
    Box=lines1[Numlist.index(int(x.strip())):(LEN-Numlist2.index(int(x.strip())))]
    list10.append(''.join(Box))
#edit for loop

for y in list10:
    outputfile.write(y)

outputfile.close()

" > $4/intermediate/scripts/targetBGCextractor


num=1
x=`wc -l $3 | awk '{print $1}'`
while read a
do
    d=`echo $a | awk '{print $1}'`
    n=`echo $a | awk '{print $2}'`
    grep ".*[_[:blank:]]$d[_[:blank:]].*" $4/intermediate/files/NA$(expr $num - 1) | awk '{print $1}' | sort -n | uniq -c > $4/intermediate/files/TargetBGCcountnum$num
    if [ ! -s $4/intermediate/files/TargetBGCcountnum$num ]
    then
        echo $d was not found!
        exit 1
    fi
    cat $4/intermediate/files/TargetBGCcountnum$num | while read y
    do
        if [ `echo $y | awk '{print $1}'` $SIG $n ]
        then
            echo $y | awk '{print $2}' >> $4/intermediate/files/TargetBGCnum$num
        fi
    done
    if [ ! -s $4/intermediate/files/TargetBGCnum$num ]
    then
        echo $d number does not full fill request!
        exit 1
    fi
    sort $4/intermediate/files/TargetBGCnum$(expr $num - 1) > $4/intermediate/files/L1; cat $4/intermediate/files/L1 > $4/intermediate/files/TargetBGCnum$(expr $num - 1)
    sort $4/intermediate/files/TargetBGCnum$num > $4/intermediate/files/L2; cat $4/intermediate/files/L2 > $4/intermediate/files/TargetBGCnum$num
    comm -12 $4/intermediate/files/TargetBGCnum$(expr $num - 1) $4/intermediate/files/TargetBGCnum$num > $4/intermediate/files/TargetBGCnumTemp
    if [ ! -s $4/intermediate/files/TargetBGCnumTemp ]
    then
        echo No BGC have all required structure number!
        exit 1
    fi
    sort -n $4/intermediate/files/TargetBGCnumTemp > $4/intermediate/files/TargetBGCnum$num
    python3 $4/intermediate/scripts/targetBGCextractor $4/intermediate/files/NA$(expr $num - 1) $4/intermediate/files/TargetBGCnum$num $4/intermediate/files/NA$num
    num=$(expr $num + 1)
done < $3

sum=`wc -l $4/intermediate/files/TargetBGCnum$x | awk '{print $1}'`
echo "$sum BGCs fullfill required standards!"
cp $4/intermediate/files/NA$x $4/BGC/BGC.txt


awk '{print $9}' $4/BGC/BGC.txt | sort -u | esl-sfetch -f $1 - > $4/BGC/BGC.fa

esl-sfetch --index $4/BGC/BGC.fa >/dev/null

#2 Extract Core PKS proteins(KS-containing protein)
echo "Extracting Core PKS proteins(KS-containing protein)..."

grep '\[.*PKS.*\]' $4/BGC/BGC.txt > $4/PKScore/PKScore.txt
awk '{print $9}' $4/PKScore/PKScore.txt | sort -u | esl-sfetch -f $4/BGC/BGC.fa - > $4/PKScore/PKScore.fa

#3 Extract KS domains and requested domains

cp $4/PKScore/PKScore.txt $4/Domains/KS/KSdomain.txt
echo "Extracting KS domains..."

#Extract KS domains.
hmmfetch $5 ketoacyl-synt > $4/intermediate/files/KSN.hmm
hmmfetch $5 Ketoacyl-synt_C > $4/intermediate/files/KSC.hmm

echo "
inputfile = open('$4/intermediate/files/Testfinal.txt')
outputfile = open('$4/intermediate/files/newdata.txt', 'w')
lines = inputfile.readlines()
inputfile.close()
LENGTH=len(lines)
Result=[]
SPACE='\t'
Num = 0
while Num < (LENGTH-1):
  ID, DOMAIN, START, STOP = lines[Num].strip().split()
  Num += 1
  ID1, DOMAIN1, START1, STOP1 = lines[Num].strip().split()
  if DOMAIN == 'N' and DOMAIN1 == 'C' and ID == ID1 and ( int(START1) - int(STOP) >= 0 ) and ( int(STOP1) - int(START) <= 600): #KSsize 600
    ADD = ID1+SPACE+START+SPACE+STOP1+'\n'
    Result.append(ADD)
  else:
    pass
for STR in Result:
  outputfile.write(STR)
outputfile.close()
" > $4/intermediate/scripts/KS.py

esl-sfetch --index $4/PKScore/PKScore.fa >/dev/null
hmmsearch --domE 1e-5 --domtblout $4/intermediate/files/KSN.dtbl $4/intermediate/files/KSN.hmm $4/PKScore/PKScore.fa >/dev/null
hmmsearch --domE 1e-5 --domtblout $4/intermediate/files/KSC.dtbl $4/intermediate/files/KSC.hmm $4/PKScore/PKScore.fa >/dev/null 
grep -v "^#" $4/intermediate/files/KSN.dtbl | awk '{print $1,"N",$20,$21}' > $4/intermediate/files/Testdata.txt
grep -v "^#" $4/intermediate/files/KSC.dtbl | awk '{print $1,"C",$20,$21}' >> $4/intermediate/files/Testdata.txt
sort -k 1,1 -k 3,3n -k 2,2 $4/intermediate/files/Testdata.txt > $4/intermediate/files/Testfinal.txt
python3 $4/intermediate/scripts/KS.py
awk '{print $1"/"$2"-"$3, $2, $3, $1}' $4/intermediate/files/newdata.txt | esl-sfetch -Cf $4/PKScore/PKScore.fa - > $4/Domains/KS/KSdomain.fa
sed -i $4/Domains/KS/KSdomain.fa -e 's/>/>KS_/g'

#Extract Labled domains.

cat $3 | while read a
do
    if [ `echo $a | awk '{print $3}'` ] && [ `echo $a | awk '{print $3}'` = E ]
    then
        b=`echo $a | awk '{print $1}'`
        echo "Extracting $b domains..."
        hmmfetch $5 $b > $4/intermediate/files/${b}.hmm
        hmmsearch --domE 1e-5 --domtblout $4/intermediate/files/${b}.dtbl $4/intermediate/files/${b}.hmm $4/BGC/BGC.fa >/dev/null
        grep -v "^#" $4/intermediate/files/${b}.dtbl | awk '{print $1"/"$20"-"$21, $20, $21, $1}' | esl-sfetch -Cf $4/BGC/BGC.fa - > $4/Domains/$b/${b}domain.fa
        sed -i $4/Domains/$b/${b}domain.fa -e "s/>/>'$b'_/g"
        grep ".*[_[:blank:]]$b[_[:blank:]].*" $4/BGC/BGC.txt > $4/Domains/$b/${b}domain.txt
    fi
done

#Clear unnecessary files
rm $4/PKScore/PKScore.fa.ssi $4/BGC/BGC.fa.ssi $4/intermediate/files/NA0
echo "All done!"