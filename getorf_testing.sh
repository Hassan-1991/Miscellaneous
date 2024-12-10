#This code extracts all ORFs from a sequence
#Then, starting from those stop codons, attempts to find the furthest start codon until another upstream stop is reached
#Helps evaluate the mismatch between actual longest ORF, vs. what getorf considers the longest ORF.

#Extract all ORFs from the sequence REL606.faa
#Use alternative initiation codons: TTG, GTG, CTG (-table 1)
#Minimum ORF length: 30bp (-minsize 30)
#Return nucleotide sequences (find 3)
getorf -sequence REL606.faa -outseq REL606.getorf -table 1 -minsize 30 -find 3
#Since the CTG start codon is barely used in bacteria (https://academic.oup.com/nar/article/45/7/3615/2990259), we remove all ORFs starting with a CTG
seqkit fx2tab REL606.getorf | grep -v -P "\tCTG" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" > REL606.getorf.noCTG

#Convert the ORF start, stop and strand information returned by getorf into a gtf file:
grep "^>" REL606.getorf.noCTG | grep -v "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$2,"\t",$4+3,"\t.\t+\t0\tgene_id \"",$1,"\";transcript_id \"",$1,"\";"}' | sed "s/_/\t/2" | cut -f1,3- > REL606.getorf.noCTG.gtf
grep "^>" REL606.getorf.noCTG | grep "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$4-3,"\t",$2,"\t.\t-\t0\tgene_id \"",$1,"\";transcript_id \"",$1,"\";"}' | sed "s/_/\t/2" | cut -f1,3- >> REL606.getorf.noCTG.gtf

#Scan the upstream 500bp of the start codon to identify a further start before a stop codon is reached
#Return the distance between the furthest start codon and the ORF stop codon
awk -F '\t' '($7=="+")' REL606.getorf.noCTG.gtf | #Starting with ORFs on the plus strand
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-500,$5,$6,$7,$8,$9}' | #Increase start coordinate by 500
awk -F '\t' '{OFS=FS}{if ($4 < 1) $4 = 1; print}' | #In case the start coordinate is less than 1, replace it by 1
gtf2bed | bedtools getfasta -s -name -fi REL606.faa -bed - | #get corresponding sequence
sed '/^>/!{ s/.*/echo "&" | rev/e }' | #Reverse all sequences
awk '{if ($0 ~ /^>/) print; else {gsub(/.{3}/, "& "); print}}' | #cut each sequence into 3-bp codons, delimited by space
cut -f2- -d " " | #Ignore the current stop codon for now
sed -E 's/AGT|AAT|GAT/%/g' | #replace other stop codons with a uniform identifier
seqkit fx2tab | cut -f1 -d "%" | #Extract sequence until an upstream stop codon is reached, returning the entire sequence delimited by two stops
sed -E 's/GTA|GTG|GTT/@/g' | #replace start codons with a uniform identifier
rev | cut -f2- -d "@" | rev | #reverse the sequence again, retain entire sequence to the 3' of the first instance of a start codon
sed "s/@/NNN/g" | sed "s/ //g" | awk -F '\t' '{print $1,length($2)+6}' | sed "s/ /\t/g" > REL606.getorf.noCTG.longestORFlengths.tsv #massage and return length of the longest ORF

awk -F '\t' '($7=="-")' REL606.getorf.noCTG.gtf | #Similarly for ORFs on minus strand
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5+500,$6,$7,$8,$9}' | #Here, the start codon is in column #5
awk -F '\t' '{OFS=FS}{if ($5 > 4629812) $5 = 4629812; print}' | #If start coordinate goes beyond total genome length, replace it with the genome length
gtf2bed | bedtools getfasta -s -name -fi REL606.faa -bed - |
sed '/^>/!{ s/.*/echo "&" | rev/e }' |
awk '{if ($0 ~ /^>/) print; else {gsub(/.{3}/, "& "); print}}' |
cut -f2- -d " " |
sed -E 's/AGT|AAT|GAT/%/g' |
seqkit fx2tab | cut -f1 -d "%" |
sed -E 's/GTA|GTG|GTT/@/g' |
rev | cut -f2- -d "@" | rev |
sed "s/@/NNN/g" | sed "s/ //g" | awk -F '\t' '{print $1,length($2)+6}' | sed "s/ /\t/g" >> REL606.getorf.noCTG.longestORFlengths.tsv

sed "s/(+)//g" REL606.getorf.noCTG.longestORFlengths.tsv | sed "s/(-)//g" | sort -k1 > temp #Store the sorted and massaged version of this file for subsequent steps

#Join the tsv with longest ORF lengths with current gtf
#Subtract length of longest ORF to stop codon to return longest ORF coordinates
#Return output to a new gtf
cut -f -2 -d "\"" REL606.getorf.noCTG.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="+")' | awk '{print $2"\t.\tCDS\t"$6-$NF+1"\t"$6"\t.\t+\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$5}' > REL606.getorf.noCTG.longest.gtf
cut -f -2 -d "\"" REL606.getorf.noCTG.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="-")' | awk '{print $2"\t.\tCDS\t"$5"\t"$5+$NF-1"\t.\t-\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$6}' >> REL606.getorf.noCTG.longest.gtf

#This version of the gtf still has the start codon coordinate of the original gtf as a 10th column
#Cases where there's a mismatch between the two:
awk -F '\t' '(($7=="+"&&$4!=$NF)||($7=="-"&&$5!=$NF))' REL606.getorf.noCTG.longest.gtf | wc -l
#8584 such cases, i.e. 12% of the total ORFs (n=72,395)

#To examine codon-delimited forms of the top five sequences from the longest ORFs:
awk -F '\t' '(($7=="+"&&$4!=$NF)||($7=="-"&&$5!=$NF))' REL606.getorf.noCTG.longest.gtf | sort -nk4 | gtf2bed | bedtools getfasta -s -name -fi REL606.faa -bed - | head |
awk '{if ($0 ~ /^>/) print; else {gsub(/.{3}/, "& "); print}}' | seqkit fx2tab | egrep "ATG|GTG|TTG"

#Same thing, but from the getorf-produced gtf:
awk -F '\t' '(($7=="+"&&$4!=$NF)||($7=="-"&&$5!=$NF))' REL606.getorf.noCTG.longest.gtf | cut -f -9 | cut -f2 -d "\"" | sed "s/.*/\"&\"/g" | #extract mismatched IDs
grep -F -f - REL606.getorf.noCTG.gtf | sort -nk 4 | gtf2bed | head | bedtools getfasta -s -name -fi REL606.faa -bed - | head | #search against getorf gtf
awk '{if ($0 ~ /^>/) print; else {gsub(/.{3}/, "& "); print}}' | seqkit fx2tab | egrep "ATG|TTG|GTG"

#egrep is used to highlight positions of all start codons
