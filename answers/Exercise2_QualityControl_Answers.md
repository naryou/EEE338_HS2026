# Answers of Excercises: Variant calling #1

Make sure you are in `VarCall` directory.

## Question 1: How many sequences are there in the reference sequence file (MedtrChr2.fa)?

```bash
$ grep -c '^>' 00_input/MedtrChr2.fa
1
$
```

`^` means the beginning of the line.

* * *

### Question 2: How many reads does each fastq file have (\*\_R1.fastq.gz)?

Check the sequence ID first.

```bash
$ zless 00_input/516950_chr2_R1.fastq.gz
@700523F:115:CAV90ANXX:6:2109:12367:86449/1
CAATACCAGTCTGCCTGTTTCAAAAAGCTACAAATAATTAAGTGATTCTAATCAAACTACTACTGATAAGGAAGGATTCTGCTATTCAGAATCTTCACGATAAAGAAATAAAACTACTGCTGATGA
+
CCCBCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGEGGGGGGGGGGGGGGGGGGGGGGDGGGGGGGGGGGGEGGGGGGGGGCGGB
@700523F:115:CAV90ANXX:6:2311:1825:75411/1
AATATACACTACTACTTATATAGGCACATAGCACCAGAGGAACAAGAAGACAACAGCTTAAATCTTGAATTACTGATTCATTAGGATTATTTTATCTTTTTAAACTTCTGTTAATTAGGAATTGTT
+
CCCCCGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGEGGGGGGGGGGGGGGGGGGGGGGGGGGEGGGGGGGGGGGGGGFGGGGGGGGGEGFGGGGGEGGGGGEEGGGGGCFGGGGGEGG>FGGGGF
@700523F:115:CAV90ANXX:6:1208:2916:99820/1
AAGAAGACAACAGCTTAAATCTTGAATTACTGATTCATTAGGATTATTTTATCTTTTTAAACTTCTGTTAATTAGGAATTGTTGATTAGGATTATTTGTCTCTGATAAATAATGCTTTTGAAGTCT
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGEBGGGGGGGGGCGEGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGEGGG?7@FGGGEGGGGBCBG>F
@700523F:115:CAV90ANXX:6:2307:11848:93404/1
ACAAGAAGACAACAGCTTAAATCTTGAATTACTGATTCATTAGGATTATTTTATCTTTTTAAACTTCTGTTAATTAGGAATTGTTGATTAGGATTATTTGTCTCTGATA
・
・
・
```

You can see all sequence IDs start with `@700`.
So, you can count the number of the sequences with counting `@700`.

```bash
$ zless 00_input/516950_chr2_R1.fastq.gz | grep "^@700" | wc -l
165033
$ zless 00_input/516950_chr2_R2.fastq.gz | grep "^@700" | wc -l
165033
$ zless 00_input/660389_chr2_R1.fastq.gz | grep "^@700" | wc -l
157539
$ zless 00_input/660389_chr2_R2.fastq.gz | grep "^@700" | wc -l
157539
$
```


Another option is counting lines and divide it by 4.  
Each sequence in fastq file has 4 lines of information.  
So, you can count the number of sequences by counting the lines and dividing them by 4.  
Before that, you need to enter the input directory. And you have to use some bash script to divide them.

```bash
$ cd 00_input
$ echo $(( $(zcat 516950_chr2_R1.fastq.gz | wc -l) /4 ))
165033
$ echo $(( $(zcat 516950_chr2_R2.fastq.gz | wc -l) /4 ))
165033
$ echo $(( $(zcat 660389_chr2_R1.fastq.gz | wc -l) /4 ))
157539
$ echo $(( $(zcat 660389_chr2_R2.fastq.gz | wc -l) /4 ))
157539
$
```


* * *

### Question 3: Does each sample have the same number of R1 and R2 reads? (Caution: Q scores can be + or @)

Yes — by design, they must match. R1 and R2 are the two ends of the same physical DNA fragment, sequenced from opposite directions in the same paired-end run. Every cluster on the flow cell that gets sequenced produces exactly one R1 read and one R2 read together, as a pair. (Claude)

* * *

### Question 4: How many bases are in the reference sequence?

```bash
$ grep -v "^>" 00_input/MedtrChr2.fa | tr -d '\n' | wc -c
2000001
$
```

* * *

### Question 5: How many missing bases (N)?

You can use `grep` command with `-o` option.  
`-o` option returns only the matched (non-empty) parts of a matching line, with each such part on a separate output line.  
So if the sequence looks like below,

```text
ATGATTNATAGCNNNATGATNGTACCATCAT
```
`grep -o "N"` returns

```text
N
N
N
N
N
```

So, you can count the line after grepping to count the number of missing bases (N).

```shell
grep -v '^>' 00_input/MedtrChr2.fa | tr -cd 'N' | wc -c
46942
$
```
