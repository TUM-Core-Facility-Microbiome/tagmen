# Tagmen Suite
## Steps
### 0) Install dependencies
Install python dependencies:
```shell script
pip install -r requirements.txt
```
Also make sure cd-hit is installed.

### 1) Extract 15N sequences from RT files
Move all RT sample files you want to analyze to a folder.

Extract sequences:
```shell script
./readtag.py <INPUT_FOLDER_RT> <OUTPUT_FOLDER_RT>
```

This will generate the output folder.
### 2) Run CD-HIT-EST
#### 2a) Generate `FASTA` file from `15N.tsv` file
```shell script
cd <OUTPUT_FOLDER_RT>
grep -v '#' 341-RT-T1_S341_L001.15N.tsv | awk '{OFS="\t"; print ">"$1"\n"$2}' > 341-RT-T1_S341_L001.15N.fasta
```

#### 2b) Run CD-Hit-EST
```shell script
cd <OUTPUT_FOLDER_RT>
cd-hit-est -d 0 -i 341-RT-T1_S341_L001.15N.fasta -o 341-RT-T1_S341_L001.15N.cdhit
```

#### 2c) Run parse clusters 
```shell script
./parse_cdhit_clusters.py <OUTPUT_FOLDER_RT>/341-RT-T1_S341_L001.15N.tsv <OUTPUT_FOLDER_RT>/341-RT-T1_S341_L001.15N.cdhit <OUTPUT_FOLDER_RT>/341-RT-T1_S341_L001.15N.cdhit.clstr
```

### 3) Extract 15N pairs from LT files
Move all LT sample files you want to analyze to a folder.

```shell script
./linktag.py <INPUT_FOLDER_LT> <OUTPUT_FOLDER_LT>
```

### 4) Aggregate results
```shell script
./aggregate.py <OUTPUT_FOLDER_LT>/337-LT-T1_S337_L001.15Npairs.tsv <INPUT_FOLDER_RT> <OUTPUT_FOLDER_RT>/341-RT-T1_S341_L001.15N.clusters.tsv <OUTPUT_FOLDER_AGGREGATE>
```
