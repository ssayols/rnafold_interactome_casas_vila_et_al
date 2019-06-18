Code assocaited with XXXXX.

### Contents

* __analysis:__ analysis of the MS SILAC experiment, domains, localization, disorder prediction, genetic interactions (Cellmap)
* __constructs:__ positional and functional annotation of the constructs
* __cytometry:__ analysis of the experiment with different stress conditions
* __DMS-seq\_data:__ analysis of the structural information from the DMS-seq data
* __functional\_follow\_up:__ analysis of modifications in known/putative binding domains
* __primers:__ primer design for the construct sequences

### Data
Data files larger tahn 10MB (uncompressed) are not included in this repo, and available upon request.
Still, scripts to generate most data are included.

```sh
find . -size -10M \( -name \*.sh -o -name \*.R -o -name \*.txt -o -name \*.csv \) -exec tar rvf x.tar --exclude "./src/*" --exclude "./bin/*" --exclude "./tmp/*" '{}' \;
```
