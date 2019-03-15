# Haplotype-phasing

input: chr22.matrix(snp binary matrix) format see tests

output: chr22_post_haplotype


python2 /yourpath/dchap.py  --reads chr22.matrix --k 2 --parsed-reads chr22_parsed_fragments --phase chr22_haplotype --assignments chr22_assignments

python2 /yourpath/dchap-postprocess.py --parsed-reads chr22_parsed_fragments --assignments chr22_assignments --blocks chr22_haplotype --corrected-blocks chr22_post_haplotype
