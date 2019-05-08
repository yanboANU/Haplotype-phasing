# DCHap: A divide-and-conquer haplotype phasing algorithm for third-generation sequencing data

DCHap is a fast and accurate haplotype phasing tool for third-generation sequencing data.

# Dependencies


# Downloading DCHap

To download DCHap, you have to clone the DCHap repository to your machine.
<pre><code> git clone https://github.com/yanboANU/Haplotype-phasing.git </code></pre>

# Input Format
input: chr22.matrix(snp binary matrix) format see tests


# Example Usage

<pre><code> cd tests
python2 /yourpath/dchap.py  --reads chr22.matrix --k 2 --parsed-reads chr22_parsed_fragments --phase chr22_haplotype --assignments chr22_assignments
python2 /yourpath/dchap-postprocess.py --parsed-reads chr22_parsed_fragments --assignments chr22_assignments --blocks chr22_haplotype --corrected-blocks chr22_post_haplotype </code></pre>

Phasing result is in chr22_post_haplotype.

# References
[1] Bansal, Vikas, and Vineet Bafna. "HapCUT: an efficient and accurate algorithm for the haplotype assembly problem." Bioinformatics 24.16 (2008): i153-i159.

[2] Duitama, Jorge, et al. "ReFHap: a reliable and fast algorithm for single individual haplotyping." Proceedings of the First ACM International Conference on Bioinformatics and Computational Biology. ACM, 2010.

[3] Edge, Peter, Vineet Bafna, and Vikas Bansal. "HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies." Genome research 27.5 (2017): 801-812.

[4] Kuleshov, Volodymyr. "Probabilistic single-individual haplotyping." Bioinformatics 30.17 (2014): i379-i385.
