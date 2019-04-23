# PopPAnTe: Population and Pedigree Association Testing for Quantitative Data

PopPAnTe is a user-friendly framework enable pairwise association testing of quantitative *-omics* variables in family-based study. Relationships between individuals can be either described by known family structures of any size and complexity, or by genetic similarity matrices (GSM) inferred from genome-wide genetic data. This approach is particularly useful when some degree of hidden relatedness (including population stratification) is expected, but extensive genealogical information is missing or incomplete. For instance, genealogical information going back more than three or four generations may be difficult to be retrieved for individuals recruited in large-scale biobank started in genetic isolates.

PopPAnTe models the data in a variance component framework to keep into account the resemblance among individuals, supports region-based testing, assesses the significance of the association through a formal likelihood ratio testing as well as through an adaptive permutation procedure, and performs basic data pre- and post-processing.

PopPAnTe is now at version 1.0.2 (released on March 20th, 2018).

## Table of contents

- [Citation](#Citation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [License](#license)
- [Acknowledgements](#acknowledgements)


## Citation

If you use PopPAnTe for research purpose, please cite:

> Visconti, Alessia, et al. *“PopPAnTe: population and pedigree association testing for quantitative data.”*, BMC genomics 18.1 (2017): 150, [DOI:10.1186/s12864-017-3527-7](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3527-7)

## Dependencies

To run PopPAnTe you will need to install [Java 7](http://www.oracle.com/technetwork/java/javase/downloads/index.html) or higher. Please note that Java should be already available in most of the POSIX compatible systems (Linux, Solaris, OS X, etc).
If you would like to compile PopPAnTe from source, you will also need [Apache Maven](https://maven.apache.org/index.html).

## Usage 

This repository includes a pre-compiled Java ARchive (JAR) file, which can be use to run PopPAnTe without any need to compile the code. You can download it [here](./bin/poppante-1.0.2.jar). 
More details on how to use PopPAnTe are available in the [PopPAnTe Users' guide](./doc/PopPAnTe-Users_guide.pdf), while the examples used in the documentation, and which can be used to familiarise oneself with the tool, are available in the [Example](./Example) folder in this repository.

## License

PopPAnTe is licensed under GNU GPL v3.

## Acknowledgements

PopPAnTe includes and modifies procedures developed by Hariklia Eleftherohorinou and Mario Falchi (variance component core, probability distribution functions), Sam Halliday ([MTJ libraries](https://github.com/fommil/matrix-toolkits-java)), Carlos Morcillo Suarez ([Q-Q and Manhattan plots](https://github.com/CarlosMorcilloSuarez/ManhattanPlotter)), Ivan Akimov & *fanweixiao* (Hashids), and Mateusz Kobos ([PCA](https://github.com/mkobos/pca_transform)). When using their code, I kept the authorship in the file. I would like to thank them for allowing the free download, modification, and distribution of their work accordingly to the terms of public licenses. 




