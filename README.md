# Metagenome annotation mapper

This is mainly a little script for the KBase Knowledge Engine project that maps metagenome annotations from eggNOG to a feature table, with the option of mapping from Gene Ontology to other annotation namespaces.

### Usage
1. Download this repo.
2. Example run with python:

```
eggnog_mapper/eggnog_to_feature_table.py \
-i <annotation_file> \
--input_cov <contig coverage file> \
--input_contig_fasta <contigs file> \
--use_cov Yes \
--binary No \
--go_xref Yes \
--eggnog_category GO
```
