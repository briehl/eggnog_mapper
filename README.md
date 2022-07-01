# Metagenome annotation mapper

This is mainly a little script for the KBase Knowledge Engine project that maps metagenome annotations from eggNOG to a feature table, with the option of mapping from Gene Ontology to other annotation namespaces.

## Usage
1. Download this repo.
2. Example run with python:
```
eggnog_mapper/eggnog_to_feature_table.py \
-i <annotation_file> \
--input_cov <contig coverage file> \
--input_contig_fasta <contigs file> \
--use_cov \
--binary No \
--go_xref Yes \
--eggnog_category GO
```
3. This creates the following output files (formats below):
   1. a contig length summary file
   2. an overall contig summary file including:
      1. contig id
      2. contig length
      3. feature hit frequency - calculated as the number of annotations of the requested eggNOG category divided by the number of annotations of that contig. E.g. if you're looking for GO annotations on a contig with 10 annotations, and there are only 3 contig genes with at least one GO hit, then that frequency is 0.3 for that contig.
      4. average contig coverage (optional, if `--use_cov` is used)
      5. consensus taxonomy for that contig (most detailed consensus)
      6. the frequency that taxonomy shows up in the contig
   3. a count table for the given eggnog_category (GO above)
      1. This table can be weighted or unweighted based on contig coverage
   4. if `--go_xref` is used, then a series of files - weighted or unweighted - is created, one for each GO cross reference. The rows have two columns - the annotation id, and the (weighted) number of times it appears in the metagenome

## All input parameters and flags
* `-i` required - Path to the annotation file
* `--input_cov` - path to the input contig coverage file
* `--input_contig_fasta` - path to the input FASTA file
* `--contig_filter` - str (optional) - contigs to retain for the output (i.e. contigs not listed will be filtered before count tables produced). If not used, all contigs (that meet other criteria below) will be used.
* `--use_cov` - is present, use the contig coverage to create weighted output files. These weights are a multiplication of the number of appearances of an annotation by the contig coverage
* `--binary` - if "Yes", make a binary representation of the annotation counts - just 1 or 0 if that annotation is present
* `--go_xref` - if "Yes", create cross-reference tables from the Gene Ontology values in the annotation to various other annotations (see below for a list)
  * EC
  * HAMAP
  * InterPro
  * KEGG
  * MetaCyc
  * Pfam
  * PIRSF
  * PROSITE
  * Reactome
  * RFam
  * Rhea
  * SMART
  * UM-BBD EnzymeId
  * UM-BBD PathwayId
  * UM-BBD ReactionId
  * UniProt KeyWord
  * UniProt Subcellular Location
  * UniRule
* `--min_contig_length` - int, default 2000 - minimum contig length used for generating outputs
* `--min_contig_coverage` - int, default 5 - minimum contig coverage used for generating outputs
* `--eggnog_category` - string, default "GO" - the eggNOG mapper ontology to use for mapping to the output files. Allowed options are:
  * ALL_CATEGORIES
  * GO
  * EC
  * KEGG_ko
  * COG
  * KEGG_Module
  * KEGG_Reaction
  * KEGG_rclass
  * BRITE
  * KEGG_TC
  * CAZy
  * BiGG_Reaction
  * eggNOG_OGs
* `--contig_taxa_threshold` - float, default 0.5 - the frequency threshold for determining contig consensus taxonomy


## Input files
### Annotation File
This is expected to be a TSV with the following columns, in order. Note that not all of these are expected to be present.
* `query_name` - string - these are expected to be contig query ids. In our datasets, these strings have the format `<unique_id>_contig_<contig_index>_<annotation>`. For example: `1085605_contig_1_32` has a unique id 1085605, this is contig index 1, annotation 32
* `seed_eggNOG_ortholog` - string - eggNOG ortholog ids
* `seed_ortholog_evalue` - float - ortholog E value 
* `seed_ortholog_score` - float - ortholog score 
* `best_tax_level` - string - name of the deepest taxonomic level identified e.g. "Actinobacteria" not "phylum"
* `Preferred_name` - string - preferred gene or protein name
* `GO` - comma-separated set of GO (Gene Ontology) ids 
* `EC` - EC (Enzyme Commission) number 
* `KEGG_ko` - KEGG (Kyoto Encyclopedia of Genes and Genomes) orthology id 
* `KEGG_Pathway` - comma-separated set of KEGG pathway ids 
* `KEGG_Module` - comma-separated set of KEGG module ids 
* `KEGG_Reaction` - comma-separated set of KEGG reaction ids where this gene/protein plays a role 
* `KEGG_rclass` - comma-separated set of KEGG reaction classes where this gene/protein plays a role
* `BRITE` - comma-separated set of KEGG BRITE hierarchy identifiers
* `KEGG_TC` - comma-separated set of KEGG TC identifiers (## TODO), 
* `CAZy` - CAZy (Carbohydrate-Active Enzyme) DB family id, 
* `BiGG_Reaction` - reaction from the BiGG database where this gene/protein plays a role, 
* `taxonomic_scope` - high-level scope of the source of the annotation (TODO get list of possible values), including:
  * Archaea
  * Bacteria
  * Eukaryota
  * Fungi
  * Mammalia
  * Metazoa
  * Opisthokonta
  * Nematoda
  * Viridiplantae
  * Viruses
* `eggNOG_OGs` - comma-separated list of eggNOG orthologous group ids 
* `best_eggNOG_OG` - 
* `COG` - COG category 
* `eggNOG_free_text_desc` - free text description from eggNOg

### Contig coverage file
This is a simple TSV where each row has a contig identifier (similar to `query_name` from the annotation file) and an integer value for the read coverage.

E.g.:
```
1085605_contig_1	28
1085605_contig_2	122
1085605_contig_3	31
1085605_contig_4	81
1085605_contig_5	24
```

### Contig FASTA file
A multi-contig FASTA file where each description line is a contig id (matching those from the coverage file and annotation file). This gets parsed for the contig lengths using BioPython.
