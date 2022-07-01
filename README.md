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
--use_cov Yes \
--binary No \
--go_xref Yes \
--eggnog_category GO
```

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
