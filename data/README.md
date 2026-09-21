# Data

`genome_manifest.tsv` lists every genome that was screened, one row per genome: the screen it
belonged to, the genus and species, the accession, the NCBI organism name, the RefSeq category,
the assembly level, the number of protein-coding genes and the NCBI directory the protein file came
from.

20,435 prokaryotic genomes were selected, one per named species, and 20,434 were screened; one
*Paracoccus* genome carried no protein annotation. The 88 eukaryote genomes are listed after them.

The protein file for a genome is `<ftp_path>/<basename>_protein.faa.gz`, where `<basename>` is the
last part of `ftp_path`. Its NCBI page is `https://www.ncbi.nlm.nih.gov/datasets/genome/<accession>/`.
Nothing here is redistributed: the manifest points at NCBI, so anyone can fetch exactly the same
genomes. `survey_scripts/fetch_proteomes.py` does that.
