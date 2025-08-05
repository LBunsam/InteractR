[INTERACTR]


InteractR is a .R bassed script meant for multi dimesional analisis of Nanopore reads especifically tailored for Trypanozoma brucei (Tb). Meant to be used for analysis of reads from Yeast Double Hybrids (Y2H) protein-protein interaction essays. The pipeline is designed to take raw Nanopore reads, do quality control, remove of overrepresented sequences and alginment of the processed reads to a Tb. genome. The script requires the utilization of several files containing expresion, subcellular localization, and domain data for the visualization of graphs thorugh the Shiny App.




Bioinformatic Analysis Workflow Documentation

1. Computational Environment
All computational analyses were performed using R (version 4.4.1) within the RStudio (version 2023.12.1+402) environment on a Windows 10 operating system. The workflow leveraged several R packages and external tools:

R Packages:

-ShortRead (v1.62.0): For reading FASTQ files.
-Biostrings (v2.72.1): For sequence manipulation.
-dplyr (v1.1.4) and tidyr (v1.3.1): For data manipulation.
-readr (v2.1.5): For file input/output operations.
-ggplot2 (v3.5.1): For data visualization.
-UpSetR (v1.4.0): For generating UpSet plots.
-DT (v0.33): For interactive tables.
-gridExtra (v2.3): For PDF table generation.
-shiny (v1.8.1.1): For building the interactive visualization application.
External Tools:
-BLAST+ (version 2.15.0): For sequence alignment, specifically the blastn and makeblastdb commands.



2. Data Acquisition
The sequencing data from Polymerase Chain Reaction (PCR) products sequenced using Oxford Nanopore Technologies (ONT).
Raw sequence reads should be stored in FASTQ format, which includes nucleotide sequences and their associated quality scores.

-Sequencing Data: .FASTQ files generated from nanopore sequencing of PCR products.

-Reference Database: "unique_list.fasta", a curated database of nucleotide sequences.

-Localization Data: "localization2.csv" ([32], Nat Commun 14, 4401, 2023). Format: (Gene ID,Localization,Product Description,Category,CtermTryptag,NtermTryptag).

-Protein Expression Data: "protein_level.csv" ([43], Wellcome Open Res, 2023). Fromat : (GeneID,bsf,pcf).

-Functional Annotations: "filtered_results.csv" with InterProScan domain annotations. Format: (Gene ID,InterproScanID,Definition,start,end,E_value).



3. Sequence Processing Workflow

3.1. FASTQ to FASTA Conversion
The script automatically detects all FASTQ files in the working directory with the .fastq pattern. Each FASTQ file should contain the reads sequences from a genes that interact in the essay. ShortRead package is used to read, extract, reformat the headers from FASTQ (@) to FASTA (>) format. The resulting sequences were saved as FASTA files using the Biostrings package.

3.2. BLAST Database Creation
A reference FASTA file (unique_list.fasta) contains a curated set of sequences from Tb. used to create a BLAST nucleotide database named Unique-CC-DB. The makeblastdb command from BLAST+ generates the database, with a check in place to prevent redundant creation if the database already exists.

3.3. BLASTN Analysis
For each converted FASTA file, a BLASTN search was performed against the Unique-CC-DB database using the blastn command. The output is formatted as a tabular file (-outfmt 6) containing the query sequence ID, subject (prey) sequence ID, and E-value. Results are saved as text files (e.g., <gene_id>_blast_results.txt).



4. BLAST Result Processing

4.1. Analysis and Filtering
BLAST results are parsed using the readr package, explicitly defining columns as Query, Prey, and E_value. For each queried sequence, only the prey with the lowest E-value is retained using the dplyr package to ensure high-confidence matches.

4.2. Contaminant Removal
Sequences listed as contaminants (in contaminants.txt) are filtered out from the BLAST results.

4.3. Data Consolidation
Results from all FASTQ files are combined into a single dataframe, ensuring no duplicate bait-prey pairs remain. The final interactors are saved to a CSV file (all_final_interactors.csv).



5. Integration of Additional Data

5.1. Subcellular Localization
Localization data from localization2.csv ([32] Nat Commun 14, 4401, 2023) is merged with the BLAST results to annotate prey sequences with their subcellular localization. Missing or unknown localizations were labeled as "Unknown". 

5.2. Protein Levels
Protein expression data from protein_level.csv ([43], Wellcome Open Res, 2023) is used for visualization in scatter plots of protein levels expression in bloodstream form (BSF) and procyclic form (PCF).

5.3. InterProScan Annotations
Functional annotations from filtered_results.csv are integrated to identify protein domains associated with the prey sequences, filtering out empty or invalid entries.



6. Data Visualization
The uses an interactive Shiny application developed to visualize the processed data. The application includes the following functionalities:

6.1. UpSet Plot
An UpSet plot is generated using the UpSetR package to display the overlap of prey sequences among selected baits, requiring at least two baits for meaningful visualization.

6.2. Subcellular Localization Bar Plot
A bar plot is created with ggplot2 showing the percentage distribution of prey localizations for selected baits, calculated in relation to the total localization counts in Tb from (Nat Commun 14, 4401, 2023).

6.3. Protein Level Scatter Plot
A scatter plot visualing BSF versus PCF protein levels for a selected bait, its prey, and the rest of the population, with distinct colors for each group (red for bait, blue for specific prey, grey for others).

6.4. InterProScan Analysis
A bar plot displaying the top 20 most enriched InterProScan domains based on unique prey counts, and an interactive table (using the DT package) lists bait-prey pairs with their E-values, localizations, and domain annotations.

6.5. Downloadable Outputs
Users can download UpSet plots in PDF format, scatter plots in PDF, and tables of unique preys (with localization and InterProScan data) in CSV format for each bait.

