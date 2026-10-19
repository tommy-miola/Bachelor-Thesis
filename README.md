# Bachelor-Thesis
Code, data, and results derived from the replication and analysis of the scientific paper "The diffusive nature of housing prices" by J. P. Bouchaud et al. (https://arxiv.org/pdf/2412.14624).

## Data
The analysis is done with 2 different datasets:
1) Cagé J., & Piketty T. (2023). *Une histoire du conflit politique. Élections et inégalités sociales en France, 1789-2022*. Link to download: https://www.unehistoireduconflitpolitique.fr/download.html.<br>
I already cleaned the dataset and uploaded the two resulting zip files in the directory "1970-2022".
The file "data" contains the original data, "data_no_inflation" contains the same data but adjusted for inflation, according to the CPI values provided by INSEE, which are also reported in the same directory.
2) "Demandes de valeurs foncière (DVF)", published by the French government in the official site: https://www.data.gouv.fr/datasets/demandes-de-valeurs-foncieres.<br>
Since only the data of the last 5 years are visible, you can find the previous files here: https://files.opendatarchives.fr/cadastre.data.gouv.fr/data/etalab-dvf/latest/csv/.<br>
Given the size of the DVF datasets, you have to download them from the websites above.

## Code
The directory "1970-2022" contains 3 scripts, that take the first dataset as an input:
- spatial_variogram.py
- temporal_variogram.py
- diffusion_tgv.py

The directory "DVF" contains only the script "dvf.py", that does everything the paper does on DVF data.

Every script has specific CLI commands to execute multiple tasks, presented at the beginning of the code through a comment.<br>
Notice that the initial comments explain what the code does and how, and give specific instructions about the input and output files.

## Results
The directories contain the main charts obtained from the scripts, that match the work carried out in the original paper.