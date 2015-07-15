# pythonGeneSeekr

The geneSeekr has very few requirements:
 
 1. Python
 2. BioPython
 3. Sequence (genome) files
 4. Target files
 5. ?
 
 There are three flags that must be provided when running the program from a system other than my own. These flag override hard-coded paths.
 
 * -p - the location you want the report folder created and populated
 * -s - the location of your query sequence files (must have a .fa* extension)
 * -t - the location of your target sequence files (must have a .fa* extension)
 
 There is one additional flag, which allows you to customise the analysis
 
  * -c - the percent identity cutoff value to determine presence/absence of gene targets
