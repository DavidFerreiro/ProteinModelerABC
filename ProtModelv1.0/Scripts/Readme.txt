Before starting a ProtModel analysis the user must do some extra work. 

	1. First, although the framework uses phy format, the user needs his/her protein alignment in fasta format. 
	2. Then, we recommend to use SWISS-MODEL to find the best template (https://swissmodel.expasy.org). If the MSA are longer than 10 sequences SWISS-MODEL will fail so, in this case, we recommend to use the “Find-WT.py” python script	, which return wild-type (a sequence composed by the most common amino acid per site). 
	
		python FindWT.py --input MSA.fasta

	3. Next, the WT sequences will be used to find the template in SWISS-MODEL that better represents the alignment. Usually, the user must download the first template of the results screen, but it is recommended to check the X-ray to avoid high values when it’s possible. 
	4. Work with template structures demands that the MSA has to have the same length as template chain. So, sequences have to be aligned with the template chain sequence and remove the positions that in the template sequence has a gap in every sequence. The “Align.py” script placed in the folder “Scripts” will make the alignment, remove the gaps positions and change the file into phylip format. Align.py” works with an input fasta MSA (--input), a template structure (--temp), the chain of the template (--chain) and the desired output file name, which must have the .phy extension (--output).
	
			python Align.py --input MSA.fasta --temp structure.pdb --chain A to Z --output MSA.phy
	
	5. Another option could be do all of this using other programs as “Muscle” (BIBLIO) for sequences alignment or Phylogeny.fr (Dereeper et al., 2008) to convert into phylip format. “Align.py” works with an input fasta MSA (--input), a template structure (--temp), the chain of the template (--chain) and the desired output file name, which must have the .phy extension (--output).
	6. To perform the simulations the user must choose a substitution rate per site prior distribution. This could be an unfamiliar measure, so we also include a script which for a desired theta (θ) value will return the sequences identity of the MSA and substitution rate per site, asking if the organism is haploid or diploid and the population size. Sequence identity will give us an idea about θ value, if it is very high, we must not use a high theta (>500).
			
			Theta.py --input NS1.phy
