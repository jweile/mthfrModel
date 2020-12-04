#Load Human dimer structure
fetch 6fcx, async=0
extract humanDimerPart1, chain A
extract humanDimerPart2, chain B
delete 6fcx

#Load E.coli monomer bound to anti-folate
fetch 2fmn, async=0
extract ecoliMonomer, chain A & 2fmn
delete 2fmn

#visualize in cartoon style and ligands as spheres
hide all
show cartoon, all
show spheres, organic

#run alignment
align  ecoliMonomer, humanDimerPart1, object=alobj
save strc_align.aln, alobj

#save structure model
save mthfr_crudeModel.pdb, (humanDimerPart1 & !segi E | (ecoliMonomer & segi E))

#find residues contacting (within 5A) the anti-folate
select activesite, byres all within 5 of (ecoliMonomer & segi E)
iterate activesite, print (model,resi,resn)

#find residues contacting the FAD
select flavo, byres humanDimerPart1 within 5 of (humanDimerPart1 & segi C)
iterate flavo, print(model, resi, resn)

#find residues contacting SAH
select inhibitor, byres humanDimerPart1 within 5 of (humanDimerPart1 & segi D)
iterate inhibitor, print(model, resi, resn)