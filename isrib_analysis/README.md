# SLURM job 7187289
### Command used from XPRESSpipe https://github.com/XPRESSyourself/XPRESSpipe/commit/552b5351e9e1f78ceabf66426cdd855b2aa349be
```
xpresspipe riboseq -i $SCRDIR/input -o $SCRDIR/output -r $REF --gtf $REF/transcripts_LCT.gtf -e isrib_riboprof -a CTGTAGGCACCATCAAT --method RPKM --sjdbOverhang 49 --quantification_method htseq
```
- New truncator that is fixed and parsed CDS
- Longest transcripts are Ensembl canonical transcripts
- Use htseq and riboseq module quantifies along CDS records (for both ribo and rna)
- No de-duplication for analysis
- Failed quality control execution due to OOM error


### Fig 3A:
- ISRIB appears to constrict translation
- Discuss how validated good run
- Why and how I set cut-off used
- Explain in validation section how the pipeline worked and what I did post-processing

### Neuro enrichment
- Run KS test
- Get random number of genes and see how many are neuro-annotated
- State these may be interesting targets to follow up on
- Explain caveats
- A lot of sodium and calcium related stuff, may play a role in neural signaling
- Scoring paradigm:
  - **GeneName** = overexpressed, or overexpressed and embryonic evidence and/or neurological function
  - *GeneName* = Only one of the above, or association   
- Need to automate this, and bootstrap
- Maybe make a better classifier to pick up pattern and not be limited to threshold

```
data_te_zero.loc[(data_te_zero['tm'] <= -1) & (data_te_zero['tmisrib'] >= 0)].index.tolist()

['AC097634.4',
 'ADAP1',
 'ADRA2C',
 'AL109811.3',
 'AL358472.7',
 'CMYA5',
 'CXorf40B',
 'DCDC2B',
 'ECHDC2',
 'EFCAB13',
 'GNB3',
 'GPR183',
 'HTRA1',
 'KIFC2',
 'LRP5L',
 'MYO5B',
 'RAD51',
 'RGPD8',
 'RPS27',
 'SCART1',
 'SLC1A1',
 'SLC9A3',
 'TLL1',
 'TMEM110-MUSTN1',
 'WNK4']
 ```

- AC097634.4
- **ADAP1**
  - GTPase activator activity and inositol 1,3,4,5 tetrakisphosphate binding
  - gene is overexpressed in Fetal Brain (15.1), Frontal cortex (12.8), Spinal cord (11.8), and Stomach (8.4)
- **ADRA2C**
  - G protein-coupled receptor superfamily
  - These receptors have a critical role in regulating neurotransmitter release from sympathetic nerves and from adrenergic neurons in the central nervous system
  - The mouse studies revealed that both the alpha2A and alpha2C subtypes were required for normal presynaptic control of transmitter release from sympathetic nerves in the heart and from central noradrenergic neurons
- AL109811.3
- AL358472.7
- CMYA5
  - Cardiomyopathy Associated 5
  - Diseases associated with CMYA5 include Tibial Muscular Dystrophy and Muscular Dystrophy
- *CXorf40B*
  - Embryonic expression in neural tube (may just be since HEK cells are progenitor-like)
- *DCDC2B*
  - bind tubulin and enhance microtubule polymerization
  - The doublecortin gene family and disorders of neuronal structure. (PMID: 20236041)
- *ECHDC2*
  - related pathways are Fatty Acid Biosynthesis
  - lyase activity
  - Embryonic expression in brain and neural tube
- EFCAB13
  - calcium ion binding
- **GNB3**
  - G Protein Subunit Beta 3
  - Heterotrimeric guanine nucleotide-binding proteins (G proteins), which integrate signals between receptors and effector proteins
  - Transmembrane signaling
  - Embryonic expression in brain, neurons, eye
  - Overexpressed in Cerebellum
  - GO term -- neuron projection (GO:0043005)
- **GPR183**
  - G Protein-Coupled Receptor 183
  - The function of this gene is unknown (Entrez)
  - G-protein coupled receptor expressed in lymphocytes that acts as a chemotactic receptor for B-cells, T-cells, splenic dendritic cells, monocytes/macrophages and astrocytes (By similarity)
  - Promotes follicular helper T (Tfh) cells differentiation by positioning activated T-cells at the follicle-T-zone interface, promoting contact of newly activated CD4 T-cells with activated dendritic cells and exposing them to Tfh-cell-promoting inducible costimulator (ICOS) ligand (By similarity)
  - Expression in splenic dendritic cells is required for their homeostasis, localization and ability to induce B- and T-cell responses: GPR183 acts as a chemotactic receptor in dendritic cells that mediates the accumulation of CD4(+) dendritic cells in bridging channels (By similarity)
  - Regulates migration of astrocytes and is involved in communication between astrocytes and macrophages (PubMed:25297897)
  - weakly expressed in brain (PubMed:16540462, PubMed:8383238). Expressed in astrocytes (PubMed:25297897)
- **HTRA1**
  - secreted enzyme that is proposed to regulate the availability of insulin-like growth factors (IGFs) by cleaving IGF-binding proteins
  - may also regulate cell growth
  - Variations in the promoter region of this gene are the cause of susceptibility to age-related macular degeneration type 7
  - Diseases associated with HTRA1 include Cerebral Arteriopathy, Autosomal Recessive, With Subcortical Infarcts And Leukoencephalopathy and Cerebral Arteriopathy, Autosomal Dominant, With Subcortical Infarcts And Leukoencephalopathy, Type 2.
  - Among its related pathways are Degradation of the extracellular matrix
  - Arteriosclerosis (increased thickness, increased stiffness, loss of elasticity) of the small arteries of the brain.
  - may regulate many physiological processes, including retinal angiogenesis and neuronal survival and maturation during development
- **KIFC2**
  - ATPase activity and microtubule motor activity
  - May play a role in microtubule-dependent retrograde axonal transport. May function as the motor for the transport of multivesicular body (MVB)-like organelles in dendrites (By similarity).
  - Embryonic expression in brain and neural tube
  - Overexpressed in brain
- *LRP5L*
  - Diseases associated with LRP5L include Epilepsy
  - annotations related to this gene include Wnt-protein binding and Wnt-activated receptor activity
- MYO5B
  - plasma membrane recycling
  - transferrin receptor recycling
  - Together with RAB11A and RAB8A participates in epithelial cell polarization
  - Together with RAB25 regulates transcytosis
- **RAD51**
  - involved in the homologous recombination and repair of DNA
  - Loss of these controls following BRCA2 inactivation may be a key event leading to genomic instability and tumorigenesis
  - Embryonic expression in brain and neural tube
  - overexpressed in Fetal Brain
- *RGPD8*
  - Ran GTPase binding
  - Embryonic expression in brain and neural tube
  - Behavior/neurological phenotype
- RPS27
  - Small Ribosomal Subunit Protein ES27, 40S ribosome
  - Mutations in this gene have been identified in numerous melanoma patients and in at least one patient with Diamond-Blackfan anemia (DBA)
- SCART1
  - Pseudogene
- **SLC1A1**
  - glutamate transporters
  -  In brain, these transporters are crucial in terminating the postsynaptic action of the neurotransmitter glutamate, and in maintaining extracellular glutamate concentrations below neurotoxic levels
  - Glutamate transporters, also known as excitatory amino acid transporters (EAATs), are sodium- and potassium-dependent members of the solute carrier family 6 (SLC1), widely distributed throughout the brain
  - Embryonic expression in brain
  - Expression in neurons and brain (in which there was dense expression in substantia nigra, red nucleus, hippocampus and in cerebral cortical layers)
- SLC9A3
  - uses an inward sodium ion gradient to expel acids from the cell
  - Involved in pH regulation to eliminate acids generated by active metabolism or to counter adverse environmental conditions
  - Plays an important role in signal transduction
- **TLL1**
  - astacin-like, zinc-dependent, metalloprotease
  - Studies in mice suggest that this gene plays multiple roles in the development of mammalian heart, and is essential for the formation of the interventricular septum
  - Embryonic expression in brain, neural tube, and neural ectoderm
  - Overexpressed in brain, cerebellar hemisphere, and Cerebellum
  - Strong heart annotation and overexpression
- *TMEM110-MUSTN1*
  - read-through transcript
  - schizophrenia-associated genetic loci. (PMID: 25056061)
- WNK4
  - regulates the balance between NaCl reabsorption and K(+) secretion
  - electrolyte homeostasis, cell signaling, survival and proliferation
  - Acts as an activator and inhibitor of sodium-coupled chloride cotransporters and potassium-coupled chloride cotransporters respectively
  - WNK4 appears to act as a molecular switch that can vary the balance between NaCl reabsorption and K(+) secretion to maintain integrated homeostasis



```
import random
random.sample(data_te_zero.index.tolist(), k=25)

['SLCO4A1',
 'KREMEN1',
 'LMO4',
 'ACAP3',
 'NCBP2',
 'VCPIP1',
 'TMEM267',
 'R3HDM4',
 'ATXN2L',
 'RAB5C',
 'TANGO2',
 'ADSL',
 'EFR3A',
 'MBNL1',
 'RMND5B',
 'BUB1B',
 'HIRA',
 'EIF4E',
 'VPS25',
 'MMAA',
 'SFMBT1',
 'AEBP2',
 'POLR2E',
 'SLC25A24',
 'RAB5B']
```




- *SLCO4A1*
  - Sodium-Independent Organic Anion Transporter E
  - Embryonic expression in brain and neural tube
- *KREMEN1*
  - modulates canonical WNT signaling
  - Embryonic expression in brain and neural tube
- **LMO4**
  - may play a role as a transcriptional regulator or as an oncogene
  - Overexpressed in brain and cortex
  - Embryonic expression in brain
  - Behavior/neurological phenotype
  - GO Terms: neural tube closure (GO:0001843), ventral spinal cord interneuron differentiation (GO:0021514), spinal cord motor neuron differentiation (GO:0021522)
- **ACAP3**
  - GTPase activator activity, endocytosis
  - overexpressed in Brain - Cerebellum
  - regulates neuronal migration in the developing cerebral cortex. (PMID: 28919417)
  - neuron migration (GO:0001764), regulation of neuron projection development (GO:0010975)
- *NCBP2*
  - component of the nuclear cap-binding protein complex (CBC)
  - behavior/neurological phenotype
- **VCPIP1**
  - Acts as a deubiquitinating enzyme
  - Overexpressed in brain
- TMEM267
  - Transmembrane Protein
- R3HDM4
  - nucleic acid binding
- *ATXN2L*
  - ataxin type 2 related protein of unknown function
  - associated with a complex group of neurodegenerative disorders
- *RAB5C*
  - ensure fidelity in the process of docking and/or fusion of vesicles with their correct acceptor compartment
  - Embryonic expression in brain
- *TANGO2*
  - secretory protein loading in the endoplasmic reticulum
  - Allelic variants of this gene are associated with rhabdomyolysis, metabolic crises with encephalopathy, and cardiac arrhythmia
  - Diseases associated with TANGO2 include Metabolic Encephalomyopathic Crises, Recurrent, With Rhabdomyolysis, Cardiac Arrhythmias, And Neurodegeneration and Cardiac Arrhythmia
- *ADSL*
  - essential enzyme involved in purine metabolism
  - Mutations in this gene are associated with adenylosuccinase deficiency (ADSLD), a disorder marked with psychomotor retardation, epilepsy or autistic features
  - Embryonic expression in brain and neural tube
- *EFR3A*
  - Whole exome sequencing studies have implicated mutations in this gene with autism spectrum disorders
  - Embryonic expression in brain
- *MBNL1*
  - The encoded protein is a C3H-type zinc finger protein that modulates alternative splicing of pre-mRNAs
  - Mice lacking this gene exhibited muscle abnormalities and cataracts
  - Behavior/neurological phenotype
- *RMND5B*
  - Required For Meiotic Nuclear Division 5 Homolog B
  - RMND5 from Xenopus laevis is an E3 ubiquitin-ligase and functions in early embryonic forebrain development. (PMID: 25793641)
- *BUB1B*
  - kinase involved in spindle checkpoint function
  - Embryonic expression in brain and neural tube
  - Behavior/neurological phenotype
- HIRA
  - histone chaperone that preferentially places the variant histone H3.3 in nucleosomes
- EIF4E
  - Eukaryotic Translation Initiation Factor 4E subunit
  - behavioral fear response (GO:0001662)
- VPS25
  - subunit of the endosomal sorting complex required for transport II (ESCRT-II)
- *MMAA*
  - translocation of cobalamin into the mitochondrion, where it is used in the final steps of adenosylcobalamin synthesis
  - Defects in this gene are a cause of methylmalonic aciduria
  - Behavior/neurological phenotype
- *SFMBT1*
  - It encodes a protein which contains four malignant brain tumor repeat (mbt) domains and may be involved in antigen recognition
- *AEBP2*
  - Among its related pathways are Chromatin Regulation / Acetylation and Chromatin organization
  - Behavior/neurological phenotype
  - Embryonic expression in neural tube
- *POLR2E*
  - fifth largest subunit of RNA polymerase II
  - Embryonic expression in brain and neural tube and eye
- *SLC25A24*
  - This gene encodes a carrier protein that transports ATP-Mg exchanging it for phosphate
  - Mediates the reversible, electroneutral exchange of Mg-ATP or Mg-ADP against phosphate ions, catalyzing the net uptake or efflux of adenine nucleotides across the mitochondrial inner membrane
  - promoting the formation of calcium-phosphate precipitates in the mitochondrial matrix, and thereby buffering calcium levels in the mitochondrial matrix
  - Embryonic expression in brain
- *RAB5B*
  - Member RAS Oncogene
  - GTP binding and GDP binding
  - Embryonic expression in brain
  - Behavior/neurological phenotype
