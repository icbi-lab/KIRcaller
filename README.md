# KIRcaller

KIRcaller is Python tool that implements our KIR3DS1/L1 genotyping method (publication under review).
It works with BAM files that were created by aligning ExomeSeq or WGS reads to genome version hg19. If a different genome assembly was used for read alignment, it is mandatory to provide a custom GFF file that describes the relevant KIR3DS1/L1 regions (see KIR3DS1_L1.gff GFF file in the "data/" subfolder) for the corresponding genome assembly.



    Requirements:
    
    - Linux
    - subread (http://subread.sourceforge.net/)
    - Python 2.7
    - Pandas (python module)
    
    Installation:
    
    callKIR3D.py and data/ including its contents need to be located in the same directory:
    
    e.g. /home/$USER/KIRCaller
    KIRcaller/
    ├── callKIR3D.py
    ├── data
        └── KIR3DS1_L1.gff
    
    featureCounts form the subread package should be located in your $PATH,
    alternatively it can be specified via the -f or --featureCounts option.
    
    Running the tool:
    
    A simple run would be:
    ./callKIR3D.py -b /data/KIR/BAM -o KIR_geontyping.tsv
    
    (tested in on CentOS Linux 7.5 using python 2.7.14, subread 1.6.0 and 1.6.1)
    
