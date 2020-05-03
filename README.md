)# SPRKKR_TOOLS
Collection of useful and not-so-useful scripts and tools for the SPR-KKR DFT package

The goal of this repository is to collect useful script, tools and source code.

Please feel free to submitt pull request or send me (Jan Minar) your contribution but please follow these guidelines:

-- Place your tool into given directory (kkrscf, kkrgen, kkrspec, kkrchi, auxilary) depending on the functionality of your script
   --AUX
   --CHI: scripts connected with linear responce calculations
   --GEN: scripts connected with all kkrgen tasks (e.g. JXC, DOS etc..)
   --SCF: scripts connected with self consistency
   --SPEC: scripts for kkrspec ARPES, AIPES, SPLEED etc.
   --VIS: Visualisation tools

-- Please create directories (with a clear unique name) which ALWAYS contain these descriptive text files:
   -- README: detailed description of your tool 
   -- REQUIREMENTS: detailed list of dependencies: libraries, compilers etc.
   -- VERSION: version of your script, and please list also version of SPRKKR package for which it was tested.
   
   
Have a fun,
Jan Minar
