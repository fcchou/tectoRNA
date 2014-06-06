####################################
Simulation of tecto-RNA with HelixMC
####################################

This piece of code uses the HelixMC package to simulate the stability of
tecto-RNA.

Setup
=====
You will need to first install python (v2.7) and following packages:

- scikit-learn (http://scikit-learn.org/stable/)

- HelixMC (https://github.com/fcchou/HelixMC)

HelixMC/scikit-learn also depends on numpy and scipy.

How to use
==========
The main simulation code is ``simu_tecto.py``. To get a list to brief description
of avaliable options, simply run ``python simu_tecto.py -h``.

Here is a simple command line example::

    python simu_tecto.py -seq1 GGGGGGGGGGGG CCCCCCCCCCCC -n_cycles 1000

After the run you will get output that looks like::

    Log likelyhood score: -11.8090588021
    Free energy: 6.96734469322

Also two files ``helix1_bp.par`` and ``helix2_bp.par`` will be generated,
containing the parameters for the best helical conformer for forming tecto-RNA
in the simulation. These files can be then used to generate pdb files.

Generating PDB files
====================
The pdb files can be generated using Web3DNA (http://w3dna.rutgers.edu/).

In the Web3DNA homepage, click *Reconstruction* on the top stripe. Then select
**Customized model**. In the dropdown menu, choose
*Customized base-pair-step/nucleotide parameter*, then click **Continue**.

Now in the **Backbone conformation** dropdown menu, select *RNA-form*. For the
**Base-pair parameter file**, upload one of the previously generated file
(e.g. ``helix1_bp.par``). Then click **Build**. Downloas the generated
pdb file, and repeat the process for the other helix.
