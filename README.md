# PureSU3
This is the final project of the computational physics course held in Bonn and
was handed in instead of an exam.

The aim was to develop software capable of doing lattice QCD calculations.
Full calculations are at the frontier of physics research at supercomputers
all around the world. The simplification needed to fit the scope of this project
was to exclude fermionic degrees of freedom and only treat the gauged SU3  
behavior.

Results obtained here still bring fascinating insights into the unperturbative
nature of lattice QCD. For example the mass of the lightest glueball can
even with this simplification be determined. Even though the physical
glueballs recieve corrections due to fermions, its still a strong hint at
dimensional transmutation, i.e. the appearance of a mass scale even though the
theory is manifestly massless.

Some parts were calculated in parallel, while others remain serial to estimate
runtime.
