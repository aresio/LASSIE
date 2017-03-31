# About LASSIE's GUI

LASSIE's Graphical User Interface was specifically designed to simplify the execution of simulations. The GUI represents a user-friendly tool that does not require any GPU programming skill or expertise in modeling biological systems with ODEs.

To launch the gui run the `lassie-gui.py` python script.

# DEPENDENCIES

The following Python libraries are mandatory to use the GUI: PyQT4; numpy; matplotlib.
In addition, to use the beta version of the SBML import module, the python-libsml library is required.

# SUPPORTED FUNCTIONALITIES

At the present time, the GUI supports the following functionalities:
* it provides a visual way to load and simulate a reaction-based model based on mass-action kinetics;
* it provides a visual way to specify the total simulation time and the number of sampling time instants of the dynamics to be saved;
* it shows a summary of the main information about the loaded model;
* it plots the output dynamics of the selected species and allows to save the output dynamics to file. 

# FUNCTIONALITIES UNDER DEVELOPMENT

We plan to extend the GUI in the next future, by introducing some additional functionalities:
* support for PySB and BioNetGen rule-based models (they will be transparently converted into reaction-based models and simulated by LASSIE);
* complete support for _de novo_ model  editing;
* support for species, reactions and parameters editing;
* advanced analysis tools (e.g., sensitivity analysis, parameter estimation);
* select and show the units of measure of, e.g., time, concentration, kinetic parameters.

# SBML support 

The possibility to import and simulate models given in the SBML standard language is currently under development.
In particular, we highlight that a well-formatted SBML file can be imported and simulated by LASSIE only if the model
is entirely based on mass-action kinetics. When a SBML file is imported, the GUI automatically creates a new directory
with the input files corresponding to the species and reactions described in the SBML file. To import a SBML file click File
> Import SBML or press CTRL+I. Note that the python-libsml library is required to use this beta version of the SBML
import module.
