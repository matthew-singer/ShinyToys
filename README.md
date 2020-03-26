# ShinyToys
Playground for making Shiny apps for interactive data visualizations

At the moment, these apps will be related to classic models of population dynamics in ecology and evolutionary biology. They are meant to be used for educational purposes as an alternative to the program Populus.


## sirapp:

An SIR model of epidemiological spread of infectious disease.
Parameters to alter are:

* beta - transmissibility (probability of interaction * probability of transfer of disease)
* gamma - recovery rate (probability of being removed from the infectious population, either by death or recovery)
* N - Total popuation size
* I - Number of initial infections at Time=0
* Days - number of days to run the simulation
