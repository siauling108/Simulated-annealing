# Simulated annealing
## The travelling salesman problem
This project consists in solving the travelling salesman problem **(TSP)** by following the simulated annealing procedure. Different moves and move-acceptance criteria can be chosen, which will be listed later on.

The travelling salesman problem is an optimization problem which tries to answer the question *"Given a list of cities and the distances between each pair of cities, what is the shortest possible route that visits each city and returns to the origin city?"*. This problem presents many variants and applications, making it very well-known in the optimization field.\
Simulated annealing **(SA)** is a probabilistic technique for approximating the global optimum of a given function, that in our case is the total lenght of the travel, which we want to minimize. In SA we introduce a "temperature" parameter T by which we imagine to cool down the system until a certain temeprature is reached. For each temperature we manipulate the visiting order of the cities by making a given number of moves; each move can be accepted or rejected depending on the fact that it satisfies a certain criterion or not. This simulation specifically consists in solving the TSP my means of different moves (swap, block-reverse) and different move-acceptance criteria (distance and Metropolis).\
Different moves, temperature profiles, and optimizable paths can be implemented in this simulation, as explained below.

### List of possible moves
Consider a set of **N** cities and a given closed path that connects them in a certain way. We can choose two different moves:
* *Swap*: two random cities get swapped in position;
* *Block-reverse*: invert the travel direction of the cities between two random cities of the set;

### List of move-acceptance criteria
* *Distance*: given an initial length *S<sub>old* , we make a move and calculate a new distance *S<sub>new*; the move will be accepted if <img src="https://render.githubusercontent.com/render/math?math=e^{-(S_{new}-S_{old})/T} \geq 1">, which automatically implies *S<sub>new* ≤ *S<sub>old*.
* *Metropolis*: <img src="https://render.githubusercontent.com/render/math?math=e^{-(S_{new}-S_{old})/T} > rand()">, where *rand()* is a random number between 0 and 1.

## Structure of the project
### Contents
Let's describe the content of each one of the main files of this project.
* The [inputs](https://github.com/MarcoCrr/Simulated-annealing/blob/master/inputs.txt) file contains the main parameters of the simulation and the folder path where data and plots are saved. List of the parameters:
  * *N*: number of cities to consider;
  * *iterations*: total number of annealing procedures;
  * *T_min*: minimum themperature reached;
  * *T*: starting temperature;
  * *alpha*: temperature scaling parameter;
  * *nsteps*: number of moves for each temperature;
  * *method*: method selected for the simulation;
  * *T_profile*: name of the chosen temperature profile;
  * *travel*: name of the selected travel.
* The [functions](https://github.com/MarcoCrr/Simulated-annealing/blob/master/functions.py) file contains all the main functions of the simulation. These functions will be imported in the [main](https://github.com/MarcoCrr/Simulated-annealing/blob/master/main.py) file;
* The [main](https://github.com/MarcoCrr/Simulated-annealing/blob/master/main.py) file is the main part of the simulation. It imports the functions from the [functions](https://github.com/MarcoCrr/Simulated-annealing/blob/master/functions.py) file and reads the [inputs](https://github.com/MarcoCrr/Simulated-annealing/blob/master/inputs.txt) file. The simulation will be executed according to the parameters in the inputs file. The data are saved into different files, which will be imported in the [plot](https://github.com/MarcoCrr/Simulated-annealing/blob/master/plot.py) file; the data being saved (in the corresponding [data]( ... ) folder) are:
  * *distances*: list containing the total distance of the path for each iteration;
  * *path*: array containing the path;
  * *Tem*: list containing the temperature profile;
  * *tot_acceptance*: array containing the move acceptance rate at each temperatue and iteration.
* The [plot](https://github.com/MarcoCrr/Simulated-annealing/blob/master/plot.py) file contains the functions for the plots. The data in the [data]( ... ) folder corresponding to that simulation are imported, and different plots can be obtained: the total travel length as a function of the iteration, the followed path, the move acceptance rate as a function of the temperature for each iteration, and the temperature profile. The results will be saved in the corresponding [plots]( ... ) folder;
* Lastly, the [tests](https://github.com/MarcoCrr/Simulated-annealing/blob/master/tests.py) file contains the tests for the functions in the [functions](https://github.com/MarcoCrr/Simulated-annealing/blob/master/functions.py) file.

### Simulation usage
1. The simulation parameters and the folder paths for the data and plots have to be set in the [inputs](https://github.com/MarcoCrr/Simulated-annealing/blob/master/inputs.txt) file. 
2. From the terminal, we launch the simulation by typing ***python main.py***;
3. If the user wants to get the plots, he has to type ***python plot.py*** from the terminal: the plots will be automatically saved in the specified folder and will also be shown in the terminal.

### Adding new options to the simulation
If needed, the user can implement new options (e.g. new methods, temperature profiles, or paths) for this simulation. The following lines suggest the main steps to do so:
* In the *functions* file, a proper function that describes the new option must be added;
* associate a name to that option and specify it in the *inputs* file in the "parameters";
* import the name of the parameter(s) in the *main* file, and add an option to select them (e.g. with an *if* statement).
