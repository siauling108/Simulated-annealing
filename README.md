# Simulated annealing
## The travelling salesman problem
This project consists in solving the travelling salesman problem **(TSP)** by following the simulated annealing procedure. Different moves and move-acceptance criteria can be chosen, which will be listed later on.

The travelling salesman problem is an optimization problem which tries to answer the question *"Given a list of cities and the distances between each pair of cities, what is the shortest possible route that visits each city and returns to the origin city?"*. This problem presents many variants and applications, making it very well-known in the optimization field.\
Simulated annealing **(SA)** is a probabilistic technique for approximating the global optimum of a given function, that in our case is the total lenght of the travel. In SA we introduce a "temperature" parameter T by which we imagine to cool down the system until a certain temeprature is reached. For each temperature we manipulate the visiting order of the cities by making a given number of moves; each move can be accepted or rejected depending on the fact that it satisfies a certain criterion or not. This project specifically consists in solving the TSP my means of different moves (swap, block-reverse, prune and graft) and different move-acceptance criteria (distance and Metropolis).\
Different moves and different acceptance criteria can be easily implemented in this simulation.

### List of possible moves
Consider a set of **N** cities and a given closed path that connects them in a certain way. We can consider three different moves:
* *Swap*: two random cities get swapped in position;
* *Block-reverse*: invert the travel direction of the cities between two random cities of the set;
* *Prune and graft*: take ("prune") a travel segment between two random sites, open a breach inside the cities list and graft it in there.

### List of move-acceptance criteria
* *Distance*: given an initial length *S<sub>old* , we make a move and calculate a new distance *S<sub>new*; the move will be accepted if <img src="https://render.githubusercontent.com/render/math?math=e^{-(S_{new}-S_{old})/T} \geq 1">, which automatically implies *S<sub>new* ≤ *S<sub>old*.
* *Metropolis*: <img src="https://render.githubusercontent.com/render/math?math=e^{-(S_{new}-S_{old})/T} > rand()">, where *rand()* is a random number between 0 and 1.

## Structure of the project
### Contents
Let's describe the content of each one of the main files of this project.
* The [inputs](https://github.com/MarcoCrr/Simulated-annealing/blob/master/inputs.txt) file contains the main parameters of the simulation (such as the number of the cities to consider, the number of iterations...), and the folder path where data and plots are saved;
* The [functions](https://github.com/MarcoCrr/Simulated-annealing/blob/master/functions.py) file contains all the main functions of the simulation: the function that generates the cities, the one that calculates the path, the one that calculates its length, all the possible moves and annealing procedures, other minor functions. These functions will be imported in the [main](https://github.com/MarcoCrr/Simulated-annealing/blob/master/main.py) file;
* The [main](https://github.com/MarcoCrr/Simulated-annealing/blob/master/main.py) file is the main part of the simulation. It imports the functions from the [functions](https://github.com/MarcoCrr/Simulated-annealing/blob/master/functions.py) file and reads the [inputs](https://github.com/MarcoCrr/Simulated-annealing/blob/master/inputs.txt) file. Then it asks the user what combination of move+acceptance criterion wants to use, and executes the calcultions accordingly. The simulation data are saved into different files, which will be imported in the [plot](https://github.com/MarcoCrr/Simulated-annealing/blob/master/plot.py) file;
* The [plot](https://github.com/MarcoCrr/Simulated-annealing/blob/master/plot.py) file contains the functions for the plots. The data obtained from the [main](https://github.com/MarcoCrr/Simulated-annealing/blob/master/main.py) file are imported, and different plots can be obtained: the total length of the travel as a function of the iteration, the followed path, and the move acceptance rate as a function of the temperature for each iteration.
* Lastly, the [tests]() file contains the tests for the functions in the [functions](https://github.com/MarcoCrr/Simulated-annealing/blob/master/functions.py) file.

### Simulation usage
1. The simulation parameters and the folder paths for the data and plots have to be set in the [inputs](https://github.com/MarcoCrr/Simulated-annealing/blob/master/inputs.txt) file. More specifically, the parameters are: the number of cities to consider, the number of iterations /annealing procedures to be executed in the simulation, the minimum and maximum temperatures, the temperature scaling parameter, the number of moves to be done for each temperature;
2. From the terminal, we launch the simulation by typing ***python main.py***: once the user chooses what combination of move+acceptance criterion wants to use (as explained in the terminal), the calculations will start;
3. If the user wants to get the plots, he has to type ***python plot.py*** from the terminal: the plots will be automatically saved in the specified folder and will also be shown in the terminal.
