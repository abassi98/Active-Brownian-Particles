# Active-Brownian-Particles
Final project for the course "Simulation in Statistical Mechanics" at Innsbruck University

The source code is written in C++, which is then binded in python though "Pybind11" <br>
To use the code download the repo with <br>
`git clone https://github.com/lupoalberto98/Active-Brownian-Particles.git` <br>
Compile the module with `make`<br>
Run the notebook to do the simulations (using the binded version) <br>
The makefile is supposed to work as it is on Mac and with python 3.9. If the working OS is Linux, then substitute the flag `-undefined dynamic_lookup` with `-fPIC` and `python3.9` with the current release that has been used.