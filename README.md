# SketPH

SketPH is a smoothed particle hydrodynamics code with a purpose of sketching


# run code 
add to the PYTHONPATH export PYTHONPATH=/path/to/the/folder/which/contains/sketph/folder/:${PYTHONPATH}

python run-AMW.py


The purpose of this program is to be used to verify the simulation techniques. 

The structure of the program assumes that it is loaded as a set of modules 
whose path, as it does in python, added to the local environment,
 and writing a particular script uses these modules as blocks.


The principle of constructing a program is somewhat different from what
 is usually in such computational programs now.
 The program is used to prototyping methods of mathematical modeling,
 so its quality is to get a minimum good program in a relatively
 short period of time (student course), with the ability to quickly
 modify it to verify and understand the properties of the particular
 model or numerical method. In contrast to the conventional approach,
 where modules assume independent operation. 


Here, modules are interconnetcted somehow in sence of global "field" for good reason.

storage - a field with a dictionary whose fields are numpy-arrays of a given size with
elements of the given type from the fields module.

field - this class introduces variables with the names of the data fields containing values,
 which are strings identical to the names of their variables.
 Separately list of fields is stored in fields_list for iterations over fields.
 The variables fields are written in the dictionary {variable : its data type}. 

neibs - contains the hash-map to make a neighbours lists 
solveSpace - a right-hand size of evolution equations to compute gradients 
solveTime - makes a time stepping model - implements a material equation of state 
run-AMW.py - run adaptive moving window simulation to obtain stationary shock waves using Hugoniot eos, perfectily plastic copper flow and feedback algorithm of a coordinate system velocity adaptation. 
