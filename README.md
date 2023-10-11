# SketPH

SketPH is a smoothed particle hydrodynamics code with a purpose of sketching


# run 
add to the PYTHONPATH by

```
export PYTHONPATH=/path/to/the/folder/which/contains/sketph/folder/:${PYTHONPATH}
```

run script
```
python run-*.py
```

# description
The purpose of this program is to be used to verify the simulation techniques. 

The structure of the program assumes that it is loaded as a set of modules 
whose path, as it does in python, added to the local environment,
 and writing a particular script uses these modules as blocks.


The principle of constructing a program is somewhat different from what
 is usually in such computational programs now.
 The program is used to prototyping methods of mathematical modeling,
 so its quality is to get a minimum good program, with a flexible 
 modifocation opportunities for simple verification and 
 understanding the properties of particular
 models or numerical methods. This approach differs from the conventional approach,
 where modules assume independent operation and general api conquers its architecture. 


All modules are interconnetcted somehow in sence of global "field"s for good reason.

"storage" is a dictionary whose fields are numpy-arrays of a given size with
elements of the given type from the fields module.

"field" is a class introducing variables with names of a data fields which contain the data.
 These are string variables which values are identical to the names of their variables.
 Separately, there is a list of fields is stored in fields_list for iterations over fields.
 The variables fields are written in the dictionary {variable : its data type}. 

"neibs" contains the hash-map to build-up a neighbours lists 

"solveSpace" is a right-hand side of evolution equations to compute gradients 

"solveTime" makes a time stepping model and implements some of material closure relations 

"run-AMW.py"/"run-MW" runs adaptive moving window AMW/ moving window (MW) simulation to obtain stationary shock wave

"models" contains Hugoniot eos, linear elastic perfectily plastic copper flow.
"stepperAMW" makes a time by using inflow and outflow boundary conditions and employs the ordinary SolveTime stepper
"AMW"/"MW" are classes containing the feedback algorithm of a coordinate system velocity adaptation in AMW and outflow velocity correction in MW.