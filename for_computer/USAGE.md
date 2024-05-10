These files will run on the lab machines.

To run the files, type one of the following into the command line while in the correct directory:
  ``julia mcmullen.jl [int::ceiling]``
  ``julia circle_counting.jl [int::ceiling]``
  ``julia bai_finch.jl [int::k0]``

For example, ``julia circle_counting.jl 50000`` runs the circle counting program at a maximum curvature of 50000.

For the driver file, run ``chmod +x driver.sh`` (replacing the file name if changed) to compile, then type ``./driver.sh`` to run the file.
An arbitrary driver file will take the form:

```#!/bin/bash```
```for ((i = [int::start] ; i < [int::end] ; i+=[int::increase] )); do julia [filename] $i; done```

The first line is necessary for it to work. 

Use ``./driver.sh > [file].tx`` to write the output to a text file.
