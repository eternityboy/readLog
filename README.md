# readLog for OpenFOAM based on python3
add alias to your .bashrc
```bash
alias readLog="python3 /path/to/readLog.py"
```

## Version 1.5
* support for both OF branches
* command line execution, example:
```bash
$ python3 readLog.py openfoam.log 1 20 50
```
```
> 1 goes for Forces
> 2 goes for Moments
> 20 start point for reference range
> 50 end point for reference range
```
## To-Do
suggestions? QT? FFT analyze?

## Version 1.4
reference point changed to reference range
### To-Do
support for both OF versions

## Version 1.2
added roll motion analize
### To-Do
system argument - reference point 

## Version 1.1
added possibilty to specify range of values to calculate mean

added mean value to plot

re-aranged part of the code with numpy

SixDoF support

Volume Phase
