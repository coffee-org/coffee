# Design steps: Polychromatic design {#page_PIAACMCsimul_DesignSteps_Poly}

## Compute complex amplitude response of each focal plane mask zone

~~~
./run optsingle 200
~~~

This command will run for a few hrs if the number of focal plane mask zones is large.

The complex amplitude response of each mask zone is computed. Several sub-processes are launched, each computing a subset of the total number of zones.

You can follow the progress of each subprocess in the corresponding tmux sessions: the tmux session names are PID followed by FPMt<index>n<NBindex>

As each subprocess computes zone responses, they are stored in FITS files FPMresp..._thread<index>.fits.tmp you can open/view these files to follow progress - each zone appears as a line, and files should fill from the bottom to the top when all threads complete, the files are merged into a single FPMresp file.

## Search for optimal solution

Once the FPMresp... file is created, it is used by the search algorithm to quickly evaluate focal plane mask designs.

The search algorithm proceeds as two nested loops (outer loop steps 1-4, inner loop steps 2-4):

1. OUTER LOOP START ITERATION: starting point is chosen, either randomly or close to the best solution. 
2. INNER LOOP START ITERATION: Derivatives around the current point.
3. For each predefined regularization coefficient :
	a. An SVD-based inversion used to define a search direction
	b. The solution is evaluated along the search direction, and the best point along this line is recorded
4. The best solution encountered in the above step is recorded. If it is a significant improvement over the last step 3 iteration, goto step 2, if not, go to step 1. If the inner loop (steps 2-3) has been running for more than 50 iterations, go to step 1.


In the C code, this optimization corresponds to execution mode 13. To stop the loop, create file "stoploop13.txt" in the working optimization directory:

	touch piaacmcconf_i000/stoploop13.txt

Note that this command will stop the loop on entering step 1, so the inner loop (steps 2-4) will still keep running until the solution stops improving.




