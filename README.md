# SFND_Radar_Detection


Implementation steps for the 2D CFAR process.
1 - Training cells and guard cells are selected,
2 - Offset is set to a small value
3 - The number of training cells is calculated
4 - Nested loops are implemented to slide the Trainig/CFAR window across the array of RDM
5 - In a new array of zeros of same dimesions of RDM, if the CUT cell value is above threshold, the corresponding cell value is set to 1, otherwise it is set to zero

Selection of Training, Guard cells and offset:
1 - Training cells were set to 5 and 5 as a starting value before tunning 
2 - Guarding cells were selected to be 2 and 2 as a starting value
3 - offset is set to 0
4 - offset is then increased to 10 -> results are suitable


Steps taken to suppress the non-thresholded cells at the edges.
The result is already inserted into an empty array of zeros so there was no need for suppression.
