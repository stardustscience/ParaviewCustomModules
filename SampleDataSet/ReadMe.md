## SampleDataSet
SampleDataSet is a Python custom module that randomly samples datasets (VTI or tetrahedral VTU) based on a scalar variable.    It uses the scalar vaariable in two ways: one, to distribute the sample candidates, and second, to determine a crowding factor.

For the first step, it assigns each cell in the dataset a probability based on the cell's average vertex (scaled from 0 to 1) value times the volume of the cell.  These probabilities are rank ordered and cells to sample are chosen using an interval search on a running sum of the probabilities.

Once a cell is chosen, a random candidate point is chosen inside the cell.   The scaled scalar value at that point is interpolated, and a local crowding radius is determined by a linear table defined by the 'spacing at max probability' and 'spacing at min probability' properties.   If no previous sample ppint is within that radius of teh candidate point, the candidate point is accepted.   This continues until the specified number of samples is reached or 1000 successive candidates have failed acceptance.

## Example
For a good time, start with the Wavelet source.   Hook this to a Programmable Filter, set the CopyArrays checkbox, and use the following:

    import numpy as np
    rtdata = output.PointData['RTData']
    midpoint = (np.max(rtdata) + np.min(rtdata)) / 2
    prob = -np.abs(rtdata - midpoint)
    prob = (prob - np.min(prob)) / (np.max(prob) - np.min(prob))
    prob = np.power(prob, 10)
    output.PointData.append(prob, 'PDF')

Then pass this result to SampleDataSet, using 0.001 as the 'spacing at max probability', 0.01 as the 'spacing at min probability', ask for 10000 samples, and set the 'String' to 'PDF' (that selects the variable to operate on).   You'll get samples clustering around an 'isosurface' and the midpoint value.

## Note
This filter retains samples from execution to execution and begins each run using the retained samples as candidates to minimize popping in time varying examples.

## Installation
Run:

	cmake -DCMAKE_INSTALL_PREFIX=$HOME
	make install
	
The python file and the shared-library it uses will be installed at $HOME/ParaviewPythonModules and you can load it into your Paraview run from there.
