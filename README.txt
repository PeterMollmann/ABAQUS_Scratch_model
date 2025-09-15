Repository to run scratch tests in Abaqus.

In Abaqus Command: "abaqus cae noGUI=SubmissionFile.py"

This runs the submission file with the abaqus kernel.

Currently only progressive load scratch tests are available.


Data generation:
As abaqus runs on an old version of python, the .csv file is made to a .py file containing a dictionary of the material parameters

Post processing:
The reaction forces on the indenter tip as well as the coordinates of the substrate surface in contact with the indenter are saved to seperate files.