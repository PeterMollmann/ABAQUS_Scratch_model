from odbAccess import *
import numpy as np
import os


def PostProcess(jobName, fileName):
    odbName = jobName + ".odb"
    odb = openOdb(path=odbName, readOnly=True)
    outputFileName = fileName + "_RFs.txt"
    outputFilePath = "SimDataOutputs/" + outputFileName
    if not os.path.exists(outputFilePath):
        os.makedirs(outputFilePath)
    steps = odb.steps.keys()
    rfs = ["RF1", "RF2", "RF3"]
    with open(outputFilePath, "w") as file:
        file.write("Step,Time,RF1,RF2,RF3\n")
        for step in steps:
            historyRegion = odb.steps[step].historyRegions.values()[-1]
            time = []
            rf1Data = []
            rf2Data = []
            rf3Data = []
            for rf in rfs:
                historyOutput = historyRegion.historyOutputs[rf].data
                if rf == "RF1":
                    for dataPoint in historyOutput:
                        time.append(dataPoint[0])
                        rf1Data.append(dataPoint[1])
                elif rf == "RF2":
                    for dataPoint in historyOutput:
                        rf2Data.append(dataPoint[1])
                elif rf == "RF3":
                    for dataPoint in historyOutput:
                        rf3Data.append(dataPoint[1])
                with open(outputFilePath, "a") as file:
                    for t, rf1, rf2, rf3 in zip(time, rf1Data, rf2Data, rf3Data):
                        file.write(
                            "%s,%5.5f,%5.5f,%5.5f,%5.5f\n" % (step, t, rf1, rf2, rf3)
                        )
    # with open(outputFilePath, "w") as file:
    #     file.write("Step,Time,RF1,RF2,RF3\n")
    #     for step in odb.steps.keys():
    #         historyRegion = odb.steps[step].historyRegions.values()[-1]
    #         time_data = []
    #         rf_data = {rf: [] for rf in historyRegion.historyOutputs.keys()}
    #         for rf in historyRegion.historyOutputs[-1].keys():
    #             historyOutput = historyRegion.historyOutputs[rf].data
    #             if rf == "RF1":
    #                 time_data = [dataPoint[0] for dataPoint in historyOutput]
    #             rf_data[rf] = [dataPoint[1] for dataPoint in historyOutput]
    #         for t, rf1, rf2, rf3 in zip(
    #             time_data, rf_data["RF1"], rf_data["RF2"], rf_data["RF3"]
    #         ):
    #             file.write(
    #                 "{},{:.5f},{:.5f},{:.5f},{:.5f}\n".format(step, t, rf1, rf2, rf3)
    #             )
    # Undeformed coordinates of contact surface
    contactSurfaceSet = "contactSurfaceSet"
    contactNodes = odb.rootAssembly.nodeSets[contactSurfaceSet.upper()].nodes[0]
    # Store undeformed coordinates in a list [(nodeLabel, x, y, z)]
    undeformedCoords = [
        (node.label, node.coordinates[0], node.coordinates[1], node.coordinates[2])
        for node in contactNodes
    ]
    # Sort nodes by x, then y, then z (before deformation)
    sortedUndeformedCoords = sorted(
        undeformedCoords, key=lambda coord: (coord[3], -coord[1], coord[2])
    )
    # Extract displacement field from the last frame
    lastStep = odb.steps.values()[-1]
    lastFrame = lastStep.frames[-1]
    displacementField = lastFrame.fieldOutputs["U"]
    # Get unique displacement values for each node
    dispSubset = displacementField.getSubset(
        region=odb.rootAssembly.nodeSets[contactSurfaceSet.upper()]
    )
    displacements = {
        value.nodeLabel: np.array(value.data) for value in dispSubset.values
    }
    # Compute deformed coordinates in the **sorted order**
    deformedCoords = []
    for nodeLabel, x, y, z in sortedUndeformedCoords:
        disp = displacements.get(
            nodeLabel, np.array([0.0, 0.0, 0.0])
        )  # Default to zero displacement if not found
        deformedCoords.append((nodeLabel, x + disp[0], y + disp[1], z + disp[2]))
    # Write sorted undeformed and deformed coordinates to file
    outputFileName = fileName + "_SurfCoords.txt"
    outputFilePath = "SimDataOutputs/" + outputFileName
    if not os.path.exists(outputFilePath):
        os.makedirs(outputFilePath)
    with open(outputFilePath, "w") as file:
        file.write(
            "NodeLabel, x_undeformed, y_undeformed, z_undeformed, x_deformed, y_deformed, z_deformed\n"
        )
        for (nodeLabel, x_und, y_und, z_und), (_, x_def, y_def, z_def) in zip(
            sortedUndeformedCoords, deformedCoords
        ):
            file.write(
                "%d, %5.5f, %5.5f, %5.5f, %5.5f, %5.5f, %5.5f\n"
                % (nodeLabel, x_und, y_und, z_und, x_def, y_def, z_def)
            )
    odb.close()
