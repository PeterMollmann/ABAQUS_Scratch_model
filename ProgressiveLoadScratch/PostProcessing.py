from odbAccess import *
import numpy as np
import os
from . import Constants as C


def PostProcess(jobName, fileName):
    odbName = jobName + ".odb"
    odb = openOdb(path=odbName, readOnly=True)

    #### ---------------- ####
    #    Get displacements
    #### ---------------- ####
    # contactRegionName = "s_Surf-1"
    # surface_names = odb.rootAssembly.surfaces.keys()
    # surface_region_name = next((s for s in surface_names if "S_" in s), None)

    allNodesInContactRegion = odb.rootAssembly.surfaces[
        C.slave_surface_name.upper()
    ].nodes[0]
    unique_nodes = {node.label: node for node in allNodesInContactRegion}.values()

    undeformedCoords = [
        (node.label, node.coordinates[0], node.coordinates[1], node.coordinates[2])
        for node in unique_nodes
    ]

    sortedUndeformedCoords = sorted(
        undeformedCoords, key=lambda coord: (coord[3], -coord[1], coord[2])
    )

    displacementField = odb.steps.values()[-1].frames[-1].fieldOutputs["U"]

    dispSubset = displacementField.getSubset(
        region=odb.rootAssembly.nodeSets[C.contact_region_nodes_name.upper()]
    )
    displacements = {
        value.nodeLabel: np.array(value.data) for value in dispSubset.values
    }

    deformedCoords = []
    for nodeLabel, x, y, z in sortedUndeformedCoords:
        disp = displacements.get(nodeLabel, np.array([0.0, 0.0, 0.0]))
        deformedCoords.append((nodeLabel, x + disp[0], y + disp[1], z + disp[2]))

    outputFileName = fileName + "_SurfCoords.txt"
    outputFilePath = "SimDataOutputs/" + outputFileName
    if not os.path.exists("SimDataOutputs/"):
        os.makedirs("SimDataOutputs/")
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

    #### ---------------- ####
    # Get forces and energies
    #### ---------------- ####
    outputFileName = fileName + "_RFsAndEnergies.txt"
    outputFilePath = "SimDataOutputs/" + outputFileName

    steps = odb.steps.keys()
    history_regions = odb.steps[steps[-1]].historyRegions.keys()
    history_indenter_region_name = next(
        (s for s in history_regions if "Node" in s), None
    )
    history_substrate_region_name = next(
        (s for s in history_regions if "Element" in s), None
    )

    def get_history_data(step_name, region_name):
        history_region = odb.steps[step_name].historyRegions[region_name]
        data = []
        time = []
        for key in history_region.historyOutputs.keys():
            output_data = np.array(history_region.historyOutputs[key].data).T
            if key == "RF1":
                time.append(output_data[0, :])
            data.append(output_data[1, :])
        return np.array(time), np.array(data)

    time1, [rf1a, rf2a, rf3a] = get_history_data(steps[0], history_indenter_region_name)
    time2, [rf1b, rf2b, rf3b] = get_history_data(steps[1], history_indenter_region_name)
    time2 = time2 + np.max(time1)

    time_all = np.append(time1, time2)
    rf1_all = np.append(rf1a, rf1b)
    rf2_all = np.append(rf2a, rf2b)
    rf3_all = np.append(rf3a, rf3b)

    _, [IEa, KEa] = get_history_data(steps[0], history_substrate_region_name)
    _, [IEb, KEb] = get_history_data(steps[1], history_substrate_region_name)

    IE_all = np.append(IEa, IEb)
    KE_all = np.append(KEa, KEb)

    # write forces and energies to file
    with open(outputFilePath, "w") as file:
        file.write("Step,Time,RF1,RF2,RF3,IE,KE \n")
        with open(outputFilePath, "a") as file:
            for t, rf1, rf2, rf3, IE, KE in zip(
                time_all, rf1_all, rf2_all, rf3_all, IE_all, KE_all
            ):
                file.write(
                    "%5.5f,%5.5f,%5.5f,%5.5f,%5.5f,%5.5f\n" % (t, rf1, rf2, rf3, IE, KE)
                )

    odb.close()
