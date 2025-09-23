import sys
import os

sys.path.append("../common")
from atomx_simulation_script_writter import atomx_simulation_script_writter

if __name__ == "__main__":
    #print(atomx_simulation_script_writter.__doc__)

    a = atomx_simulation_script_writter("ATOMX_34Ar")
    a.set_path("/mnt/nptool/atomx/data", "input", "/mnt/nptool/atomx/figures")
    a.make_project_configuration()
    a.make_geant4_vis_file()

    detector = a.make_detector_file(reaction_z=[150,160], step_limit=10)
    reaction = a.make_eventgen_file_beam(beam_particle="34Ar", energy=81.6)
    a.add_eventgen_elastic(reaction, beam_particle="34Ar", target="4He", heavy="37K", light="1H")

    a.make_script(detector.name, reaction.name)
