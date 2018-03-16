"""

Short demonstration of how to use P3CAN module of the tribology package.

"""


from tribology.p3can.p3can import p3can, generate_input_file
import os


def p3can_demo():
    """

    Generate a P3CAN input file, then start a simulation run. This function
    performs a number of simulation runs for different simulation types;
    see the P3CAN documentation for more information (also have a look at the
    P3CAN input file).

    A directory `results` will be created in the current working directory that
    contains the simulation output files. See terminal output for more
    information.

    """

    for sim_type in (3, 5, 6, 7, 8):
        out_file = generate_input_file(
            sim_type,
            os.getcwd() + os.sep + 'demo.txt')
        p3can(out_file)


if __name__ == "__main__":
    p3can_demo()
