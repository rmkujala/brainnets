#! /bin/python
"""
Creates the necessary file for submitting a job to slurm and submits it.

Usage::

    python slurm_submit.py myscript.py [batch]

Default behavior submits to queue ``play`` with settings::

    #SBATCH --time=0-00:15:00 --mem-per-cpu=5000
    #SBATCH -p play

if "batch" is the last parameter these settings are used::

    #SBATCH --time=5-00:00:00 --mem-per-cpu=5000
    #SBATCH -p batch

(The filename to be submitted is myscript.sh)
"""

import sys
import subprocess


def main():
    # script to be executed
    pythonscript = sys.argv[1]
    if pythonscript[-3:] != ".py":
        print "trying to submit a non-python script, exiting.."
        print "make changes if you want this to be possible"
        exit()
    basename = pythonscript[:-3]  # corresponds to "myscript"

    try:
        f = open(pythonscript, "r")
        f.close()
    except:
        print "something wrong with the filename perhaps, exiting..."

    # default goes to play
    play = True
    try:
        batchstring = sys.argv[2]
        if batchstring == "batch":
            play = False
            print "submitting to batch..."
    except:
        print "submitting to play..."

    shfname = basename + ".sh"
    outfname = basename + ".out"

    stringtowrite = ""
    stringtowrite += "#!/bin/bash\n"
    if play:
        stringtowrite += "#SBATCH --time=0-00:15:00 --mem-per-cpu=5000 \n"
        stringtowrite += "#SBATCH -p play \n"
    else:
        stringtowrite += "#SBATCH --time=5-00:00:00 --mem-per-cpu=5000 \n"
        stringtowrite += "#SBATCH -p batch \n"
    stringtowrite += "#SBATCH -o " + outfname + "\n"
    stringtowrite += "#SBATCH -n 12 \n"
    stringtowrite += "#SBATCH -N1-1 \n"
    stringtowrite += "\n"
    stringtowrite += \
        "module load igraph/0.6 python-igraph/0.6 scipy/0.11.0 matlab\n"
    stringtowrite += "\n"
    stringtowrite += "python " + pythonscript

    f = open(shfname, "w")
    print shfname
    f.write(stringtowrite)
    f.close()  # flushes the conent
    subprocess.call(["sbatch", shfname])

if __name__ == "__main__":
    main()
