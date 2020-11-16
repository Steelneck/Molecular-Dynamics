import Calculations.MSD
import os
import sys

if os.path.exists("atoms.traj"):
    os.unlink("atoms.traj")
    MSD.run_md()
else:
    sys.exit(1)
