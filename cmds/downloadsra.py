import subprocess
import pathlib
from cmds.runtable import RunTable

acc_num = 'SRR6664514'
subprocess.run(["prefetch", acc_num])
