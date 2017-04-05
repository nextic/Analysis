#!/Applications/anaconda3/bin/python

from __future__ import print_function
import sys
import subprocess
import tables as tb

def is_ok(filename):
    try:
        tb.open_file(filename).close()
        return True
    except:
        return False

def delete_file(filename):
    subprocess.call(["rm", "-f", filename])

for filename in sys.argv[1:]:
    if not is_ok(filename):
        delete_file(filename)
        print(filename, "has been deleted")
