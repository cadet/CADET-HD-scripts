#!/usr/bin/python3

import re
import os
import shutil
import subprocess
import argparse
from tempfile import NamedTemporaryFile
# from scipy.optimize import minimize
from random import randint, seed

LOG="optpack.log"

def sed_inplace(filename, pattern, repl):
    '''
    Perform the pure-Python equivalent of in-place `sed` substitution: e.g.,
    `sed -i -e 's/'${pattern}'/'${repl}' "${filename}"`.
    '''
    # For efficiency, precompile the passed regular expression.
    pattern_compiled = re.compile(pattern)

    # For portability, NamedTemporaryFile() defaults to mode "w+b" (i.e., binary
    # writing with updating). This is usually a good thing. In this case,
    # however, binary writing imposes non-trivial encoding constraints trivially
    # resolved by switching to text writing. Let's do that.
    with NamedTemporaryFile(mode='w', delete=False) as tmp_file:
        with open(filename) as src_file:
            for line in src_file:
                tmp_file.write(pattern_compiled.sub(repl, line))
    # Overwrite the original file with the munged temporary file in a
    # manner preserving file attributes (e.g., permissions).
    shutil.copystat(filename, tmp_file.name)
    shutil.move(tmp_file.name, filename)

# def reprint(string):
#     thread = current_process().name
#     print(thread, ": " , string)

def run_command_in_shell(cmdstring):
    # print(cmdstring)
    p = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
    stdout, stderr = p.communicate()
    # LOG = Path.cwd().joinpath(LOG_FILE)
    with open(LOG, 'a') as logfile:
        logfile.write(stdout)
        logfile.write(stderr)
    if stderr.strip() != "":
        print('---')
        print(stderr)
        print('---')
    else:
        return stdout

def genpack(rate, target_porosity):

    try:
        os.remove('packing.nfo')
    except OSError:
        pass

    # r = r[0]
    seed(a=None, version=2)
    myseed = randint(0,2147483647)

    shutil.copy('packing.fba.xyzd', 'packing.xyzd')
    sed_inplace('generation.conf', 'Contraction rate.*', 'Contraction rate: ' + str(rate) )
    sed_inplace('generation.conf', 'Seed.*', 'Seed: ' + str(myseed) )
    out = run_command_in_shell('PackingGeneration.exe -ls')
    porosity = 0
    # print(out)
    for line in out.split('\n'):
        if 'Calc. porosity' in line:
            # print(line.split(' '))
            porosity = float(line.strip().split()[-1])
    error = target_porosity - porosity
    print("r: {r}, seed: {seed},  por: {porosity}, err: {error}".format(r=rate, seed=myseed,porosity=porosity, error=error))
    return (porosity, error)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--porosity", type=float, required=True, help="Target porosity value")
    ap.add_argument("-r", "--rate", type=float, default=1, help="Contraction rate")
    ap.add_argument("-s", "--seed", type=int, default=1, help="Seed")
    # ap.add_argument("--seedmax", type=int, help="Max seed value")
    ap.add_argument("-t","--tolerance", type=float, default=1e-4, help="Tolerance of absolute error.")
    args = vars(ap.parse_args())

    error = 1;
    porosities = []
    while (abs(error) > args['tolerance']):
        por, error = genpack(args['rate'], args['porosity'])
        porosities.append(por)
        print ("[{0:04d}]: Avg: {1}, Max: {2}, Min: {3}\n".format(len(porosities), sum(porosities)/len(porosities), max(porosities), min(porosities)))

    ## Ideally we'd use something like scipy.minimize. But seeing how random seeds affect the final porosity, this is a decent preliminary method to get results.
    ## TODO: I should talk to Bill about methods like MCMC and monte carlo

if __name__ == "__main__":
    main()
