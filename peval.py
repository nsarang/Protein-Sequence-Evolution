import os
import argparse
import subprocess
import xgboost as xgb
import numpy as np

EXE_PATH = "bin/gendata"
XGB_MODEL = "data/xgbreg.bin"


def get_args():
    parser = argparse.ArgumentParser(
        description="Protein Quality Assesment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(dest="pdb", help="Path to protein pdb file")
    parser.add_argument(
        "-d",
        dest="profdir",
        default="/tmp/profdir",
        help="Path to directory to load/save profiles",
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.7")
    args = parser.parse_args()
    return args


def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)
     

def get_features(params):
    for line in execute(["./bin/gendata", *params]):
        if line == "\n":
            pass
        elif line.startswith("OUTPUT"):
            return list(map(float, line.split()[3:-1]))
        else:
            print(line, end="")


def get_prediction(X):
    X = xgb.DMatrix(np.array(X)[None, :])
    bst = xgb.Booster({'nthread': 4})  # init model
    bst.load_model(XGB_MODEL)  # load data

    return bst.predict(X)[0]

if __name__ == "__main__":
    args = get_args()

    pdb_path = os.path.normpath(args.pdb)
    if not os.path.exists(pdb_path):
        raise ValueError("PDB file not found.")

    profile_dir = os.path.join(os.path.normpath(args.profdir), "")
    os.makedirs(profile_dir, exist_ok=True)

    features = get_features((pdb_path, profile_dir))
    print("Protein Score:", get_prediction(features))