from tecto_rna import TectoRNA, PARAMS_TLR
import numpy as np
from scipy.spatial.distance import mahalanobis
from sklearn.covariance import EmpiricalCovariance
import argparse

###### Load input arguments from command line ######
parser = argparse.ArgumentParser(
    description='Perform simulation to compute statbility of tecto RNA.')
parser.add_argument(
    '-seq1', type=str, nargs=2, help="Helix seqeuence of first molecule")
parser.add_argument(
    '-seq2', type=str, nargs=2, default=None,
    help="Helix seqeuence of second molecule. "
    "Assumed to be same as seq1 is not input.")
parser.add_argument(
    '-n_cycles', type=int, default=10000,
    help="Number of cycles for the simulation.")
parser.add_argument(
    '-bp_par_prefix', type=str, default="helix",
    help="Prefix for output base-pair parameters.")

args = parser.parse_args()

n_cycles = args.n_cycles
seq1 = args.seq1
seq2 = args.seq2
if seq2 is None:
    seq2 = seq1

###### Actual simulation #######
tecto = TectoRNA(seq1, seq2)
params = np.empty((n_cycles, 6))

precision = np.load('precision_test.npy')  # Needed for finding best conformer
lowest_dist = float('inf')
lowest_pose = []  # Best pose, saved for latter output
for i in xrange(n_cycles):
    params_new = tecto.update_get_connection()
    params[i] = params_new

    # Find the distance to ideal tetraloop-receptor params
    diff = PARAMS_TLR - params_new
    # Make sure the diff angles are within (-pi, pi]
    for i in xrange(3, 6):  # index 3,4,5 are angles, others are distances
        if diff[i] > np.pi:
            diff[i] -= 2 * np.pi
        elif diff[i] <= -np.pi:
            diff[i] += 2 * np.pi
    dist = mahalanobis(diff, np.zeros(6), precision)
    if dist < lowest_dist:
        lowest_dist = dist
        lowest_pose = [tecto.pose1.copy(), tecto.pose2.copy()]

###### Likelyhood Computation ######
# Fold the angles in params into proper range, such that
# they centered at the mean.
N_CYCLE_FOLD_ANGLE = 10
for j in xrange(N_CYCLE_FOLD_ANGLE):
    mean = np.mean(params, axis=0)
    for i in xrange(3, 6):  # index 3,4,5 are angles, others are distances
        params[:, i][params[:, i] > mean[i] + np.pi] -= 2 * np.pi
        params[:, i][params[:, i] < mean[i] - np.pi] += 2 * np.pi
        if PARAMS_TLR[i] > mean[i] + np.pi:
            PARAMS_TLR[i] += 2 * np.pi
        if PARAMS_TLR[i] < mean[i] - np.pi:
            PARAMS_TLR[i] -= 2 * np.pi

est = EmpiricalCovariance(True, False)
est.fit(params)
log_likelyhood = est.score(PARAMS_TLR[None, :])
KT = 0.59
free_e = -log_likelyhood * KT

print 'Log likelyhood score:', log_likelyhood
print 'Free energy:', free_e


###### Output the best conformer to pdb ######
def generate_bp_par_file(params, bps, out_name):
    assert(len(params) == len(bps))
    n_bp = len(params)
    # convert from radians to degrees
    params[:, 3:] = np.degrees(params[:, 3:])

    out_str = "%4d # base-pairs\n" % n_bp
    out_str += "   0 # ***local base-pair & step parameters***\n"
    out_str += (
        "#        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening     "
        "Shift     Slide     Rise      Tilt      Roll      Twist\n")
    for i in xrange(n_bp):
        out_str += "%s-%s " % tuple(bps[i])
        # First 6 entries are 0
        for j in xrange(6):
            out_str += "%10.3f" % 0
        for j in xrange(6):
            out_str += "%10.3f" % params[i, j]
        out_str += '\n'
    with open(out_name, 'w') as out:
        out.write(out_str)

# Process the parameters and sequence
bps1 = []
for base1, base2 in zip(seq1[1], reversed(seq1[0])):
    bps1.append(base1 + base2)

bps2 = []
for base1, base2 in zip(seq2[0], reversed(seq2[1])):
    bps2.append(base1 + base2)
params1 = np.vstack([
    np.zeros(6),
    lowest_pose[0].params])
params2 = np.vstack([
    PARAMS_TLR,
    lowest_pose[1].params])

# Save bp params files to disk
FILE1 = "%s1_bp.par" % args.bp_par_prefix
FILE2 = "%s2_bp.par" % args.bp_par_prefix
generate_bp_par_file(params1, bps1, FILE1)
generate_bp_par_file(params2, bps2, FILE2)
