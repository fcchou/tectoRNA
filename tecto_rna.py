import numpy as np
from helixmc import util, pose, random_step


# Frame for Tetraloop-Receptor, compute from P4-P6 ()
ORIG1 = np.array([51.4917, 76.1789, 99.1397])
FRAME1 = np.array([
    [0.2324, 0.8820, -0.4101],
    [0.6346, -0.4570, -0.6233],
    [-0.7371, -0.1154, -0.6659]])

ORIG2 = np.array([60.2665, 65.1550, 99.1096])
FRAME2 = np.array([
    [0.5198, -0.8456, 0.1216],
    [-0.4753, -0.1679, 0.8637],
    [-0.7099, -0.5067, -0.4892]])

# Parameters for tetraloop to receptor
# Order: Shift, Slide, Rise, Tilt, Roll, Twist
PARAMS_TLR = util.frames2params_3dna(
    ORIG1, ORIG2, FRAME1, FRAME2)


class TectoRNA(object):
    def __init__(self, seq1=None, seq2=None):
        """
        Tecto RNA object for simulating tecto-RNA conformations.

        Illustration:

        Receptor  ---- Tetraloop
          A--U            A--U
          G--C            U--A
          C--G            C--G
          A--U            C--G
        Tetraloop ---- Receptor

        =>
        seq1 = ('ACGA', 'UCGU')
        seq2 = ('UAGG', 'CCUA')

        Parameters
        ----------
        seq1 : (str, str)
            Sequence of the first helix (extend out from the tetra-loop).
        seq2 : (str, str)
            Sequence of the second helix (extend out from the tetra-loop).
        """
        self.params_tlr = PARAMS_TLR
        self.o_tlr, self.r_tlr = util.params2coords(PARAMS_TLR)
        self.rand_step = random_step.RandomStepAgg('RNA_2.8_all_wgu.npz')
        self._symmerize_rand_step()

        self._validate_input_seq(seq1)
        self._seq1 = seq1
        self.pose1 = self._get_pose(self._seq1)

        # For seq2, reverse the order so that it extends out from the
        # receptor.
        self._validate_input_seq(seq2)
        self._seq2 = (seq2[0], seq2[1])
        self.pose2 = self._get_pose(self._seq2)

    def update(self):
        """
        Update the system, generate a new snapshot of
        the two-heliex construct.
        """
        length = len(self._seq1[0])
        for i in xrange(length-1):
            seq = (
                self._seq1[0][i:(i+2)] +
                self._seq1[1][(length-2-i):(length-i)])
            p, o, R = self.rand_step(seq)
            self.pose1.update(i, p, o, R)

        length = len(self._seq2[0])
        for i in xrange(length-1):
            seq = (
                self._seq2[0][i:(i+2)] +
                self._seq2[1][(length-2-i):(length-i)])
            p, o, R = self.rand_step(seq)
            self.pose2.update(i, p, o, R)

    def get_connection(self):
        "Get the 6D parameters connecting the two helix ends on the top."
        o1 = self.pose1.coord_terminal
        f1 = self.pose1.frame_terminal
        o2 = self.pose2.coord_terminal
        f2 = self.pose2.frame_terminal

        o2 = self.r_tlr.dot(o2) + self.o_tlr
        f2 = self.r_tlr.dot(f2)

        f1[:, 1:] *= -1
        f2[:, 1:] *= -1

        params = util.frames2params(o2, o1, f2, f1)
        return params

    def update_get_connection(self):
        "Update the system then get the connection parameters."
        self.update()
        return self.get_connection()

    def _symmerize_rand_step(self):
        "Symmetrize the randome step generator."
        name_corr = dict()
        for name in self.rand_step._names:
            name_corr[name] = name[2:] + name[:2]
        new_params = {}
        for name in self.rand_step._names:
            name_corr = name[2:] + name[:2]
            params1 = self.rand_step.get_rand_step(name).params
            params2 = self.rand_step.get_rand_step(name_corr).params
            params2[:, 0] *= -1
            params2[:, 3] *= -1
            new_params[name] = np.vstack((params1, params2))

        for key in new_params:
            self.rand_step.get_rand_step(key).params = new_params[key]

    def _get_pose(self, seq):
        self._validate_input_seq(seq)
        params = np.tile(self.rand_step.params_avg, (len(seq[0])-1, 1))
        return pose.HelixPose(params)

    def _validate_input_seq(self, seq):
        if len(seq) != 2:
            raise ValueError(
                "Input sequence is not a tuple of length 2! "
                "E.g. ('AAGG', 'CCUU').")
        if len(seq[0]) != len(seq[1]):
            raise ValueError(
                "The two sequence strings in the sequence tuple "
                "is not of equal length!")
