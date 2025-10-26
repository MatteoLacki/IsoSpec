"""Here we store formulas for approximating numbers of isotopologues and subisotopologues."""
from math import exp, log, lgamma, pi

from .IsoSpecPy import IsoParamsFromFormula


def log_multinomial_confs_cnt(n, i):
	"""Return the natural logarithm of the number of configurations of the multinomial distribution.

	There are n independent trials of an experiment that results in one of the i outcomes.
	For each outcome, count how many times it occured.
	These counts then follow the multinomial distribution.

	Args:
		n (int): The number of trials.
		i (int): The number of possible outcomes.
	"""
	return lgamma(n + i) - \
		   lgamma(n + 1) - \
		   lgamma(i)

def test_log_multinomial_confs_cnt():
	x = log_multinomial_confs_cnt(10,3)
	assert abs(x - 4.189655) < 10**(-5)


def multinomial_confs_cnt(n, i):
	"""Return the approximate number of configurations of the multinomial distribution.

	Args:
		n (int): The number of trials.
		i (int): The number of possible outcomes.
	"""
	return exp(log_subiso_cnt_simplex(n, i))


def log_V_simplex(n, i):
	"""Get the natural logarithm of the volume of a simplex {(x_1,..,x_{i-1}): \sum_{j=1}^i = n}.

	Args:
		n (int): The number of atoms of the element.
		i (int): The number of isotopes of the element.
	"""
	return (i-1)*log(n) - lgamma(i)


def test_log_V_simplex():
	x = log_V_simplex(10,3)
	assert abs(x - 3.912023) < 10**(-5)


def V_simplex(n, i):
	"""Get the volume of a simplex {(x_1,..,x_{i-1}): \sum_{j=1}^i = n}.

	Args:
		n (int): The number of atoms of the element.
		i (int): The number of isotopes of the element.
	"""
	return exp(log_V_simplex(n,i))

def log_V_ellipsoid(n, R2, probs):
	"""Get the natural logarithm of the volume of the ellipsoid.

	The ellipsoid is defined by x' W^{-1} x <= R2,
	where W = diag(probs[0:-1]) - probs[0:-1] * probs[0:-1]'
	and R2 is the square of the radius.
	diag(probs) is a matrix with values probs[0:-1] on the diagonal,
	and probs[0:-1] * probs[0:-1]' is a projection on probs[0:-1].
	Since probabilities are nonzero and sum to one, then det W != 0.
	Also, the expression does not really depend upon the choice of one of the
	ommited probability term, e.g. the outcome would stay the same if we remover p[1].

	Args:
		n (int): The number of atoms of the element.
		R2 (float): The square of the radius of the ellipsoid.
		probs (list): List with the natural frequencies of isotopes.
	"""
	i = len(probs)
	log_det = sum(log(p) for p in probs) 
	return (log_det + (i-1)*(log(n) + log(R2) + log(pi)))/2.0 - lgamma((i+1)/2.0)


def test_log_V_ellipsoid():
	assert abs(log_V_ellipsoid(100, 10, [.2,.3,.5]) - 6.299206) < 10**(-5)


def V_ellipsoid(n, R2, probs):
	"""Get the volume of the ellipsoid."""
	return exp(log_V_ellipsoid(n, R2, probs))

def log_subisotopologue_cnt(atoms_cnt, isotope_frequencies, ellipsoid_R2):
	"""Get the natural logarithm of the approximate number of subisotopologues.

	Args:
		atoms_cnt (int): The number of atoms of the given element.
		isotope_frequencies (list): The natural frequencies of isotopes (sum to one).
		ellipsoid_R2 (float): The radius of the ellipsoid used to approximate the optimal P-set.
	"""
	isotopes_cnt = len(isotope_frequencies)
	return log_multinomial_confs_cnt(atoms_cnt, isotopes_cnt) + \
		   log_V_ellipsoid(atoms_cnt, ellipsoid_R2, isotope_frequencies) - \
		   log_V_simplex(atoms_cnt, isotopes_cnt)

def test_log_subisotopologue_cnt():
	assert abs(log_subisotopologue_cnt(100, [.2,.3,.5], 10) - 6.328959) < 10**(-5)


def subisotopologue_cnt(atoms_cnt, isotope_frequencies, ellipsoid_R2):
	"""Get the approximate number of subisotopologues."""
	return exp(log_subisotopologue_cnt(atoms_cnt, isotope_frequencies, ellipsoid_R2))


def approximate_subisotopologues(molecule, P):
	"""Approximate the number of subisotopologues.

	Args:
		molecule (str): A string with molecule, e.g. 'C100H202'.
		P (float): The joint probability threshold.
	"""
	from scipy.stats import chi2
	assert P >= 0 and P <= 1, 'That is not a probability.'
	mol = IsoParamsFromFormula(molecule)
	chi2_df = sum(len(p) for p in mol.probs) - len(mol.probs)
	R2 = chi2.ppf(q=P, df=chi2_df)
	return {e: subisotopologue_cnt(n, p, R2) for n, p, e in \
			zip(mol.atomCounts, mol.probs, mol.elems)}


if __name__ == '__main__':
	print(approximate_subisotopologues('C100H202', .999))
	print(approximate_subisotopologues("C100H100", .999))
