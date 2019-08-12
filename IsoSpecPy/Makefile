# This Makefile is not meant for end-users, only as a convenience for developers. Do not touch.
reinstall:
	rm -rf IsoSpecPy.egg-info build dist
	pip uninstall -y IsoSpecPy ; pip uninstall -y IsoSpecPy ; pip uninstall -y IsoSpecPy || true
	python setup.py install --user
