# This Makefile is not meant for end-users, only as a convenience for developers. Do not touch.
reinstall: clean
	pip uninstall -y IsoSpecPy ; pip uninstall -y IsoSpecPy ; pip uninstall -y IsoSpecPy || true
	python setup.py install --user

reinstallall: reinstall2 reinstall3
rea: reinstallall

reinstall2: clean
	pip2 uninstall -y IsoSpecPy ; pip2 uninstall -y IsoSpecPy ; pip2 uninstall -y IsoSpecPy || true
	python2 setup.py install --user
reinstall3: clean
	pip3 uninstall -y IsoSpecPy ; pip3 uninstall -y IsoSpecPy ; pip3 uninstall -y IsoSpecPy || true
	python3 setup.py install --user
clean:
	rm -rf IsoSpecPy.egg-info build dist
twine: clean
	python setup.py sdist
	twine upload dist/*
