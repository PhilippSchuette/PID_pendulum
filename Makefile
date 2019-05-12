# Make a local package for import in api/app.py.
all:
	python3 setup.py bdist_wheel
	sudo pip3 install dist/PID_pendulum-0.0.1-py3-none-any.whl
