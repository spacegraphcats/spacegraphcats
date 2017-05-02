all:

lint:
	mypy spacegraphcats/*.py --ignore-missing-imports

test:
	py.test spacegraphcats
