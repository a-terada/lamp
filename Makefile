all: lcm53/lcm
test: 
	python -m unittest discover
lcm53/lcm: 
	cd lcm53 && make && cd ..
clean:
	rm lcm53/lcm *.pyc frepattern/*.pyc functions/*.pyc
