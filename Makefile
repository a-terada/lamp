all: lcm53
test: 
	python -m unittest discover
lcm53: lcm53.zip
	unzip lcm53.zip -d lcm53
	cd lcm53 && make && cd ..
clean:
	rm -rf lcm53 *.pyc frepattern/*.pyc functions/*.pyc
