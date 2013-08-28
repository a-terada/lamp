all: lcm53
test: 
	python -m unittest discover
lcm53: lcm53.zip
	unzip lcm53.zip -d lcm53
	pushd ./lcm53 > /dev/null && make && popd > /dev/null
clean:
	rm -rf lcm53 *.pyc frepattern/*.pyc functions/*.pyc
