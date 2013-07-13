lcm53: lcm53.zip
	mkdir lcm53
	tar zxvf lcm53.zip -C lcm53
	pushd ./lcm53 > /dev/null && make && popd > /dev/null
clean:
	rm -rf lcm53 *.pyc frepattern/*.pyc functions/*.pyc
