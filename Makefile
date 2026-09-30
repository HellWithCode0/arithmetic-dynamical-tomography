.PHONY: test reproduce

test:
	python -m pytest -q

reproduce:
	adt-reproduce collision-65-119 --output results/reproduced/collision_65_119.json
	adt-reproduce prime-separation --bound 100000 --output results/reproduced/prime_separation_100000.json
	adt-reproduce squarefree-search --bound 100000 --output results/reproduced/squarefree_100000.json
	adt-reproduce prime-power-depths --prime-bound 100 --probes -3 -2 -1 0 1 2 3 --output results/reproduced/prime_power_depths.csv
